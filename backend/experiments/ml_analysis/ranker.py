"""One ranker, trained one way, so every experiment in this suite is
comparing the same thing.

The architecture is the benchmark's `InvariantRanker` (conv encoder ->
pairwise ranking head, optionally chemistry-conditioned or
chemistry-adversarial). What this module adds is a per-pair WEIGHT, which
is what the class-imbalance experiment needs and what the benchmark's
trainer has no way to express.

Why weight pairs rather than rows: the loss is pairwise within an
experiment group, so a row has no loss of its own. The unit that carries
gradient is the pair, so that is the unit a weight has to attach to.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parents[2].parent))

from backend.experiments.benchmark.invariant_ranker import (  # noqa: E402
    InvariantRanker, seq_to_onehot, MAX_LEN,
)

MARGIN = 0.1
PAIRS_PER_EXP = 64
BATCH_PAIRS = 2048
MIN_EXP_ROWS = 5

CACHE_DIR = Path(__file__).resolve().parent / "results" / "models"

# One epoch over the full benchmark is ~2.6 minutes on this machine, and
# several experiments below want the same trained model. Cache on a hash
# of (config, training rows) so a rerun of the suite is cheap and so two
# experiments quoting "the same model" are quoting the same weights.



@dataclass
class TrainConfig:
    mode: str = "conditioned"          # seqonly | conditioned | invariant
    epochs: int = 12
    lr: float = 1e-3
    pairs_per_exp: int = PAIRS_PER_EXP
    seed: int = 0
    modality_weights: dict | None = None   # modality -> pair weight
    verbose: bool = False


class Encoded:
    """One-hot + chemistry ids for a frame, computed once and reused."""

    def __init__(self, df: pd.DataFrame, chem_vocab: dict[str, int] | None = None):
        self.df = df.reset_index(drop=True)
        self.oh = seq_to_onehot(self.df["seq"].tolist(), MAX_LEN)
        if chem_vocab is None:
            chem_vocab = {c: i for i, c in enumerate(sorted(self.df["chemistry"].unique()))}
        self.chem_vocab = chem_vocab
        self.chem_id = self.df["chemistry"].map(
            lambda c: chem_vocab.get(c, 0)).to_numpy(dtype=np.int64)


def _group_indices(df: pd.DataFrame, min_rows: int = MIN_EXP_ROWS):
    idx = df.groupby("experiment_id").indices
    return {e: i for e, i in idx.items() if len(i) >= min_rows}


def train_ranker(train_df: pd.DataFrame, cfg: TrainConfig,
                 chem_vocab: dict[str, int] | None = None):
    """Train and return (model, encoded_train). Deterministic given cfg.seed."""
    torch.manual_seed(cfg.seed)
    rng = np.random.default_rng(cfg.seed)

    enc = Encoded(train_df, chem_vocab)
    df = enc.df
    groups = _group_indices(df)
    if not groups:
        raise ValueError("no experiment group has enough rows to form pairs")

    rank = df["rank_label"].to_numpy(dtype=np.float32)
    modality = df["modality"].to_numpy()
    weights = cfg.modality_weights or {}

    model = InvariantRanker(cfg.mode, n_chem=max(len(enc.chem_vocab), 1))
    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr)
    loss_fn = nn.MarginRankingLoss(margin=MARGIN, reduction="none")
    ce = nn.CrossEntropyLoss()

    exp_ids = list(groups)
    for epoch in range(cfg.epochs):
        # Build the whole epoch's pairs first, then step on large batches.
        # One optimizer step per experiment (1,941 steps of 64 pairs) spends
        # nearly all its time in Python and autograd overhead rather than in
        # the convolution; the same pairs in batches of BATCH_PAIRS train the
        # same model roughly 30x faster.
        rng.shuffle(exp_ids)
        ai_all, bi_all = [], []
        for e in exp_ids:
            idx = groups[e]
            m = len(idx)
            n_pairs = min(cfg.pairs_per_exp, m * (m - 1) // 2)
            a = rng.integers(0, m, n_pairs)
            b = (a + rng.integers(1, m, n_pairs)) % m
            ai_all.append(idx[a])
            bi_all.append(idx[b])
        ai_all = np.concatenate(ai_all)
        bi_all = np.concatenate(bi_all)
        perm = rng.permutation(len(ai_all))
        ai_all, bi_all = ai_all[perm], bi_all[perm]

        total, seen = 0.0, 0
        for s0 in range(0, len(ai_all), BATCH_PAIRS):
            ai = ai_all[s0:s0 + BATCH_PAIRS]
            bi = bi_all[s0:s0 + BATCH_PAIRS]
            sign = np.where(rank[ai] >= rank[bi], 1.0, -1.0).astype(np.float32)
            w = np.array([weights.get(mm, 1.0) for mm in modality[ai]],
                         dtype=np.float32)

            oh_a = torch.from_numpy(enc.oh[ai])
            oh_b = torch.from_numpy(enc.oh[bi])
            ch_a = torch.from_numpy(enc.chem_id[ai])
            ch_b = torch.from_numpy(enc.chem_id[bi])
            sg = torch.from_numpy(sign)
            wt = torch.from_numpy(w)

            sa = model.score(oh_a, ch_a)
            sb = model.score(oh_b, ch_b)
            per_pair = loss_fn(sa, sb, sg)
            loss = (per_pair * wt).sum() / wt.sum()

            if cfg.mode == "invariant":
                logits = model.chem_logits(torch.cat([oh_a, oh_b]))
                loss = loss + ce(logits, torch.cat([ch_a, ch_b]))

            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss.detach()) * len(ai)
            seen += len(ai)
        if cfg.verbose:
            print(f"  epoch {epoch + 1}/{cfg.epochs} loss={total / max(seen,1):.4f}",
                  flush=True)
    model.eval()
    return model, enc


@torch.no_grad()
def predict(model, df: pd.DataFrame, chem_vocab: dict[str, int],
            batch: int = 4096) -> np.ndarray:
    enc = Encoded(df, chem_vocab)
    out = []
    for s in range(0, len(enc.df), batch):
        oh = torch.from_numpy(enc.oh[s:s + batch])
        ch = torch.from_numpy(enc.chem_id[s:s + batch])
        out.append(model.score(oh, ch).numpy())
    return np.concatenate(out) if out else np.zeros(0, dtype=np.float32)


def _cache_path(train_df: pd.DataFrame, cfg: TrainConfig,
                chem_vocab: dict, tag: str) -> Path:
    h = hashlib.sha256()
    h.update(json.dumps(asdict(cfg), sort_keys=True, default=str).encode())
    h.update(str(len(train_df)).encode())
    h.update(pd.util.hash_pandas_object(train_df["seq"], index=False)
             .values.tobytes())
    h.update(pd.util.hash_pandas_object(train_df["rank_label"], index=False)
             .values.tobytes())
    h.update(json.dumps(sorted(chem_vocab.items())).encode())
    return CACHE_DIR / f"{tag}-{h.hexdigest()[:16]}.pt"


def train_or_load(train_df: pd.DataFrame, cfg: TrainConfig,
                  chem_vocab: dict[str, int] | None = None,
                  tag: str = "ranker"):
    """train_ranker, memoised on (config, training data, vocabulary)."""
    if chem_vocab is None:
        chem_vocab = {c: i for i, c in enumerate(sorted(train_df["chemistry"].unique()))}
    path = _cache_path(train_df, cfg, chem_vocab, tag)
    if path.exists():
        model = InvariantRanker(cfg.mode, n_chem=max(len(chem_vocab), 1))
        model.load_state_dict(torch.load(path, map_location="cpu"))
        model.eval()
        if cfg.verbose:
            print(f"  [cache hit] {path.name}", flush=True)
        return model, chem_vocab
    model, enc = train_ranker(train_df, cfg, chem_vocab)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), path)
    return model, enc.chem_vocab
