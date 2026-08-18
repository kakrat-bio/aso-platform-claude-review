"""E11 — An MLP ranker on RNA-FM embeddings.

RNA-FM (Chen et al. 2022) is an RNA foundation model: a 12-layer, 640-wide
transformer pretrained on non-coding RNA. The question here is whether its
representation carries anything about oligo activity that a small
supervised conv net trained on this benchmark does not already have.

Protocol. Sequences are embedded once (mean-pooled final layer over real
nucleotide tokens, `<cls>`/`<eos>` excluded), then a 2-layer MLP is trained
with the SAME pairwise margin loss, the SAME grouped split and the SAME
metrics as every other ranker in this suite. Three references on identical
data: the conv ranker trained from scratch, a ridge regression on 4-mer
counts, and drawn-random scores.

**RNA-FM is frozen.** Only the MLP head is trained. Fine-tuning the encoder
on 30k oligos of 12-28 nt would be a different experiment and a much larger
one; this measures the pretrained representation as it ships.

**The data is subsampled, by gene.** Embedding all 165,237 sequences takes
about 70 minutes of CPU; the RNase-H arm is cut to whole genes totalling
~30k rows while siRNA and splice-switching are kept whole. Every model in
this experiment sees exactly the same subsample, so the comparison is fair
even though the absolute numbers are not comparable to E1-E7, which use the
full benchmark. That is stated on every row rather than left to be noticed.

One caveat about the encoder itself: RNA-FM was pretrained on non-coding RNA
of typical length in the hundreds of nucleotides. These oligos are 12-28 nt.
A representation can be excellent on its training distribution and carry
little for inputs an order of magnitude shorter, and a null result here is
evidence about the transfer, not about RNA-FM.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C, rnafm  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    BATCH_PAIRS, MIN_EXP_ROWS, PAIRS_PER_EXP, TrainConfig, predict,
    train_or_load,
)

NAME = "exp11_rnafm_mlp"
RNASE_H_ROW_BUDGET = 30_000
EPOCHS_MLP = 30
EMB_CACHE = C.OUT_DIR / "rnafm_embeddings.npz"
HEADLINE = ["ndcg@10", "top10", "pearson", "mrr"]


def subsample(df: pd.DataFrame, seed: int = C.SEED) -> pd.DataFrame:
    """Whole genes from RNase-H up to the row budget; other arms kept whole."""
    rng = np.random.default_rng(seed)
    keep = [df[df["modality"] != "rnase_h"]]
    rn = df[df["modality"] == "rnase_h"]
    genes = rn["target_gene"].unique()
    genes = genes[rng.permutation(len(genes))]
    sizes = rn.groupby("target_gene").size()
    chosen, total = [], 0
    for g in genes:
        if total >= RNASE_H_ROW_BUDGET:
            break
        chosen.append(g)
        total += int(sizes[g])
    keep.append(rn[rn["target_gene"].isin(chosen)])
    return pd.concat(keep).reset_index(drop=True)


class MLPRanker(nn.Module):
    def __init__(self, in_d: int, hidden: int = 256, p_drop: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(in_d),
            nn.Linear(in_d, hidden), nn.ReLU(), nn.Dropout(p_drop),
            nn.Linear(hidden, hidden // 2), nn.ReLU(), nn.Dropout(p_drop),
            nn.Linear(hidden // 2, 1),
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_mlp(X: np.ndarray, df: pd.DataFrame, epochs: int = EPOCHS_MLP,
              lr: float = 1e-3, seed: int = 0) -> MLPRanker:
    """Same pairwise margin objective and pair sampling as the conv ranker."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    groups = {e: i for e, i in df.groupby("experiment_id").indices.items()
              if len(i) >= MIN_EXP_ROWS}
    rank = df["rank_label"].to_numpy(dtype=np.float32)
    Xt = torch.from_numpy(X)

    model = MLPRanker(X.shape[1])
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    loss_fn = nn.MarginRankingLoss(margin=0.1)
    exp_ids = list(groups)
    for ep in range(epochs):
        model.train()
        rng.shuffle(exp_ids)
        ai, bi = [], []
        for e in exp_ids:
            idx = groups[e]
            m = len(idx)
            n_pairs = min(PAIRS_PER_EXP, m * (m - 1) // 2)
            a = rng.integers(0, m, n_pairs)
            b = (a + rng.integers(1, m, n_pairs)) % m
            ai.append(idx[a])
            bi.append(idx[b])
        ai, bi = np.concatenate(ai), np.concatenate(bi)
        perm = rng.permutation(len(ai))
        ai, bi = ai[perm], bi[perm]
        total = 0.0
        for s in range(0, len(ai), BATCH_PAIRS):
            a, b = ai[s:s + BATCH_PAIRS], bi[s:s + BATCH_PAIRS]
            sign = torch.from_numpy(
                np.where(rank[a] >= rank[b], 1.0, -1.0).astype(np.float32))
            loss = loss_fn(model(Xt[a]), model(Xt[b]), sign)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss.detach()) * len(a)
        if (ep + 1) % 10 == 0:
            print(f"  mlp epoch {ep + 1}/{epochs} loss={total / len(ai):.4f}",
                  flush=True)
    model.eval()
    return model


@torch.no_grad()
def mlp_predict(model, X: np.ndarray) -> np.ndarray:
    return model(torch.from_numpy(X)).numpy()


def ridge_scores(Xtr, ytr, Xte, alpha: float = 1.0) -> np.ndarray:
    """Closed-form ridge on centred features — no extra dependency."""
    mu = Xtr.mean(0)
    A = Xtr - mu
    y = ytr - ytr.mean()
    w = np.linalg.solve(A.T @ A + alpha * np.eye(A.shape[1]), A.T @ y)
    return (Xte - mu) @ w


def main(epochs_conv: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    sub = subsample(df)
    split = C.grouped_split(sub)
    train = sub[split.train].reset_index(drop=True)
    test = sub[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    seqs = sub["seq"].tolist()
    if EMB_CACHE.exists():
        z = np.load(EMB_CACHE, allow_pickle=True)
        if list(z["seqs"]) == seqs:
            X = z["X"]
            print(f"[{NAME}] embeddings from cache", flush=True)
        else:
            X = None
    else:
        X = None
    if X is None:
        print(f"[{NAME}] embedding {len(seqs)} sequences with RNA-FM "
              f"(frozen)", flush=True)
        X = rnafm.embed(seqs)
        C.OUT_DIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(EMB_CACHE, X=X, seqs=np.array(seqs, dtype=object))

    Xtr, Xte = X[split.train], X[split.test]

    scorers: dict[str, np.ndarray] = {}

    print(f"[{NAME}] MLP on RNA-FM embeddings", flush=True)
    mlp = train_mlp(Xtr, train)
    scorers["rnafm_mlp"] = mlp_predict(mlp, Xte)

    print(f"[{NAME}] conv ranker on the same subsample", flush=True)
    conv, _ = train_or_load(train, TrainConfig(mode="conditioned",
                                               epochs=epochs_conv, verbose=True),
                            chem_vocab, tag="e11-conv-subsample")
    scorers["conv_ranker"] = predict(conv, test, chem_vocab)

    print(f"[{NAME}] ridge on 4-mer counts", flush=True)
    Ktr = C.build_features(train, k=4, with_chemistry=False)
    Kte = C.build_features(test, k=4, with_chemistry=False)
    scorers["kmer_ridge"] = ridge_scores(
        Ktr, train["rank_label"].to_numpy(dtype=float), Kte)

    print(f"[{NAME}] ridge on RNA-FM embeddings (linear probe)", flush=True)
    scorers["rnafm_ridge"] = ridge_scores(
        Xtr, train["rank_label"].to_numpy(dtype=float), Xte)

    scorers["random"] = np.random.default_rng(C.SEED).random(len(test))

    results, table = {}, []
    for name, s in scorers.items():
        pg = C.per_group_metrics(test, s)
        overall = C.summarise(pg)
        overall["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(pg, "ndcg@10")
        by_mod = {m: C.summarise(g) for m, g in pg.groupby("modality")}
        results[name] = {"overall": overall, "by_modality": by_mod,
                         "_pg": pg.set_index("experiment_id")}
        table.append({"scorer": name, "n_groups": overall["n_groups"],
                      **{m: round(overall.get(m, np.nan), 4) for m in HEADLINE}})

    # Paired against the conv ranker on identical held-out groups.
    ref = results["conv_ranker"]["_pg"]
    paired = {}
    for name, r in results.items():
        if name == "conv_ranker":
            continue
        cur = r["_pg"]
        shared = ref.index.intersection(cur.index)
        paired[f"{name}_vs_conv_ranker"] = {
            m: C.paired_bootstrap_p(cur.loc[shared, m].to_numpy(),
                                    ref.loc[shared, m].to_numpy())
            for m in HEADLINE
        }
    for r in results.values():
        r.pop("_pg")

    payload = {
        "experiment": NAME,
        "question": "Does a frozen RNA-FM representation carry oligo-activity "
                    "signal a small supervised conv net does not?",
        "protocol": {
            "encoder": "multimolecule/rnafm, 12 layers x 640, FROZEN; "
                       "mean-pooled final layer over nucleotide tokens",
            "weight_loading": "checkpoint loaded explicitly with the `model.` "
                              "prefix stripped and verified; from_pretrained "
                              "silently returns a randomly initialised model "
                              "on this transformers version",
            "head": "2-layer MLP, pairwise margin loss, same pair sampling as "
                    "the conv ranker",
            "mlp_epochs": EPOCHS_MLP,
            "conv_epochs": epochs_conv,
            "split": split.kind, "split_note": split.note,
            "subsample": {
                "rnase_h_row_budget": RNASE_H_ROW_BUDGET,
                "selection": "whole genes, so the gene split stays honest",
                "rows": int(len(sub)),
                "by_modality": sub.groupby("modality").size().to_dict(),
                "warning": "absolute numbers here are NOT comparable to E1-E7, "
                           "which use the full benchmark; only the comparison "
                           "between scorers in this table is meaningful",
            },
        },
        "caveats": [
            "RNA-FM was pretrained on non-coding RNA hundreds of nucleotides "
            "long; these oligos are 12-28 nt. A null result is evidence about "
            "transfer to very short sequences, not about RNA-FM in general.",
            "The encoder is frozen. Fine-tuning it is a different experiment.",
        ],
        "table": table,
        "results": results,
        "paired_vs_conv_ranker": paired,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs_conv=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
