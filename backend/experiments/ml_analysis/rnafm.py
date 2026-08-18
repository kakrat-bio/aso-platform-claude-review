"""RNA-FM embeddings, loaded without the broken tokenizer path.

Two things go wrong with the obvious `from_pretrained` route on this
machine, and both fail QUIETLY, which is why this module exists:

1. `RnaTokenizer.from_pretrained` raises under transformers 4.57 —
   `extra_special_tokens` is a list in the published tokenizer config and
   the current transformers expects a mapping. Loud failure, easy to spot.

2. `RnaFmModel.from_pretrained` "succeeds" and returns a model with every
   single weight randomly initialised. The published checkpoint stores the
   base model under a `model.` prefix (it is an `RnaFmForPreTraining`
   checkpoint) while `RnaFmModel.base_model_prefix` is `rnafm`, so no key
   matches and transformers silently initialises all 205 tensors. The only
   sign is a warning in a wall of text. Embedding anything with that model
   and calling the result "RNA-FM" would be reporting a randomly
   initialised transformer as a foundation model.

So: the tokenizer is reimplemented from the published `vocab.txt` (it is
character-level, 28 tokens), and the checkpoint is loaded explicitly with
the prefix stripped, with `assert_loaded()` verifying that the encoder
tensors actually arrived. The pooler is absent from the checkpoint — it is
a pre-training head, not part of the encoder — so embeddings are mean-pooled
over real tokens rather than read from `pooler_output`.
"""

from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import torch

MODEL_ID = "multimolecule/rnafm"
SPECIALS = ["<pad>", "<cls>", "<eos>", "<unk>", "<mask>", "<null>"]
PAD, CLS, EOS, UNK = 0, 1, 2, 3


def _snapshot() -> Path:
    hits = glob.glob(str(Path.home() / ".cache/huggingface/hub/"
                         "models--multimolecule--rnafm/snapshots/*/vocab.txt"))
    if not hits:
        raise FileNotFoundError(
            "RNA-FM snapshot not in the HuggingFace cache. Fetch it with "
            "`python -c \"from multimolecule import RnaFmModel; "
            "RnaFmModel.from_pretrained('multimolecule/rnafm')\"` (the call "
            "warns about uninitialised weights — that is the bug this module "
            "works around — but it does download the files)."
        )
    return Path(hits[0]).parent


def load_vocab() -> dict[str, int]:
    lines = (_snapshot() / "vocab.txt").read_text().splitlines()
    return {tok: i for i, tok in enumerate(lines)}


def encode(seqs, vocab: dict[str, int], max_len: int | None = None):
    """Character-level tokenisation: <cls> seq <eos>, T mapped to U."""
    toks = [[CLS] + [vocab.get(c, UNK) for c in s.upper().replace("T", "U")]
            + [EOS] for s in seqs]
    L = max_len or max(len(t) for t in toks)
    ids = np.full((len(toks), L), PAD, dtype=np.int64)
    mask = np.zeros((len(toks), L), dtype=np.int64)
    for i, t in enumerate(toks):
        t = t[:L]
        ids[i, :len(t)] = t
        mask[i, :len(t)] = 1
    return torch.from_numpy(ids), torch.from_numpy(mask)


def load_model():
    """RnaFmModel with the published encoder weights actually in it."""
    from multimolecule import RnaFmConfig, RnaFmModel
    from safetensors.torch import load_file

    cfg = RnaFmConfig.from_pretrained(MODEL_ID)
    model = RnaFmModel(cfg)
    sd = load_file(str(_snapshot() / "model.safetensors"))
    stripped = {k[len("model."):]: v for k, v in sd.items()
                if k.startswith("model.")}
    report = model.load_state_dict(stripped, strict=False)
    assert_loaded(report, stripped)
    model.eval()
    return model, report


def assert_loaded(report, stripped: dict) -> None:
    """Fail loudly if the encoder did not actually receive pretrained weights.

    Only the pooler may be missing: it is a pre-training head and is not in
    the published checkpoint. Anything else missing means the prefix
    handling broke again and the embeddings would be noise.
    """
    unexpected = list(report.unexpected_keys)
    missing = [k for k in report.missing_keys if not k.startswith("pooler.")]
    if missing or unexpected:
        raise RuntimeError(
            f"RNA-FM weights did not load cleanly: {len(missing)} missing "
            f"(non-pooler) {missing[:5]}, {len(unexpected)} unexpected "
            f"{unexpected[:5]}. Refusing to emit embeddings from a partly "
            f"random model."
        )
    n_enc = sum(1 for k in stripped if k.startswith("encoder.layer."))
    if n_enc < 100:
        raise RuntimeError(
            f"only {n_enc} encoder tensors found in the checkpoint; expected "
            f"~192 for a 12-layer model")


@torch.no_grad()
def embed(seqs, batch: int = 256, verbose: bool = True) -> np.ndarray:
    """Mean-pooled final-layer embeddings over real tokens. (n, 640).

    Deduplicated first — the benchmark repeats sequences across experiments,
    and embedding each copy would multiply the cost for nothing.
    """
    model, _ = load_model()
    vocab = load_vocab()
    uniq = list(dict.fromkeys(seqs))
    index = {s: i for i, s in enumerate(uniq)}
    out = np.zeros((len(uniq), model.config.hidden_size), dtype=np.float32)
    for start in range(0, len(uniq), batch):
        chunk = uniq[start:start + batch]
        ids, mask = encode(chunk, vocab)
        h = model(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).float()
        # <cls> and <eos> are excluded from the mean: they carry no
        # nucleotide and would dilute short oligos more than long ones.
        m[:, 0, :] = 0
        for i, s in enumerate(chunk):
            m[i, len(s) + 1, :] = 0
        pooled = (h * m).sum(1) / m.sum(1).clamp(min=1)
        out[start:start + len(chunk)] = pooled.numpy()
        if verbose and (start // batch) % 50 == 0:
            print(f"  embedded {start + len(chunk)}/{len(uniq)}", flush=True)
    return out[[index[s] for s in seqs]]
