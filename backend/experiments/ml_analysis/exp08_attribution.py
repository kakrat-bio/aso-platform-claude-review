"""E8 — Which nucleotide positions does the ranker actually use?

Two attribution methods, deliberately different in what they assume:

**Integrated gradients** (Sundararajan et al. 2017). Path integral of the
gradient from a baseline to the input. The baseline is the all-zero one-hot
vector, which for this encoding means "no nucleotide" rather than an
average nucleotide — the natural absent state for a padded one-hot input,
and the padding the model already sees beyond each sequence's length.
Completeness (attributions sum to score minus baseline score) is checked
and reported, because a violated completeness check means the number of
integration steps was too small and the attributions should not be read.

**Occlusion.** Zero out one position, re-score, take the drop. Model-
agnostic and assumption-free, but it moves the input off the data manifold
and cannot see interactions.

They are reported side by side rather than averaged. Agreement between them
is evidence; disagreement is a result about the model, not an error.

Positions are reported two ways: absolute (5' index) and relative (fraction
along the oligo), because the benchmark mixes 12-28 nt oligos and an
absolute-position profile silently mixes a gapmer's 5' wing with another
oligo's centre. The gapmer-relevant profile — a 5-10-5 architecture — is
only meaningful within a fixed length, so the 20-mers are also profiled
alone.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    Encoded, TrainConfig, train_or_load,
)
from backend.experiments.benchmark.invariant_ranker import MAX_LEN  # noqa: E402

NAME = "exp08_attribution"
NUCLEOTIDES = "ACGU"
IG_STEPS = 64
SAMPLE_PER_MODALITY = 400
COMPLETENESS_TOL = 0.05     # fraction of |score - baseline_score|


def integrated_gradients(model, oh: np.ndarray, chem: np.ndarray,
                         steps: int = IG_STEPS):
    """IG from an all-zero baseline. Returns (attributions, completeness)."""
    x = torch.from_numpy(oh)
    ch = torch.from_numpy(chem)
    baseline = torch.zeros_like(x)
    total = torch.zeros_like(x)
    # Riemann midpoint rule: unbiased at a given step count, unlike the left
    # rule, which systematically under-attributes.
    for i in range(steps):
        alpha = (i + 0.5) / steps
        xi = (baseline + alpha * (x - baseline)).requires_grad_(True)
        s = model.score(xi, ch).sum()
        g, = torch.autograd.grad(s, xi)
        total += g
    attr = (x - baseline) * total / steps
    with torch.no_grad():
        f_x = model.score(x, ch)
        f_0 = model.score(baseline, ch)
    delta = (f_x - f_0).detach().numpy()
    summed = attr.sum(dim=(1, 2)).detach().numpy()
    denom = np.abs(delta) + 1e-9
    completeness = np.abs(summed - delta) / denom
    return attr.detach().numpy(), completeness


@torch.no_grad()
def occlusion(model, oh: np.ndarray, chem: np.ndarray) -> np.ndarray:
    """Drop in score when each position is zeroed. (n, MAX_LEN)."""
    x = torch.from_numpy(oh)
    ch = torch.from_numpy(chem)
    base = model.score(x, ch).numpy()
    out = np.zeros((len(oh), MAX_LEN), dtype=np.float32)
    for p in range(MAX_LEN):
        xp = x.clone()
        xp[:, p, :] = 0.0
        out[:, p] = base - model.score(xp, ch).numpy()
    return out


def profile(values: np.ndarray, lengths: np.ndarray, n_bins: int = 20):
    """Mean |attribution| per relative position, length-normalised."""
    bins = np.zeros(n_bins)
    counts = np.zeros(n_bins)
    for row, L in zip(values, lengths):
        L = int(L)
        if L < 2:
            continue
        for p in range(L):
            b = min(n_bins - 1, int(p / L * n_bins))
            bins[b] += abs(row[p])
            counts[b] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, bins / np.maximum(counts, 1), np.nan)


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    model, _ = train_or_load(train, TrainConfig(mode="conditioned", epochs=epochs,
                                                modality_weights={}, verbose=True),
                             chem_vocab, tag="e1-none")

    rng = np.random.default_rng(C.SEED)
    results = {}
    for mod in sorted(test["modality"].unique()):
        sub = test[test["modality"] == mod]
        if len(sub) > SAMPLE_PER_MODALITY:
            sub = sub.iloc[rng.permutation(len(sub))[:SAMPLE_PER_MODALITY]]
        sub = sub.reset_index(drop=True)
        enc = Encoded(sub, chem_vocab)
        lengths = sub["seq_len"].to_numpy()

        print(f"[{NAME}] {mod}: IG over {len(sub)} sequences", flush=True)
        ig, completeness = integrated_gradients(model, enc.oh, enc.chem_id)
        ig_pos = np.abs(ig).sum(axis=2)          # (n, MAX_LEN)
        print(f"[{NAME}] {mod}: occlusion", flush=True)
        occ = occlusion(model, enc.oh, enc.chem_id)

        bad = float((completeness > COMPLETENESS_TOL).mean())
        # Per-nucleotide signed attribution, over real positions only.
        nt_attr = {}
        for j, nt in enumerate(NUCLEOTIDES):
            vals = []
            for i, L in enumerate(lengths):
                m = enc.oh[i, :int(L), j] > 0
                if m.any():
                    vals.append(ig[i, :int(L), j][m].mean())
            nt_attr[nt] = (round(float(np.mean(vals)), 6) if vals else None)

        entry = {
            "n_sequences": int(len(sub)),
            "ig_completeness": {
                "tolerance": COMPLETENESS_TOL,
                "fraction_violating": round(bad, 4),
                "median_relative_error": round(float(np.median(completeness)), 5),
                "usable": bool(bad < 0.2),
                "note": ("attributions with a violated completeness check "
                         "should not be read; raise IG_STEPS" if bad >= 0.2
                         else "completeness holds at the reported step count"),
            },
            "ig_relative_position_profile": [
                None if np.isnan(v) else round(float(v), 6)
                for v in profile(ig_pos, lengths)],
            "occlusion_relative_position_profile": [
                None if np.isnan(v) else round(float(v), 6)
                for v in profile(np.abs(occ), lengths)],
            "per_nucleotide_mean_ig": nt_attr,
        }

        # Fixed-length profile: only here does "position 3 of a 5-10-5" mean
        # the same thing across sequences.
        for L in sorted(pd.Series(lengths).value_counts().head(2).index):
            m = lengths == L
            if m.sum() < 20:
                continue
            ig_L = ig_pos[m][:, :int(L)].mean(axis=0)
            occ_L = np.abs(occ[m][:, :int(L)]).mean(axis=0)
            entry[f"fixed_length_{int(L)}"] = {
                "n": int(m.sum()),
                "ig_by_position": [round(float(v), 6) for v in ig_L],
                "occlusion_by_position": [round(float(v), 6) for v in occ_L],
                "ig_argmax_position": int(np.argmax(ig_L)),
                "occlusion_argmax_position": int(np.argmax(occ_L)),
            }

        # Do the two methods agree on which positions matter?
        agree = []
        for i, L in enumerate(lengths):
            L = int(L)
            if L < 4:
                continue
            a, b = ig_pos[i, :L], np.abs(occ[i, :L])
            if np.std(a) > 0 and np.std(b) > 0:
                agree.append(np.corrcoef(a, b)[0, 1])
        entry["ig_vs_occlusion_agreement"] = {
            "n": len(agree),
            "mean_pearson_over_positions": (round(float(np.mean(agree)), 4)
                                            if agree else None),
        }
        results[mod] = entry

    payload = {
        "experiment": NAME,
        "question": "Which nucleotide positions drive the ranker's score?",
        "protocol": {
            "methods": ["integrated gradients (midpoint rule, %d steps, "
                        "all-zero baseline)" % IG_STEPS,
                        "occlusion (zero one position, measure score drop)"],
            "split": split.kind, "epochs": epochs,
            "sample_per_modality": SAMPLE_PER_MODALITY,
            "positions": "reported relative (20 bins along the oligo) and, for "
                         "the two commonest lengths, absolute",
        },
        "results": results,
    }
    C.write_result(NAME, payload)
    for mod, r in results.items():
        print(mod, "completeness ok:", r["ig_completeness"]["usable"],
              "| IG/occlusion agreement:",
              r["ig_vs_occlusion_agreement"]["mean_pearson_over_positions"])
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
