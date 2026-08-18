"""E6 — Is the conditioned model significantly better than the alternatives?

Three model classes, identical data, identical split, identical training
budget:

    seqonly       sequence only; chemistry is not an input
    conditioned   chemistry embedding concatenated into the ranking head
    invariant     domain-adversarial; chemistry is removed from the
                  representation via gradient reversal, and is not an input
                  at inference

Every pair is compared on the SAME held-out experiment groups, paired by
`experiment_id`, so a group's intrinsic difficulty cancels. Two tests are
reported because neither alone is convincing on this data: a Wilcoxon
signed-rank test (distribution-free, but assumes symmetric differences) and
a paired bootstrap on the mean difference (no symmetry assumption).

Multiplicity: three model classes give three pairwise comparisons per
metric. Holm-Bonferroni adjusted p-values are reported alongside the raw
ones. A raw p below 0.05 that does not survive Holm is reported as not
significant.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

try:
    from scipy.stats import wilcoxon
except ImportError:                                  # reported, never faked
    wilcoxon = None

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    TrainConfig, predict, train_or_load,
)

NAME = "exp06_significance"
MODES = ["seqonly", "conditioned", "invariant"]
METRICS = ["ndcg@10", "top10", "pearson", "mrr"]


def holm(pvals: dict[str, float]) -> dict[str, float]:
    """Holm-Bonferroni step-down adjustment, monotonicity enforced."""
    items = [(k, v) for k, v in pvals.items() if not np.isnan(v)]
    items.sort(key=lambda kv: kv[1])
    n, out, running = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        adj = min(1.0, (n - i) * p)
        running = max(running, adj)
        out[k] = running
    for k, v in pvals.items():
        out.setdefault(k, float("nan"))
    return out


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    per_group, summaries = {}, {}
    for mode in MODES:
        print(f"[{NAME}] training mode={mode}", flush=True)
        tag = "e1-none" if mode == "conditioned" else f"e6-{mode}"
        cfg = (TrainConfig(mode=mode, epochs=epochs, modality_weights={},
                           verbose=True)
               if mode == "conditioned"
               else TrainConfig(mode=mode, epochs=epochs, verbose=True))
        model, _ = train_or_load(train, cfg, chem_vocab, tag=tag)
        pg = C.per_group_metrics(test, predict(model, test, chem_vocab))
        per_group[mode] = pg.set_index("experiment_id")
        s = C.summarise(pg)
        s["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(pg, "ndcg@10")
        s["by_modality"] = {m: C.summarise(sub)
                            for m, sub in pg.groupby("modality")}
        summaries[mode] = s

    comparisons, raw_p = {}, {}
    for i, a in enumerate(MODES):
        for b in MODES[i + 1:]:
            A, B = per_group[a], per_group[b]
            shared = A.index.intersection(B.index)
            key = f"{a}_vs_{b}"
            comparisons[key] = {"n_paired_groups": int(len(shared))}
            for m in METRICS:
                x = A.loc[shared, m].to_numpy(dtype=float)
                y = B.loc[shared, m].to_numpy(dtype=float)
                ok = ~(np.isnan(x) | np.isnan(y))
                boot = C.paired_bootstrap_p(x, y)
                if wilcoxon is not None and ok.sum() >= 10 and np.any(x[ok] != y[ok]):
                    stat, p = wilcoxon(x[ok], y[ok])
                    w = {"statistic": float(stat), "p": float(p)}
                else:
                    w = {"statistic": None, "p": float("nan"),
                         "unavailable": ("scipy not installed" if wilcoxon is None
                                         else "fewer than 10 usable pairs or "
                                              "all differences zero")}
                comparisons[key][m] = {"paired_bootstrap": boot, "wilcoxon": w}
                raw_p[f"{key}|{m}"] = w["p"]

    adjusted = holm(raw_p)
    for key in comparisons:
        for m in METRICS:
            p_raw = raw_p.get(f"{key}|{m}", float("nan"))
            p_adj = adjusted.get(f"{key}|{m}", float("nan"))
            comparisons[key][m]["wilcoxon"]["p_holm"] = p_adj
            comparisons[key][m]["significant_at_0.05_after_holm"] = bool(
                not np.isnan(p_adj) and p_adj < 0.05)

    table = []
    for mode in MODES:
        s = summaries[mode]
        table.append({"model": mode, "n_groups": s["n_groups"],
                      **{m: round(s[m], 4) for m in METRICS}})

    payload = {
        "experiment": NAME,
        "question": "Is chemistry conditioning significantly better than "
                    "sequence-only or chemistry-invariant ranking?",
        "protocol": {"split": split.kind, "split_note": split.note,
                     "epochs": epochs,
                     "pairing": "by experiment_id; identical held-out groups",
                     "multiplicity": "Holm-Bonferroni over all "
                                     f"{len(raw_p)} comparisons",
                     "scipy_available": wilcoxon is not None},
        "table": table,
        "summaries": summaries,
        "comparisons": comparisons,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
