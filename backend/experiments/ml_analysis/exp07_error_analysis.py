"""E7 — Where the ranker fails.

Per held-out experiment group, the model's NDCG@10 is the unit of analysis.
The output ranks the hardest and easiest genes, mechanisms, chemistry
classes and cell lines, and then asks what distinguishes them: group size,
label spread, sequence-length spread, GC spread.

The failure mode worth naming ahead of time is a group with no label
spread. If every oligo in an experiment has nearly the same measured
activity, the *true* order is noise, and no scorer can recover it. Those
groups are reported separately rather than counted as model errors, because
"the model failed" and "there was nothing to find" are different findings.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    TrainConfig, predict, train_or_load,
)

NAME = "exp07_error_analysis"
METRIC = "ndcg@10"
FLAT_LABEL_SD = 5.0        # percentile points; below this the order is noise


def _round_or_none(v, nd: int = 4):
    """NaN is reported as null, not as the string "NaN" in JSON."""
    v = float(v)
    return None if v != v else round(v, nd)


def gc(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / max(len(s), 1)


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    model, _ = train_or_load(train, TrainConfig(mode="conditioned", epochs=epochs,
                                                modality_weights={}, verbose=True),
                             chem_vocab, tag="e1-none")
    pred = predict(model, test, chem_vocab)
    pg = C.per_group_metrics(test, pred)

    # Per-group descriptors, joined onto the metric table.
    test = test.assign(_gc=test["seq"].map(gc))
    desc = test.groupby("experiment_id").agg(
        label_sd=("rank_label", "std"),
        label_range=("rank_label", lambda s: float(s.max() - s.min())),
        len_sd=("seq_len", "std"),
        gc_sd=("_gc", "std"),
        gc_mean=("_gc", "mean"),
        cell_line=("cell_line", "first"),
        source=("source", "first"),
    ).reset_index()
    pg = pg.merge(desc, on="experiment_id", how="left")

    flat = pg[pg["label_sd"].fillna(0) < FLAT_LABEL_SD]
    live = pg[pg["label_sd"].fillna(0) >= FLAT_LABEL_SD]

    def rank_by(col, frame, min_groups=3):
        agg = frame.groupby(col).agg(
            n_groups=(METRIC, "size"),
            mean=(METRIC, "mean"),
            median=(METRIC, "median"),
            mean_top10=("top10", "mean"),
            mean_pearson=("pearson", "mean"),
        ).reset_index()
        agg = agg[agg["n_groups"] >= min_groups].sort_values("mean")
        return agg

    by_gene = rank_by("gene", live)
    by_mod = rank_by("modality", live, min_groups=1)
    by_chem = rank_by("chem_class", live, min_groups=1)
    by_cell = rank_by("cell_line", live)

    # What correlates with a group being hard?
    corr = {}
    for col in ["n", "label_sd", "label_range", "len_sd", "gc_sd", "gc_mean"]:
        v = live[col].to_numpy(dtype=float)
        m = live[METRIC].to_numpy(dtype=float)
        ok = ~(np.isnan(v) | np.isnan(m))
        corr[col] = (round(float(np.corrcoef(v[ok], m[ok])[0, 1]), 4)
                     if ok.sum() > 3 and np.std(v[ok]) > 0 else None)

    worst = live.nsmallest(15, METRIC)[
        ["experiment_id", "gene", "modality", "chem_class", "n",
         METRIC, "top10", "pearson", "label_sd"]]
    best = live.nlargest(15, METRIC)[
        ["experiment_id", "gene", "modality", "chem_class", "n",
         METRIC, "top10", "pearson", "label_sd"]]

    # Sequence-level: what do the misses look like? For each group, compare
    # the oligos the model put in its top 10 to the ones that belonged there.
    miss_gc, hit_gc, miss_len, hit_len = [], [], [], []
    d = test.assign(_pred=pred)
    for _, g in d.groupby("experiment_id"):
        if len(g) < 10:
            continue
        yt = g["rank_label"].to_numpy(dtype=float)
        yp = g["_pred"].to_numpy(dtype=float)
        true_top = set(np.argsort(-yt)[:10])
        pred_top = list(np.argsort(-yp)[:10])
        gcv = g["_gc"].to_numpy()
        lv = g["seq_len"].to_numpy()
        for i in pred_top:
            (hit_gc if i in true_top else miss_gc).append(gcv[i])
            (hit_len if i in true_top else miss_len).append(lv[i])

    def dist(v):
        v = np.asarray(v, dtype=float)
        if not len(v):
            return None
        return {"n": int(len(v)), "mean": round(float(v.mean()), 4),
                "sd": round(float(v.std()), 4)}

    payload = {
        "experiment": NAME,
        "question": "Which genes, mechanisms and chemistries does the ranker "
                    "fail on, and what do the failures look like?",
        "protocol": {"split": split.kind, "split_note": split.note,
                     "epochs": epochs, "metric": METRIC,
                     "flat_label_sd_threshold": FLAT_LABEL_SD},
        "group_counts": {
            "total": int(len(pg)),
            "flat_label_excluded": int(len(flat)),
            "analysed": int(len(live)),
        },
        "flat_label_groups": {
            "note": "label sd below the threshold: the true order is within "
                    "assay noise, so a low score here is not a model error",
            "mean_metric": _round_or_none(flat[METRIC].mean()) if len(flat) else None,
            "examples": flat.nsmallest(10, METRIC)[
                ["experiment_id", "gene", "modality", "n", METRIC, "label_sd"]
            ].to_dict("records") if len(flat) else [],
        },
        "hardest_genes": by_gene.head(15).to_dict("records"),
        "easiest_genes": by_gene.tail(15).iloc[::-1].to_dict("records"),
        "by_modality": by_mod.to_dict("records"),
        "by_chemistry_class": by_chem.to_dict("records"),
        "by_cell_line": by_cell.head(15).to_dict("records"),
        "difficulty_correlations": {
            "note": f"Pearson r between the group descriptor and {METRIC}",
            "values": corr,
        },
        "worst_groups": worst.to_dict("records"),
        "best_groups": best.to_dict("records"),
        "what_misses_look_like": {
            "note": "oligos the model placed in its top 10; 'hit' also belongs "
                    "in the true top 10, 'miss' does not",
            "gc_fraction": {"hit": dist(hit_gc), "miss": dist(miss_gc)},
            "length": {"hit": dist(hit_len), "miss": dist(miss_len)},
        },
    }
    C.write_result(NAME, payload)
    print("hardest genes:\n", by_gene.head(10).to_string(index=False))
    print("\nby modality:\n", by_mod.to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
