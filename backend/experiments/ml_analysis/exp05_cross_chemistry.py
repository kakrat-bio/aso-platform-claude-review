"""E5 — Cross-chemistry transfer matrix: train on X, test on Y.

The 228 raw chemistry strings in the benchmark are structured annotations
("L13 MOE|sugar|1,2,3,11,12,13 PS|backbone|...") describing which positions
carry which modification. Treating them as 228 categories makes a transfer
matrix that is mostly empty cells. They collapse to a handful of chemistry
CLASSES by the modifications present, which is the level the question is
actually about.

Restricted to RNase-H on purpose. Across the whole benchmark, chemistry and
mechanism are almost perfectly confounded — every splice-switching row is
MOE, every siRNA row is unmodified — so a pooled transfer matrix would be
measuring mechanism transfer and calling it chemistry. Within RNase-H,
chemistry varies while mechanism is fixed, which is the only place the
question is answerable.

Diagonal cells are same-chemistry generalisation (still a gene split);
off-diagonal cells are transfer. Reading the matrix by ROW tells you what a
model trained on that chemistry can do elsewhere; reading by COLUMN tells
you how hard that chemistry is to predict.

The models here are SEQUENCE-ONLY, unlike everywhere else in the suite. A
chemistry-conditioned model trained on MOE alone has never updated the
embedding for cEt, so scoring cEt with it would measure a randomly
initialised embedding vector and report the result as failed transfer. Once
chemistry is removed as an input, an off-diagonal cell measures the one
thing the question is about: whether sequence preferences learned under one
chemistry hold under another.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    TrainConfig, predict, train_or_load,
)

NAME = "exp05_cross_chemistry"
MIN_ROWS = 2000          # below this a chemistry cannot train its own model
MIN_TEST_GROUPS = 5
METRIC = "ndcg@10"


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    rn = df[df["modality"] == "rnase_h"].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    counts = rn.groupby("chem_class").size().sort_values(ascending=False)
    usable = [c for c, n in counts.items() if n >= MIN_ROWS]
    excluded = {c: int(n) for c, n in counts.items() if n < MIN_ROWS}

    # One gene split for the whole arm, so every cell of the matrix holds out
    # the SAME genes. Splitting per chemistry would let a gene sit in one
    # cell's training set and another cell's test set.
    split = C.gene_split(rn)
    train_all = rn[split.train].reset_index(drop=True)
    test_all = rn[split.test].reset_index(drop=True)

    models = {}
    for c in usable:
        tr = train_all[train_all["chem_class"] == c].reset_index(drop=True)
        if tr["experiment_id"].nunique() < 2:
            continue
        print(f"[{NAME}] training on chem_class={c} ({len(tr)} rows)", flush=True)
        models[c], _ = train_or_load(
            tr, TrainConfig(mode="seqonly", epochs=epochs, verbose=True),
            chem_vocab, tag=f"e5-seqonly-{c.replace('+', '_')}")

    matrix, cells = {}, []
    for tr_c, model in models.items():
        matrix[tr_c] = {}
        for te_c in usable:
            te = test_all[test_all["chem_class"] == te_c].reset_index(drop=True)
            if te.empty:
                matrix[tr_c][te_c] = None
                continue
            pg = C.per_group_metrics(te, predict(model, te, chem_vocab))
            if len(pg) < MIN_TEST_GROUPS:
                matrix[tr_c][te_c] = {
                    "value": None,
                    "unavailable": f"{len(pg)} held-out groups, need "
                                   f"{MIN_TEST_GROUPS}",
                }
                continue
            s = C.summarise(pg)
            matrix[tr_c][te_c] = {
                "value": round(s[METRIC], 4),
                "top10": round(s["top10"], 4),
                "pearson": round(s["pearson"], 4),
                "mrr": round(s["mrr"], 4),
                "n_groups": s["n_groups"],
                "n_genes": s["n_genes"],
                "ci95": C.bootstrap_ci_over_genes(pg, METRIC),
            }
            cells.append({
                "train_chem": tr_c, "test_chem": te_c,
                "diagonal": tr_c == te_c,
                METRIC: matrix[tr_c][te_c]["value"],
                "top10": matrix[tr_c][te_c]["top10"],
                "n_groups": s["n_groups"],
            })

    cell_df = pd.DataFrame(cells)
    summary = {}
    if not cell_df.empty:
        diag = cell_df[cell_df["diagonal"]][METRIC]
        off = cell_df[~cell_df["diagonal"]][METRIC]
        summary = {
            "mean_diagonal": round(float(diag.mean()), 4) if len(diag) else None,
            "mean_off_diagonal": round(float(off.mean()), 4) if len(off) else None,
            "transfer_penalty": (round(float(diag.mean() - off.mean()), 4)
                                 if len(diag) and len(off) else None),
        }

    payload = {
        "experiment": NAME,
        "question": "Does a ranker trained on one chemistry class rank a "
                    "different chemistry?",
        "protocol": {
            "restricted_to": "rnase_h",
            "why_restricted": "chemistry and mechanism are confounded across "
                              "the benchmark (all splice_switching rows are "
                              "MOE, all siRNA rows unmodified); only inside "
                              "rnase_h does chemistry vary at fixed mechanism",
            "split": "one gene split shared by every cell",
            "model": "SEQUENCE-ONLY conv ranker; chemistry is not an input, so an off-diagonal cell is not measuring an untrained chemistry embedding",
            "epochs": epochs,
            "metric": METRIC,
            "min_rows_to_train": MIN_ROWS,
        },
        "chemistry_class_counts": {c: int(n) for c, n in counts.items()},
        "excluded_too_small": excluded,
        "matrix": matrix,
        "cells": cells,
        "summary": summary,
    }
    C.write_result(NAME, payload)
    if not cell_df.empty:
        print(cell_df.pivot(index="train_chem", columns="test_chem",
                            values=METRIC).to_string())
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
