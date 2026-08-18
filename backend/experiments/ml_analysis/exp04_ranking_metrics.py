"""E4 — The standard ranking metrics, at every k, with intervals.

NDCG@k, MAP@k, MRR and top-k overlap for k in {1, 3, 5, 10, 20}, computed
within each held-out experiment group and averaged over groups, for the
pooled ranker and for two reference points that make the numbers readable:

    random     scores drawn uniformly at random, fresh per group. This is
               the floor. Permuting a constant vector is NOT a random
               baseline — it leaves every tie intact and scores far too
               well — so the floor is drawn, not shuffled.
    oracle     scores equal to the labels. This is the ceiling, and it is
               not 1.0 for every metric: with ties in `rank_label`, top-k
               overlap and MAP@k cannot reach 1.

Reading NDCG here needs one caution. `rank_label` is a within-experiment
percentile in [0, 100], so it is dense and positive, and NDCG's gain is
never zero for any item. That inflates NDCG relative to a
binary-relevance setting — hence NDCG@1 near 0.68 while top-1 overlap is
near 0.05. Both are correct; they answer different questions. Top-k
overlap and MRR are the honest headline for "did we surface the good
ones", and the oracle/random columns bracket everything.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    TrainConfig, predict, train_or_load,
)

NAME = "exp04_ranking_metrics"


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

    rng = np.random.default_rng(C.SEED)
    scored = {
        "model": pred,
        "random": rng.random(len(test)),
        "oracle": test["rank_label"].to_numpy(dtype=float),
    }

    results, table = {}, []
    for name, s in scored.items():
        pg = C.per_group_metrics(test, s)
        overall = C.summarise(pg)
        by_mod = {}
        for mod, sub in pg.groupby("modality"):
            m = C.summarise(sub)
            m["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(sub, "ndcg@10")
            m["ci95_top10"] = C.bootstrap_ci_over_genes(sub, "top10")
            by_mod[mod] = m
        overall["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(pg, "ndcg@10")
        overall["ci95_top10"] = C.bootstrap_ci_over_genes(pg, "top10")
        overall["ci95_mrr"] = C.bootstrap_ci_over_genes(pg, "mrr")
        results[name] = {"overall": overall, "by_modality": by_mod}
        for k in C.K_VALUES:
            table.append({
                "scorer": name, "k": k,
                "ndcg@k": round(overall[f"ndcg@{k}"], 4),
                "map@k": round(overall[f"map@{k}"], 4),
                "top-k overlap": round(overall[f"top{k}"], 4),
                "mrr": round(overall["mrr"], 4),
            })

    payload = {
        "experiment": NAME,
        "question": "Standard ranking metrics at multiple k, bracketed by a "
                    "drawn-random floor and a label oracle ceiling.",
        "protocol": {
            "split": split.kind, "split_note": split.note, "epochs": epochs,
            "k_values": list(C.K_VALUES),
            "grouping": "metrics computed within experiment_id, then averaged "
                        "over groups with >= %d rows" % C.MIN_GROUP,
            "ci": "percentile bootstrap resampling GENES, not rows",
        },
        "interpretation": {
            "ndcg_is_inflated": "rank_label is a dense percentile in [0,100], "
                                "so every item has positive gain and NDCG "
                                "cannot fall near zero; compare against the "
                                "random row of this table, not against 0.",
            "oracle_below_one": "ties in rank_label cap top-k overlap and "
                                "MAP@k below 1.0 even for a perfect scorer.",
        },
        "table": table,
        "results": results,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
