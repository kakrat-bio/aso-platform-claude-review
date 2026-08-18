"""E1 — Does upweighting the scarce modalities change the ranker?

The unified benchmark is 96.3% RNase-H by row (159,109 / 165,237) and 98.3%
by experiment group (1,941 / 1,975). If the pooled ranker is really an
RNase-H ranker wearing a pooled label, then paying the scarce modalities
more attention during training should move their held-out numbers.

Four weightings, one training protocol, one held-out set:

    none          every pair weight 1.0 (the current benchmark behaviour)
    inverse_rows  weight_m proportional to 1 / rows_m
    inverse_exps  weight_m proportional to 1 / experiment_groups_m
    sqrt_inverse  square root of inverse_rows, the usual compromise

The weight attaches to the PAIR, because the loss is pairwise within an
experiment group and a single row carries no loss of its own.

Read the output as a difference, not as a leaderboard: the question is
whether the scarce-modality columns move, and whether RNase-H pays for it.
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

NAME = "exp01_class_imbalance"
HEADLINE = ["ndcg@10", "top10", "pearson", "mrr"]


def weight_schemes(train: pd.DataFrame) -> dict[str, dict]:
    rows = train.groupby("modality").size()
    exps = train.groupby("modality")["experiment_id"].nunique()

    def norm(s):
        w = s / s.min()          # scarcest modality gets weight 1 before invert
        w = 1.0 / w
        return (w / w.min()).to_dict()

    inv_rows = {k: float(v) for k, v in norm(rows).items()}
    inv_exps = {k: float(v) for k, v in norm(exps).items()}
    sqrt_rows = {k: float(np.sqrt(v)) for k, v in inv_rows.items()}
    return {
        "none": {},
        "inverse_rows": inv_rows,
        "inverse_exps": inv_exps,
        "sqrt_inverse": sqrt_rows,
    }


def main(epochs: int = 10) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)

    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}
    schemes = weight_schemes(train)

    composition = {
        "rows": train.groupby("modality").size().to_dict(),
        "experiment_groups": train.groupby("modality")["experiment_id"]
                                  .nunique().to_dict(),
        "row_share": (train.groupby("modality").size() / len(train))
                     .round(4).to_dict(),
    }

    per_scheme, per_group_store = {}, {}
    for scheme, weights in schemes.items():
        print(f"[{NAME}] training scheme={scheme} weights={weights}", flush=True)
        model, _ = train_or_load(
            train,
            TrainConfig(mode="conditioned", epochs=epochs,
                        modality_weights=weights, verbose=True),
            chem_vocab=chem_vocab, tag=f"e1-{scheme}",
        )
        pred = predict(model, test, chem_vocab)
        pg = C.per_group_metrics(test, pred)
        per_group_store[scheme] = pg
        by_mod = {}
        for mod, sub in pg.groupby("modality"):
            s = C.summarise(sub)
            s["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(sub, "ndcg@10")
            by_mod[mod] = s
        per_scheme[scheme] = {
            "weights": weights,
            "overall": C.summarise(pg),
            "by_modality": by_mod,
        }

    # Paired against the unweighted baseline, matched on experiment_id so a
    # group's intrinsic difficulty cancels.
    base = per_group_store["none"].set_index("experiment_id")
    tests = {}
    for scheme, pg in per_group_store.items():
        if scheme == "none":
            continue
        cur = pg.set_index("experiment_id")
        shared = base.index.intersection(cur.index)
        per_mod = {}
        for mod in sorted(base.loc[shared, "modality"].unique()):
            m = base.loc[shared, "modality"] == mod
            ids = shared[m]
            per_mod[mod] = {
                metric: C.paired_bootstrap_p(cur.loc[ids, metric].to_numpy(),
                                             base.loc[ids, metric].to_numpy())
                for metric in HEADLINE
            }
        tests[scheme] = per_mod

    table = []
    for scheme, res in per_scheme.items():
        for mod, s in res["by_modality"].items():
            table.append({
                "scheme": scheme, "modality": mod,
                "n_groups": s["n_groups"],
                **{m: round(s.get(m, float("nan")), 4) for m in HEADLINE},
            })

    # How much of a change could this design have detected at all?
    #
    # The held-out set is built by holding out whole groups, and the scarce
    # arms have very few groups to hold out: siRNA has 11 experiments in
    # total and splice-switching has 5 genes. So the arms this experiment is
    # ABOUT are the arms with the least evaluation power, and a "no change"
    # result on them has to be read as "no change large enough to see with
    # this many groups", not as "no change".
    #
    # The half-width of the paired-bootstrap CI on the mean difference is a
    # rough resolution limit: an effect much smaller than it is invisible
    # here. It is NOT a formal minimum detectable effect, and on an arm with
    # a handful of groups the interval is itself estimated from a handful of
    # numbers — so a narrow half-width there is not evidence of precision.
    power = {}
    base_pg = per_group_store["none"]
    for mod in sorted(test["modality"].unique()):
        sub = base_pg[base_pg["modality"] == mod]
        widths = []
        for scheme, per_mod in tests.items():
            b = per_mod.get(mod, {}).get("ndcg@10", {})
            ci = b.get("ci95")
            if ci and not any(np.isnan(ci)):
                widths.append((ci[1] - ci[0]) / 2)
        spread = float(sub["ndcg@10"].std()) if len(sub) > 1 else float("nan")
        power[mod] = {
            "held_out_experiment_groups": int(len(sub)),
            "held_out_genes": int(sub["gene"].nunique()) if len(sub) else 0,
            "total_groups_in_benchmark": int(
                df.loc[df["modality"] == mod, "experiment_id"].nunique()),
            "between_group_sd_of_ndcg@10": (round(spread, 4)
                                            if spread == spread else None),
            "median_ci_half_width_ndcg@10": (
                round(float(np.median(widths)), 4) if widths else None),
            "reading": ("adequately powered" if len(sub) >= 20 else
                        f"UNDERPOWERED: {len(sub)} held-out groups. A null "
                        f"result on this arm means no effect large enough to "
                        f"see, not no effect. The interval itself is also "
                        f"unreliable here: a paired bootstrap over "
                        f"{len(sub)} values resamples {len(sub)} numbers, so "
                        f"it can come out narrow by accident and must not be "
                        f"read as a tight bound. Treat every p-value on this "
                        f"arm as indicative only."),
        }

    payload = {
        "experiment": NAME,
        "question": "Does upweighting siRNA / splice-switching pairs change "
                    "held-out ranking quality?",
        "protocol": {
            "split": split.kind,
            "split_note": split.note,
            "epochs": epochs,
            "model": "conditioned conv ranker, pairwise margin loss",
            "weight_unit": "pair (the loss has no per-row term)",
        },
        "training_composition": composition,
        "table": table,
        "evaluation_power": power,
        "results": per_scheme,
        "paired_vs_none": tests,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 10)
