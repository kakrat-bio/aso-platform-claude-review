"""E3 — What a random row split buys you that a gene split does not.

Same mechanism, same model, same metrics; the only thing that changes is
whether the held-out rows come from genes the model has never seen or from
genes it has already been trained on.

The gap between the two is the amount of the reported number that is
memorisation of a target rather than sequence knowledge. It is reported per
mechanism because the three arms have very different structure:

    rnase_h           339 genes, 469 rows/gene   -> a real gene split exists
    splice_switching  5 genes                    -> a gene split is 1 gene wide
    sirna             one "gene" per row         -> NO gene split is possible

For siRNA the two protocols are the same protocol, and the row labelled
`gene` is a random split under another name. That is reported as the finding
rather than hidden by running it anyway.
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

NAME = "exp03_within_mechanism_split"
HEADLINE = ["ndcg@10", "top10", "pearson", "mrr"]


def run_one(sub: pd.DataFrame, split: C.Split, chem_vocab: dict,
            tag: str, epochs: int) -> tuple[dict, pd.DataFrame]:
    tr = sub[split.train].reset_index(drop=True)
    te = sub[split.test].reset_index(drop=True)
    model, _ = train_or_load(tr, TrainConfig(mode="conditioned", epochs=epochs,
                                             verbose=True),
                             chem_vocab, tag=tag)
    pg = C.per_group_metrics(te, predict(model, te, chem_vocab))
    s = C.summarise(pg)
    s["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(pg, "ndcg@10")
    s["ci95_pearson"] = C.bootstrap_ci_over_genes(pg, "pearson")
    s["split_kind"] = split.kind
    s["split_note"] = split.note
    s["train_rows"] = int(split.train.sum())
    s["test_rows"] = int(split.test.sum())
    return s, pg


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    out, rows = {}, []
    for mod in sorted(df["modality"].unique()):
        sub = df[df["modality"] == mod].reset_index(drop=True)
        entry: dict = {
            "n_rows": len(sub),
            "n_genes": int(sub["target_gene"].nunique()),
            "rows_per_gene": round(len(sub) / max(sub["target_gene"].nunique(), 1), 2),
            "n_experiments": int(sub["experiment_id"].nunique()),
        }

        print(f"[{NAME}] {mod}: random row split", flush=True)
        rnd_s, rnd_pg = run_one(sub, C.random_row_split(sub), chem_vocab,
                                f"e3-{mod}-random", epochs)
        entry["random_row_split"] = rnd_s

        try:
            gs = C.gene_split(sub)
        except ValueError as exc:
            entry["gene_split"] = None
            entry["gene_split_blocked"] = str(exc)
            entry["gap"] = None
        else:
            print(f"[{NAME}] {mod}: gene split", flush=True)
            gen_s, gen_pg = run_one(sub, gs, chem_vocab, f"e3-{mod}-gene", epochs)
            entry["gene_split"] = gen_s
            entry["gene_split_blocked"] = None
            entry["gap"] = {
                m: round(rnd_s.get(m, np.nan) - gen_s.get(m, np.nan), 4)
                for m in HEADLINE
            }
            # Unpaired: the two protocols hold out different rows, so there is
            # no matched group to pair on. Reported as CIs, not as a p-value.
            entry["note"] = ("random and gene splits hold out different rows, "
                             "so the difference is reported with CIs rather "
                             "than a paired test")
        out[mod] = entry

        rows.append({
            "modality": mod,
            "rows_per_gene": entry["rows_per_gene"],
            **{f"random_{m}": round(rnd_s.get(m, np.nan), 4) for m in HEADLINE},
            **{f"gene_{m}": (round(entry["gene_split"].get(m, np.nan), 4)
                             if entry.get("gene_split") else None)
               for m in HEADLINE},
        })

    payload = {
        "experiment": NAME,
        "question": "Within one mechanism, how much of the reported ranking "
                    "quality is memorising the target gene?",
        "protocol": {"epochs": epochs, "test_frac": 0.25,
                     "model": "conditioned conv ranker, pairwise margin loss",
                     "trained_separately_per_mechanism": True},
        "table": rows,
        "results": out,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(rows).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
