"""E2 — Per-mechanism breakdown of the pooled ranker.

Two things reported side by side.

**Dominance.** How much of the training signal each mechanism actually
contributes: rows, experiment groups, and — the number that matters for a
pairwise loss — training PAIRS per epoch.

**Specialist vs pooled.** For each mechanism, a ranker trained only on that
mechanism, evaluated on the same held-out groups as the pooled ranker. If
pooling helps the scarce mechanisms, the pooled model wins there. If pooling
is RNase-H drowning them out, the specialists win.

A caveat that survives into every conclusion: the splice-switching arm has
5 genes in total, so its held-out set is ONE gene. Its numbers are a single
gene's numbers and are labelled that way rather than averaged into silence.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.ranker import (  # noqa: E402
    MIN_EXP_ROWS, PAIRS_PER_EXP, TrainConfig, predict, train_or_load,
)

NAME = "exp02_per_mechanism"
HEADLINE = ["ndcg@10", "top10", "pearson", "mrr"]


def training_pairs(df: pd.DataFrame) -> dict:
    """Pairs each mechanism actually contributes per epoch."""
    out = {}
    for mod, sub in df.groupby("modality"):
        n = 0
        for _, g in sub.groupby("experiment_id"):
            m = len(g)
            if m >= MIN_EXP_ROWS:
                n += min(PAIRS_PER_EXP, m * (m - 1) // 2)
        out[mod] = int(n)
    return out


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    pairs = training_pairs(train)
    total_pairs = sum(pairs.values())
    dominance = {
        mod: {
            "rows": int((train["modality"] == mod).sum()),
            "row_share": round(float((train["modality"] == mod).mean()), 4),
            "experiment_groups": int(train.loc[train["modality"] == mod,
                                               "experiment_id"].nunique()),
            "training_pairs_per_epoch": pairs.get(mod, 0),
            "pair_share": round(pairs.get(mod, 0) / max(total_pairs, 1), 4),
        }
        for mod in sorted(train["modality"].unique())
    }

    # Same config and tag as E1's unweighted arm, so "the pooled model" means
    # the same weights in both experiments rather than two similar ones.
    cfg = TrainConfig(mode="conditioned", epochs=epochs,
                      modality_weights={}, verbose=True)
    print(f"[{NAME}] pooled model", flush=True)
    pooled, _ = train_or_load(train, cfg, chem_vocab, tag="e1-none")
    pooled_pred = predict(pooled, test, chem_vocab)
    pooled_pg = C.per_group_metrics(test, pooled_pred)

    spec_cfg = TrainConfig(mode="conditioned", epochs=epochs, verbose=True)
    rows, detail = [], {}
    for mod in sorted(train["modality"].unique()):
        tr_m = train[train["modality"] == mod].reset_index(drop=True)
        te_m = test[test["modality"] == mod].reset_index(drop=True)
        if te_m.empty:
            continue

        p_sub = pooled_pg[pooled_pg["modality"] == mod]
        pooled_s = C.summarise(p_sub)
        pooled_s["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(p_sub, "ndcg@10")
        pooled_s["ci95_pearson"] = C.bootstrap_ci_over_genes(p_sub, "pearson")

        try:
            print(f"[{NAME}] specialist model: {mod}", flush=True)
            spec, _ = train_or_load(tr_m, spec_cfg, chem_vocab,
                                    tag=f"e2-spec-{mod}")
            s_pred = predict(spec, te_m, chem_vocab)
            s_pg = C.per_group_metrics(te_m, s_pred)
            spec_s = C.summarise(s_pg)
            spec_s["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(s_pg, "ndcg@10")
            a = s_pg.set_index("experiment_id")
            b = p_sub.set_index("experiment_id")
            shared = a.index.intersection(b.index)
            paired = {
                m: C.paired_bootstrap_p(a.loc[shared, m].to_numpy(),
                                        b.loc[shared, m].to_numpy())
                for m in HEADLINE
            }
            spec_err = None
        except Exception as exc:
            spec_s, paired, spec_err = None, None, f"{type(exc).__name__}: {exc}"

        detail[mod] = {
            "dominance": dominance.get(mod),
            "pooled": pooled_s,
            "specialist": spec_s,
            "specialist_error": spec_err,
            "specialist_minus_pooled": paired,
        }
        rows.append({
            "modality": mod,
            "pair_share": dominance[mod]["pair_share"],
            "test_groups": pooled_s["n_groups"],
            "test_genes": pooled_s["n_genes"],
            **{f"pooled_{m}": round(pooled_s.get(m, np.nan), 4) for m in HEADLINE},
            **{f"spec_{m}": (round(spec_s.get(m, np.nan), 4) if spec_s else None)
               for m in HEADLINE},
        })

    payload = {
        "experiment": NAME,
        "question": "Is the pooled ranker an RNase-H ranker, and do the scarce "
                    "mechanisms do better trained alone?",
        "protocol": {"split": split.kind, "split_note": split.note,
                     "epochs": epochs,
                     "model": "conditioned conv ranker, pairwise margin loss"},
        "dominance": dominance,
        "table": rows,
        "detail": detail,
        "caveats": [
            "splice_switching has 5 genes in total, so its held-out set is a "
            "single gene; read its row as one gene's result, not a mean.",
            "sirna is grouped by experiment_id because its target_gene column "
            "holds one value per row.",
        ],
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(rows).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
