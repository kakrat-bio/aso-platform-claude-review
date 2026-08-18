"""E12 — The ML ranker against the platform's own heuristic, head to head.

The platform already ranks candidates. `gene_silencing_service` tiles a
target transcript and sorts the windows by `_composite_score`, which is
0.65 * (normalised ViennaRNA target-duplex dG) + 0.35 * (chemistry-adjusted
Tm fit). The question is whether the learned ranker beats it.

Comparing a de-novo generator against a tiling designer directly has no
ground truth — two sequence sets, no measurements. So the comparison is run
where ground truth exists: on held-out benchmark experiments, where dozens
of real oligos against one target have measured activities. Those oligos ARE
tiling candidates. Scoring them with the platform's own functions, and with
the learned ranker, and comparing both against the measurements, is a fair
head-to-head of two *selection policies* over the same candidate pool.

Reported in two parts:

**Ranking quality.** The same NDCG/top-k/Pearson/MRR the rest of the suite
uses, for the heuristic score, the learned score, and drawn-random.

**Selection outcome.** The number a user actually experiences: if you take
the top 5 (or 10) candidates each policy proposes, what is the mean measured
activity of what you got, and how often does the true best oligo appear?

Two facts about the heuristic that come out of running it, and which are
findings in their own right:

* `CHEM_TM_BOOST` has entries for `gapmer`, `lna_gapmer`, `pmo` and `2ome`.
  It has none for MOE or cEt — the two chemistries that make up essentially
  the whole RNase-H arm and every splice-switching row. Those get a boost of
  0, so the Tm-fit term is computed from the unmodified-DNA Tm.
* `OPTIMAL_TM_RANGES` is keyed by platform mechanism id (A1, A2, A12, A15),
  which has no entry for siRNA. The siRNA arm therefore falls through to the
  default (50, 70) window and its heuristic column should be read as
  "the heuristic applied outside its intended scope", not as a fair test of
  it.
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
from backend.services.gene_silencing_service import (  # noqa: E402
    CHEM_TM_BOOST, OPTIMAL_TM_RANGES, _calc_tm, _composite_score,
    _target_duplex_energy, _tm_fit_score,
)

NAME = "exp12_ml_vs_heuristic"
HEADLINE = ["ndcg@10", "top10", "pearson", "mrr"]
SELECT_K = (5, 10)

# The platform mechanism whose design window each benchmark arm corresponds
# to. siRNA has none — the platform's ASO windows do not cover it.
MECHANISM_FOR = {
    "rnase_h": "A1",                 # RNase H1 gapmer, 50-65 C
    "splice_switching": "A2",        # steric block, 62-75 C
    "sirna": None,                   # falls through to the (50, 70) default
}
COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


def revcomp_rna(seq: str) -> str:
    return "".join(COMPLEMENT.get(c, "N") for c in reversed(seq.upper()))


def heuristic_score(seq: str, mechanism_id: str | None) -> float:
    """The platform's composite score for a perfect-match candidate.

    A tiled candidate's target region is its own reverse complement, which is
    what `generate_candidates` passes to `_target_duplex_energy`; the same
    thing is reconstructed here because the benchmark stores the oligo but
    not the transcript it came from.
    """
    dna = seq.upper().replace("U", "T")
    tm = _calc_tm(dna)
    dg = _target_duplex_energy(seq, revcomp_rna(seq))
    fit = _tm_fit_score(tm, "gapmer", [], mechanism_id or "")
    return _composite_score(dg, fit)


def selection_outcome(df_test: pd.DataFrame, score: np.ndarray,
                      k: int) -> dict:
    """What you get if you take this policy's top k, per experiment group."""
    d = df_test.assign(_s=score)
    rng = np.random.default_rng(C.SEED)
    means, best_hits, sizes = [], [], []
    for _, g in d.groupby("experiment_id"):
        if len(g) < max(k, C.MIN_GROUP):
            continue
        # Same shuffle-then-sort as `common.per_group_metrics`: the heuristic
        # ties 27% of oligos within an experiment, and a stable sort would
        # resolve those ties by transcript position rather than arbitrarily.
        perm = rng.permutation(len(g))
        yt = g["rank_label"].to_numpy(dtype=float)[perm]
        order = np.argsort(-g["_s"].to_numpy(dtype=float)[perm])[:k]
        means.append(float(yt[order].mean()))
        best_hits.append(1.0 if int(np.argmax(yt)) in set(order) else 0.0)
        sizes.append(len(g))
    if not means:
        return {"n_groups": 0}
    return {
        "n_groups": len(means),
        "mean_percentile_of_selected": round(float(np.mean(means)), 3),
        "best_oligo_recovered_rate": round(float(np.mean(best_hits)), 4),
        "mean_group_size": round(float(np.mean(sizes)), 1),
        "note": "mean_percentile_of_selected is the mean within-experiment "
                "activity percentile of the k chosen oligos; picking at "
                "random gives ~50",
    }


def main(epochs: int = 8) -> dict:
    df = C.add_chemistry_class(C.load_benchmark())
    split = C.grouped_split(df)
    train = df[split.train].reset_index(drop=True)
    test = df[split.test].reset_index(drop=True)
    chem_vocab = {c: i for i, c in enumerate(sorted(df["chemistry"].unique()))}

    model, _ = train_or_load(train, TrainConfig(mode="conditioned", epochs=epochs,
                                                modality_weights={}, verbose=True),
                             chem_vocab, tag="e1-none")
    ml = predict(model, test, chem_vocab)

    print(f"[{NAME}] scoring {len(test)} held-out oligos with the platform "
          f"heuristic", flush=True)
    heur = np.array([heuristic_score(s, MECHANISM_FOR.get(m))
                     for s, m in zip(test["seq"], test["modality"])],
                    dtype=float)

    scorers = {
        "platform_heuristic": heur,
        "ml_ranker": ml,
        "random": np.random.default_rng(C.SEED).random(len(test)),
    }

    results, table, pgs = {}, [], {}
    for name, s in scorers.items():
        pg = C.per_group_metrics(test, s)
        pgs[name] = pg.set_index("experiment_id")
        overall = C.summarise(pg)
        overall["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(pg, "ndcg@10")
        by_mod = {}
        for mod, sub in pg.groupby("modality"):
            m = C.summarise(sub)
            m["ci95_ndcg@10"] = C.bootstrap_ci_over_genes(sub, "ndcg@10")
            by_mod[mod] = m
        selection = {}
        for k in SELECT_K:
            selection[f"top{k}"] = selection_outcome(test, s, k)
            for mod in sorted(test["modality"].unique()):
                sub = test[test["modality"] == mod]
                selection[f"top{k}_{mod}"] = selection_outcome(
                    sub, s[np.asarray(test["modality"] == mod)], k)
        results[name] = {"overall": overall, "by_modality": by_mod,
                         "selection": selection}
        table.append({"scorer": name, "n_groups": overall["n_groups"],
                      **{m: round(overall.get(m, np.nan), 4) for m in HEADLINE},
                      "top5_mean_percentile":
                          selection["top5"].get("mean_percentile_of_selected"),
                      "top5_best_recovered":
                          selection["top5"].get("best_oligo_recovered_rate")})

    ref = pgs["platform_heuristic"]
    paired = {}
    for name in ("ml_ranker", "random"):
        cur = pgs[name]
        shared = ref.index.intersection(cur.index)
        paired[f"{name}_vs_platform_heuristic"] = {
            "n_paired_groups": int(len(shared)),
            **{m: C.paired_bootstrap_p(cur.loc[shared, m].to_numpy(),
                                       ref.loc[shared, m].to_numpy())
               for m in HEADLINE},
        }

    # Do the two policies even agree on which oligos are good?
    agree = []
    d = test.assign(_h=heur, _m=ml)
    for _, g in d.groupby("experiment_id"):
        if len(g) < C.MIN_GROUP:
            continue
        h, m = g["_h"].to_numpy(), g["_m"].to_numpy()
        if np.std(h) > 0 and np.std(m) > 0:
            agree.append(np.corrcoef(h, m)[0, 1])

    payload = {
        "experiment": NAME,
        "question": "Does the learned ranker beat the platform's existing "
                    "composite heuristic on real held-out measurements?",
        "protocol": {
            "heuristic": "backend.services.gene_silencing_service."
                         "_composite_score, called with the platform's own "
                         "_calc_tm / _target_duplex_energy / _tm_fit_score",
            "candidate_pool": "held-out benchmark oligos; the same pool is "
                              "ranked by both policies, so this compares "
                              "SELECTION, not sequence generation",
            "tie_breaking": "random, seeded. The heuristic ties 27% of oligos "
                            "within an experiment (it rounds to one decimal "
                            "and 54% of oligos sit on the tm_fit plateau at "
                            "100), and a stable sort would have broken those "
                            "ties by transcript position.",
            "split": split.kind, "split_note": split.note, "epochs": epochs,
            "mechanism_window_used": MECHANISM_FOR,
        },
        "heuristic_coverage_gaps": {
            "chem_tm_boost_keys": sorted(CHEM_TM_BOOST),
            "optimal_tm_range_keys": sorted(OPTIMAL_TM_RANGES),
            "finding": "CHEM_TM_BOOST has no MOE or cEt entry, and those two "
                       "chemistries are essentially the entire RNase-H arm and "
                       "all of splice-switching. Every benchmark oligo "
                       "therefore gets a Tm boost of 0 and its Tm-fit term is "
                       "computed from the unmodified-DNA Tm.",
            "sirna_out_of_scope": "OPTIMAL_TM_RANGES has no siRNA key; that arm "
                                  "falls through to the default (50, 70) "
                                  "window, so its heuristic row is the "
                                  "heuristic applied outside its intended "
                                  "scope.",
        },
        "table": table,
        "results": results,
        "paired_vs_heuristic": paired,
        "policy_agreement": {
            "n_groups": len(agree),
            "mean_within_group_pearson": (round(float(np.mean(agree)), 4)
                                          if agree else None),
            "note": "correlation between the two policies' scores inside each "
                    "experiment; near zero means they are ranking on "
                    "essentially unrelated grounds",
        },
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).to_string(index=False))
    return payload


if __name__ == "__main__":
    main(epochs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
