"""Shared harness for the ML analysis suite.

Data loading, splitting, chemistry parsing and ranking metrics, in one place
so every experiment measures the same way. Where an experiment departs from
these defaults it says so.

TWO CONVENTIONS THAT MATTER THROUGHOUT
--------------------------------------
**Effective n is unique genes, not rows.** 159,215 RNase-H rows come from 339
genes. Every confidence interval here resamples genes.

**The siRNA arm cannot be gene-split.** Its `target_gene` column holds the
mRNA target site, one per row, so a "gene split" there is a random row split
wearing a label. Experiments that need a real gene split restrict to
RNase-H; experiments that deliberately use a random split say so in their
name and output.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

BACKEND = Path(__file__).resolve().parents[2]
DATA_PATH = BACKEND / "data" / "benchmark" / "unified_benchmark.parquet"
OUT_DIR = Path(__file__).resolve().parent / "results"

SEED = 42
K_VALUES = (1, 3, 5, 10, 20)
MIN_GROUP = 5


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_benchmark(drop_cross_modality: bool = True) -> pd.DataFrame:
    """The unified benchmark, with the known leak removed by default.

    106 sequences appear in both `rnase_h` and `splice_switching`. Any
    comparison across modality is measuring that overlap unless they go.
    """
    df = pd.read_parquet(DATA_PATH)
    if drop_cross_modality:
        n_mod = df.groupby("seq")["modality"].nunique()
        df = df[~df["seq"].isin(set(n_mod[n_mod > 1].index))]
    return df.reset_index(drop=True)


# Chemistry strings are structured, e.g.
#   "L13 MOE|sugar|1,2,3,11,12,13 PS|backbone|1,2,..."
# 228 distinct strings collapse to a handful of chemistry CLASSES, which is
# the level a transfer experiment is actually about.
_CHEM_TOKENS = ("MOE", "cEt", "LNA", "PS", "PMO", "OMe", "F")


def chemistry_class(raw: str) -> str:
    """Collapse a chemistry string to the modifications it actually uses."""
    if not isinstance(raw, str) or not raw.strip():
        return "unspecified"
    found = [t for t in _CHEM_TOKENS if re.search(rf"\b{re.escape(t)}\b", raw)]
    # PS alone is a backbone; it only distinguishes a chemistry when it is
    # the sole modification.
    sugar = [t for t in found if t != "PS"]
    if sugar:
        return "+".join(sorted(sugar))
    return "PS" if "PS" in found else "unmodified"


def add_chemistry_class(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["chem_class"] = df["chemistry"].map(chemistry_class)
    return df


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------

ALPH = "ACGU"


def kmer_features(seqs, k: int = 4) -> np.ndarray:
    n = len(seqs)
    X = np.zeros((n, 4 ** k), dtype=np.float32)
    for i, s in enumerate(seqs):
        row = X[i]
        for j in range(len(s) - k + 1):
            idx, ok = 0, True
            for ch in s[j:j + k]:
                pos = ALPH.find(ch)
                if pos < 0:
                    ok = False
                    break
                idx = idx * 4 + pos
            if ok:
                row[idx] += 1
        row /= row.sum() + 1e-9
    return X


def build_features(df: pd.DataFrame, k: int = 4,
                   with_chemistry: bool = True) -> np.ndarray:
    Xk = kmer_features(df["seq"].tolist(), k)
    if not with_chemistry:
        return Xk.astype(np.float32)
    codes = pd.Categorical(df["chemistry"]).codes
    Xc = np.zeros((len(df), int(codes.max()) + 1), dtype=np.float32)
    Xc[np.arange(len(df)), codes] = 1.0
    return np.hstack([Xk, Xc]).astype(np.float32)


# ---------------------------------------------------------------------------
# Splits
# ---------------------------------------------------------------------------

@dataclass
class Split:
    train: np.ndarray
    test: np.ndarray
    kind: str
    note: str = ""
    n_train_genes: int = 0
    n_test_genes: int = 0


def gene_split(df: pd.DataFrame, test_frac: float = 0.25,
               seed: int = SEED) -> Split:
    """Hold out whole genes. Refuses when the column does not group."""
    genes = df["target_gene"].values
    uniq = np.unique(genes)
    rows_per_gene = len(df) / max(len(uniq), 1)
    if rows_per_gene < 2.0:
        raise ValueError(
            f"gene split refused: {len(uniq)} genes for {len(df)} rows "
            f"({rows_per_gene:.2f}/gene) — the column is not grouping anything"
        )
    rng = np.random.default_rng(seed)
    test_genes = uniq[rng.permutation(len(uniq))[: int(len(uniq) * test_frac)]]
    te = np.isin(genes, test_genes)
    return Split(~te, te, "gene", n_train_genes=len(uniq) - len(test_genes),
                 n_test_genes=len(test_genes))


def random_row_split(df: pd.DataFrame, test_frac: float = 0.25,
                     seed: int = SEED) -> Split:
    """Split rows at random, ignoring gene membership.

    Deliberately optimistic: the same gene appears on both sides, so the
    model can memorise it. Reported as an upper bound, never as the headline.
    """
    rng = np.random.default_rng(seed)
    mask = rng.random(len(df)) < test_frac
    return Split(~mask, mask, "random_row",
                 note="same gene on both sides; upper bound, not generalisation")


# ---------------------------------------------------------------------------
# Ranking metrics
# ---------------------------------------------------------------------------

def _dcg(rels: np.ndarray) -> float:
    return float(np.sum(rels / np.log2(np.arange(2, len(rels) + 2))))


def ndcg_at_k(y_true: np.ndarray, y_pred: np.ndarray, k: int) -> float:
    order = np.argsort(-y_pred)[:k]
    ideal = np.sort(y_true)[::-1][:k]
    idcg = _dcg(ideal)
    return _dcg(y_true[order]) / idcg if idcg > 0 else float("nan")


def topk_overlap(y_true: np.ndarray, y_pred: np.ndarray, k: int) -> float:
    kk = min(k, len(y_true))
    true_top = set(np.argsort(-y_true)[:kk])
    pred_top = set(np.argsort(-y_pred)[:kk])
    return len(true_top & pred_top) / kk


def average_precision_at_k(y_true: np.ndarray, y_pred: np.ndarray,
                           k: int) -> float:
    """AP over the true top-k treated as the relevant set."""
    kk = min(k, len(y_true))
    relevant = set(np.argsort(-y_true)[:kk])
    order = np.argsort(-y_pred)[:kk]
    hits, precisions = 0, []
    for rank, idx in enumerate(order, start=1):
        if idx in relevant:
            hits += 1
            precisions.append(hits / rank)
    return float(np.mean(precisions)) if precisions else 0.0


def reciprocal_rank(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """1/rank of the single best true item in the predicted order."""
    best = int(np.argmax(y_true))
    order = np.argsort(-y_pred)
    rank = int(np.where(order == best)[0][0]) + 1
    return 1.0 / rank


def per_group_metrics(df_test: pd.DataFrame, pred: np.ndarray,
                      k_values=K_VALUES) -> pd.DataFrame:
    """Every ranking metric, computed WITHIN each experiment group.

    Ranking quality is a within-experiment question — comparing an oligo in
    one assay against one in another is not what the product asks.
    """
    d = df_test.assign(_pred=pred)
    rng = np.random.default_rng(SEED)
    rows = []
    for (gid, gdf) in d.groupby("experiment_id", sort=False):
        if len(gdf) < MIN_GROUP:
            continue
        yt = gdf["rank_label"].to_numpy(dtype=float)
        yp = gdf["_pred"].to_numpy(dtype=float)
        # TIES ARE BROKEN AT RANDOM, not by row order.
        #
        # `np.argsort` is stable, so equal scores keep their original order —
        # and the original order in this benchmark is the order the oligos
        # were tiled along the transcript. Any scorer with ties would then be
        # silently ranking by transcript position. That is not hypothetical:
        # the platform's own composite score rounds to one decimal and ties
        # 27% of oligos within an experiment, with a tied top-5 cutoff in 58%
        # of groups. Shuffling each group first makes tie-breaking arbitrary,
        # which is what a tie means.
        perm = rng.permutation(len(gdf))
        yt, yp = yt[perm], yp[perm]
        rec = {
            "experiment_id": gid,
            "gene": gdf["target_gene"].iloc[0],
            "modality": gdf["modality"].iloc[0],
            "chem_class": gdf.get("chem_class", pd.Series(["?"])).iloc[0],
            "n": len(gdf),
            "mrr": reciprocal_rank(yt, yp),
        }
        if len(np.unique(yt)) > 1 and len(np.unique(yp)) > 1:
            rec["pearson"] = float(np.corrcoef(yt, yp)[0, 1])
        else:
            rec["pearson"] = np.nan
        for k in k_values:
            # A top-k metric on a group of k or fewer items is degenerate:
            # every scorer selects the whole group, so top-k overlap is 1.0
            # and MAP@k is 1.0 regardless of the ranking. Those cells are NaN
            # and drop out of the mean rather than inflating it. The
            # consequence is that different k are averaged over different
            # numbers of groups; `n_groups@k` records how many.
            if len(gdf) <= k:
                rec[f"ndcg@{k}"] = np.nan
                rec[f"top{k}"] = np.nan
                rec[f"map@{k}"] = np.nan
                continue
            rec[f"ndcg@{k}"] = ndcg_at_k(yt, yp, k)
            rec[f"top{k}"] = topk_overlap(yt, yp, k)
            rec[f"map@{k}"] = average_precision_at_k(yt, yp, k)
        rows.append(rec)
    return pd.DataFrame(rows)


def summarise(per_group: pd.DataFrame, k_values=K_VALUES) -> dict:
    if per_group.empty:
        return {"n_groups": 0}
    out = {"n_groups": int(len(per_group)),
           "n_genes": int(per_group["gene"].nunique())}
    for col in ["mrr", "pearson"] + [f"{m}@{k}" if m != "top" else f"top{k}"
                                     for k in k_values
                                     for m in ("ndcg", "top", "map")]:
        if col in per_group:
            vals = per_group[col].to_numpy(dtype=float)
            usable = int(np.sum(~np.isnan(vals)))
            out[col] = float(np.nanmean(vals)) if usable else float("nan")
            out[f"n_groups_{col}"] = usable
    return out


def bootstrap_ci_over_genes(per_group: pd.DataFrame, column: str,
                            n_boot: int = 1000, seed: int = SEED):
    """Percentile CI resampling GENES. Rows within a gene are not independent."""
    if per_group.empty or column not in per_group:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    by_gene = {g: sub[column].to_numpy(dtype=float)
               for g, sub in per_group.groupby("gene")}
    keys = list(by_gene)
    means = []
    for _ in range(n_boot):
        pick = rng.choice(len(keys), size=len(keys), replace=True)
        pool = np.concatenate([by_gene[keys[i]] for i in pick])
        pool = pool[~np.isnan(pool)]
        if len(pool):
            means.append(pool.mean())
    if not means:
        return (float("nan"), float("nan"))
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


def write_result(name: str, payload: dict) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


# The strongest grouping key each modality can actually support.
#   rnase_h          339 genes over 159k rows      -> gene
#   splice_switching   6 genes over 2.2k rows      -> gene (coarse but real)
#   sirna            one "gene" per row (broken)   -> experiment_id
# Splitting on a key the column cannot support produces a random split
# reported as a gene split, which is the exact failure this suite exists
# to measure rather than commit.
GROUP_KEY = {
    "rnase_h": "target_gene",
    "splice_switching": "target_gene",
    "sirna": "experiment_id",
}


def grouped_split(df: pd.DataFrame, test_frac: float = 0.25,
                  seed: int = SEED) -> Split:
    """Hold out whole groups, per modality, using each modality's best key.

    Applied within modality so every arm keeps its share of the test set;
    a global split would put whole minority modalities on one side.
    """
    rng = np.random.default_rng(seed)
    test = np.zeros(len(df), dtype=bool)
    notes = []
    n_tr = n_te = 0
    for mod, sub in df.groupby("modality"):
        key = GROUP_KEY.get(mod, "experiment_id")
        uniq = np.unique(sub[key].values)
        n_hold = max(1, int(round(len(uniq) * test_frac)))
        hold = set(uniq[rng.permutation(len(uniq))[:n_hold]])
        test[np.asarray(df["modality"] == mod) &
             np.asarray(df[key].isin(hold))] = True
        notes.append(f"{mod}: {key}, {len(uniq) - n_hold} train / {n_hold} test")
        n_tr += len(uniq) - n_hold
        n_te += n_hold
    return Split(~test, test, "grouped", note="; ".join(notes),
                 n_train_genes=n_tr, n_test_genes=n_te)


def paired_bootstrap_p(a: np.ndarray, b: np.ndarray, n_boot: int = 10000,
                       seed: int = SEED) -> dict:
    """Two-sided paired bootstrap on the mean difference a - b.

    Paired by group, so the shared difficulty of a group cancels. Reported
    alongside Wilcoxon because the per-group metric distributions are
    bounded and skewed, and neither test alone is convincing on its own.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = ~(np.isnan(a) | np.isnan(b))
    d = a[ok] - b[ok]
    if len(d) < 2:
        return {"n": int(len(d)), "mean_diff": float("nan"), "p": float("nan")}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    draws = d[idx].mean(axis=1)
    obs = d.mean()
    # p under H0: shift the resampled distribution to zero mean.
    centred = draws - obs
    p = float((np.abs(centred) >= abs(obs)).mean())
    return {
        "n": int(len(d)),
        "mean_diff": float(obs),
        "ci95": [float(np.percentile(draws, 2.5)),
                 float(np.percentile(draws, 97.5))],
        "p": p,
    }
