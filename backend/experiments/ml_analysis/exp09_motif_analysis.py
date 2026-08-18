"""E9 — Motif and architecture analysis.

Three questions, three sections.

**1. Gapmer architecture (the 5-10-5 question).** The benchmark's chemistry
strings are position-level annotations, e.g.
`L20 MOE|sugar|1,2,3,4,5,16,17,18,19,20 PS|backbone|...`, which is literally
a 5-10-5 MOE gapmer. So the architecture does not have to be inferred — it
can be parsed and counted. What the generator emits, however, is a bare
sequence: it is conditioned on chemistry but does not produce a
modification pattern. So the honest question for the DESIGNS is the weaker
one it can answer — whether their lengths are compatible with the
architectures the training data actually uses.

**2. Reynolds rules.** The eight sequence criteria of Reynolds et al. 2004
(Nat Biotechnol 22:326), scored on 19-mers. The siRNA arm is 3,947 19-mers,
so the rules apply exactly. They are scored on the training data too, which
is the part that makes it a finding: if high-activity training siRNAs do
not satisfy the rules any more often than low-activity ones, then a
generator matching the rules is matching a convention rather than activity.

**3. Recurrent subsequences.** Enriched k-mers in the generated designs,
per mechanism, against a dinucleotide-shuffled control — which holds
length, base composition AND dinucleotide frequency fixed, so an enriched
k-mer is not just a GC artefact. Shared and mechanism-specific k-mers are
listed separately.
"""

from __future__ import annotations

import collections
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.generate_sequences import (  # noqa: E402
    load_generated,
)

NAME = "exp09_motif_analysis"
KMER_K = 4
TOP_KMERS = 15


# --------------------------------------------------------------------------
# 1. Architecture, parsed from the chemistry annotation
# --------------------------------------------------------------------------

_SUGAR = re.compile(r"(\w+)\|sugar\|([\d,]+)")


def parse_architecture(chem: str) -> dict | None:
    """Return {'wings': (5', 3'), 'gap': n, 'sugar': 'MOE', 'length': L}.

    Only meaningful for a gapmer layout: modified positions at both termini
    and an unmodified stretch in the middle. A uniformly modified oligo (a
    steric blocker) has no gap and is labelled `uniform`.
    """
    if not isinstance(chem, str) or "|sugar|" not in chem:
        return None
    m = re.match(r"L(\d+)\s", chem)
    length = int(m.group(1)) if m else None
    sm = _SUGAR.search(chem)
    if not sm or length is None:
        return None
    sugar = sm.group(1)
    pos = sorted(int(p) for p in sm.group(2).split(",") if p)
    if not pos:
        return None
    if len(pos) == length:
        return {"sugar": sugar, "length": length, "layout": "uniform",
                "wing5": length, "gap": 0, "wing3": 0}
    # Leading run from position 1, trailing run ending at position L.
    w5 = 0
    while w5 < len(pos) and pos[w5] == w5 + 1:
        w5 += 1
    w3 = 0
    while w3 < len(pos) and pos[len(pos) - 1 - w3] == length - w3:
        w3 += 1
    gap = length - w5 - w3
    if w5 == 0 or w3 == 0 or gap <= 0 or w5 + w3 != len(pos):
        return {"sugar": sugar, "length": length, "layout": "irregular",
                "wing5": w5, "gap": gap, "wing3": w3}
    return {"sugar": sugar, "length": length, "layout": f"{w5}-{gap}-{w3}",
            "wing5": w5, "gap": gap, "wing3": w3}


def architecture_section(df: pd.DataFrame, gen: pd.DataFrame) -> dict:
    rows = []
    for chem in df["chemistry"].drop_duplicates():
        a = parse_architecture(chem)
        if a:
            rows.append({"chemistry": chem, **a})
    adf = pd.DataFrame(rows)
    counts = df["chemistry"].value_counts()
    if not adf.empty:
        adf["rows"] = adf["chemistry"].map(counts).fillna(0).astype(int)
        adf["modality"] = adf["chemistry"].map(
            df.drop_duplicates("chemistry").set_index("chemistry")["modality"])

    by_layout = (adf.groupby(["modality", "sugar", "layout"])["rows"].sum()
                 .sort_values(ascending=False).reset_index()
                 if not adf.empty else pd.DataFrame())

    gen_len = (gen[gen["kind"] == "generated"]
               .assign(L=lambda d: d["seq"].str.len())
               .groupby("mechanism")["L"]
               .agg(["mean", "std", "min", "max", "count"]).round(3))

    # An architecture has a fixed length: a 5-10-5 gapmer is 20 nt and a
    # 3-10-3 is 16 nt. So the question a bare generated sequence CAN answer is
    # whether its length matches the length the dominant architectures need.
    compatible = {}
    if not adf.empty:
        for mech, sub in gen[gen["kind"] == "generated"].groupby("mechanism"):
            lens = sub["seq"].str.len()
            hist = lens.value_counts().sort_index()
            arch = adf[adf["modality"] == mech]
            # Training rows by the length their architecture requires.
            arch_rows = arch.groupby("length")["rows"].sum().sort_values(
                ascending=False)
            layout_by_len = (arch.groupby(["length", "layout"])["rows"].sum()
                             .sort_values(ascending=False))
            top_len = int(arch_rows.index[0]) if len(arch_rows) else None
            compatible[mech] = {
                "generated_length_histogram": {int(k): int(v)
                                               for k, v in hist.items()},
                "generated_modal_length": int(lens.mode().iloc[0]),
                "training_rows_by_architecture_length": {
                    int(k): int(v) for k, v in arch_rows.items()},
                "dominant_architecture": (
                    {"length": top_len,
                     "layout": str(layout_by_len.index[0][1]),
                     "rows": int(layout_by_len.iloc[0])}
                    if len(layout_by_len) else None),
                "generated_share_at_dominant_length": (
                    round(float((lens == top_len).mean()), 4)
                    if top_len is not None else None),
                "reading": "the generator emits length but not a modification "
                           "pattern; a design whose length is not the "
                           "dominant architecture's length cannot be built as "
                           "that architecture at all",
            }

    return {
        "note": "architecture is PARSED from the benchmark's chemistry "
                "annotation, not inferred from sequence; the generator emits "
                "bare sequences and does not produce a modification pattern, "
                "so only length compatibility can be checked for the designs",
        "training_layouts_by_rows": by_layout.to_dict("records")
                                    if not by_layout.empty else [],
        "is_5_10_5_dominant": (
            bool(not by_layout.empty
                 and by_layout.iloc[0]["layout"] == "5-10-5")
            if not by_layout.empty else None),
        "generated_length_distribution": gen_len.to_dict("index"),
        "generated_length_compatibility": compatible,
    }


# --------------------------------------------------------------------------
# 2. Reynolds rules
# --------------------------------------------------------------------------

def reynolds(seq: str) -> dict | None:
    """The eight Reynolds et al. 2004 criteria, on a 19-mer sense strand."""
    s = seq.upper().replace("T", "U")
    if len(s) != 19 or set(s) - set("ACGU"):
        return None
    gc = (s.count("G") + s.count("C")) / 19
    au_15_19 = sum(1 for c in s[14:19] if c in "AU")
    longest_run = max((len(r) for r in re.findall(r"(.)\1*", s)), default=1)
    crit = {
        "gc_30_52": 0.30 <= gc <= 0.52,
        "au_rich_15_19": au_15_19 >= 3,
        "no_internal_repeat": longest_run <= 4,
        "A_at_19": s[18] == "A",
        "A_at_3": s[2] == "A",
        "U_at_10": s[9] == "U",
        "not_GC_at_19": s[18] not in "GC",
        "not_G_at_13": s[12] != "G",
    }
    crit["score"] = int(sum(bool(v) for k, v in crit.items() if k != "score"))
    crit["gc"] = round(gc, 4)
    return crit


def _unpaired_bootstrap(x: np.ndarray, y: np.ndarray,
                        n_boot: int = 10000, seed: int = C.SEED) -> dict:
    """Two-sample bootstrap on the difference in means.

    UNPAIRED on purpose. The top-20% and bottom-20% siRNAs are different
    sequences with no correspondence between them, so pairing them by
    position — which is what truncating both arrays to a common length and
    subtracting would do — invents a relationship that is not there.
    """
    rng = np.random.default_rng(seed)
    diffs = (rng.choice(x, (n_boot, len(x))).mean(axis=1)
             - rng.choice(y, (n_boot, len(y))).mean(axis=1))
    obs = float(x.mean() - y.mean())
    centred = diffs - obs
    return {
        "n_x": int(len(x)), "n_y": int(len(y)),
        "mean_diff": round(obs, 4),
        "ci95": [round(float(np.percentile(diffs, 2.5)), 4),
                 round(float(np.percentile(diffs, 97.5)), 4)],
        "p": float((np.abs(centred) >= abs(obs)).mean()),
    }


def reynolds_section(df: pd.DataFrame, gen: pd.DataFrame) -> dict:
    out = {"criteria": "Reynolds et al. 2004, Nat Biotechnol 22:326; scored "
                       "on 19-mer sense strands only"}

    def score_set(seqs):
        rs = [reynolds(s) for s in seqs]
        rs = [r for r in rs if r]
        if not rs:
            return None
        keys = [k for k in rs[0] if k not in ("score", "gc")]
        return {
            "n_scored": len(rs),
            "n_skipped_not_19mer": int(len(seqs) - len(rs)),
            "mean_score_out_of_8": round(float(np.mean([r["score"] for r in rs])), 3),
            "criterion_pass_rate": {k: round(float(np.mean([r[k] for r in rs])), 4)
                                    for k in keys},
        }

    si = df[df["modality"] == "sirna"]
    out["training_sirna_all"] = score_set(si["seq"].tolist())

    # The part that makes this a finding rather than a checklist: do the
    # ACTIVE training siRNAs satisfy the rules more often than the inactive?
    if len(si) > 40:
        hi = si[si["rank_label"] >= si["rank_label"].quantile(0.8)]
        lo = si[si["rank_label"] <= si["rank_label"].quantile(0.2)]
        h, l = score_set(hi["seq"].tolist()), score_set(lo["seq"].tolist())
        if h and l:
            hs = [reynolds(s)["score"] for s in hi["seq"] if reynolds(s)]
            ls = [reynolds(s)["score"] for s in lo["seq"] if reynolds(s)]
            out["training_sirna_top20pct"] = h
            out["training_sirna_bottom20pct"] = l
            out["active_vs_inactive"] = {
                "n_top20pct": len(hs), "n_bottom20pct": len(ls),
                "mean_score_top20pct": round(float(np.mean(hs)), 4),
                "mean_score_bottom20pct": round(float(np.mean(ls)), 4),
                "mean_score_difference": round(float(np.mean(hs) - np.mean(ls)), 4),
                "bootstrap": _unpaired_bootstrap(np.array(hs), np.array(ls)),
                "reading": "a difference near zero means the Reynolds rules do "
                           "not separate active from inactive siRNAs in THIS "
                           "benchmark, so matching them is matching a "
                           "convention, not predicted activity",
            }

    for kind in ("generated", "shuffled", "random"):
        sub = gen[(gen["mechanism"] == "sirna") & (gen["kind"] == kind)]
        out[f"designs_{kind}"] = score_set(sub["seq"].tolist())
    return out


# --------------------------------------------------------------------------
# 3. Recurrent subsequences
# --------------------------------------------------------------------------

def kmer_counts(seqs, k: int) -> collections.Counter:
    c = collections.Counter()
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - k + 1):
            sub = s[i:i + k]
            if set(sub) <= set("ACGU"):
                c[sub] += 1
    return c


def enrichment(fg: collections.Counter, bg: collections.Counter) -> pd.DataFrame:
    keys = set(fg) | set(bg)
    nf, nb = max(sum(fg.values()), 1), max(sum(bg.values()), 1)
    rows = []
    for kmer in keys:
        f, b = fg.get(kmer, 0), bg.get(kmer, 0)
        pf, pb = (f + 1) / (nf + len(keys)), (b + 1) / (nb + len(keys))
        rows.append({"kmer": kmer, "count": f, "control_count": b,
                     "freq": round(pf, 6), "control_freq": round(pb, 6),
                     "log2_enrichment": round(float(np.log2(pf / pb)), 4)})
    return pd.DataFrame(rows).sort_values("log2_enrichment", ascending=False)


def kmer_section(gen: pd.DataFrame, k: int = KMER_K) -> dict:
    out, top_sets = {}, {}
    for mech, sub in gen.groupby("mechanism"):
        fg = kmer_counts(sub[sub["kind"] == "generated"]["seq"], k)
        bg = kmer_counts(sub[sub["kind"] == "shuffled"]["seq"], k)
        e = enrichment(fg, bg)
        out[mech] = {
            "control": "dinucleotide-shuffled designs (length, composition "
                       "and dinucleotide frequency held fixed)",
            "distinct_kmers": int(len(fg)),
            "enriched": e.head(TOP_KMERS).to_dict("records"),
            "depleted": e.tail(TOP_KMERS).iloc[::-1].to_dict("records"),
        }
        top_sets[mech] = set(e.head(TOP_KMERS)["kmer"])

    mechs = sorted(top_sets)
    shared = set.intersection(*top_sets.values()) if len(top_sets) > 1 else set()
    out["_across_mechanisms"] = {
        "k": k,
        "shared_top_kmers": sorted(shared),
        "mechanism_specific": {
            m: sorted(top_sets[m] - set.union(*[top_sets[o] for o in mechs
                                                if o != m]))
            for m in mechs} if len(mechs) > 1 else {},
        "pairwise_jaccard": {
            f"{a}|{b}": round(len(top_sets[a] & top_sets[b])
                              / max(len(top_sets[a] | top_sets[b]), 1), 4)
            for i, a in enumerate(mechs) for b in mechs[i + 1:]},
    }
    return out


def main() -> dict:
    df = C.load_benchmark()
    gen = load_generated()

    payload = {
        "experiment": NAME,
        "question": "What structure do the generated designs have — gapmer "
                    "architecture, Reynolds compliance, recurrent motifs?",
        "design_sets": {
            "source": "backend/results/benchmark/generative_v3/generator.pt",
            "counts": gen.groupby(["mechanism", "kind"]).size()
                         .unstack(fill_value=0).to_dict("index"),
        },
        "architecture": architecture_section(df, gen),
        "reynolds": reynolds_section(df, gen),
        "recurrent_kmers": kmer_section(gen),
    }
    C.write_result(NAME, payload)
    a = payload["architecture"]
    print("5-10-5 dominant:", a["is_5_10_5_dominant"])
    for r in a["training_layouts_by_rows"][:6]:
        print(" ", r)
    print("shared top k-mers:",
          payload["recurrent_kmers"]["_across_mechanisms"]["shared_top_kmers"])
    return payload


if __name__ == "__main__":
    main()
