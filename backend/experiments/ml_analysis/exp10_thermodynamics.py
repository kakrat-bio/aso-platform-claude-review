"""E10 — Thermodynamic validation of the generated designs.

Four quantities per oligo, each from a tool that is actually appropriate
for it:

    Tm            nearest-neighbour melting temperature of the oligo against
                  its perfect complement (primer3, SantaLucia parameters)
    hairpin dG    self-structure of the oligo alone (primer3)
    MFE           minimum free energy fold of the oligo as RNA (ViennaRNA
                  RNA.fold)
    duplex dG     hybridisation free energy against the perfect complement
                  (ViennaRNA RNA.duplexfold)

Two caveats stated up front rather than buried.

**These are unmodified-backbone numbers.** Every oligo in this benchmark
that matters clinically carries 2'-MOE, cEt or PS modifications, and those
shift Tm substantially (MOE roughly +1 to +2 degrees per modified residue).
Neither primer3 nor ViennaRNA models them. So these values are comparable
BETWEEN sequence sets — which is exactly what a distribution-matching
question needs — and are not predictions of the real duplex stability of a
modified oligo.

**Matching a distribution is not evidence of activity.** A generator that
reproduces the training Tm distribution has learned the composition of the
training set. That is worth confirming and is not the same as designing
something that works. The test here is two-sample (Kolmogorov-Smirnov plus
a difference in means with a bootstrap interval) against the real training
sequences of the same mechanism, with the random and shuffled controls
alongside so the reader can see what "not matching" looks like.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

try:
    import RNA
except ImportError:
    RNA = None
try:
    import primer3
except ImportError:
    primer3 = None
try:
    from scipy.stats import ks_2samp
except ImportError:
    ks_2samp = None

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.ml_analysis.generate_sequences import (  # noqa: E402
    load_generated,
)

NAME = "exp10_thermodynamics"
N_TRAINING_SAMPLE = 1500
COMPLEMENT_RNA = {"A": "U", "U": "A", "G": "C", "C": "G"}
COMPLEMENT_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}


def to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def revcomp(seq: str, table: dict) -> str:
    return "".join(table.get(c, "N") for c in reversed(seq.upper()))


def thermo(seq: str) -> dict:
    """All four quantities for one oligo. Missing tools yield None, not zero."""
    out = {"tm": None, "hairpin_dg": None, "mfe": None, "duplex_dg": None}
    dna = to_dna(seq)
    rna = to_rna(seq)
    if set(dna) - set("ACGT"):
        return out
    if primer3 is not None:
        out["tm"] = float(primer3.calc_tm(dna))
        # primer3 reports dG in cal/mol; kcal/mol keeps it on the same scale
        # as the ViennaRNA numbers next to it.
        out["hairpin_dg"] = float(primer3.calc_hairpin(dna).dg) / 1000.0
    if RNA is not None:
        _, mfe = RNA.fold(rna)
        out["mfe"] = float(mfe)
        out["duplex_dg"] = float(
            RNA.duplexfold(rna, revcomp(rna, COMPLEMENT_RNA)).energy)
    return out


def describe(values) -> dict | None:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    v = v[~np.isnan(v)]
    if not len(v):
        return None
    return {
        "n": int(len(v)),
        "mean": round(float(v.mean()), 4),
        "sd": round(float(v.std(ddof=1)) if len(v) > 1 else 0.0, 4),
        "p05": round(float(np.percentile(v, 5)), 4),
        "median": round(float(np.median(v)), 4),
        "p95": round(float(np.percentile(v, 95)), 4),
    }


def compare(a, b, seed: int = C.SEED) -> dict:
    """Two-sample comparison. Unpaired: different sequences, no matching."""
    x = np.asarray([v for v in a if v is not None], dtype=float)
    y = np.asarray([v for v in b if v is not None], dtype=float)
    x, y = x[~np.isnan(x)], y[~np.isnan(y)]
    if len(x) < 5 or len(y) < 5:
        return {"unavailable": f"{len(x)} vs {len(y)} usable values"}
    rng = np.random.default_rng(seed)
    diffs = (rng.choice(x, (2000, len(x))).mean(axis=1)
             - rng.choice(y, (2000, len(y))).mean(axis=1))
    out = {
        "mean_difference": round(float(x.mean() - y.mean()), 4),
        "ci95_difference": [round(float(np.percentile(diffs, 2.5)), 4),
                            round(float(np.percentile(diffs, 97.5)), 4)],
        # Standardised so Tm (degrees) and dG (kcal/mol) are readable together.
        "cohens_d": round(float((x.mean() - y.mean())
                                / np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2)), 4),
    }
    if ks_2samp is not None:
        ks = ks_2samp(x, y)
        out["ks_statistic"] = round(float(ks.statistic), 4)
        out["ks_p"] = float(ks.pvalue)
        out["reading"] = (
            "KS p below 0.05 means the two distributions differ detectably; "
            "with n in the thousands this happens for differences too small "
            "to matter, so read Cohen's d alongside it")
    else:
        out["ks_statistic"] = None
        out["ks_unavailable"] = "scipy not installed"
    return out


def main() -> dict:
    if RNA is None and primer3 is None:
        payload = {"experiment": NAME, "status": "blocked",
                   "blocker": "neither ViennaRNA nor primer3 is importable; "
                              "no thermodynamic quantity can be computed"}
        C.write_result(NAME, payload)
        print(payload)
        return payload

    df = C.load_benchmark()
    gen = load_generated()
    rng = np.random.default_rng(C.SEED)

    metrics = ["tm", "hairpin_dg", "mfe", "duplex_dg"]
    sets: dict[tuple[str, str], list[dict]] = {}

    for mech in sorted(gen["mechanism"].unique()):
        real = df[df["modality"] == mech]
        if len(real) > N_TRAINING_SAMPLE:
            real = real.iloc[rng.permutation(len(real))[:N_TRAINING_SAMPLE]]
        print(f"[{NAME}] {mech}: training n={len(real)}", flush=True)
        sets[(mech, "training")] = [thermo(s) for s in real["seq"]]
        # High-activity training oligos: the distribution a design set would
        # ideally match, rather than the mixed pool.
        hi = df[(df["modality"] == mech)
                & (df["rank_label"] >= df[df["modality"] == mech]
                   ["rank_label"].quantile(0.8))]
        if len(hi) > N_TRAINING_SAMPLE:
            hi = hi.iloc[rng.permutation(len(hi))[:N_TRAINING_SAMPLE]]
        sets[(mech, "training_top20pct")] = [thermo(s) for s in hi["seq"]]
        for kind in ("generated", "shuffled", "random"):
            sub = gen[(gen["mechanism"] == mech) & (gen["kind"] == kind)]
            print(f"[{NAME}] {mech}/{kind}: n={len(sub)}", flush=True)
            sets[(mech, kind)] = [thermo(s) for s in sub["seq"]]

    distributions, comparisons, table = {}, {}, []
    for (mech, kind), vals in sets.items():
        d = {m: describe([v[m] for v in vals]) for m in metrics}
        distributions.setdefault(mech, {})[kind] = d
        for m in metrics:
            if d[m]:
                table.append({"mechanism": mech, "set": kind, "metric": m,
                              "mean": d[m]["mean"], "sd": d[m]["sd"],
                              "n": d[m]["n"]})

    for mech in sorted(gen["mechanism"].unique()):
        comparisons[mech] = {}
        for ref in ("training", "training_top20pct"):
            for kind in ("generated", "shuffled", "random"):
                comparisons[mech][f"{kind}_vs_{ref}"] = {
                    m: compare([v[m] for v in sets[(mech, kind)]],
                               [v[m] for v in sets[(mech, ref)]])
                    for m in metrics
                }

    payload = {
        "experiment": NAME,
        "question": "Do the generated designs land in the same thermodynamic "
                    "distribution as the real training oligos?",
        "tools": {
            "primer3": getattr(primer3, "__version__", "installed")
                       if primer3 else "MISSING",
            "ViennaRNA": getattr(RNA, "__version__", "installed")
                         if RNA else "MISSING",
            "scipy": "installed" if ks_2samp else "MISSING",
        },
        "caveats": [
            "Tm, hairpin dG, MFE and duplex dG are computed on the UNMODIFIED "
            "backbone. 2'-MOE, cEt and PS shift real duplex stability and "
            "neither tool models them; these numbers compare sequence sets to "
            "each other, they are not predictions of a modified oligo's Tm.",
            "Matching the training distribution shows the generator learned "
            "the training composition. It is not evidence of activity.",
            "primer3 hairpin dG converted from cal/mol to kcal/mol to sit on "
            "the same scale as the ViennaRNA energies.",
        ],
        "distributions": distributions,
        "comparisons": comparisons,
        "table": table,
    }
    C.write_result(NAME, payload)
    print(pd.DataFrame(table).pivot_table(index=["mechanism", "metric"],
                                          columns="set", values="mean")
          .round(3).to_string())
    return payload


if __name__ == "__main__":
    main()
