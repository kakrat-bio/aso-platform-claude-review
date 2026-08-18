"""TG06 candidate generation — translational regulation.

Designs steric-blocking ASOs against the elements that control how much
protein a transcript makes: the 5' UTR, upstream ORFs, 3' UTR miRNA sites,
structured elements, IRES domains, the Kozak consensus and the poly(A) site.

WHAT IS REAL HERE, AND WHAT IS NOT
----------------------------------
Every number in a candidate is one of two kinds and they are kept apart in
the response.

  realMetrics         computed from the actual sequence — ViennaRNA duplex
                      free energy, self-structure MFE, nearest-neighbour Tm,
                      GC content, element overlap in nucleotides.

  heuristicEstimates  chemistry and length rules of thumb carried over from
                      the silencing designer. Ordinal, not measured.

There is deliberately **no predicted fold-change**. The spec for this module
asked for `_translational_change_score()` returning a log2 fold-change from a
"mechanism-specific sensitivity coefficient", but no such coefficient has
been fitted, there is no calibration set for translational output, and the
mechanism-recovery benchmark cannot score one. A number like "2.3x
upregulation" would be indistinguishable from a measurement once rendered.

What replaces it is `elementEngagement`: a transparent 0-1 ranking signal
built only from quantities that were actually computed — how much of the
regulatory element the oligo covers, and how favourable the duplex is. It
orders candidates, which is what the page needs. It does not claim to
predict the effect size, and `interpretation` on every candidate says so.

See docs/planning/therapeutic_goal_scope_plan_v3.md and the TG06 spec.
"""

from __future__ import annotations

import logging

from services.gene_silencing_service import (
    _calc_gc,
    _calc_tm,
    _cellular_uptake_score,
    _ensembl_get,
    _nuclease_resistance_score,
    _reverse_complement,
    _self_complement_mfe,
    _target_duplex_energy,
    get_target_analysis,
    ENSEMBL_REST,
)
from services.gene_feature_service import _scan_uorfs
from services import real_data_cache as RDC

logger = logging.getLogger(__name__)

# Which regions each TG06 mechanism acts on, and what chemistry it allows.
# All seven are RNase H-independent: the point is to occupy a site, not to
# destroy the transcript, so a gapmer would defeat the mechanism.
TRANSLATIONAL_MECHANISM_CHEMISTRY: dict[str, dict] = {
    "A2":  {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "5p_utr"},
    "A5":  {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "5p_uorf"},
    "A6":  {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "3p_utr_mirna"},
    "A27": {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "structured_element"},
    "A29": {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "ires_element"},
    "A30": {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "kozak_consensus"},
    "A31": {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"], "target_region": "polya_site"},
}

TRANSLATION_CHEMISTRY_OPTIONS = [
    {"id": "pmo", "label": "PMO (Phosphorodiamidate Morpholino)",
     "description": "Non-ionic steric blocker; no RNase H recruitment."},
    {"id": "moe_full_ps", "label": "2'-MOE Full Phosphorothioate",
     "description": "Fully modified steric blocker; the standard non-cleaving backbone."},
    {"id": "lna_dna_mixmer", "label": "LNA/DNA Mixmer",
     "description": "High-affinity steric blocker; LNA content raises Tm."},
]

TRANSLATION_LENGTH_RANGE = {"min": 18, "max": 30, "default": 20, "step": 1}

# The Kozak consensus, gcc(A/G)ccAUGG. Position -3 and +4 are the two that
# dominate initiation efficiency.
_KOZAK_STRONG_MINUS3 = set("AG")
_KOZAK_STRONG_PLUS4 = "G"
KOZAK_UPSTREAM_NT = 6
KOZAK_DOWNSTREAM_NT = 4

# How far either side of the poly(A) signal an A31 oligo may sit. PABP binds
# across the tail, so selectivity has to come from flanking 3' UTR sequence.
POLYA_WINDOW_NT = 60

# A structured element is called where local folding is this much more stable
# than the transcript-wide mean. Not calibrated against activity; it selects
# where to look, and the reported MFE is the real quantity.
STRUCTURE_MFE_Z = -1.0
STRUCTURE_WINDOW_NT = 60


def _ires_domains(utr5: str) -> list[dict]:
    """Candidate IRES domains in the 5' UTR.

    Reports the most stably folded windows, which is where the structured
    domains an ASO would disrupt actually are. This is a STRUCTURAL
    observation, not an IRES classifier: no validated IRES predictor is
    wired, and `predicted` is False on every entry to say so. A29 is gated
    on the user's defect selection, not on this.
    """
    return [
        {**w, "domain": f"structured region {i + 1}", "predicted": False,
         "note": "Stably folded window in the 5' UTR. Not an IRES call — no "
                 "validated IRES predictor is wired."}
        for i, w in enumerate(_structured_windows(utr5)[:3])
    ]


def _structured_windows(seq: str) -> list[dict]:
    """Locally stable windows, by real ViennaRNA folding."""
    if not seq or len(seq) < STRUCTURE_WINDOW_NT:
        return []
    try:
        import RNA
    except ImportError:
        return []

    step = max(10, STRUCTURE_WINDOW_NT // 3)
    windows = []
    for start in range(0, len(seq) - STRUCTURE_WINDOW_NT + 1, step):
        sub = seq[start:start + STRUCTURE_WINDOW_NT]
        structure, mfe = RNA.fold(sub)
        windows.append({"start": start, "end": start + STRUCTURE_WINDOW_NT,
                        "mfe": round(float(mfe), 2), "structure": structure})
    if not windows:
        return []
    mfes = [w["mfe"] for w in windows]
    mean = sum(mfes) / len(mfes)
    sd = (sum((m - mean) ** 2 for m in mfes) / len(mfes)) ** 0.5 or 1.0
    stable = [w for w in windows if (w["mfe"] - mean) / sd <= STRUCTURE_MFE_Z]
    return sorted(stable, key=lambda w: w["mfe"])


def _kozak_context(mrna: str, cds_start: int | None) -> dict | None:
    """The Kozak consensus around the initiator AUG, and its strength.

    Strength is reported as which consensus positions match, not as a
    predicted initiation rate — the -3 purine and the +4 G are the two
    positions with the largest documented effect, so the call is 'strong'
    (both), 'adequate' (one) or 'weak' (neither).
    """
    if not mrna or cds_start is None or cds_start < KOZAK_UPSTREAM_NT:
        return None
    start = cds_start - KOZAK_UPSTREAM_NT
    end = min(len(mrna), cds_start + 3 + KOZAK_DOWNSTREAM_NT)
    context = mrna[start:end]
    if len(context) < KOZAK_UPSTREAM_NT + 3:
        return None

    minus3 = mrna[cds_start - 3] if cds_start >= 3 else ""
    plus4 = mrna[cds_start + 3] if cds_start + 3 < len(mrna) else ""
    matches = []
    if minus3 in _KOZAK_STRONG_MINUS3:
        matches.append("-3 purine")
    if plus4 == _KOZAK_STRONG_PLUS4:
        matches.append("+4 G")
    strength = ("strong" if len(matches) == 2
                else "adequate" if len(matches) == 1 else "weak")
    return {
        "sequence": context, "start": start, "end": end,
        "minus3": minus3, "plus4": plus4,
        "consensusMatches": matches, "strength": strength,
        "note": ("Consensus-position match against gcc(A/G)ccAUGG. Not a "
                 "predicted initiation rate."),
    }


def get_translational_target(
    ensembl_gene_id: str,
    gene_symbol: str = "",
    organism: str = "",
) -> dict:
    """Transcript structure for TG06 targeting.

    Live from Ensembl, else a real earlier fetch, else explicitly
    unavailable — the same three-outcome contract as everything else. No
    element is synthesised when the source cannot answer.
    """
    def fetch() -> dict | None:
        base = get_target_analysis(ensembl_gene_id, gene_symbol, organism)
        if not base or not base.get("mrnaSequence"):
            return None
        mrna = (base.get("mrnaSequence") or "").upper().replace("T", "U")
        utr5 = (base.get("utr5Sequence") or "").upper().replace("T", "U")
        utr3 = (base.get("utr3Sequence") or "").upper().replace("T", "U")
        cds_start = len(utr5) if utr5 else None

        return {
            "canonicalTranscript": base.get("canonicalTranscript"),
            "mrnaSequence": mrna,
            "utr5": {"sequence": utr5, "start": 0, "end": len(utr5)} if utr5 else None,
            "utr3": {"sequence": utr3,
                     "start": max(0, len(mrna) - len(utr3)),
                     "end": len(mrna)} if utr3 else None,
            "cdsStart": cds_start,
            "uorfs": _scan_uorfs(utr5.replace("U", "T")) if utr5 else [],
            "kozak": _kozak_context(mrna, cds_start),
            "structuredElements": _structured_windows(utr5 or mrna[:600]),
            "ires": _ires_domains(utr5),
            "polyASite": ({"start": max(0, len(mrna) - POLYA_WINDOW_NT),
                           "end": len(mrna),
                           "sequence": mrna[-POLYA_WINDOW_NT:]}
                          if len(mrna) > POLYA_WINDOW_NT else None),
            "exons": base.get("exons", []),
        }

    resolved = RDC.resolve(
        "translational_target", f"{organism}:{gene_symbol or ensembl_gene_id}",
        fetch, source="Ensembl REST",
    )
    return {
        "status": resolved["status"],
        "dataProvenance": {k: resolved[k] for k in
                           ("status", "source", "fetchedAt", "ageSeconds",
                            "stale", "note")},
        **(resolved["data"] or {}),
    }


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------

def _element_for(target_element: str, target: dict) -> tuple[int, int, str] | None:
    """(start, end, label) of the region a mechanism acts on."""
    mrna = target.get("mrnaSequence") or ""
    utr5 = target.get("utr5") or {}
    utr3 = target.get("utr3") or {}

    if target_element == "5p_utr" and utr5:
        return utr5["start"], utr5["end"], "5' UTR"
    if target_element == "kozak_consensus" and target.get("kozak"):
        k = target["kozak"]
        return k["start"], k["end"], "Kozak consensus"
    if target_element == "5p_uorf":
        uorfs = target.get("uorfs") or []
        if uorfs:
            u = uorfs[0]
            return int(u.get("start", 1)) - 1, int(u.get("end", 1)), "uORF"
    if target_element == "3p_utr_mirna" and utr3:
        return utr3["start"], utr3["end"], "3' UTR"
    if target_element in ("structured_element", "ires_element"):
        elements = (target.get("ires") if target_element == "ires_element"
                    else target.get("structuredElements")) or []
        if elements:
            e = elements[0]
            offset = (utr5.get("start", 0) if utr5 else 0)
            return offset + e["start"], offset + e["end"], (
                "IRES domain" if target_element == "ires_element"
                else "structured element")
    if target_element == "polya_site" and target.get("polyASite"):
        p = target["polyASite"]
        return p["start"], p["end"], "poly(A) site"
    if mrna:
        return None
    return None


def _element_engagement(overlap_nt: int, oligo_len: int,
                        duplex_dg: float) -> float:
    """A transparent 0-1 ranking signal. NOT a predicted fold-change.

    Two computed quantities only:

      coverage   what fraction of the oligo actually sits on the regulatory
                 element. An oligo half off the element engages it half as
                 much, whatever else is true.
      affinity   the real ViennaRNA duplex free energy, normalised over a
                 range that spans what these lengths produce.

    Combined as a plain mean so the contribution of each is obvious. This
    orders candidates; it does not estimate how much protein output moves.
    """
    coverage = max(0.0, min(1.0, overlap_nt / max(oligo_len, 1)))
    # -40 kcal/mol is comfortably past what a 20-30mer duplex reaches, so
    # this saturates rather than rewarding runaway affinity.
    affinity = max(0.0, min(1.0, -duplex_dg / 40.0))
    return round((coverage + affinity) / 2.0, 4)


def generate_translational_candidates(
    target_element: str,
    translational_goal: str,
    mechanism_id: str,
    aso_length: int,
    chemistry: str,
    modifications: list[str] | None,
    target: dict,
    delivery_context: str | None = None,
    target_rbp: str | None = None,
    max_candidates: int = 20,
) -> dict:
    """Tile the regulatory element and rank the oligos that cover it."""
    modifications = modifications or []
    spec = TRANSLATIONAL_MECHANISM_CHEMISTRY.get(mechanism_id)
    if spec is None:
        return {"status": "unsupported_mechanism", "candidates": [],
                "message": f"{mechanism_id} is not a TG06 mechanism."}
    if chemistry not in spec["allowed"]:
        return {
            "status": "incompatible_chemistry", "candidates": [],
            "message": (
                f"{chemistry} is not usable for {mechanism_id}. Translational "
                f"regulation needs an RNase H-independent steric blocker — a "
                f"cleaving chemistry destroys the transcript instead of "
                f"occupying the element. Allowed: {', '.join(spec['allowed'])}."
            ),
        }

    mrna = target.get("mrnaSequence") or ""
    if not mrna:
        return {
            "status": RDC.UNAVAILABLE, "candidates": [],
            "message": ("No transcript sequence is available for this target, "
                        "so no oligo can be designed against it."),
        }

    region = _element_for(target_element, target)
    if region is None:
        return {
            "status": "element_not_found", "candidates": [],
            "message": (
                f"No {target_element.replace('_', ' ')} was located in this "
                f"transcript, so there is nothing for {mechanism_id} to "
                f"occupy. This is a property of the transcript, not a "
                f"failure to compute."
            ),
        }

    el_start, el_end, el_label = region
    # Tile a window that extends one oligo length past the element on each
    # side, so partial-overlap designs are considered and then ranked down
    # by coverage rather than silently excluded.
    scan_start = max(0, el_start - aso_length)
    scan_end = min(len(mrna), el_end + aso_length)

    candidates = []
    for pos in range(scan_start, scan_end - aso_length + 1):
        site = mrna[pos:pos + aso_length]
        if len(site) < aso_length or "N" in site:
            continue
        overlap = max(0, min(pos + aso_length, el_end) - max(pos, el_start))
        if overlap == 0:
            continue

        # The shared _reverse_complement translates ATGC only, so a U in an
        # RNA target passes through untouched and yields a mixed-alphabet
        # oligo that primer3's Tm rejects. The silencing designer never hits
        # this because it works on DNA-alphabet cDNA. Normalise here, and
        # keep the ASO in the DNA alphabet, which is how oligos are written
        # and what every downstream helper expects.
        aso = _reverse_complement(site.replace("U", "T"))
        gc = _calc_gc(aso)
        tm = _calc_tm(aso)
        duplex_dg = _target_duplex_energy(aso, site)
        self_mfe = _self_complement_mfe(aso)
        engagement = _element_engagement(overlap, aso_length, duplex_dg)

        candidates.append({
            "sequence": aso,
            "sequenceAlphabet": "DNA",
            "targetSite": site,
            "targetSiteAlphabet": "RNA",
            "targetStart": pos,
            "targetEnd": pos + aso_length,
            "targetRegion": el_label,
            "targetElement": target_element,
            "mechanismId": mechanism_id,
            "chemistry": chemistry,
            "modifications": modifications,
            "deliveryContext": delivery_context,
            "elementOverlapNt": overlap,
            "realMetrics": {
                "targetDuplexEnergy": round(duplex_dg, 2),
                "meltingTempC": round(tm, 1),
                "selfStructureMfe": round(self_mfe, 2),
                "gcContent": round(gc, 3),
                "lengthNt": aso_length,
                "elementOverlapNt": overlap,
                "provenance": "ViennaRNA duplexfold / fold; nearest-neighbour Tm",
            },
            "heuristicEstimates": {
                "nucleaseResistance": _nuclease_resistance_score(chemistry, modifications),
                "cellularUptake": _cellular_uptake_score(chemistry, aso_length),
                "provenance": ("Ordinal chemistry and length heuristics, not "
                               "measurements."),
            },
            "elementEngagement": engagement,
            "interpretation": (
                f"Covers {overlap}/{aso_length} nt of the {el_label}. "
                f"elementEngagement is a ranking signal built from that "
                f"coverage and the computed duplex energy — it is NOT a "
                f"predicted change in protein output, and no calibrated "
                f"model for that exists here."
            ),
        })

    if not candidates:
        return {"status": "no_candidates", "candidates": [],
                "message": f"No {aso_length} nt window overlaps the {el_label}."}

    candidates.sort(key=lambda c: (-c["elementEngagement"],
                                   c["realMetrics"]["targetDuplexEnergy"]))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1

    rbp_note = None
    if target_rbp and target_rbp.strip():
        rbp_note = (
            f"You named {target_rbp.strip()} as the competing RNA-binding "
            f"protein. No binding-energy model for it is wired, so no "
            f"displacement score is reported — a competition number without "
            f"the RBP's own affinity would be invented. The duplex energies "
            f"above are the real quantity to compare against once you have it."
        )

    return {
        "status": "ok",
        "mechanismId": mechanism_id,
        "targetElement": target_element,
        "translationalGoal": translational_goal,
        "elementRegion": {"start": el_start, "end": el_end, "label": el_label},
        "candidates": candidates[:max_candidates],
        "totalScanned": len(candidates),
        "rbpNote": rbp_note,
        "scoringNote": (
            "Candidates are ordered by elementEngagement, then by duplex free "
            "energy. Every value under realMetrics is computed from the "
            "sequence; every value under heuristicEstimates is an ordinal "
            "rule of thumb. No predicted fold-change is emitted."
        ),
    }


def get_translational_design_options() -> dict:
    """Form options for the TG06 page."""
    return {
        "chemistries": TRANSLATION_CHEMISTRY_OPTIONS,
        "lengthRange": TRANSLATION_LENGTH_RANGE,
        "mechanismRegions": {
            mid: spec["target_region"]
            for mid, spec in TRANSLATIONAL_MECHANISM_CHEMISTRY.items()
        },
    }
