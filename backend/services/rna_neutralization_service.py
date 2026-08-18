"""TG05 candidate generation — RNA neutralization.

Designs oligos that occupy a toxic RNA rather than degrade it. Three modes:

  steric_repeat_masking (A14)  mask an expanded repeat tract so it stops
                               sequestering RNA-binding proteins into foci
  microrna_antagomir (A12)     sequester a pathogenic miRNA
  aptamer_decoy (A25)          flagged, not designed here

WHERE THE TARGET SEQUENCE COMES FROM
------------------------------------
A repeat tract is the one design target this platform can construct
legitimately: an expanded (CAG)n or (CUG)n tract *is* the repeat unit
repeated, so building it from a verified unit and count is reproducing a
known sequence, not inventing one.

That only holds if the unit is real. The order is:

  1. the curated repeat-expansion catalogue (F12), when the gene is in it —
     provenance CONFIRMED, and the catalogue also gives the transcript
     region, which decides whether the tract is targetable at all
  2. the unit the user typed, at user_asserted provenance

and `tractProvenance` on every response says which fired. A tract built
from a user-typed unit is a hypothesis about the target; one from the
catalogue is a lookup. They are not the same evidence and the response
does not pretend otherwise.

WHAT IS NOT COMPUTED
--------------------
No off-target scan runs, so no off-target repeat count is reported — and
that number matters here more than anywhere else, because a (CAG)n oligo is
complementary to *every* CAG-repeat transcript in the transcriptome by
construction. Saying "3 off-targets" without a scan would be actively
misleading. The response says so instead.

No RBP binding-affinity model is wired, so no displacement score is
emitted. The real duplex energies are reported for comparison once a
measured RBP affinity is available.
"""

from __future__ import annotations

import logging

from services.gene_silencing_service import (
    _calc_gc,
    _calc_tm,
    _cellular_uptake_score,
    _nuclease_resistance_score,
    _reverse_complement,
    _self_complement_mfe,
    _target_duplex_energy,
)
from services import reference_tables as RT
from services.feature_service import (
    KNOWN_REPEAT_UNITS,
    PATHOGENIC_REPEAT_THRESHOLD,
    _extract_repeat_count,
    _normalize_repeat_unit,
)

logger = logging.getLogger(__name__)

NEUTRALIZATION_MECHANISM_CHEMISTRY: dict[str, dict] = {
    # Every TG05 mechanism is occupancy-only: destroying the transcript is a
    # different therapy (that is A1), so a cleaving chemistry is refused.
    "A14": {"allowed": ["pmo", "moe_full_ps", "lna_dna_mixmer"],
            "target": "repeat_tract"},
    "A12": {"allowed": ["moe_full_ps", "lna_dna_mixmer"],
            "target": "mirna"},
}

NEUTRALIZATION_CHEMISTRY_OPTIONS = [
    {"id": "pmo", "label": "PMO (Morpholino)",
     "description": "Non-ionic steric blocker; RNase H-independent."},
    {"id": "moe_full_ps", "label": "2'-O-MOE Full Phosphorothioate",
     "description": "Fully modified steric blocker; standard anti-miR backbone."},
    {"id": "lna_dna_mixmer", "label": "LNA/DNA Mixmer",
     "description": "High-affinity steric blocker."},
]

NEUTRALIZATION_LENGTH_RANGE = {"min": 12, "max": 25, "default": 17, "step": 1}

# How many repeat units of context to build around the oligo when computing
# the duplex, so the energy reflects binding inside a tract rather than to a
# bare oligo-length fragment.
TRACT_CONTEXT_UNITS = 6


# Decorations a repeat unit is commonly written with: "(CAG)n", "CAG x 50",
# "CAG-CAG". These are stripped before validation. Anything else that is a
# letter must be a nucleotide.
_UNIT_DECORATION = str.maketrans("", "", "()[]{}<> \t\r\n-_*x×,.0123456789")


def _strict_repeat_unit(raw: str | None) -> str | None:
    """Accept a repeat unit only if it is genuinely nucleotides.

    The shared _normalize_repeat_unit strips every non-ACGTU character, which
    is right for "(CAG)n" but means arbitrary prose yields a motif:
    "not a motif!" keeps its t, a and t and comes back as "TAT". A tract
    would then be built from a typo and designed against.

    Here the decorations above are removed, a trailing repeat-count "n" is
    dropped, and whatever remains must be entirely nucleotides or the input
    is refused.
    """
    if not raw or not raw.strip():
        return None
    cleaned = raw.strip().translate(_UNIT_DECORATION)
    if cleaned[-1:] in ("n", "N") and len(cleaned) > 1:
        cleaned = cleaned[:-1]          # the "n" of "(CAG)n"
    cleaned = cleaned.upper()
    if not cleaned or not set(cleaned) <= set("ACGTU"):
        return None
    return cleaned.replace("U", "T")


def resolve_repeat_tract(gene_symbol: str, repeat_unit: str | None,
                         estimated_repeat_count: str | None) -> dict:
    """The repeat unit and tract for this target, with its provenance."""
    row = RT.row_for("repeat_expansion_loci", gene_symbol)
    if row and row.get("repeat_unit"):
        unit = _normalize_repeat_unit(row["repeat_unit"])
        if unit:
            return {
                "unit": unit,
                "region": (row.get("transcript_region") or "").strip() or None,
                "pathogenicMin": row.get("pathogenic_min_repeats"),
                "provenance": "confirmed",
                "source": RT.provenance_of(row) or "repeat-expansion catalogue",
                "note": (
                    f"{gene_symbol} is a curated repeat-expansion locus; the "
                    f"unit and its transcript region come from the catalogue."
                ),
            }

    unit = _strict_repeat_unit(repeat_unit)
    if not unit:
        return {
            "unit": None, "region": None, "pathogenicMin": None,
            "provenance": "unavailable", "source": None,
            "note": (
                f"No repeat unit is available for {gene_symbol}: it is not in "
                f"the repeat-expansion catalogue and no valid nucleotide "
                f"motif was supplied. A tract cannot be constructed, so no "
                f"oligo can be designed against it."
            ),
        }

    known = KNOWN_REPEAT_UNITS.get(unit)
    count = _extract_repeat_count(estimated_repeat_count)
    return {
        "unit": unit,
        "region": None,
        "pathogenicMin": None,
        "provenance": "user_asserted",
        "source": "repeat unit supplied on the form",
        "note": (
            f"Tract built from the unit you supplied ({unit}"
            + (f", {known}" if known else ", not in the reference list")
            + "). This is a hypothesis about the target, not a catalogue "
              "lookup — the catalogue has no entry for this gene."
            + (f" Stated expansion ~{count} copies is below the pathogenic "
               f"threshold (~{PATHOGENIC_REPEAT_THRESHOLD}); repeat masking "
               f"targets expanded tracts."
               if count is not None and count < PATHOGENIC_REPEAT_THRESHOLD
               else "")
        ),
    }


def generate_neutralization_candidates(
    gene_symbol: str,
    mechanism_id: str,
    neutralization_mode: str,
    repeat_unit: str | None,
    estimated_repeat_count: str | None,
    oligo_length: int,
    chemistry: str,
    modifications: list[str] | None = None,
    delivery_context: str | None = None,
    target_rbp: str | None = None,
    mirna_sequence: str | None = None,
    max_candidates: int = 12,
) -> dict:
    """Phase-shifted oligos across the target, ranked by real duplex energy."""
    modifications = modifications or []

    if mechanism_id == "A25" or neutralization_mode == "aptamer_decoy":
        return {
            "status": "flagged_not_designed", "candidates": [],
            "message": (
                "A25 (RNA aptamer) is surfaced as a flag, not designed here. "
                "Aptamer selection is structure-based (SELEX) against a "
                "protein surface, not antisense complementarity."
            ),
        }

    spec = NEUTRALIZATION_MECHANISM_CHEMISTRY.get(mechanism_id)
    if spec is None:
        return {"status": "unsupported_mechanism", "candidates": [],
                "message": f"{mechanism_id} is not a designable TG05 mechanism."}
    if chemistry not in spec["allowed"]:
        return {
            "status": "incompatible_chemistry", "candidates": [],
            "message": (
                f"{chemistry} is not usable for {mechanism_id}. TG05 "
                f"neutralization is occupancy-only — a cleaving chemistry "
                f"destroys the transcript, which is a different therapy. "
                f"Allowed: {', '.join(spec['allowed'])}."
            ),
        }

    # --- A12: needs the actual miRNA sequence --------------------------------
    if mechanism_id == "A12":
        seq = (mirna_sequence or "").upper().replace("T", "U")
        seq = "".join(b for b in seq if b in "ACGU")
        if not seq:
            return {
                "status": "target_unavailable", "candidates": [],
                "message": (
                    "An anti-miR is the reverse complement of a specific "
                    "mature miRNA. No miRNA sequence was supplied and no "
                    "miRBase lookup is wired, so there is nothing to design "
                    "against. Supply the mature sequence."
                ),
            }
        target_seq, tract = seq, {
            "unit": None, "region": "mature miRNA", "provenance": "user_asserted",
            "source": "miRNA sequence supplied on the form",
            "note": "Anti-miR designed against the supplied mature sequence.",
        }
    else:
        # --- A14: build the tract from a resolved unit -----------------------
        tract = resolve_repeat_tract(gene_symbol, repeat_unit,
                                     estimated_repeat_count)
        if not tract["unit"]:
            return {"status": "target_unavailable", "candidates": [],
                    "tractProvenance": tract, "message": tract["note"]}
        unit = tract["unit"]
        reps = max(TRACT_CONTEXT_UNITS,
                   (oligo_length * 3) // len(unit) + TRACT_CONTEXT_UNITS)
        target_seq = (unit * reps).replace("T", "U")

    if len(target_seq) < oligo_length:
        return {"status": "target_too_short", "candidates": [],
                "tractProvenance": tract,
                "message": (f"The target is {len(target_seq)} nt, shorter than "
                            f"the {oligo_length} nt oligo requested.")}

    # A repeat tract is periodic, so only the phases within one unit are
    # distinct designs — tiling further just repeats them.
    unit_len = len(tract["unit"]) if tract.get("unit") else 1
    n_phases = (unit_len if mechanism_id == "A14"
                else max(1, len(target_seq) - oligo_length + 1))

    seen: set[str] = set()
    candidates = []
    for phase in range(min(n_phases, len(target_seq) - oligo_length + 1)):
        site = target_seq[phase:phase + oligo_length]
        if len(site) < oligo_length or site in seen:
            continue
        seen.add(site)

        aso = _reverse_complement(site.replace("U", "T"))
        gc = _calc_gc(aso)
        tm = _calc_tm(aso)
        duplex_dg = _target_duplex_energy(aso, site)
        self_mfe = _self_complement_mfe(aso)

        candidates.append({
            "sequence": aso,
            "sequenceAlphabet": "DNA",
            "targetSite": site,
            "targetSiteAlphabet": "RNA",
            "phase": phase,
            "tilingPattern": (
                f"{oligo_length}-mer, phase {phase} of {unit_len}"
                if mechanism_id == "A14"
                else f"{oligo_length}-mer anti-miR"
            ),
            "mechanismId": mechanism_id,
            "chemistry": chemistry,
            "modifications": modifications,
            "deliveryContext": delivery_context,
            "realMetrics": {
                "targetDuplexEnergy": round(duplex_dg, 2),
                "meltingTempC": round(tm, 1),
                "selfStructureMfe": round(self_mfe, 2),
                "gcContent": round(gc, 3),
                "lengthNt": oligo_length,
                "provenance": "ViennaRNA duplexfold / fold; nearest-neighbour Tm",
            },
            "heuristicEstimates": {
                "nucleaseResistance": _nuclease_resistance_score(chemistry, modifications),
                "cellularUptake": _cellular_uptake_score(chemistry, oligo_length),
                "provenance": "Ordinal chemistry and length heuristics, not measurements.",
            },
        })

    if not candidates:
        return {"status": "no_candidates", "candidates": [],
                "tractProvenance": tract,
                "message": "No distinct oligo phase could be built."}

    candidates.sort(key=lambda c: c["realMetrics"]["targetDuplexEnergy"])
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1

    notes = []
    if mechanism_id == "A14":
        notes.append(
            "OFF-TARGET RISK IS NOT QUANTIFIED. An oligo complementary to a "
            f"({tract['unit']})n tract is complementary to every transcript "
            "carrying that repeat, by construction — this is the central "
            "selectivity problem for repeat masking. No transcriptome scan "
            "is wired, so no off-target count is reported. Treat every "
            "candidate as requiring an experimental specificity check."
        )
    if target_rbp and target_rbp.strip():
        notes.append(
            f"You named {target_rbp.strip()} as the sequestered protein. No "
            f"binding-affinity model for it is wired, so no displacement "
            f"score is reported — a competition number without the RBP's own "
            f"affinity would be invented. Compare against the duplex "
            f"energies above once a measured affinity is available."
        )

    return {
        "status": "ok",
        "mechanismId": mechanism_id,
        "neutralizationMode": neutralization_mode,
        "tractProvenance": tract,
        "targetSequence": target_seq[:120],
        "candidates": candidates[:max_candidates],
        "notes": notes,
        "scoringNote": (
            "Ordered by ViennaRNA duplex free energy. realMetrics are "
            "computed from the sequence; heuristicEstimates are ordinal rules "
            "of thumb. No efficacy or displacement score is emitted."
        ),
    }


def get_neutralization_design_options() -> dict:
    return {
        "chemistries": NEUTRALIZATION_CHEMISTRY_OPTIONS,
        "lengthRange": NEUTRALIZATION_LENGTH_RANGE,
        "knownRepeatUnits": [
            {"unit": u, "disease": d} for u, d in sorted(KNOWN_REPEAT_UNITS.items())
        ],
        "pathogenicThreshold": PATHOGENIC_REPEAT_THRESHOLD,
    }
