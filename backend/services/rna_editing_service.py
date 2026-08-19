"""
RNA Editing / Correction — guide RNA candidate generation for TG03.

Generates guide RNA candidates for ADAR (A-to-I), APOBEC (C-to-U), and
SMaRT (trans-splicing) mechanisms. Reuses target analysis from
gene_silencing_service and biophysical helpers.
"""

from __future__ import annotations

import re
from typing import Any

import RNA
import primer3

from services.gene_silencing_service import (
    get_target_analysis,
    _calc_gc,
    _calc_tm,
    _self_complement_mfe,
    _polyg_score,
    _cpg_count,
    _longest_homopolymer,
    _purine_content,
    _sequence_complexity,
    _gc_skew,
    _molecular_weight,
    _extinction_coefficient,
    _reverse_complement,
)


# ---------------------------------------------------------------------------
# Chemistry & modification options for TG03
# ---------------------------------------------------------------------------

EDITING_CHEMISTRY_OPTIONS = [
    {
        "id": "2ome_ps",
        "label": "2'-O-Methyl + PS",
        "description": "Fully modified 2'-O-Me with phosphorothioate backbone. High nuclease stability; non-cleaving backbone protects guide RNA.",
        "detail": "The gold standard for guide RNAs in ADAR recruitment. 2'-O-Me modifications increase binding affinity and resist nuclease degradation without triggering RNase H-mediated cleavage of the target mRNA.",
    },
    {
        "id": "moe_ps",
        "label": "2'-O-MOE PS Oligo",
        "description": "High thermodynamic affinity (Tm), ideal for long LEAPER/ADAR recruiting guide RNAs.",
        "detail": "2'-O-Methoxyethyl (MOE) provides the highest binding affinity per modification, making it ideal for longer guide RNAs (70-120 nt) where maintaining duplex stability is critical.",
    },
    {
        "id": "stereopure_ps",
        "label": "Stereopure Phosphorothioate",
        "description": "Precise stereochemistry optimization for recruiting endogenous ADAR without cytotoxicity.",
        "detail": "Stereopure PS linkages with defined chirality (Sp/Rp) offer improved pharmacokinetics and reduced off-target immune activation compared to mixed-chirality PS backbones.",
    },
]

EDITING_MODIFICATION_OPTIONS = [
    {
        "id": "phosphorothioate",
        "label": "Phosphorothioate Backbone",
        "description": "PS linkages at guide RNA termini for nuclease resistance and protein binding.",
    },
    {
        "id": "2ome_wings",
        "label": "2'-O-Me Wing Modifications",
        "description": "Additional 2'-O-Me at guide RNA ends for enhanced stability.",
    },
    {
        "id": "lna_wings",
        "label": "LNA Wing Modifications",
        "description": "Locked nucleic acid at guide RNA ends for maximum binding affinity.",
    },
]

EDITING_LENGTH_RANGE = {"min": 30, "max": 120, "default": 71, "step": 1}


# ---------------------------------------------------------------------------
# HGVS parsing helpers
# ---------------------------------------------------------------------------

def _parse_hgvs_position(hgvs: str) -> dict[str, Any] | None:
    """Extract position and base change from a simple HGVS substitution.

    Handles both exonic variants (c.82G>A) and intronic/splice-site
    variants (c.10921+2T>C, c.10921-1G>T).

    Returns {"position": int, "ref": str, "alt": str, "type": str, "intron_offset": int|None} or None.
    """
    hgvs = hgvs.strip()
    # Match: c.82G>A, g.12345C>T, n.100A>G (with optional transcript prefix)
    # Also match intronic variants: c.10921+2T>C, c.10921-1G>T
    m = re.search(r"[cgn]\.(\d+)([+-]\d+)?([ATGC])>([ATGC])", hgvs, re.IGNORECASE)
    if not m:
        return None
    intron_offset = None
    if m.group(2):
        intron_offset = int(m.group(2))
    return {
        "position": int(m.group(1)),
        "ref": m.group(3).upper().replace("T", "U"),
        "alt": m.group(4).upper().replace("T", "U"),
        "type": hgvs[0].lower() if hgvs[0] in "cgnCNG" else "c",
        "intronOffset": intron_offset,
    }


def _editing_base_compatible(alt_base: str, edit_type: str) -> bool:
    """Check if a variant's mutant base is compatible with the editing modality.

    RNA editing acts on the CURRENT (mutant) base in the transcript.
    ADAR deaminates an adenosine (A->I, read as G) so the mutant base must be A.
    APOBEC deaminates a cytidine (C->U) so the mutant base must be C.
    """
    alt = alt_base.upper().replace("T", "U")
    if edit_type == "a_to_i":
        return alt == "A"  # ADAR requires adenosine at mutant position
    if edit_type == "c_to_u":
        return alt == "C"  # APOBEC requires cytidine at mutant position
    return True  # trans-splicing accepts any


# ---------------------------------------------------------------------------
# Exon boundary / splice junction mapping (for SMaRT targeting)
# ---------------------------------------------------------------------------

def _compute_exon_cds_ranges(exons: list[dict]) -> list[tuple[int, int, int, int]]:
    """Compute exon boundaries in CDS / mRNA coordinates (0-indexed).

    Assumes the mRNA sequence is the concatenation of exon sequences.
    Returns ``(exon_index, cds_start, cds_end, exon_length)`` for each exon.
    """
    ranges: list[tuple[int, int, int, int]] = []
    cumulative = 0
    for exon in exons:
        exon_len = exon.get("length") or 0
        idx = exon.get("index") or 0
        ranges.append((idx, cumulative, cumulative + exon_len, exon_len))
        cumulative += exon_len
    return ranges


def _find_variant_exon(
    exons_cds: list[tuple[int, int, int, int]],
    position: int,
) -> tuple[int, int, int, int] | None:
    """Find the exon that contains the given CDS position.

    Returns ``(exon_index, cds_start, cds_end, exon_length)`` or ``None``.
    """
    for idx, cds_start, cds_end, length in exons_cds:
        if cds_start <= position < cds_end:
            return (idx, cds_start, cds_end, length)
    # Fallback: nearest exon
    if exons_cds:
        return exons_cds[-1]
    return None


def _resolve_splice_junction_position(
    exons: list[dict],
    var_position: int,
    splicing_direction: str,
    guide_length: int,
    mrna_len: int,
) -> tuple[int, int, str]:
    """Determine the splice-junction target position and search base offset
    for SMaRT trans-splicing guides.

    Returns ``(junction_pos, base_offset, region_label)``.
    """
    exons_cds = _compute_exon_cds_ranges(exons)
    variant_exon = _find_variant_exon(exons_cds, var_position)

    if not variant_exon:
        # No exon data — fall back to variant position
        junction_pos = var_position
        base_offset = max(0, min(var_position, mrna_len - guide_length))
        return junction_pos, base_offset, "trans-splicing junction"

    _exon_idx, exon_start, exon_end, _exon_len = variant_exon

    if splicing_direction == "five_prime":
        # 5' exon replacement: the guide ABD spans the 5' boundary of the
        # variant's exon (where the upstream 3'SS acceptor lives).
        junction_pos = exon_start
        base_offset = junction_pos
        region_label = "5' Splice Junction (Acceptor → Donor)"
    else:
        # three_prime (default): 3' exon replacement — the guide ABD spans
        # the 3' boundary of the variant's exon (where the downstream 5'SS
        # donor lives).
        junction_pos = exon_end
        base_offset = max(0, junction_pos - guide_length)
        region_label = "3' Splice Junction (Donor → Acceptor)"

    return junction_pos, base_offset, region_label


# ---------------------------------------------------------------------------
# Editing efficiency scoring
# ---------------------------------------------------------------------------

def _score_editing_efficiency(
    candidate_seq: str,
    target_seq: str,
    edit_position: int,
    edit_type: str,
    mismatch_pocket: str = "c",
) -> float:
    """Score predicted on-target editing efficiency (0-100).

    Based on:
    - Duplex stability (ΔG of guide-target binding)
    - Mismatch pocket quality (A-C mismatch boosts ADAR)
    - Sequence context (flanking bases affect deamination rate)
    """
    if not candidate_seq or not target_seq:
        return 0.0

    # Base score from duplex stability
    duplex = RNA.duplexfold(candidate_seq.upper(), target_seq.upper())
    duplex_score = max(0, min(50, -duplex.energy * 2))

    # Mismatch pocket bonus (for A-to-I: C at opposing position boosts efficiency)
    pocket_bonus = 0
    if edit_type == "a_to_i" and mismatch_pocket.upper() == "C":
        pocket_bonus = 20  # A-C mismatch is optimal for ADAR
    elif edit_type == "c_to_u":
        pocket_bonus = 15  # C-to-U is generally efficient

    # Context bonus (bases flanking the edit site)
    ctx_start = max(0, edit_position - 2)
    ctx_end = min(len(target_seq), edit_position + 3)
    context = target_seq[ctx_start:ctx_end].upper()
    context_bonus = 0
    if "AG" in context:  # ADAR prefers AG context
        context_bonus = 10
    elif "GA" in context:
        context_bonus = 5

    score = duplex_score + pocket_bonus + context_bonus
    return round(max(0, min(100, score)), 1)


def _count_bystander_adenosines(
    candidate_seq: str,
    target_seq: str,
    edit_position: int,
    window_radius: int = 20,
) -> list[dict]:
    """Find adenosines within the binding window that could be bystander-edited.

    Returns list of {position, risk, context} dicts.
    """
    bystanders = []
    if not target_seq or edit_position is None:
        return bystanders

    window_start = max(0, edit_position - window_radius)
    window_end = min(len(target_seq), edit_position + window_radius + 1)

    for i in range(window_start, window_end):
        if i == edit_position:
            continue  # Skip the target position
        if target_seq[i].upper() == "A":
            # Risk based on distance from target and sequence context
            distance = abs(i - edit_position)
            if distance <= 5:
                risk = "high"
            elif distance <= 12:
                risk = "medium"
            else:
                risk = "low"

            ctx_start = max(0, i - 2)
            ctx_end = min(len(target_seq), i + 3)
            context = target_seq[ctx_start:ctx_end]

            bystanders.append({
                "position": i,
                "risk": risk,
                "context": context,
            })

    return bystanders


def _score_adar_recruitment(candidate_seq: str, target_seq: str) -> float:
    """Score endogenous ADAR recruitment potential (0-100).

    Based on dsRNA duplex length, stability, and structure.
    """
    if not candidate_seq or not target_seq:
        return 0.0

    # Duplex energy
    duplex = RNA.duplexfold(candidate_seq.upper(), target_seq.upper())
    energy_score = max(0, min(40, -duplex.energy * 1.5))

    # Length bonus (longer duplexes recruit ADAR better)
    length = min(len(candidate_seq), len(target_seq))
    length_score = max(0, min(30, (length - 30) * 0.5))

    # Self-structure penalty (if guide folds on itself, it can't bind target)
    self_mfe = _self_complement_mfe(candidate_seq)
    structure_penalty = max(0, min(30, -self_mfe * 3))

    score = energy_score + length_score + structure_penalty
    return round(max(0, min(100, score)), 1)


# ---------------------------------------------------------------------------
# SMaRT trans-splicing scoring helpers
# ---------------------------------------------------------------------------

# Consensus splice donor (5'SS) and acceptor (3'SS) motifs
_SPLICE_DONOR = "GTAAGT"   #Canonical 5' splice site donor
_SPLICE_ACCEPTOR = "YNCAGG"  # Simplified 3' splice site acceptor (Y=C/T, N=any)
_SPLICE_DONOR_INTRON = "GT"  # Minimal invariant dinucleotide at donor
_SPLICE_ACCEPTOR_INTRON = "AG"  # Minimal invariant dinucleotide at acceptor


def _score_splice_site_strength(
    seq: str,
    direction: str,
) -> float:
    """Score how well the target region matches splice site consensus (0-100).

    For trans-splicing, the guide must span or be adjacent to a splice junction.
    Stronger splice site consensus = more efficient trans-splicing.
    """
    if not seq:
        return 0.0

    score = 0.0
    seq_upper = seq.upper()

    if direction == "five_prime":
        # For 5' exon replacement, the guide targets the 3' end of the upstream exon
        # and should be near the 5' splice site (donor) of the intron
        # Check for GT dinucleotide (intron start) in the region
        if _SPLICE_DONOR_INTRON in seq_upper:
            score += 30
        # Check for extended donor motif
        if _SPLICE_DONOR[:4] in seq_upper:
            score += 20
        # Check for pyrimidine-rich region upstream (typical of acceptor context)
        pyrimidines = sum(1 for c in seq_upper if c in "CT")
        pyrimidine_fraction = pyrimidines / max(1, len(seq_upper))
        score += max(0, min(20, pyrimidine_fraction * 30))
    else:
        # For 3' exon replacement (default), the guide targets the 5' end of the
        # downstream exon and should be near the 3' splice site (acceptor) of the intron
        # Check for AG dinucleotide (intron end) in the region
        if _SPLICE_ACCEPTOR_INTRON in seq_upper:
            score += 30
        # Check for extended acceptor motif
        if "AGG" in seq_upper or "CAG" in seq_upper:
            score += 20
        # Check for polypyrimidine tract context
        pyrimidines = sum(1 for c in seq_upper if c in "CT")
        pyrimidine_fraction = pyrimidines / max(1, len(seq_upper))
        score += max(0, min(20, pyrimidine_fraction * 30))

    # Bonus for GC content in functional range (40-60%)
    gc = _calc_gc(seq_upper)
    if 0.40 <= gc <= 0.60:
        score += 15
    elif 0.35 <= gc <= 0.65:
        score += 8

    # Bonus for moderate self-structure (some structure helps spliceosome recognition)
    self_mfe = _self_complement_mfe(seq_upper)
    if -30 <= self_mfe <= -5:
        score += 15  # Moderate structure is good
    elif self_mfe < -30:
        score += 5   # Too much structure may hinder binding
    else:
        score += 10  # Little structure is acceptable

    return round(max(0, min(100, score)), 1)


def _score_binding_domain(
    candidate_seq: str,
    target_seq: str,
    abd_length: int,
) -> float:
    """Score the antisense binding domain quality for trans-splicing (0-100).

    The ABD must form a stable duplex with the pre-mRNA to recruit the
    spliceosome. Scoring considers duplex stability, length adequacy, and
    thermodynamic properties.
    """
    if not candidate_seq or not target_seq:
        return 0.0

    # Duplex stability (ΔG)
    duplex = RNA.duplexfold(candidate_seq.upper(), target_seq.upper())
    duplex_score = max(0, min(40, -duplex.energy * 1.2))

    # Length adequacy: guide should cover the ABD length
    effective_length = min(len(candidate_seq), len(target_seq))
    length_score = max(0, min(25, (effective_length / max(1, abd_length)) * 25))

    # Tm in optimal range for spliceosomal recognition (55-75°C)
    tm = _calc_tm(candidate_seq)
    if 55 <= tm <= 75:
        tm_score = 20
    elif 45 <= tm <= 85:
        tm_score = 10
    else:
        tm_score = 0

    # Self-structure penalty (guide must be free to bind target)
    self_mfe = _self_complement_mfe(candidate_seq)
    structure_score = max(0, min(15, -self_mfe * 2)) if self_mfe < -10 else 15

    score = duplex_score + length_score + tm_score + structure_score
    return round(max(0, min(100, score)), 1)


def _score_splice_compatibility(
    candidate_seq: str,
    target_seq: str,
    splicing_direction: str,
) -> float:
    """Score compatibility with spliceosomal machinery (0-100).

    Evaluates whether the guide-target duplex would be recognized by the
    spliceosome for trans-splicing. Factors include:
    - Exon definition signals (ESE/ESS)
    - Splice site proximity
    - RNA secondary structure at junction
    """
    if not candidate_seq or not target_seq:
        return 0.0

    score = 0.0
    target_upper = target_seq.upper()
    guide_upper = candidate_seq.upper()

    # Check for exonic splicing enhancer (ESE) motifs in target
    ese_motifs = ["GAAGAA", "GACAGA", "TTCGAG", "GAGGAG", "CACGTC"]
    for motif in ese_motifs:
        if motif in target_upper:
            score += 8
            break

    # Check for exonic splicing silencer (ESS) motifs — penalize
    ess_motifs = ["TGCATG", "TGACTG", "TCATTC"]
    for motif in ess_motifs:
        if motif in target_upper:
            score -= 5

    # Bonus for direction-specific junction proximity
    if splicing_direction == "five_prime":
        # 5' exon replacement: guide should end near exon-intron boundary
        if guide_upper[-3:] in ("GTA", "GTC", "GTG", "GTT"):
            score += 15  # Ends near GT donor
    else:
        # 3' exon replacement: guide should start near intron-exon boundary
        if guide_upper[:3] in ("AGG", "AGC", "AGA", "AGT"):
            score += 15  # Starts near AG acceptor

    # Duplex energy contribution
    duplex = RNA.duplexfold(guide_upper, target_upper)
    score += max(0, min(30, -duplex.energy * 1.0))

    # Complexity bonus (low complexity sequences may form R-loops)
    complexity = _sequence_complexity(candidate_seq)
    score += max(0, min(15, complexity * 20))

    return round(max(0, min(100, score)), 1)


def _chemistry_modification_bonus(
    chemistry: str,
    modifications: list[str],
    guide_length: int,
) -> float:
    """Calculate score bonus/penalty based on chemistry and modifications (0-100).

    Different chemistries have different thermodynamic properties, nuclease
    resistance, and spliceosome compatibility.
    """
    score = 50.0  # Baseline

    # Chemistry-specific adjustments
    chemistry_effects = {
        "2ome_ps": {
            "tm_bonus": 2,       # 2'-O-Me increases Tm ~2°C per modification
            "stability": 8,      # Good nuclease resistance
            "splice_compat": 5,  # Compatible with spliceosome
        },
        "moe_ps": {
            "tm_bonus": 4,       # MOE increases Tm ~4°C per modification
            "stability": 10,     # Excellent nuclease resistance
            "splice_compat": 3,  # Slightly lower spliceosome compatibility (higher affinity may resist remodeling)
        },
        "stereopure_ps": {
            "tm_bonus": 2,       # Similar to 2'-O-Me
            "stability": 9,      # Good with defined stereochemistry
            "splice_compat": 7,  # Best spliceosome compatibility (precise chirality)
        },
    }

    chem = chemistry_effects.get(chemistry, chemistry_effects["2ome_ps"])
    score += chem["stability"]
    score += chem["splice_compat"]

    # Length-dependent chemistry bonus
    if guide_length > 100:
        # Longer guides benefit more from MOE (higher Tm)
        if chemistry == "moe_ps":
            score += 8
        elif chemistry == "stereopure_ps":
            score += 5

    # Modification-specific adjustments
    mod_effects = {
        "phosphorothioate": 3,   # Standard nuclease resistance
        "2ome_wings": 5,         # Enhanced stability at termini
        "lna_wings": 8,          # Maximum binding affinity at termini
    }
    for mod in (modifications or []):
        score += mod_effects.get(mod, 0)

    # Cap at 100
    return round(max(0, min(100, score)), 1)


# ---------------------------------------------------------------------------
# Chemistry / modification metadata helpers
# ---------------------------------------------------------------------------

# Per-chemistry estimated Tm lift (°C) over the unmodified duplex
_CHEM_TM_LIFTS = {
    "2ome_ps": 8.0,
    "moe_ps": 12.0,
    "stereopure_ps": 5.0,
}

# Per-modification estimated Tm lift (°C)
_MOD_TM_LIFTS = {
    "phosphorothioate": 0.0,
    "2ome_wings": 3.0,
    "lna_wings": 5.0,
}

# Nuclease-resistance baseline per chemistry (0-100)
_CHEM_NUCLEASE_RESISTANCE = {
    "2ome_ps": 85,
    "moe_ps": 90,
    "stereopure_ps": 80,
}

# Nuclease-resistance bonus per modification (0-100)
_MOD_NUCLEASE_BONUS = {
    "phosphorothioate": 5,
    "2ome_wings": 8,
    "lna_wings": 10,
}


def _adjust_tm_for_chemistry(
    tm: float,
    chemistry: str,
    modifications: list[str],
) -> tuple[float, float]:
    """Return ``(adjusted_tm_celsius, binding_affinity_adjustment_celsius)``.

    The adjustment approximates the cumulative Tm lift contributed by the
    chosen chemistry (full-backbone ribose modifications) and wing
    modifications.  Values are estimates consistent with the per-modification
    Tm bonuses documented in :func:`_chemistry_modification_bonus`.
    """
    adjustment = _CHEM_TM_LIFTS.get(chemistry, 5.0)
    for mod in (modifications or []):
        adjustment += _MOD_TM_LIFTS.get(mod, 0.0)
    adjusted = round(tm + adjustment, 1)
    return adjusted, round(adjustment, 1)


def _calculate_nuclease_resistance(
    chemistry: str,
    modifications: list[str],
) -> int:
    """Score (0-100) how well the guide resists exonuclease/nuclease degradation."""
    score = _CHEM_NUCLEASE_RESISTANCE.get(chemistry, 50)
    for mod in (modifications or []):
        score += _MOD_NUCLEASE_BONUS.get(mod, 0)
    return min(100, max(0, score))


def _splicing_efficiency_score(
    splice_site_score: float,
    splice_compat_score: float,
) -> float:
    """Overall SMaRT splicing efficiency (0-100).

    Combines splice-site strength with spliceosome compatibility to predict
    how efficiently trans-splicing will occur at the target junction.
    """
    return round(max(0, min(100, splice_site_score * 0.55 + splice_compat_score * 0.45)), 1)


def _smart_mechanism_notes(
    chemistry: str,
    modifications: list[str],
    splicing_direction: str,
    abd_length: int,
) -> str:
    """Human-readable mechanism notes for SMaRT trans-splicing candidates."""
    dir_label = (
        "5' exon replacement (targets 5' splice donor junction)"
        if splicing_direction == "five_prime"
        else "3' exon replacement (targets 3' splice acceptor junction)"
    )
    chem_label = next(
        (c["label"] for c in EDITING_CHEMISTRY_OPTIONS if c["id"] == chemistry),
        chemistry,
    )
    mod_label = ", ".join(m.replace("_", " ") for m in (modifications or [])) or "none"
    return (
        f"SMaRT {dir_label}. ABD={abd_length} nt • chemistry={chem_label} • "
        f"modifications: {mod_label}."
    )


def _editing_mechanism_notes(edit_type: str, chemistry: str) -> str:
    """Human-readable mechanism notes for ADAR / APOBEC candidates."""
    chem_label = next(
        (c["label"] for c in EDITING_CHEMISTRY_OPTIONS if c["id"] == chemistry),
        chemistry,
    )
    if edit_type == "a_to_i":
        return f"ADAR A→I editing via guide RNA recruitment. Chemistry: {chem_label}."
    if edit_type == "c_to_u":
        return f"APOBEC C→U editing via guide RNA recruitment. Chemistry: {chem_label}."
    return ""


# ---------------------------------------------------------------------------
# Main candidate generation
# ---------------------------------------------------------------------------

# Which mechanism this service can actually build a guide for, and under
# which edit type.
#
# `mechanism_id` used to be a passenger: it was echoed into every candidate
# and changed nothing. Asking for A19 returned an ADAR-recruiting guide
# labelled "A19", which is a REPAIR crRNA in name only — REPAIR guides carry
# a Cas13b direct-repeat scaffold this service never emits. The design is
# driven by `edit_type`, which is correct; what was missing was a check that
# the requested mechanism is one this architecture produces.
#
# A13 and A17 both recruit ENDOGENOUS ADAR with a linear antisense guide and
# differ mainly in guide length (SDRE ~25-40 nt, LEAPER arRNA ~71-111 nt),
# which the `guide_length` parameter already covers. A16 is the C-to-U
# counterpart. A20 is built by the trans-splicing path, which emits a
# pre-trans-splicing molecule rather than a deaminase guide.
EDITING_MECHANISM_EDIT_TYPES: dict[str, set[str]] = {
    "A13": {"a_to_i"},
    "A16": {"c_to_u"},
    "A17": {"a_to_i"},
    "A20": {"trans_splicing"},
}


def _validate_editing_mechanism(mechanism_id: str, edit_type: str) -> None:
    """Refuse to label a guide with a mechanism this service does not build."""
    allowed = EDITING_MECHANISM_EDIT_TYPES.get(mechanism_id)
    if allowed is None:
        raise ValueError(
            f"{mechanism_id} is not designable by this service. It builds "
            f"linear guides that recruit an endogenous deaminase, plus "
            f"trans-splicing molecules; mechanisms needing a delivered "
            f"protein and a scaffolded crRNA (CIRTS, REPAIR) are out of "
            f"scope. Designable here: "
            f"{', '.join(sorted(EDITING_MECHANISM_EDIT_TYPES))}."
        )
    if edit_type not in allowed:
        raise ValueError(
            f"{mechanism_id} does not perform '{edit_type}' editing. It "
            f"supports: {', '.join(sorted(allowed))}."
        )


def generate_guide_rna_candidates(
    ensembl_gene_id: str,
    variant_hgvs: str,
    edit_type: str,
    guide_length: int = 71,
    chemistry: str = "2ome_ps",
    modifications: list[str] | None = None,
    mismatch_pocket: str = "c",
    max_bystander_edits: int = 0,
    splicing_direction: str | None = None,
    abd_length: int = 150,
    delivery_context: str | None = None,
    mechanism_id: str = "A13",
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
) -> list[dict]:
    """Generate guide RNA candidates for RNA editing mechanisms.

    Scans the target mRNA region around the variant position and generates
    complementary guide RNA sequences that recruit editing machinery.

    For SMaRT trans-splicing (``edit_type == "trans_splicing"``) the
    ``splicing_direction`` parameter controls *where* each guide is anchored
    relative to the exon-intron boundary, and ``abd_length`` controls the
    actual length of the antisense binding domain (i.e. the guide length).
    """
    _validate_editing_mechanism(mechanism_id, edit_type)

    if modifications is None:
        modifications = ["phosphorothioate"]

    candidates = []

    # Get target mRNA sequence and exon data
    target = get_target_analysis(ensembl_gene_id, gene_symbol, organism)
    mrna_seq = target.get("mrnaSequence", "")
    exons = target.get("exons", [])
    if not mrna_seq:
        return candidates

    mrna_seq = mrna_seq.upper()
    mrna_len = len(mrna_seq)

    # Parse the variant to get position
    var_info = _parse_hgvs_position(variant_hgvs)
    if not var_info:
        return candidates

    var_position = var_info["position"]
    var_ref = var_info["ref"]
    var_alt = var_info["alt"]

    # Check if the variant base is compatible with the editing modality
    if not _editing_base_compatible(var_alt, edit_type):
        return candidates  # Incompatible: e.g., APOBEC needs C at mutant position

    # --- Trans-splicing setup ---
    is_trans_splicing = edit_type == "trans_splicing"
    if is_trans_splicing:
        # Normalise direction (default to 3' exon replacement)
        if not splicing_direction or splicing_direction == "":
            splicing_direction = "three_prime"
        # For SMaRT, the ABD length *is* the guide length — the entire guide
        # acts as the antisense binding domain.
        effective_guide_length = abd_length
    else:
        effective_guide_length = guide_length

    # --- Determine search window / junction anchor ---
    if is_trans_splicing and exons:
        junction_pos, base_offset, junction_label = _resolve_splice_junction_position(
            exons, var_position, splicing_direction, effective_guide_length, mrna_len
        )
    else:
        junction_pos = var_position
        target_center = min(var_position, mrna_len - effective_guide_length)
        base_offset = target_center
        junction_label = ""

    search_radius = 30  # Search ±30 nt around the anchor position
    search_start = max(0, base_offset - search_radius)
    search_end = min(mrna_len - effective_guide_length, base_offset + search_radius)

    # Chemistry-adjusted metadata (computed once)
    nuclease_resistance = _calculate_nuclease_resistance(chemistry, modifications)

    seen = set()
    step = max(1, effective_guide_length // 5)

    for offset in range(search_start, search_end + 1, step):
        # The guide RNA is complementary to the target region
        target_region = mrna_seq[offset:offset + effective_guide_length]
        if len(target_region) < effective_guide_length:
            continue

        guide_seq = _reverse_complement(target_region)
        if guide_seq in seen:
            continue
        seen.add(guide_seq)

        # Calculate biophysical properties
        gc = _calc_gc(guide_seq)
        tm_raw = _calc_tm(guide_seq)
        self_mfe = _self_complement_mfe(guide_seq)
        pg = _polyg_score(guide_seq)
        cpg = _cpg_count(guide_seq)

        # Chemistry-adjusted Tm and metadata
        adjusted_tm, tm_adjustment = _adjust_tm_for_chemistry(
            tm_raw, chemistry, modifications
        )
        duplex_energy = round(
            RNA.duplexfold(guide_seq.upper(), target_region.upper()).energy, 2
        )

        # Quality score components
        gc_score = max(0, 100 - abs(gc - 0.50) * 400)
        tm_score = max(0, 100 - abs(tm_raw - 72) * 2)  # Optimal Tm ~72°C for guide RNA
        mfe_penalty = max(0, min(100, -self_mfe * 10))

        junction_offset = offset - base_offset

        if is_trans_splicing:
            # ---- SMaRT trans-splicing specific scoring ----
            splice_site_score = _score_splice_site_strength(
                target_region, splicing_direction
            )
            binding_domain_score = _score_binding_domain(
                guide_seq, target_region, abd_length
            )
            splice_compat_score = _score_splice_compatibility(
                guide_seq, target_region, splicing_direction
            )
            chem_mod_score = _chemistry_modification_bonus(
                chemistry, modifications, effective_guide_length
            )

            # SMaRT-specific composite scores (replace ADAR metrics)
            splicing_efficiency = _splicing_efficiency_score(
                splice_site_score, splice_compat_score
            )
            spliceosome_recruitment = round(
                max(0, min(100, binding_domain_score * 0.60 + chem_mod_score * 0.40)), 1
            )

            quality = max(0, min(100, (
                splice_site_score * 0.30
                + binding_domain_score * 0.25
                + splice_compat_score * 0.20
                + chem_mod_score * 0.10
                + gc_score * 0.10
                - mfe_penalty * 0.05
            )))

            if splicing_direction == "five_prime":
                region_label = (
                    f"5' Exon Replacement · Junction offset {junction_offset:+d}"
                )
                target_bp = "5' Splice Donor Junction"
            else:
                region_label = (
                    f"3' Exon Replacement · Junction offset {junction_offset:+d}"
                )
                target_bp = "3' Splice Acceptor Junction"

            candidate = {
                "sequence": guide_seq,
                "guideLength": effective_guide_length,
                "targetBasePair": target_bp,
                "meltingTempC": tm_raw,
                "adjustedTmC": adjusted_tm,
                "bindingAffinityAdjustment": tm_adjustment,
                "gcContent": round(gc * 100, 1),
                "gcScore": round(gc_score, 1),
                "tmScore": round(tm_score, 1),
                "selfStructureMfe": self_mfe,
                "mfePenalty": round(mfe_penalty, 1),
                "targetDuplexEnergy": duplex_energy,
                "nucleaseResistance": nuclease_resistance,
                "polygTracts": pg,
                "cpgCount": cpg,
                "longestHomopolymer": _longest_homopolymer(guide_seq),
                "purineContent": _purine_content(guide_seq),
                "sequenceComplexity": _sequence_complexity(guide_seq),
                "gcSkew": _gc_skew(guide_seq),
                "molecularWeight": _molecular_weight(guide_seq),
                "extinctionCoefficient": _extinction_coefficient(guide_seq),
                "qualityScore": round(quality, 1),
                "targetRegion": region_label,
                "mechanismId": mechanism_id,
                "chemistry": chemistry,
                "modifications": modifications,
                "editType": edit_type,
                "variantHgvs": variant_hgvs,
                "mechanismNotes": _smart_mechanism_notes(
                    chemistry, modifications, splicing_direction, abd_length
                ),
                # ADAR metrics are null for trans-splicing
                "onTargetEditScore": None,
                "bystanderRiskCount": 0,
                "bystanderRiskDetails": [],
                "adarRecruitmentScore": None,
                # SMaRT trans-splicing specific scores
                "spliceSiteScore": round(splice_site_score, 1),
                "bindingDomainScore": round(binding_domain_score, 1),
                "spliceCompatibilityScore": round(splice_compat_score, 1),
                "chemistryModificationScore": round(chem_mod_score, 1),
                "splicingEfficiencyScore": splicing_efficiency,
                "spliceosomeRecruitmentScore": spliceosome_recruitment,
                "splicingDirection": splicing_direction,
                "abdLength": abd_length,
                "spliceJunctionPosition": junction_pos,
                "junctionOffset": junction_offset,
                "spliceJunctionLabel": junction_label,
            }
        else:
            # ---- Standard ADAR / APOBEC scoring ----
            edit_rel_position = var_position - offset
            if 0 <= edit_rel_position < effective_guide_length:
                # The variant falls within this guide's binding window
                edit_score = _score_editing_efficiency(
                    guide_seq, target_region, edit_rel_position, edit_type, mismatch_pocket
                )
                bystanders = _count_bystander_adenosines(
                    guide_seq, target_region, edit_rel_position, window_radius=20
                )
                adar_score = _score_adar_recruitment(guide_seq, target_region)

                # Determine target base pair type
                if edit_type == "a_to_i":
                    target_bp = f"A:{mismatch_pocket.upper()} Mismatch"
                elif edit_type == "c_to_u":
                    target_bp = "C:U Deamination"
                else:
                    target_bp = "Editing Site"
            else:
                # Guide doesn't cover the variant position — lower score
                edit_score = 20.0
                bystanders = []
                adar_score = 30.0
                target_bp = "No Variant Coverage"

            # Count high-risk bystanders
            high_risk_count = sum(1 for b in bystanders if b["risk"] == "high")

            quality = max(0, min(100, (
                edit_score * 0.40
                + adar_score * 0.20
                + gc_score * 0.15
                + tm_score * 0.15
                - mfe_penalty * 0.10
                - high_risk_count * 5
            )))
            relative_pos = offset - base_offset
            region_label = f"Target offset {relative_pos:+d}"

            candidate = {
                "sequence": guide_seq,
                "guideLength": effective_guide_length,
                "targetBasePair": target_bp,
                "meltingTempC": tm_raw,
                "adjustedTmC": adjusted_tm,
                "bindingAffinityAdjustment": tm_adjustment,
                "gcContent": round(gc * 100, 1),
                "gcScore": round(gc_score, 1),
                "tmScore": round(tm_score, 1),
                "selfStructureMfe": self_mfe,
                "mfePenalty": round(mfe_penalty, 1),
                "targetDuplexEnergy": duplex_energy,
                "nucleaseResistance": nuclease_resistance,
                "onTargetEditScore": edit_score,
                "bystanderRiskCount": high_risk_count,
                "bystanderRiskDetails": bystanders[:10],  # Limit to top 10
                "adarRecruitmentScore": adar_score,
                "polygTracts": pg,
                "cpgCount": cpg,
                "longestHomopolymer": _longest_homopolymer(guide_seq),
                "purineContent": _purine_content(guide_seq),
                "sequenceComplexity": _sequence_complexity(guide_seq),
                "gcSkew": _gc_skew(guide_seq),
                "molecularWeight": _molecular_weight(guide_seq),
                "extinctionCoefficient": _extinction_coefficient(guide_seq),
                "qualityScore": round(quality, 1),
                "targetRegion": region_label,
                "mechanismId": mechanism_id,
                "chemistry": chemistry,
                "modifications": modifications,
                "editType": edit_type,
                "variantHgvs": variant_hgvs,
                "mechanismNotes": _editing_mechanism_notes(edit_type, chemistry),
                # SMaRT fields are null for ADAR/APOBEC
                "spliceSiteScore": None,
                "bindingDomainScore": None,
                "spliceCompatibilityScore": None,
                "chemistryModificationScore": None,
                "splicingEfficiencyScore": None,
                "spliceosomeRecruitmentScore": None,
                "splicingDirection": None,
                "abdLength": None,
                "spliceJunctionPosition": None,
                "junctionOffset": None,
                "spliceJunctionLabel": None,
            }

        candidates.append(candidate)

    # Sort by quality score descending
    candidates.sort(key=lambda c: c["qualityScore"], reverse=True)
    return candidates[:10]


# ---------------------------------------------------------------------------
# Design options
# ---------------------------------------------------------------------------

def get_rna_editing_design_options() -> dict:
    """Return available chemistry, modification, and length options for TG03."""
    return {
        "chemistryOptions": EDITING_CHEMISTRY_OPTIONS,
        "modificationOptions": EDITING_MODIFICATION_OPTIONS,
        "lengthRange": EDITING_LENGTH_RANGE,
    }
