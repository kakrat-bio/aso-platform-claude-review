"""Per-goal entry points, kept as thin filters over the unified arbitration.

WHAT CHANGED AND WHY
--------------------
This module used to hold nine separate scoring paths, one per therapeutic
goal, plus four central compatibility tables. All of that is gone. Scoring
now happens once, in `mechanism_arbitration`, over every designable
mechanism; the therapeutic goal is an OUTPUT of that pass, not an input to
it.

The functions below survive because the frontend still routes users by goal.
They are now **filters applied after scoring**, never separate scorers. Two
callers get the same answer for the same transcript regardless of which door
they came through, which is the property that was impossible when nine
mechanism sets overlapping in nine ways each had their own ranking code.

Deleted here, and why:

  DEFECT_COMPATIBILITY, SCOPE_COMPATIBILITY,
  UPREGULATION_DEFECT_COMPATIBILITY, SPLICE_DEFECT_COMPATIBILITY
      These restated each mechanism's own `suitableVariantTypes` and
      `molecularDefect` in a second place, where they could silently diverge
      from the rulebook they were derived from. They now live in each
      mechanism's rule.json under `arbitration`.

  rank_isoform_engineering_mechanisms       (TG07)
      Scored a mechanism set that was a strict subset of TG04's, so two code
      paths over an identical set could only diverge. TG07 survives as a
      display tag; all its mechanisms remain available.

TG06 was retired here too and has since been RESTORED as a scored goal with
seven mechanisms (A2, A5, A6, A27 plus A29 IRES, A30 Kozak, A31 PABP). It is
still served by a filter over the single shared pass rather than by a scorer
of its own — that property is what the retirement was protecting, and it is
kept.

  generate_rna_engineering_candidates and its helpers  (TG09)
      Generated sequences, melting temperatures, folding free energies and
      binding affinities from `hash()` of the input strings. Those were not
      predictions, they were deterministic noise formatted to look like
      measurements. TG09 is now a rulebook lookup with no score and no rank.

See docs/planning/therapeutic_goal_scope_plan_v3.md.
"""

from __future__ import annotations

from .mechanism_arbitration import (  # noqa: F401  (re-exported)
    DELIVERY_CONTEXTS,
    DELIVERY_PRECEDENT,
    DELIVERY_TIER_WEIGHT,
    EVIDENCE_WEIGHT,
    MOLECULAR_DEFECTS,
    RETIRED_AS_SCORING_PARTITION,
    ArbitrationContext,
    arbitrate,
    canonical_defect,
    load_rule,
    parse_hgvs_substitution,
)
from .feature_service import (  # noqa: F401  (re-exported)
    KNOWN_REPEAT_UNITS,
    PATHOGENIC_REPEAT_THRESHOLD,
)

# ---------------------------------------------------------------------------
# Input vocabularies
#
# These are the form options each goal's page renders. They are now VIEWS onto
# the single unified defect vocabulary rather than four independent lists, so
# a term cannot mean one thing on the silencing page and another on the
# neutralization page.
# ---------------------------------------------------------------------------


def _defects(*keys: str) -> dict[str, str]:
    return {k: MOLECULAR_DEFECTS[k] for k in keys}


DEFECT_TYPES = _defects(
    "gain_of_function",
    "overexpression",
    "mirna_dysregulation",
    "viral_toxic_rna",
    "therapeutic_reduction",
)

SILENCING_SCOPES = {
    "total_knockdown": "Total transcript knockdown",
    "allele_specific": "Allele-specific silencing (spare wild-type)",
}

GENE_UPREGULATION_DEFECT_TYPES = _defects(
    "haploinsufficiency",
    "poison_exon_inclusion",
    "nat_mediated_repression",
    "uorf_mediated_repression",
    "mirna_mediated_repression",
    "rbp_mediated_repression",
    "epigenetic_promoter_silencing",
)

SPLICE_DEFECT_TYPES = _defects(
    "exon_skipping_mutation",
    "exon_inclusion_defect",
    "cryptic_splice_site",
    "pseudoexon_activation",
    "apa_dysregulation",
)

NEUTRALIZATION_DEFECT_TYPES = {
    **_defects("toxic_rna_gain_of_function", "rbp_sequestration"),
    # Kept as selectable input for backward compatibility. Both alias onto
    # unified terms; see DEFECT_ALIASES for why each is a genuine synonym
    # rather than a convenience merge.
    "pathogenic_mirna": "Pathogenic microRNA / ncRNA overexpression",
    "loss_of_function": "Pure loss-of-function (haploinsufficiency / null)",
}

EDIT_TYPES = {
    "a_to_i": "A-to-I Editing (ADAR Recruitment)",
    "c_to_u": "C-to-U Editing (APOBEC / RESCUE)",
    "trans_splicing": "Trans-Splicing / Pre-mRNA Repair (SMaRT)",
}

ENZYME_RECRUITMENT = {
    "adar1": "Endogenous ADAR1 (p110/p150)",
    "adar2": "Endogenous ADAR2",
    "exogenous_deaminase": "Exogenous Deaminase (engineered)",
}

MISMATCH_POCKET = {
    "c": "C (A-C Mismatch — High Efficiency)",
    "g": "G",
    "u": "U",
}

SPLICING_DIRECTIONS = {
    "three_prime": "3' Exon Replacement",
    "five_prime": "5' Exon Replacement",
}

INTRON_SITES = {
    "acceptor_junction": "Acceptor Junction",
    "donor_junction": "Donor Junction",
}

NEUTRALIZATION_MODES = {
    "steric_repeat_masking": "Steric Repeat Masking (RNase H-Independent)",
    "microrna_antagomir": "MicroRNA / ncRNA Antagomir",
    "aptamer_decoy": "Aptamer Decoy Sequestration",
}

NEUTRALIZATION_MODE_MECHANISMS = {
    "steric_repeat_masking": ["A14"],
    "microrna_antagomir": ["A12"],
    "aptamer_decoy": ["A25"],
}

STERIC_CHEMISTRIES = {
    "pmo": "PMO (Morpholino)",
    "moe_full_ps": "2'-O-MOE Full Phosphorothioate",
    "lna_dna_mixmer": "LNA / DNA Mixmer",
}

TRANSLATIONAL_GOALS = {
    "enhance": "Enhance Translation (Upregulate Protein)",
    "suppress": "Suppress Translation (Downregulate Protein)",
}

TRANSLATIONAL_TARGET_ELEMENTS = {
    "5p_utr": "5' UTR / Kozak Sequence",
    "3p_utr_mirna": "3' UTR miRNA Seed Site",
    "uorf": "5' UTR uORF / Upstream AUG",
    "structured_element": "IRES / G-quadruplex / Riboswitch",
    # TG06 restoration — the three elements the new mechanisms target.
    "ires_element": "IRES Domain (cap-independent initiation)",
    "kozak_consensus": "Kozak Consensus (initiator AUG context)",
    "polya_site": "Poly(A) Site / PABP Interface",
}

TRANSLATIONAL_CHEMISTRIES = {
    "pmo": "PMO (Phosphorodiamidate Morpholino) — Recommended for 5' UTR blocking",
    "moe_full_ps": "2'-O-MOE Full Phosphorothioate",
    "lna_dna_mixmer": "LNA / DNA Mixmer (Steric Blockade)",
}

# TG06 is a scored partition again, but the mapping stays on the unified
# defect vocabulary rather than reverting to a private mechanism list — the
# (goal, element) pair is a form convenience, not a second source of truth
# about which mechanism serves which defect. That still lives in each
# rule.json.
TRANSLATIONAL_ELEMENT_DEFECT: dict[tuple[str, str], str] = {
    ("suppress", "5p_utr"): "gain_of_function",
    ("enhance", "3p_utr_mirna"): "mirna_mediated_repression",
    ("enhance", "uorf"): "uorf_mediated_repression",
    ("suppress", "structured_element"): "structured_element_dysregulation",
    ("enhance", "structured_element"): "structured_element_dysregulation",
    # A29 / A30 / A31. Both directions map to the same defect class: the
    # element is either driving translation or it is not, and which way the
    # user wants to push it is a design choice, not a different defect.
    ("enhance", "ires_element"): "ires_mediated_translation",
    ("suppress", "ires_element"): "ires_mediated_translation",
    ("enhance", "kozak_consensus"): "kozak_context_dysregulation",
    ("suppress", "kozak_consensus"): "kozak_context_dysregulation",
    ("enhance", "polya_site"): "pabp_competition_defect",
    ("suppress", "polya_site"): "pabp_competition_defect",
}

# All seven TG06 mechanisms, for callers that need the roster.
TRANSLATIONAL_MECHANISM_IDS = ["A2", "A5", "A6", "A27", "A29", "A30", "A31"]

# TG07 is a scored goal again, and no longer a strict subset of TG04: A11 is
# now dual-tagged and A32/A33 are unique to it.
ISO_ENGINEERING_MECHANISM_IDS = ["A7", "A8", "A9", "A10", "A11", "A32", "A33"]

# The old `intron_retention` row mapped to `apa_dysregulation`, a defect
# served only by A11 — which was not a TG07 mechanism at the time, so that
# row always returned an empty ranking. Both halves of that are now fixed:
# A11 carries the TG07 tag, and intron retention has its own defect class
# and its own mechanism (A33) rather than being aliased onto APA.
ISOFORM_GOAL_DEFECT_MAP = {
    "exon_skipping": "exon_skipping_mutation",
    "exon_inclusion": "exon_inclusion_defect",
    "intron_retention": "intron_retention_defect",
    "alternative_splice_site": "cryptic_splice_site",
    "mutually_exclusive_exon": "exon_inclusion_defect",
    "apa_modulation": "apa_dysregulation",
    "alt_promoter_switch": "alt_promoter_dysregulation",
}

# TG09 form vocabularies. Retained so the existing page still renders; the
# endpoint now returns a rulebook lookup rather than generated candidates.
RNA_ENGINEERING_STRUCTURAL_CLASSES = {
    "rna_aptamer": "RNA Aptamer (Protein / Ligand Binding)",
    "catalytic_ribozyme": "Catalytic Ribozyme (mRNA Cleavage)",
    "riboswitch": "Riboswitch / Inducible RNA Sensor",
    "multivalent_scaffold": "Multivalent / Chimeric RNA Scaffold",
}

RNA_ENGINEERING_TARGET_TYPES = {
    "protein_active_site": "Protein Active Site / Surface Domain",
    "cell_surface_receptor": "Cell Surface Receptor (Internalizing)",
    "small_molecule": "Small Molecule / Metabolite",
    "target_rna": "Target RNA Transcript (Cleavage Site)",
}

RNA_ENGINEERING_SCAFFOLDS = {
    "selex_refinement": "SELEX Motif Structural Refinement",
    "hammerhead": "Hammerhead Architecture (Type I / Type III)",
    "three_way_junction": "3-Way Junction / Stable Stem-Loop Scaffold",
}

RNA_ENGINEERING_CHEM_STABILIZATIONS = {
    "two_f_pyrimidine": "2'-Fluoro-Pyrimidine (2'-F)",
    "two_ome_ps": "2'-O-Methyl (2'-OMe) / Phosphorothioate Stems",
    "inverted_abasic": "Inverted Abasic End-Cap (3'-3' Attachment)",
}

RNA_ENGINEERING_K_D_GOALS = {
    "nanomolar": "Nanomolar (1–10 nM)",
    "sub_nanomolar": "Sub-nanomolar (< 1 nM)",
}


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

def _flagged_for_goal(ctx: ArbitrationContext, goal: str) -> list[dict]:
    """Flag-only mechanisms tagged with this goal.

    Returned alongside the ranking, never inside it. Without this a
    goal-routed page shows nothing at all for TG08 and TG09, whose every
    mechanism is flag_only — being unscorable is not a reason to be
    unlistable.
    """
    ctx.goal_filter = None
    out = arbitrate(ctx)
    return [m for m in out.get("flaggedMechanisms", [])
            if goal in (m.get("goalTags") or [])]


def _filtered(ctx: ArbitrationContext, goal: str,
              restrict_to: list[str] | None = None) -> list[dict]:
    """Run the unified arbitration, then keep the goal's mechanisms.

    The filter is applied to a finished ranking. It changes what the caller
    SEES, never what was scored or how — which is the whole point of making
    the goal an output.
    """
    ctx.goal_filter = [goal]
    results = arbitrate(ctx)["results"]
    if restrict_to is not None:
        allowed = set(restrict_to)
        results = [r for r in results if r["id"] in allowed]
    return [_legacy_shape(r) for r in results]


def _legacy_shape(result: dict) -> dict:
    """Add the two fields the pre-arbitration response carried.

    `keywordMatch` was a soft bonus for the user's free-text variant
    description overlapping a mechanism's `suitableVariantTypes`. It is not
    evidence about the transcript, so it no longer influences the ranking;
    the field is kept at False so existing consumers do not break.
    """
    return {**result, "keywordMatch": False}


# Note: there is deliberately no "drop the ones we cannot design" filter here
# any more. A21 is scored, ranked and shown, carrying `designAvailable: false`
# and a reason string, because siRNA is a genuine competitor to RNase H
# knockdown for any knockdown target and the scientist is really choosing
# between them. Hiding it because this pipeline cannot emit a duplex removed a
# real option from a real decision.
#
# A24, A25 and A26 never reach these functions: they are FLAGGED, held out of
# the ranking by `arbitrate`, and surfaced through `modalityFlags`.


# ---------------------------------------------------------------------------
# TG01 — Gene Silencing
# ---------------------------------------------------------------------------

def rank_gene_silencing_mechanisms(
    defect_type: str,
    silencing_scope: str,
    delivery_context: str | None,
    known_variant: str | None,
    transcript_sequence: str | None = None,
    cds_start: int | None = None,
) -> list[dict]:
    """TG01 view of the unified ranking.

    `transcript_sequence` / `cds_start` are optional but they are what
    separates A1 from A2: A1 needs an accessible cleavable site anywhere in
    the transcript (F10a), A2 needs one at the 5' UTR or start codon (F10b).
    Without a sequence neither query can run and the two tie on evidence
    alone, exactly as they did before.
    """
    ctx = ArbitrationContext(
        molecular_defect=defect_type,
        allele_selective=(silencing_scope == "allele_specific"),
        delivery_context=delivery_context,
        known_variant=known_variant,
        transcript_sequence=transcript_sequence,
        cds_start=cds_start,
    )
    return (_filtered(ctx, "TG01"))


# ---------------------------------------------------------------------------
# TG02 — Gene Activation / Upregulation
# ---------------------------------------------------------------------------

def rank_gene_upregulation_mechanisms(
    defect_type: str,
    delivery_context: str | None,
    known_regulatory_element: str | None,
    gene_features: dict | None = None,
) -> list[dict]:
    ctx = ArbitrationContext(
        molecular_defect=defect_type,
        delivery_context=delivery_context,
        known_variant=known_regulatory_element,
        gene_features=gene_features,
    )
    return (_filtered(ctx, "TG02"))


# ---------------------------------------------------------------------------
# TG04 — RNA Processing Modulation
# ---------------------------------------------------------------------------

def rank_rna_processing_mechanisms(
    splice_defect_type: str,
    target_exon: str | None,
    delivery_context: str | None,
    known_variant: str | None,
) -> list[dict]:
    ctx = ArbitrationContext(
        molecular_defect=splice_defect_type,
        delivery_context=delivery_context,
        known_variant=known_variant,
        extras={"targetExon": target_exon},
    )
    return (_filtered(ctx, "TG04"))


# ---------------------------------------------------------------------------
# TG03 — RNA Editing / Correction
# ---------------------------------------------------------------------------

# The plan defers TG03: mechanism choice there is near-bijective on the
# variant (A→I for a G>A, C→U for a T>C, trans-splicing for a larger lesion),
# and the hard part is guide design, not mechanism selection. The variant →
# edit-type lookup below is cheap and honest; no arbitration is claimed.
EDIT_TYPE_DEFECT = {
    "a_to_i": "correctable_point_variant",
    "c_to_u": "correctable_point_variant",
    "trans_splicing": "coding_region_lesion",
}


def rank_rna_editing_mechanisms(
    edit_type: str,
    variant_hgvs: str | None,
    enzyme_recruitment: str | None,
    delivery_context: str | None,
    guide_length: int | None,
    mismatch_pocket: str | None,
    max_bystander_edits: int | None,
    exon_count: int | None = None,
    intron_count: int | None = None,
    total_transcripts: int | None = None,
) -> list[dict]:
    ctx = ArbitrationContext(
        molecular_defect=EDIT_TYPE_DEFECT.get(edit_type),
        edit_type=edit_type,
        variant_hgvs=variant_hgvs,
        delivery_context=delivery_context,
        exon_count=exon_count,
        intron_count=intron_count,
        total_transcripts=total_transcripts,
        extras={
            "enzymeRecruitment": enzyme_recruitment,
            "guideLength": guide_length,
            "mismatchPocket": mismatch_pocket,
            "maxBystanderEdits": max_bystander_edits,
        },
    )
    return (_filtered(ctx, "TG03"))


# ---------------------------------------------------------------------------
# TG05 — RNA Neutralization
# ---------------------------------------------------------------------------

def rank_rna_neutralization_mechanisms(
    molecular_defect: str,
    neutralization_mode: str,
    repeat_unit: str | None = None,
    estimated_repeat_count: str | None = None,
    steric_chemistry: str | None = None,
    target_rbp: str | None = None,
    oligo_length: int | None = None,
    delivery_context: str | None = None,
    target_gene_type: str | None = None,
) -> list[dict]:
    ctx = ArbitrationContext(
        molecular_defect=molecular_defect,
        transcript_class=target_gene_type,
        delivery_context=delivery_context,
        repeat_unit=repeat_unit,
        repeat_count=estimated_repeat_count,
        oligo_length=oligo_length or 18,
        extras={"stericChemistry": steric_chemistry, "targetRbp": target_rbp},
    )
    return _filtered(
        ctx, "TG05", restrict_to=NEUTRALIZATION_MODE_MECHANISMS.get(
            neutralization_mode, []
        )
    )


# ---------------------------------------------------------------------------
# TG06 — Translational Regulation (retired as a scoring partition)
# ---------------------------------------------------------------------------

def rank_translational_regulation_mechanisms(
    translational_goal: str | None,
    target_element: str | None,
    delivery_context: str | None = None,
    oligo_length: int | None = None,
) -> list[dict]:
    """TG06 view of the unified ranking.

    TG06 is a scored goal again, now with seven mechanisms: A2, A5 and A6
    (shared with gene silencing and activation), A27, and the three added
    with the restoration — A29 IRES, A30 Kozak, A31 PABP competition.

    This is still a filter over the one shared pass, not a separate scorer.
    The (goal, element) pair the page collects is translated into the unified
    defect vocabulary and the shared ranking does the rest, so a mechanism
    cannot score differently here than it does anywhere else.
    """
    defect = TRANSLATIONAL_ELEMENT_DEFECT.get(
        (translational_goal or "", target_element or "")
    )
    ctx = ArbitrationContext(
        molecular_defect=defect,
        delivery_context=delivery_context,
        oligo_length=oligo_length or 18,
    )
    return (_filtered(ctx, "TG06"))


# ---------------------------------------------------------------------------
# TG07 — Isoform Engineering (retired as a scoring partition)
# ---------------------------------------------------------------------------

def rank_isoform_engineering_mechanisms(
    isoform_goal: str,
    target_exon_locus: str | None = None,
    splice_element_target: str | None = None,
    delivery_context: str | None = None,
) -> list[dict]:
    """TG07 view of the unified ranking.

    TG07 is a scored goal again and no longer a strict subset of TG04. It now
    covers seven mechanisms: A7-A10 (shared with RNA processing), A11 (APA,
    dual-tagged), and A32 / A33, which are unique to it. The framing that
    separates them is intent — TG04 fixes broken splicing, TG07 chooses
    between isoforms that are all functional.

    Still a filter over the one shared pass, so a mechanism shared with TG04
    cannot score differently depending on which page asked.
    """
    ctx = ArbitrationContext(
        molecular_defect=ISOFORM_GOAL_DEFECT_MAP.get(isoform_goal, isoform_goal),
        delivery_context=delivery_context,
        extras={
            "targetExonLocus": target_exon_locus,
            "spliceElementTarget": splice_element_target,
        },
    )
    return (_filtered(ctx, "TG07"))


# ---------------------------------------------------------------------------
# TG08 and TG09 — flagged modalities
#
# A24 (mRNA), A25 (aptamer) and A26 (circRNA) are surfaced qualitatively and
# never scored. The reason is validatability, not modality: there is no
# approved mRNA protein replacement therapy and no approved circRNA therapy,
# so a stated applicability for either could not be checked against any case
# in the world. Contrast A21, which IS scored despite being undesignable here,
# because five approved siRNA drugs make it fully validatable.
#
# Whether a flag actually fires for a given target is decided by
# mechanism_arbitration.modality_flags() from features P2, P6 and B1. These
# helpers return the rulebook content behind each flag.
# ---------------------------------------------------------------------------

def _flagged_mechanism(mechanism_id: str) -> dict:
    rule = load_rule(mechanism_id) or {}
    arb = rule.get("arbitration", {})
    return {
        "mechanismId": mechanism_id,
        "name": rule.get("name"),
        "category": rule.get("category"),
        "status": "FLAGGED",
        "scored": False,
        "flagReason": arb.get("flagReason"),
        "designUnavailableReason": arb.get("designUnavailableReason"),
        "evidenceLevel": rule.get("evidenceLevel"),
        "fdaApprovedDrugs": rule.get("fdaApprovedDrugs"),
        "clinicalTrialExamples": rule.get("clinicalTrialExamples"),
        "molecularDefect": rule.get("molecularDefect"),
        "designRules": rule.get("designRules"),
        "advantages": rule.get("advantages"),
        "limitations": rule.get("limitations"),
        "offTargetConsiderations": rule.get("offTargetConsiderations"),
        "references": rule.get("references", [])[:3],
    }


# TG09 now carries four flagged mechanisms; TG08 carries five.
PROTEIN_FUNCTION_MECHANISM_IDS = ["A25", "A37", "A38", "A39"]
PROTEIN_REPLACEMENT_MECHANISM_IDS = ["A24", "A26", "A34", "A35", "A36"]


def lookup_protein_function_modulation() -> dict:
    """A25 and the three aptamer variants, as rulebook content.

    No score, no rank. `mechanisms` carries all four; the top-level fields
    stay on A25 so existing callers keep working.
    """
    primary = _flagged_mechanism("A25")
    primary["mechanisms"] = [
        _flagged_mechanism(m) for m in PROTEIN_FUNCTION_MECHANISM_IDS
    ]
    return primary


def protein_replacement_scope_notice() -> dict:
    """A24 and A26, as rulebook content. No score, no rank."""
    return {
        "status": "FLAGGED",
        "scored": False,
        "mechanisms": [_flagged_mechanism(m)
                       for m in PROTEIN_REPLACEMENT_MECHANISM_IDS],
        "goalNotice": RETIRED_AS_SCORING_PARTITION["TG08"],
    }
