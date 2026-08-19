"""
Mechanism selection endpoints. Kept separate from main.py's gene retrieval
pipeline since this is a distinct concern (Rulebook Engine, not the
Biological Information Retrieval Engine).

`POST /api/mechanisms/arbitrate` is the primary endpoint: it scores every
designable mechanism in one pass and reports the therapeutic goal as an
OUTPUT. The per-goal endpoints below remain for the existing goal-routed
pages and are now filters over that same ranking, not separate scorers.

Goal status after the scope review (see
docs/planning/therapeutic_goal_scope_plan_v3.md):

- TG01 Gene Silencing              scored   A1, A2, A12, A15; A21 scored too
- TG02 Gene Activation             scored   A3, A4, A5, A6, A23; A28 halts
- TG03 RNA Editing                 deferred A13, A16-A20; no arbitration claimed
- TG04 RNA Processing              scored   A7, A8, A9, A10; A11 halts
- TG05 RNA Neutralization          narrow   A14 (halts pending F12)
- TG06 Translational Regulation    scored   A2, A5, A6, A27, A29, A30, A31
- TG07 Isoform Engineering         scored   A7-A11, A32, A33 (A11 dual-tagged)
- TG08 Protein Replacement         flag only; A24, A26, A34-A36 never scored
- TG09 Protein Function Modulation flag only; A25, A37-A39 never scored

Four mechanism states, and only the last is absence — nothing is in it:

  SCORED + DESIGNABLE            competes; this platform emits candidates
  SCORED + DESIGN UNAVAILABLE    competes; another pipeline required (A21)
  HALTED                         in the choice set, required feature absent
  FLAGGED                        surfaced qualitatively, never scored
  REMOVED                        nothing

Flagged mechanisms are returned in `flaggedMechanisms`, ALWAYS — not only
when a modality flag fires. A flag needs tissue expression or subcellular
localisation, inputs the pages do not collect, so gating visibility on one
made all nine invisible in the case where they matter most. They carry a
null applicability, confidence and score so nothing can rank them against a
scored mechanism.

A21 is scored despite being undesignable here because siRNA is a genuine
alternative to RNase H knockdown and five approved drugs make it fully
validatable. A24/A25/A26 are flagged instead because no approved therapy
exists in their indication class, so any number attached to them could not be
checked against anything. The distinction is validatability, not modality.

No mechanism was deleted and every rulebook is retained.
"""

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from typing import Optional

from services.mechanism_arbitration import (
    ArbitrationContext,
    load_rule,
    MOLECULAR_DEFECTS,
    RETIRED_AS_SCORING_PARTITION,
    arbitrate,
    therapeutic_goals,
)
from services.mechanism_service import (
    _flagged_for_goal,
    rank_gene_silencing_mechanisms,
    rank_gene_upregulation_mechanisms,
    rank_rna_processing_mechanisms,
    rank_rna_editing_mechanisms,
    rank_rna_neutralization_mechanisms,
    rank_translational_regulation_mechanisms,
    rank_isoform_engineering_mechanisms,
    lookup_protein_function_modulation,
    protein_replacement_scope_notice,
    DEFECT_TYPES,
    SILENCING_SCOPES,
    GENE_UPREGULATION_DEFECT_TYPES,
    DELIVERY_CONTEXTS,
    SPLICE_DEFECT_TYPES,
    EDIT_TYPES,
    ENZYME_RECRUITMENT,
    MISMATCH_POCKET,
    SPLICING_DIRECTIONS,
    INTRON_SITES,
    NEUTRALIZATION_DEFECT_TYPES,
    NEUTRALIZATION_MODES,
    STERIC_CHEMISTRIES,
    TRANSLATIONAL_GOALS,
    TRANSLATIONAL_TARGET_ELEMENTS,
    TRANSLATIONAL_CHEMISTRIES,
    RNA_ENGINEERING_STRUCTURAL_CLASSES,
    RNA_ENGINEERING_TARGET_TYPES,
    RNA_ENGINEERING_SCAFFOLDS,
    RNA_ENGINEERING_CHEM_STABILIZATIONS,
    RNA_ENGINEERING_K_D_GOALS,
    ISO_ENGINEERING_MECHANISM_IDS,
    ISOFORM_GOAL_DEFECT_MAP,
)
from services.gene_feature_service import analyze_gene_features
from services.reference_tables import status as reference_table_status
from services.spliceai_service import status as spliceai_status

router = APIRouter()


class GeneSilencingRequest(BaseModel):
    gene_symbol: str
    defect_type: str
    silencing_scope: str
    delivery_context: Optional[str] = None
    known_variant: Optional[str] = None


class GeneUpregulationRequest(BaseModel):
    gene_symbol: str
    defect_type: str
    delivery_context: Optional[str] = None
    known_regulatory_element: Optional[str] = None
    gene_features: Optional[dict] = None


class RnaProcessingRequest(BaseModel):
    gene_symbol: str
    splice_defect_type: str
    target_exon: Optional[str] = None
    delivery_context: Optional[str] = None
    known_variant: Optional[str] = None


class RnaEditingRequest(BaseModel):
    gene_symbol: str
    edit_type: str
    variant_hgvs: Optional[str] = None
    enzyme_recruitment: Optional[str] = None
    delivery_context: Optional[str] = None
    guide_length: Optional[int] = None
    mismatch_pocket: Optional[str] = None
    max_bystander_edits: Optional[int] = None
    splicing_direction: Optional[str] = None
    intron_site: Optional[str] = None
    abd_length: Optional[int] = None
    exon_count: Optional[int] = None
    intron_count: Optional[int] = None
    total_transcripts: Optional[int] = None


class RnaNeutralizationRequest(BaseModel):
    gene_symbol: str
    molecular_defect: str
    neutralization_mode: str
    repeat_unit: Optional[str] = None
    estimated_repeat_count: Optional[str] = None
    steric_chemistry: Optional[str] = None
    target_rbp: Optional[str] = None
    oligo_length: Optional[int] = None
    delivery_context: Optional[str] = None
    target_gene_type: Optional[str] = None


class TranslationalRegulationRequest(BaseModel):
    gene_symbol: str
    translational_goal: Optional[str] = None
    target_element: Optional[str] = None
    steric_chemistry: Optional[str] = None
    target_rbp: Optional[str] = None
    oligo_length: Optional[int] = None
    delivery_context: Optional[str] = None


class RnaEngineeringRequest(BaseModel):
    gene_symbol: str
    structural_class: str
    target_type: str
    scaffold: str
    chem_stabilization: str
    kd_goal: str
    delivery_context: Optional[str] = None


class IsoformEngineeringRequest(BaseModel):
    gene_symbol: str
    isoform_goal: str
    target_exon_locus: Optional[str] = None
    splice_element_target: Optional[str] = None
    steric_chemistry: Optional[str] = None
    delivery_context: Optional[str] = None


@router.get("/api/mechanisms/options")
async def mechanism_options():
    """Input options for all mechanism selection forms."""
    return {
        "geneSilencing": {
            "defectTypes": [{"id": k, "label": v} for k, v in DEFECT_TYPES.items()],
            "silencingScopes": [{"id": k, "label": v} for k, v in SILENCING_SCOPES.items()],
        },
        "geneUpregulation": {
            "defectTypes": [{"id": k, "label": v} for k, v in GENE_UPREGULATION_DEFECT_TYPES.items()],
        },
        "rnaProcessing": {
            "spliceDefectTypes": [{"id": k, "label": v} for k, v in SPLICE_DEFECT_TYPES.items()],
        },
        "rnaEditing": {
            "editTypes": [{"id": k, "label": v} for k, v in EDIT_TYPES.items()],
            "enzymeRecruitment": [{"id": k, "label": v} for k, v in ENZYME_RECRUITMENT.items()],
            "mismatchPocket": [{"id": k, "label": v} for k, v in MISMATCH_POCKET.items()],
            "splicingDirections": [{"id": k, "label": v} for k, v in SPLICING_DIRECTIONS.items()],
            "intronSites": [{"id": k, "label": v} for k, v in INTRON_SITES.items()],
        },
        "rnaNeutralization": {
            "molecularDefects": [{"id": k, "label": v} for k, v in NEUTRALIZATION_DEFECT_TYPES.items()],
            "neutralizationModes": [{"id": k, "label": v} for k, v in NEUTRALIZATION_MODES.items()],
            "stericChemistries": [{"id": k, "label": v} for k, v in STERIC_CHEMISTRIES.items()],
        },
        "translationalRegulation": {
            "translationalGoals": [{"id": k, "label": v} for k, v in TRANSLATIONAL_GOALS.items()],
            "targetElements": [{"id": k, "label": v} for k, v in TRANSLATIONAL_TARGET_ELEMENTS.items()],
            "stericChemistries": [{"id": k, "label": v} for k, v in TRANSLATIONAL_CHEMISTRIES.items()],
        },
        "rnaEngineering": {
            "structuralClasses": [{"id": k, "label": v} for k, v in RNA_ENGINEERING_STRUCTURAL_CLASSES.items()],
            "targetTypes": [{"id": k, "label": v} for k, v in RNA_ENGINEERING_TARGET_TYPES.items()],
            "scaffolds": [{"id": k, "label": v} for k, v in RNA_ENGINEERING_SCAFFOLDS.items()],
            "chemStabilizations": [{"id": k, "label": v} for k, v in RNA_ENGINEERING_CHEM_STABILIZATIONS.items()],
            "kdGoals": [{"id": k, "label": v} for k, v in RNA_ENGINEERING_K_D_GOALS.items()],
        },
        "deliveryContexts": [{"id": k, "label": v} for k, v in DELIVERY_CONTEXTS.items()],
    }


@router.post("/api/mechanisms/gene-silencing")
async def gene_silencing_mechanisms(payload: GeneSilencingRequest):
    if payload.defect_type not in DEFECT_TYPES:
        raise HTTPException(status_code=400, detail=f"Unknown defect_type: {payload.defect_type}")
    if payload.silencing_scope not in SILENCING_SCOPES:
        raise HTTPException(status_code=400, detail=f"Unknown silencing_scope: {payload.silencing_scope}")

    results = rank_gene_silencing_mechanisms(
        defect_type=payload.defect_type,
        silencing_scope=payload.silencing_scope,
        delivery_context=payload.delivery_context,
        known_variant=payload.known_variant,
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "Gene Silencing",
        "inputs": {
            "defectType": payload.defect_type,
            "silencingScope": payload.silencing_scope,
            "deliveryContext": payload.delivery_context,
            "knownVariant": payload.known_variant,
        },
        "results": results,
    }


def _resolve_gene_features(gene_symbol: str, supplied: Optional[dict],
                           organism: str = "homo_sapiens") -> Optional[dict]:
    """Gene features for the arbitration context, fetched if not supplied.

    The endpoint used to pass `payload.gene_features` straight through, so a
    caller that did not send them left F4 (poison exon) and F6 (natural
    antisense transcript) UNRESOLVED. Both are REQUIRED features, so A3 and A4
    halted for every real target and scored only when the user hand-asserted
    the matching molecular defect — the platform never checked whether the
    gene actually HAS a poison exon or a NAT.

    Note the shape: `FeatureContext.gene_feature` reads
    `gene_features["features"][key]`, so the WHOLE payload from
    `analyze_gene_features` must be passed, not its `features` sub-dict.
    Passing the inner dict silently resolves nothing.
    """
    if supplied:
        # A caller may send either shape; normalise the inner one.
        return supplied if "features" in supplied else {"features": supplied}
    symbol = (gene_symbol or "").strip().upper()
    if not symbol:
        return None
    try:
        return analyze_gene_features(symbol, organism=organism)
    except Exception as exc:  # never fail the ranking over an enrichment step
        logger.warning("Gene-feature resolution failed for %s: %s", symbol, exc)
        return None


@router.post("/api/mechanisms/gene-upregulation")
async def gene_upregulation_mechanisms(payload: GeneUpregulationRequest):
    if payload.defect_type not in GENE_UPREGULATION_DEFECT_TYPES:
        raise HTTPException(status_code=400, detail=f"Unknown defect_type: {payload.defect_type}")

    results = rank_gene_upregulation_mechanisms(
        defect_type=payload.defect_type,
        delivery_context=payload.delivery_context,
        known_regulatory_element=payload.known_regulatory_element,
        gene_features=_resolve_gene_features(
            payload.gene_symbol, payload.gene_features),
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "Gene Activation / Upregulation",
        "inputs": {
            "defectType": payload.defect_type,
            "deliveryContext": payload.delivery_context,
            "knownRegulatoryElement": payload.known_regulatory_element,
        },
        "results": results,
    }


@router.post("/api/mechanisms/rna-processing")
async def rna_processing_mechanisms(payload: RnaProcessingRequest):
    if payload.splice_defect_type not in SPLICE_DEFECT_TYPES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown splice_defect_type: {payload.splice_defect_type}",
        )

    results = rank_rna_processing_mechanisms(
        splice_defect_type=payload.splice_defect_type,
        target_exon=payload.target_exon,
        delivery_context=payload.delivery_context,
        known_variant=payload.known_variant,
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "RNA Processing Modulation",
        "inputs": {
            "spliceDefectType": payload.splice_defect_type,
            "targetExon": payload.target_exon,
            "deliveryContext": payload.delivery_context,
            "knownVariant": payload.known_variant,
        },
        "results": results,
    }


@router.post("/api/mechanisms/rna-editing")
async def rna_editing_mechanisms(payload: RnaEditingRequest):
    if payload.edit_type not in EDIT_TYPES:
        raise HTTPException(status_code=400, detail=f"Unknown edit_type: {payload.edit_type}")
    if payload.enzyme_recruitment and payload.enzyme_recruitment not in ENZYME_RECRUITMENT:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown enzyme_recruitment: {payload.enzyme_recruitment}",
        )
    if payload.mismatch_pocket and payload.mismatch_pocket not in MISMATCH_POCKET:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown mismatch_pocket: {payload.mismatch_pocket}",
        )
    if payload.splicing_direction and payload.splicing_direction not in SPLICING_DIRECTIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown splicing_direction: {payload.splicing_direction}",
        )

    results = rank_rna_editing_mechanisms(
        edit_type=payload.edit_type,
        variant_hgvs=payload.variant_hgvs,
        enzyme_recruitment=payload.enzyme_recruitment,
        delivery_context=payload.delivery_context,
        guide_length=payload.guide_length,
        mismatch_pocket=payload.mismatch_pocket,
        max_bystander_edits=payload.max_bystander_edits,
        exon_count=payload.exon_count,
        intron_count=payload.intron_count,
        total_transcripts=payload.total_transcripts,
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "RNA Editing / Correction",
        "inputs": {
            "editType": payload.edit_type,
            "variantHgvs": payload.variant_hgvs,
            "enzymeRecruitment": payload.enzyme_recruitment,
            "deliveryContext": payload.delivery_context,
            "guideLength": payload.guide_length,
            "mismatchPocket": payload.mismatch_pocket,
            "maxBystanderEdits": payload.max_bystander_edits,
            "splicingDirection": payload.splicing_direction,
            "intronSite": payload.intron_site,
            "abdLength": payload.abd_length,
        },
        "results": results,
    }


@router.post("/api/mechanisms/rna-neutralization")
async def rna_neutralization_mechanisms(payload: RnaNeutralizationRequest):
    if payload.molecular_defect not in NEUTRALIZATION_DEFECT_TYPES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown molecular_defect: {payload.molecular_defect}",
        )
    if payload.neutralization_mode not in NEUTRALIZATION_MODES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown neutralization_mode: {payload.neutralization_mode}",
        )
    if payload.steric_chemistry and payload.steric_chemistry not in STERIC_CHEMISTRIES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown steric_chemistry: {payload.steric_chemistry}",
        )

    results = rank_rna_neutralization_mechanisms(
        molecular_defect=payload.molecular_defect,
        neutralization_mode=payload.neutralization_mode,
        repeat_unit=payload.repeat_unit,
        estimated_repeat_count=payload.estimated_repeat_count,
        steric_chemistry=payload.steric_chemistry,
        target_rbp=payload.target_rbp,
        oligo_length=payload.oligo_length,
        delivery_context=payload.delivery_context,
        target_gene_type=payload.target_gene_type,
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "RNA Neutralization",
        "inputs": {
            "molecularDefect": payload.molecular_defect,
            "neutralizationMode": payload.neutralization_mode,
            "repeatUnit": payload.repeat_unit,
            "estimatedRepeatCount": payload.estimated_repeat_count,
            "stericChemistry": payload.steric_chemistry,
            "targetRbp": payload.target_rbp,
            "oligoLength": payload.oligo_length,
            "deliveryContext": payload.delivery_context,
            "targetGeneType": payload.target_gene_type,
        },
         "results": results,
     }


@router.post("/api/mechanisms/translational-regulation")
async def translational_regulation_mechanisms(payload: TranslationalRegulationRequest):
     if payload.translational_goal and payload.translational_goal not in TRANSLATIONAL_GOALS:
         raise HTTPException(
             status_code=400,
             detail=f"Unknown translational_goal: {payload.translational_goal}",
         )
     if payload.target_element and payload.target_element not in TRANSLATIONAL_TARGET_ELEMENTS:
         raise HTTPException(
             status_code=400,
             detail=f"Unknown target_element: {payload.target_element}",
         )
     if payload.steric_chemistry and payload.steric_chemistry not in TRANSLATIONAL_CHEMISTRIES:
         raise HTTPException(
             status_code=400,
             detail=f"Unknown steric_chemistry: {payload.steric_chemistry}",
         )

     # TG06 is a scored goal again (7 mechanisms). Still a filter over the
     # one shared pass, so it cannot disagree with any other route.
     results = rank_translational_regulation_mechanisms(
         translational_goal=payload.translational_goal,
         target_element=payload.target_element,
         delivery_context=payload.delivery_context,
         oligo_length=payload.oligo_length,
     )

     return {
         "geneSymbol": payload.gene_symbol.strip().upper(),
         "therapeuticGoal": "Translational Regulation",
         "inputs": {
             "translationalGoal": payload.translational_goal,
             "targetElement": payload.target_element,
             "stericChemistry": payload.steric_chemistry,
             "targetRbp": payload.target_rbp,
             "oligoLength": payload.oligo_length,
             "deliveryContext": payload.delivery_context,
         },
         "results": results,
      }


def _aptamer_candidate_from_rule(mechanism_id: str, rank: int) -> dict:
    """One flagged aptamer mechanism, shaped like a candidate row.

    Every quantitative field is None. These are NOT predictions and there is
    nothing to predict from: no aptamer has been selected, so there is no
    sequence to fold, no Tm, no dG and no measured Kd. The earlier version of
    this endpoint filled those fields from hash() of the form inputs, which
    rendered as measurements.

    The guidance strings say what a real SELEX campaign would determine,
    which is genuinely useful and is not a number.
    """
    rule = load_rule(mechanism_id) or {}
    arb = rule.get("arbitration", {})
    return {
        "rank": rank,
        "constructId": f"APT-{mechanism_id}-001",
        "mechanismId": mechanism_id,
        "mechanismName": rule.get("name"),
        "structuralMotif": "SELEX-derived",
        "length": "~20-100 nt",
        "tm": None,
        "deltaGFolding": None,
        "kdPrediction": "nanomolar range (SELEX-dependent)",
        "targetSpecificityScore": None,
        "serumStability": "requires chemical stabilization",
        "structuralRigidityFlag": "dependent on selected aptamer",
        "sequence": None,
        "dotBracket": None,
        "rationale": rule.get("designRules", ""),
        "foldingScore": None,
        "tHalfScore": None,
        "scored": False,
        "flagReason": arb.get("flagReason"),
        "designUnavailableReason": arb.get("designUnavailableReason"),
    }


@router.post("/api/mechanisms/rna-engineering")
async def rna_engineering_mechanisms(payload: RnaEngineeringRequest):
    if payload.structural_class not in RNA_ENGINEERING_STRUCTURAL_CLASSES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown structural_class: {payload.structural_class}",
        )
    if payload.target_type not in RNA_ENGINEERING_TARGET_TYPES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown target_type: {payload.target_type}",
        )
    if payload.scaffold not in RNA_ENGINEERING_SCAFFOLDS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown scaffold: {payload.scaffold}",
        )
    if payload.chem_stabilization not in RNA_ENGINEERING_CHEM_STABILIZATIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown chem_stabilization: {payload.chem_stabilization}",
        )
    if payload.kd_goal not in RNA_ENGINEERING_K_D_GOALS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown kd_goal: {payload.kd_goal}",
        )

    # TG09 is a flag, not a ranking: it contains one mechanism (A25), and a
    # ranking over a single item is not a ranking. The generated aptamer
    # candidates this endpoint used to return — sequence, Tm, ΔG, predicted
    # Kd, serum half-life — were all derived from hash() of the form inputs.
    # They were not predictions and they are not returned any more.
    #
    # Whether the aptamer flag actually fires for a given target is decided by
    # the unified pass from feature B1; call /api/mechanisms/arbitrate for
    # that. This endpoint returns the rulebook content behind the flag.
    lookup = lookup_protein_function_modulation()
    mechanisms = lookup.get("mechanisms", [])
    candidates = [
        _aptamer_candidate_from_rule(m["mechanismId"], i + 1)
        for i, m in enumerate(mechanisms)
    ]

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "Protein Function Modulation",
        "status": "FLAGGED",
        "goalNotice": RETIRED_AS_SCORING_PARTITION["TG09"],
        "flaggedMechanisms": lookup.get("mechanisms", []),
        "inputs": {
            "structuralClass": payload.structural_class,
            "targetType": payload.target_type,
            "scaffold": payload.scaffold,
            "chemStabilization": payload.chem_stabilization,
            "kdGoal": payload.kd_goal,
            "deliveryContext": payload.delivery_context,
        },
        "mechanism": lookup,
        "mechanisms": mechanisms,
        "candidates": candidates,
    }


@router.post("/api/mechanisms/isoform-engineering")
async def isoform_engineering_mechanisms(payload: IsoformEngineeringRequest):
    if payload.isoform_goal not in ISOFORM_GOAL_DEFECT_MAP:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown isoform_goal: {payload.isoform_goal}",
        )

    # TG07 is a scored goal again (7 mechanisms, no longer a TG04 subset).
    # Still a filter over the one shared pass.
    results = rank_isoform_engineering_mechanisms(
        isoform_goal=payload.isoform_goal,
        target_exon_locus=payload.target_exon_locus,
        splice_element_target=payload.splice_element_target,
        delivery_context=payload.delivery_context,
    )

    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "therapeuticGoal": "Isoform Engineering",
        "inputs": {
            "isoformGoal": payload.isoform_goal,
            "targetExonLocus": payload.target_exon_locus,
            "spliceElementTarget": payload.splice_element_target,
            "stericChemistry": payload.steric_chemistry,
            "deliveryContext": payload.delivery_context,
        },
        "results": results,
    }


class ArbitrationRequest(BaseModel):
    """Inputs to the unified pass.

    Note what is NOT here: a therapeutic goal. `goal_filter` narrows the
    finished ranking for a user who already knows what they want; it never
    decides what gets scored.
    """

    gene_symbol: str
    molecular_defect: Optional[str] = None
    allele_selective: Optional[bool] = None
    transcript_class: Optional[str] = None
    edit_type: Optional[str] = None
    variant_hgvs: Optional[str] = None
    known_variant: Optional[str] = None
    delivery_context: Optional[str] = None
    exon_count: Optional[int] = None
    intron_count: Optional[int] = None
    total_transcripts: Optional[int] = None
    gene_features: Optional[dict] = None
    transcript_sequence: Optional[str] = None
    cds_start: Optional[int] = None
    repeat_unit: Optional[str] = None
    repeat_count: Optional[str] = None
    oligo_length: int = 18
    # Modality-flag inputs. Never scored; they decide only whether a
    # qualitative signpost toward an undesigned modality is shown.
    tissue_tpm: Optional[float] = None
    protein_localisation: Optional[str] = None
    goal_filter: Optional[list[str]] = None


@router.post("/api/mechanisms/arbitrate")
async def arbitrate_mechanisms(payload: ArbitrationRequest):
    """Score every designable mechanism in one pass.

    The therapeutic goal comes back as an output label on the winning
    mechanism. This is the endpoint that makes the nusinersen class of
    failure impossible: a splice-modulating answer surfaces for an
    upregulation intent because nothing was filtered before scoring.
    """
    if payload.molecular_defect and payload.molecular_defect not in MOLECULAR_DEFECTS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown molecular_defect: {payload.molecular_defect}",
        )
    if payload.edit_type and payload.edit_type not in EDIT_TYPES:
        raise HTTPException(
            status_code=400, detail=f"Unknown edit_type: {payload.edit_type}"
        )
    goals = therapeutic_goals()
    for goal in payload.goal_filter or []:
        if goal not in goals:
            raise HTTPException(
                status_code=400, detail=f"Unknown therapeutic goal: {goal}"
            )

    ctx = ArbitrationContext(
        gene_symbol=payload.gene_symbol,
        molecular_defect=payload.molecular_defect,
        allele_selective=payload.allele_selective,
        transcript_class=payload.transcript_class,
        edit_type=payload.edit_type,
        variant_hgvs=payload.variant_hgvs,
        known_variant=payload.known_variant,
        delivery_context=payload.delivery_context,
        exon_count=payload.exon_count,
        intron_count=payload.intron_count,
        total_transcripts=payload.total_transcripts,
        gene_features=payload.gene_features,
        transcript_sequence=payload.transcript_sequence,
        cds_start=payload.cds_start,
        repeat_unit=payload.repeat_unit,
        repeat_count=payload.repeat_count,
        oligo_length=payload.oligo_length,
        tissue_tpm=payload.tissue_tpm,
        protein_localisation=payload.protein_localisation,
        goal_filter=payload.goal_filter,
    )
    return arbitrate(ctx)


@router.get("/api/mechanisms/scope")
async def mechanism_scope():
    """What this designer scores, defers, and refuses — and why."""
    return {
        "goals": therapeutic_goals(),
        "retiredAsScoringPartition": RETIRED_AS_SCORING_PARTITION,
        "proteinReplacement": protein_replacement_scope_notice(),
        "proteinFunctionModulation": lookup_protein_function_modulation(),
        "molecularDefects": [
            {"id": k, "label": v} for k, v in MOLECULAR_DEFECTS.items()
        ],
        # Which reference tables are populated. An unpopulated table is why a
        # mechanism halts or a modality flag is withheld, so the answer to
        # "why did I get nothing?" is visible rather than buried.
        "referenceTables": reference_table_status(),
        # F1/F2/F3 resolve from SpliceAI when a pre-mRNA sequence is supplied
        # and fall back to the user-asserted stand-in otherwise. This says
        # which of those is in play.
        "spliceai": spliceai_status(),
    }


@router.get("/api/mechanisms/gene-features")
async def gene_features(
    gene_symbol: str = Query(..., description="Gene symbol to analyze"),
    organism: str = Query("homo_sapiens", description="Species slug"),
    ensembl_id: Optional[str] = Query(None, description="Ensembl gene ID (optional)"),
    tissue_tpm: Optional[float] = Query(None, description="Baseline tissue TPM for overexpression warning"),
    exon_count: Optional[int] = Query(None, description="Known exon count from the gene pipeline (optional)"),
    total_transcripts: Optional[int] = Query(None, description="Known transcript count from the gene pipeline (optional)"),
    gene_type: Optional[str] = Query(None, description="Gene biotype, e.g. protein_coding (optional)"),
):
    """Analyze gene structural features to determine TG02 mechanism availability."""
    result = analyze_gene_features(
        gene_symbol=gene_symbol.strip(),
        organism=organism.strip(),
        ensembl_id=ensembl_id.strip() if ensembl_id else None,
        tissue_tpm=tissue_tpm,
        exon_count=exon_count,
        total_transcripts=total_transcripts,
        gene_type=gene_type,
    )
    return result

# ---------------------------------------------------------------------------
# Designers for mechanisms the goal-specific services cannot build
# ---------------------------------------------------------------------------

class SirnaDuplexRequest(BaseModel):
    ensembl_gene_id: str
    gene_symbol: str = ""
    organism: str = "homo_sapiens"
    max_candidates: int = 12


@router.post("/api/mechanisms/A21/sirna-duplex")
async def design_a21_duplex(payload: SirnaDuplexRequest):
    """A21 — siRNA duplex design (guide + passenger with 3' overhangs)."""
    from services.sirna_duplex_service import design_sirna_duplexes
    try:
        return design_sirna_duplexes(
            ensembl_gene_id=payload.ensembl_gene_id,
            gene_symbol=payload.gene_symbol,
            organism=payload.organism,
            max_candidates=payload.max_candidates,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


class EditorGuideRequest(BaseModel):
    mechanism_id: str
    ensembl_gene_id: str
    edit_position: int
    gene_symbol: str = ""
    organism: str = "homo_sapiens"
    max_candidates: int = 8


@router.post("/api/mechanisms/editor-guide")
async def design_editor_guide(payload: EditorGuideRequest):
    """A18 / A19 — spacers for protein-dependent RNA editors."""
    from services.programmable_editor_service import design_editor_guides
    try:
        return design_editor_guides(
            mechanism_id=payload.mechanism_id,
            ensembl_gene_id=payload.ensembl_gene_id,
            edit_position=payload.edit_position,
            gene_symbol=payload.gene_symbol,
            organism=payload.organism,
            max_candidates=payload.max_candidates,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


class IntronRetentionRequest(BaseModel):
    ensembl_gene_id: str
    intron_number: int
    gene_symbol: str = ""
    organism: str = "homo_sapiens"
    oligo_length: int = 20
    splice_element: str = "both"
    max_candidates: int = 12


@router.post("/api/mechanisms/A33/intron-retention")
async def design_a33_intron(payload: IntronRetentionRequest):
    """A33 — steric blockers across an intron's splice sites."""
    from services.transcript_architecture_service import design_intron_retention
    try:
        return design_intron_retention(**payload.dict())
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


class AlternativePromoterRequest(BaseModel):
    ensembl_gene_id: str
    gene_symbol: str = ""
    organism: str = "homo_sapiens"
    oligo_length: int = 20
    promoter_index: int = 1
    max_candidates: int = 12


@router.post("/api/mechanisms/A32/alternative-promoter")
async def design_a32_promoter(payload: AlternativePromoterRequest):
    """A32 — oligos across the TSS-proximal region of one promoter."""
    from services.transcript_architecture_service import design_alternative_promoter
    try:
        return design_alternative_promoter(**payload.dict())
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
