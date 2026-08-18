"""TG05 RNA-neutralization design endpoints.

Mechanism *ranking* for TG05 lives on `/api/mechanisms/rna-neutralization`
(a filter over the shared arbitration pass); this is the design step after
it. Every response carries `tractProvenance`, because whether the repeat
unit came from the curated catalogue or from what the user typed is the
difference between a lookup and a hypothesis.
"""

from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from services.rna_neutralization_service import (
    NEUTRALIZATION_MECHANISM_CHEMISTRY,
    generate_neutralization_candidates,
    get_neutralization_design_options,
    resolve_repeat_tract,
)

router = APIRouter()

VALID_MODES = {"steric_repeat_masking", "microrna_antagomir", "aptamer_decoy"}


class NeutralizationCandidateRequest(BaseModel):
    gene_symbol: str
    mechanism_id: str
    neutralization_mode: str
    repeat_unit: Optional[str] = None
    estimated_repeat_count: Optional[str] = None
    mirna_sequence: Optional[str] = None
    oligo_length: int = 17
    chemistry: str = "moe_full_ps"
    modifications: Optional[list[str]] = None
    delivery_context: Optional[str] = None
    target_rbp: Optional[str] = None
    max_candidates: int = 12


class RepeatTractRequest(BaseModel):
    gene_symbol: str
    repeat_unit: Optional[str] = None
    estimated_repeat_count: Optional[str] = None


@router.get("/api/rna-neutralization/options")
async def neutralization_options():
    """Chemistries, length range, and the curated repeat-unit reference."""
    return get_neutralization_design_options()


@router.post("/api/rna-neutralization/repeat-tract")
async def repeat_tract(payload: RepeatTractRequest):
    """Which repeat unit applies to this gene, and how we know."""
    tract = resolve_repeat_tract(
        payload.gene_symbol.strip().upper(),
        payload.repeat_unit,
        payload.estimated_repeat_count,
    )
    return {"geneSymbol": payload.gene_symbol.strip().upper(), **tract}


@router.post("/api/rna-neutralization/candidates")
async def neutralization_candidates(payload: NeutralizationCandidateRequest):
    """Phase-shifted oligos across the toxic RNA, ranked by duplex energy."""
    if payload.neutralization_mode not in VALID_MODES:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown neutralization_mode: {payload.neutralization_mode}",
        )
    if (payload.mechanism_id not in NEUTRALIZATION_MECHANISM_CHEMISTRY
            and payload.mechanism_id != "A25"):
        raise HTTPException(
            status_code=400,
            detail=(
                f"{payload.mechanism_id} is not a TG05 mechanism. Expected one "
                f"of {', '.join(sorted(NEUTRALIZATION_MECHANISM_CHEMISTRY))} "
                f"or A25 (flagged, not designed)."
            ),
        )

    result = generate_neutralization_candidates(
        gene_symbol=payload.gene_symbol.strip().upper(),
        mechanism_id=payload.mechanism_id,
        neutralization_mode=payload.neutralization_mode,
        repeat_unit=payload.repeat_unit,
        estimated_repeat_count=payload.estimated_repeat_count,
        oligo_length=payload.oligo_length,
        chemistry=payload.chemistry,
        modifications=payload.modifications,
        delivery_context=payload.delivery_context,
        target_rbp=payload.target_rbp,
        mirna_sequence=payload.mirna_sequence,
        max_candidates=payload.max_candidates,
    )
    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "inputs": {
            "mechanismId": payload.mechanism_id,
            "neutralizationMode": payload.neutralization_mode,
            "repeatUnit": payload.repeat_unit,
            "estimatedRepeatCount": payload.estimated_repeat_count,
            "oligoLength": payload.oligo_length,
            "chemistry": payload.chemistry,
            "deliveryContext": payload.delivery_context,
            "targetRbp": payload.target_rbp,
        },
        **result,
    }
