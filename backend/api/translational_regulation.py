"""TG06 translational-regulation design endpoints.

Mirrors backend/api/gene_silencing.py. Mechanism *ranking* for TG06 lives on
`/api/mechanisms/translational-regulation` (a filter over the shared
arbitration pass); this module is the design step that follows it — fetching
the transcript's regulatory elements and tiling oligos against one of them.

Every response carries `dataProvenance`, because the transcript may have come
from a live Ensembl call, from a real earlier fetch replayed during an
outage, or not at all. Those are different situations and the caller is told
which one it got.
"""

from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from services import real_data_cache as RDC
from services.translational_regulation_service import (
    TRANSLATIONAL_MECHANISM_CHEMISTRY,
    generate_translational_candidates,
    get_translational_design_options,
    get_translational_target,
)

router = APIRouter()

VALID_ELEMENTS = {
    "5p_utr", "5p_uorf", "3p_utr_mirna", "structured_element",
    "ires_element", "kozak_consensus", "polya_site",
}
VALID_GOALS = {"enhance", "suppress"}


class TranslationalTargetRequest(BaseModel):
    ensembl_gene_id: str
    gene_symbol: str = ""
    organism: str = "homo_sapiens"


class TranslationalCandidateRequest(BaseModel):
    ensembl_gene_id: str
    gene_symbol: str = ""
    organism: str = "homo_sapiens"
    target_element: str
    translational_goal: str
    mechanism_id: str
    aso_length: int = 20
    chemistry: str = "pmo"
    modifications: Optional[list[str]] = None
    delivery_context: Optional[str] = None
    target_rbp: Optional[str] = None
    max_candidates: int = 20


@router.get("/api/translational-regulation/options")
async def translational_options():
    """Chemistries, length range, and which region each mechanism acts on."""
    return get_translational_design_options()


@router.post("/api/translational-regulation/target")
async def translational_target(payload: TranslationalTargetRequest):
    """Transcript structure: UTRs, uORFs, Kozak context, structured elements.

    Returns status `unavailable` rather than a constructed transcript when
    Ensembl cannot answer and nothing real has been cached.
    """
    result = get_translational_target(
        payload.ensembl_gene_id, payload.gene_symbol, payload.organism,
    )
    if result.get("status") == RDC.UNAVAILABLE:
        return {
            "geneSymbol": payload.gene_symbol.strip().upper(),
            "status": RDC.UNAVAILABLE,
            "dataProvenance": result.get("dataProvenance"),
            "message": (
                f"No transcript structure is available for "
                f"{payload.gene_symbol or payload.ensembl_gene_id}: Ensembl "
                f"did not answer and nothing real has been cached for it. "
                f"Regulatory elements are read off the sequence, so none can "
                f"be reported."
            ),
        }
    return {"geneSymbol": payload.gene_symbol.strip().upper(), **result}


@router.post("/api/translational-regulation/candidates")
async def translational_candidates(payload: TranslationalCandidateRequest):
    """Tile oligos across the chosen regulatory element and rank them."""
    if payload.target_element not in VALID_ELEMENTS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown target_element: {payload.target_element}",
        )
    if payload.translational_goal not in VALID_GOALS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown translational_goal: {payload.translational_goal}",
        )
    if payload.mechanism_id not in TRANSLATIONAL_MECHANISM_CHEMISTRY:
        raise HTTPException(
            status_code=400,
            detail=(
                f"{payload.mechanism_id} is not a TG06 mechanism. Expected one "
                f"of {', '.join(sorted(TRANSLATIONAL_MECHANISM_CHEMISTRY))}."
            ),
        )

    target = get_translational_target(
        payload.ensembl_gene_id, payload.gene_symbol, payload.organism,
    )
    provenance = target.get("dataProvenance")
    if target.get("status") == RDC.UNAVAILABLE:
        return {
            "geneSymbol": payload.gene_symbol.strip().upper(),
            "status": RDC.UNAVAILABLE,
            "dataProvenance": provenance,
            "candidates": [],
            "message": (
                "No transcript sequence is available for this target, so no "
                "oligo can be designed against it."
            ),
        }

    result = generate_translational_candidates(
        target_element=payload.target_element,
        translational_goal=payload.translational_goal,
        mechanism_id=payload.mechanism_id,
        aso_length=payload.aso_length,
        chemistry=payload.chemistry,
        modifications=payload.modifications,
        target=target,
        delivery_context=payload.delivery_context,
        target_rbp=payload.target_rbp,
        max_candidates=payload.max_candidates,
    )
    return {
        "geneSymbol": payload.gene_symbol.strip().upper(),
        "dataProvenance": provenance,
        "inputs": {
            "targetElement": payload.target_element,
            "translationalGoal": payload.translational_goal,
            "mechanismId": payload.mechanism_id,
            "asoLength": payload.aso_length,
            "chemistry": payload.chemistry,
            "modifications": payload.modifications or [],
            "deliveryContext": payload.delivery_context,
            "targetRbp": payload.target_rbp,
        },
        **result,
    }
