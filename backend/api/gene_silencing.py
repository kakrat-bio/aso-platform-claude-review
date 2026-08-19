"""
Gene Silencing API — target analysis + ASO candidate generation.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional
from fastapi import Depends
from database.models import User
from services.auth_service import get_current_user
from services.email_service import send_report_email

from services.gene_silencing_service import (
    get_target_analysis,
    generate_candidates,
    parse_hgvs_c,
    CHEMISTRY_OPTIONS,
    MODIFICATION_OPTIONS,
    LENGTH_RANGE,
)
from services.rna_processing_service import (
    generate_rna_processing_candidates,
    RNA_PROCESSING_MECHANISMS,
)
from services.variant_details_service import get_clinvar_variants
from services.notification_service import add_notification
from services.sequence_liability_service import get_sequence_liabilities

router = APIRouter()


class CandidateRequest(BaseModel):
    ensembl_gene_id: str
    mechanism_id: str
    target_exon_indices: Optional[list[int]] = None
    aso_length: int = 18
    chemistry: str = "gapmer"
    modifications: list[str] = []
    delivery_context: Optional[str] = None
    defect_type: Optional[str] = None
    silencing_scope: Optional[str] = None
    known_variant: Optional[str] = None


class ReportEmailRequest(BaseModel):
    report_content: str
    filename: str


@router.get("/api/gene-silencing/target/{ensembl_gene_id}")
async def target_analysis(
    ensembl_gene_id: str,
    gene_symbol: Optional[str] = None,
    organism: Optional[str] = None,
):
    """Fetch transcript / exon structure for the confirmed gene.

    Optional query params ``gene_symbol`` and ``organism`` enable a symbol-based
    fallback when the Ensembl ID lookup returns no exon data.
    """
    result = get_target_analysis(ensembl_gene_id, gene_symbol=gene_symbol or "", organism=organism or "")
    if not result.get("exons"):
        raise HTTPException(status_code=404, detail="No exon data found for this gene.")
    return result


@router.get("/api/gene-silencing/options")
async def design_options():
    """Available chemistry, modification, and length options."""
    return {
        "chemistryOptions": CHEMISTRY_OPTIONS,
        "modificationOptions": MODIFICATION_OPTIONS,
        "lengthRange": LENGTH_RANGE,
    }


@router.post("/api/gene-silencing/generate")
async def generate_aso_candidates(payload: CandidateRequest):
    """Generate ranked candidates using the selected mechanism's design rules.

    A7-A11 (TG04, RNA Processing Modulation) dispatch to
    rna_processing_service instead of the RNase H-oriented gene-silencing
    designer — they are steric-blocking, splice-/3'-end-processing
    mechanisms with a different targeting geometry (splice junctions /
    polyadenylation signal), not CDS-region knockdown.
    """
    target = get_target_analysis(payload.ensembl_gene_id)
    if not target.get("mrnaSequence"):
        raise HTTPException(status_code=404, detail="Could not fetch mRNA sequence.")

    try:
        if payload.mechanism_id in RNA_PROCESSING_MECHANISMS:
            candidates = generate_rna_processing_candidates(
                target_exon_indices=payload.target_exon_indices,
                aso_length=payload.aso_length,
                chemistry=payload.chemistry,
                modifications=payload.modifications,
                exons=target.get("exons", []),
                canonical_transcript=target.get("canonicalTranscript"),
                mechanism_id=payload.mechanism_id,
                utr3_sequence=target.get("utr3Sequence"),
                delivery_context=payload.delivery_context,
            )
        else:
            candidates = generate_candidates(
                target_exon_indices=payload.target_exon_indices,
                aso_length=payload.aso_length,
                chemistry=payload.chemistry,
                modifications=payload.modifications,
                mrna_sequence=target["mrnaSequence"],
                exons=target.get("exons", []),
                mechanism_id=payload.mechanism_id,
                delivery_context=payload.delivery_context,
                defect_type=payload.defect_type,
                silencing_scope=payload.silencing_scope,
                known_variant=payload.known_variant,
            )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    add_notification(
        "analysis",
        f"Generated {len(candidates)} ASO candidates",
        f"Candidate design completed for {payload.ensembl_gene_id}.",
    )

    # Was an "ADMET" block per candidate. Absorption, distribution,
    # metabolism, excretion and half-life are set by the backbone chemistry and
    # any conjugate, and this call never even passed the chemistry — see the
    # module docstring in services/admet_service.py. What a sequence does
    # determine (CpG/TLR9, G-quadruplex, uridine tracts) is reported instead,
    # with the chemistry echoed alongside so the reader can see what the flags
    # do and do not account for.
    top_liabilities = {}
    if candidates:
        top_liabilities = get_sequence_liabilities(
            aso_sequence=candidates[0]["sequence"],
            chemistry=candidates[0].get("chemistry") or payload.chemistry,
        )

    enriched_candidates = []
    for candidate in candidates:
        enriched = {
            **candidate,
            "sequenceLiabilities": get_sequence_liabilities(
                aso_sequence=candidate["sequence"],
                chemistry=candidate.get("chemistry") or payload.chemistry,
            ),
        }
        enriched_candidates.append(enriched)

    return {
        "geneId": payload.ensembl_gene_id,
        "mechanismId": payload.mechanism_id,
        "targetExons": None if payload.mechanism_id == "A2" else payload.target_exon_indices,
        "chemistry": candidates[0]["chemistry"] if candidates else payload.chemistry,
        "modifications": payload.modifications,
        "asoLength": candidates[0]["length"] if candidates else payload.aso_length,
        "totalExons": len(target.get("exons", [])),
        "cdsLength": target.get("cdsLength"),
        "mechanismNotes": candidates[0].get("mechanismNotes", "") if candidates else "",
        "isAlleleSpecific": (payload.silencing_scope or "").lower().strip() == "allele_specific",
        "variantParse": parse_hgvs_c(payload.known_variant) if payload.known_variant else None,
        "sequenceLiabilities": top_liabilities,
        "candidates": enriched_candidates,
    }


@router.post("/api/gene-silencing/email-report")
async def email_aso_report(
    payload: ReportEmailRequest,
    user: User = Depends(get_current_user),
):
    """Email a generated report to the authenticated user's account email."""
    if not payload.report_content.strip():
        raise HTTPException(status_code=422, detail="Report content is required.")
    sent, message = send_report_email(
        user.email,
        user.name,
        payload.report_content,
        payload.filename or "aso-candidate-report.txt",
    )
    if not sent:
        raise HTTPException(status_code=503, detail=message)
    return {"ok": True, "message": message}


@router.get("/api/gene-silencing/variants/{ensembl_gene_id}")
async def clinvar_variants(ensembl_gene_id: str):
    """Fetch ClinVar pathogenic variants for allele-specific silencing."""
    variants = get_clinvar_variants(ensembl_gene_id)
    return {"geneId": ensembl_gene_id, "variants": variants}
