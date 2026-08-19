"""
Isoform Engineering endpoints.
"""

from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import Optional

from services.email_service import send_report_email
from services.auth_service import get_current_user
from services.isoform_engineering_service import (
    get_isoform_engineering_design_options,
    generate_isoform_candidates,
)
from database.models import User


class ReportEmailRequest(BaseModel):
    report_content: str
    filename: str = "isoform-engineering-report.txt"


router = APIRouter()


@router.get("/api/isoform-engineering/options")
async def isoform_engineering_options():
    """Input options for isoform engineering form."""
    return get_isoform_engineering_design_options()


@router.post("/api/isoform-engineering/generate")
async def generate_isoform_constructs(payload: dict):
    """Generate isoform engineering construct candidates from real exon data."""
    target_symbol = payload.get("target_symbol", "")
    isoform_goal = payload.get("isoform_goal", "")
    target_exon_locus = payload.get("target_exon_locus", "")
    splice_element_target = payload.get("splice_element_target", "")
    steric_chemistry = payload.get("steric_chemistry", "")
    enforce_in_frame = payload.get("enforce_in_frame", True)
    aso_length = int(payload.get("aso_length", 20))
    max_candidates = int(payload.get("max_candidates", 12))
    organism = payload.get("organism", "homo_sapiens")

    try:
        result = generate_isoform_candidates(
            target_symbol=target_symbol,
            isoform_goal=isoform_goal,
            target_exon_locus=target_exon_locus,
            splice_element_target=splice_element_target,
            steric_chemistry=steric_chemistry,
            enforce_in_frame=enforce_in_frame,
            aso_length=aso_length,
            max_candidates=max_candidates,
            organism=organism,
        )
        return result
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@router.post("/api/isoform-engineering/email-report")
async def email_isoform_report(
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
        payload.filename or "isoform-engineering-report.txt",
    )
    if not sent:
        raise HTTPException(status_code=503, detail=message)
    return {"ok": True, "message": message}
