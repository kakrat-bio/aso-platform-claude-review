import {
  IsoformEngineeringInputs,
  IsoformEngineeringResponse,
  DesignOptions,
  IsoformCandidate,
} from "@/types/isoformEngineering";
import { getToken } from "@/lib/auth";

const API_BASE = process.env.NEXT_PUBLIC_API_BASE_URL ?? "http://localhost:8000";

export async function fetchDesignOptions(): Promise<DesignOptions> {
  const res = await fetch(`${API_BASE}/api/isoform-engineering/options`, {
    cache: "no-store",
  });
  if (!res.ok) throw new Error("Could not load design options.");
  return res.json();
}

export async function generateConstructs(
  params: IsoformEngineeringInputs
): Promise<IsoformEngineeringResponse> {
  const res = await fetch(`${API_BASE}/api/isoform-engineering/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      target_symbol: params.targetSymbol,
      isoform_goal: params.isoformGoal,
      target_exon_locus: params.targetExonLocus,
      splice_element_target: params.spliceElementTarget,
      steric_chemistry: params.stericChemistry,
      enforce_in_frame: params.enforceInFrame,
      aso_length: params.asoLength ?? 20,
      max_candidates: params.maxCandidates ?? 12,
    }),
  });
  if (!res.ok) {
    const body = await res.json().catch(() => ({}));
    throw new Error(body.detail || "Could not generate constructs.");
  }
  return res.json();
}

export async function emailIsoformReport(reportContent: string, filename: string): Promise<{ message: string }> {
  const token = getToken();
  if (!token) throw new Error("Sign in to email a report to your registered address.");
  const res = await fetch(`${API_BASE}/api/isoform-engineering/email-report`, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      Authorization: `Bearer ${token}`,
    },
    body: JSON.stringify({ report_content: reportContent, filename }),
  });
  const body = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(body.detail || "Could not email the report.");
  return body;
}
