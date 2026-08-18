import {
  TranslationalTargetResponse,
  TranslationalCandidateResponse,
  TranslationalDesignOptions,
} from "@/types/translational";

const API_BASE = process.env.NEXT_PUBLIC_API_BASE_URL ?? "http://localhost:8000";

async function post<T>(path: string, body: unknown, fallbackError: string): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || fallbackError);
  }
  return res.json();
}

export async function fetchTranslationalDesignOptions(): Promise<TranslationalDesignOptions> {
  const res = await fetch(`${API_BASE}/api/translational-regulation/options`, {
    cache: "no-store",
  });
  if (!res.ok) throw new Error("Could not load translational design options.");
  return res.json();
}

/**
 * Transcript structure: UTRs, uORFs, Kozak context, structured elements.
 *
 * May come back with status "unavailable" — Ensembl unreachable and nothing
 * real cached. That is a legitimate answer, not an error, and the caller
 * should render the message rather than retrying blindly.
 */
export async function fetchTranslationalTarget(params: {
  ensemblGeneId: string;
  geneSymbol?: string;
  organism?: string;
}): Promise<TranslationalTargetResponse> {
  return post<TranslationalTargetResponse>(
    "/api/translational-regulation/target",
    {
      ensembl_gene_id: params.ensemblGeneId,
      gene_symbol: params.geneSymbol ?? "",
      organism: params.organism ?? "homo_sapiens",
    },
    "Could not load the translational target.",
  );
}

export async function generateTranslationalCandidates(params: {
  ensemblGeneId: string;
  geneSymbol?: string;
  organism?: string;
  targetElement: string;
  translationalGoal: string;
  mechanismId: string;
  asoLength?: number;
  chemistry?: string;
  modifications?: string[];
  deliveryContext?: string | null;
  targetRbp?: string | null;
  maxCandidates?: number;
}): Promise<TranslationalCandidateResponse> {
  return post<TranslationalCandidateResponse>(
    "/api/translational-regulation/candidates",
    {
      ensembl_gene_id: params.ensemblGeneId,
      gene_symbol: params.geneSymbol ?? "",
      organism: params.organism ?? "homo_sapiens",
      target_element: params.targetElement,
      translational_goal: params.translationalGoal,
      mechanism_id: params.mechanismId,
      aso_length: params.asoLength ?? 20,
      chemistry: params.chemistry ?? "pmo",
      modifications: params.modifications ?? [],
      delivery_context: params.deliveryContext ?? null,
      target_rbp: params.targetRbp ?? null,
      max_candidates: params.maxCandidates ?? 20,
    },
    "Could not generate translational candidates.",
  );
}
