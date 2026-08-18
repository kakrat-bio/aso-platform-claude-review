import {
  RnaNeutralizationResponse,
  RepeatTractResponse,
  NeutralizationDesignOptions,
} from "@/types/rnaNeutralization";

const API_BASE = process.env.NEXT_PUBLIC_API_BASE_URL ?? "http://localhost:8000";

async function post<T>(path: string, body: unknown, fallback: string): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || fallback);
  }
  return res.json();
}

export async function fetchNeutralizationOptions(): Promise<NeutralizationDesignOptions> {
  const res = await fetch(`${API_BASE}/api/rna-neutralization/options`, {
    cache: "no-store",
  });
  if (!res.ok) throw new Error("Could not load neutralization design options.");
  return res.json();
}

/** Which repeat unit applies to this gene, and whether it is curated or asserted. */
export async function fetchRepeatTract(params: {
  geneSymbol: string;
  repeatUnit?: string | null;
  estimatedRepeatCount?: string | null;
}): Promise<RepeatTractResponse> {
  return post<RepeatTractResponse>(
    "/api/rna-neutralization/repeat-tract",
    {
      gene_symbol: params.geneSymbol,
      repeat_unit: params.repeatUnit ?? null,
      estimated_repeat_count: params.estimatedRepeatCount ?? null,
    },
    "Could not resolve the repeat tract.",
  );
}

export async function generateNeutralizationCandidates(params: {
  geneSymbol: string;
  mechanismId: string;
  neutralizationMode: string;
  repeatUnit?: string | null;
  estimatedRepeatCount?: string | null;
  mirnaSequence?: string | null;
  oligoLength?: number;
  chemistry?: string;
  modifications?: string[];
  deliveryContext?: string | null;
  targetRbp?: string | null;
  maxCandidates?: number;
}): Promise<RnaNeutralizationResponse> {
  return post<RnaNeutralizationResponse>(
    "/api/rna-neutralization/candidates",
    {
      gene_symbol: params.geneSymbol,
      mechanism_id: params.mechanismId,
      neutralization_mode: params.neutralizationMode,
      repeat_unit: params.repeatUnit ?? null,
      estimated_repeat_count: params.estimatedRepeatCount ?? null,
      mirna_sequence: params.mirnaSequence ?? null,
      oligo_length: params.oligoLength ?? 17,
      chemistry: params.chemistry ?? "moe_full_ps",
      modifications: params.modifications ?? [],
      delivery_context: params.deliveryContext ?? null,
      target_rbp: params.targetRbp ?? null,
      max_candidates: params.maxCandidates ?? 12,
    },
    "Could not generate neutralization candidates.",
  );
}
