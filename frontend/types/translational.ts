/** Where the transcript data came from. Never silently a substitute. */
export interface DataProvenance {
  status: "live" | "cached" | "curated" | "unavailable";
  source?: string | null;
  fetchedAt?: number | null;
  ageSeconds?: number | null;
  stale?: boolean;
  note?: string | null;
}

export interface TranslationalDesignOptions {
  chemistries: { id: string; label: string; description: string }[];
  lengthRange: { min: number; max: number; default: number; step: number };
  mechanismRegions: Record<string, string>;
}

export interface KozakContext {
  sequence: string;
  start: number;
  end: number;
  minus3: string;
  plus4: string;
  consensusMatches: string[];
  strength: "strong" | "adequate" | "weak";
  note: string;
}

export interface StructuredWindow {
  start: number;
  end: number;
  mfe: number;
  structure: string;
  domain?: string;
  /** Always false for IRES entries — no validated IRES predictor is wired. */
  predicted?: boolean;
  note?: string;
}

export interface TranslationalTargetResponse {
  geneSymbol: string;
  status: string;
  dataProvenance?: DataProvenance;
  message?: string;
  canonicalTranscript?: { id?: string } | null;
  mrnaSequence?: string;
  utr5?: { sequence: string; start: number; end: number } | null;
  utr3?: { sequence: string; start: number; end: number } | null;
  cdsStart?: number | null;
  uorfs?: { start: number; end: number; sequence?: string }[];
  kozak?: KozakContext | null;
  structuredElements?: StructuredWindow[];
  ires?: StructuredWindow[];
  polyASite?: { start: number; end: number; sequence: string } | null;
}

/**
 * Computed from the sequence. Every field here is a real measurement of the
 * molecule, not an estimate.
 */
export interface TranslationalRealMetrics {
  targetDuplexEnergy: number;
  meltingTempC: number;
  selfStructureMfe: number;
  gcContent: number;
  lengthNt: number;
  elementOverlapNt: number;
  provenance: string;
}

/** Ordinal chemistry and length rules of thumb. Not measurements. */
export interface TranslationalHeuristics {
  nucleaseResistance: number;
  cellularUptake: number;
  provenance: string;
}

export interface TranslationalCandidate {
  rank?: number;
  sequence: string;
  sequenceAlphabet: "DNA";
  targetSite: string;
  targetSiteAlphabet: "RNA";
  targetStart: number;
  targetEnd: number;
  targetRegion: string;
  targetElement: string;
  mechanismId: string;
  chemistry: string;
  modifications: string[];
  deliveryContext?: string | null;
  elementOverlapNt: number;
  realMetrics: TranslationalRealMetrics;
  heuristicEstimates: TranslationalHeuristics;
  /**
   * A 0-1 ranking signal from element coverage and computed duplex energy.
   * Deliberately NOT a predicted fold-change: no coefficient has been fitted
   * and there is no calibration set for translational effect size.
   */
  elementEngagement: number;
  interpretation: string;
}

export interface TranslationalCandidateResponse {
  geneSymbol: string;
  status: string;
  dataProvenance?: DataProvenance;
  message?: string;
  mechanismId?: string;
  targetElement?: string;
  translationalGoal?: string;
  elementRegion?: { start: number; end: number; label: string };
  candidates: TranslationalCandidate[];
  totalScanned?: number;
  rbpNote?: string | null;
  scoringNote?: string;
  inputs?: Record<string, unknown>;
}
