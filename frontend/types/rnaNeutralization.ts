/**
 * How the repeat unit was established. A curated catalogue entry and a unit
 * the user typed are different evidence, and the UI must be able to say which
 * one produced the tract.
 */
export interface RepeatTractProvenance {
  unit: string | null;
  region: string | null;
  pathogenicMin?: string | null;
  provenance: "confirmed" | "user_asserted" | "unavailable";
  source: string | null;
  note: string;
}

export interface RepeatTractResponse extends RepeatTractProvenance {
  geneSymbol: string;
}

export interface NeutralizationDesignOptions {
  chemistries: { id: string; label: string; description: string }[];
  lengthRange: { min: number; max: number; default: number; step: number };
  knownRepeatUnits: { unit: string; disease: string }[];
  pathogenicThreshold: number;
}

/** Computed from the sequence. */
export interface NeutralizationRealMetrics {
  targetDuplexEnergy: number;
  meltingTempC: number;
  selfStructureMfe: number;
  gcContent: number;
  lengthNt: number;
  provenance: string;
}

/** Ordinal chemistry and length rules of thumb. Not measurements. */
export interface NeutralizationHeuristics {
  nucleaseResistance: number;
  cellularUptake: number;
  provenance: string;
}

export interface RnaNeutralizationCandidate {
  rank?: number;
  sequence: string;
  sequenceAlphabet: "DNA";
  targetSite: string;
  targetSiteAlphabet: "RNA";
  phase: number;
  tilingPattern: string;
  mechanismId: string;
  chemistry: string;
  modifications: string[];
  deliveryContext?: string | null;
  realMetrics: NeutralizationRealMetrics;
  heuristicEstimates: NeutralizationHeuristics;
}

export interface RnaNeutralizationResponse {
  geneSymbol: string;
  status: string;
  message?: string;
  mechanismId?: string;
  neutralizationMode?: string;
  tractProvenance?: RepeatTractProvenance;
  targetSequence?: string;
  candidates: RnaNeutralizationCandidate[];
  /**
   * Caveats that must be shown, not buried. For repeat masking this carries
   * the off-target warning: a (CAG)n oligo is complementary to every
   * CAG-repeat transcript by construction, and no scan is wired.
   */
  notes?: string[];
  scoringNote?: string;
  inputs?: Record<string, unknown>;
}

export interface DesignPipelineState {
  step: number;
  selectedCandidate: RnaNeutralizationCandidate | null;
  chemicalModifications: string[];
  deliveryConjugation: string;
  secondaryStructurePassed: boolean;
  selfDimerPassed: boolean;
}
