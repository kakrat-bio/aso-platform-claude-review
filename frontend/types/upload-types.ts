export interface OrfHit {
  strand: "+" | "-";
  frame: number;
  start: number;
  end: number;
  length: number;
  proteinLength: number;
}

export interface ValidateResponse {
  valid: boolean;
  error?: string;
  sequence?: string;
  sequenceType?: "dna" | "rna" | "protein" | "unknown";
  length?: number;
  gcContent?: number | null;
  invalidChars?: string[];
  features?: string[];
  orfs?: OrfHit[];
  filename?: string | null;
  hasPolyA?: boolean;
  hasPolyG?: boolean;
}

export interface SpecificityHeuristic {
  lengthBasedRiskEstimate: "High" | "Medium" | "Low";
  note: string;
  internalRepetitiveness: number;
  recommendedMinLength: number;
  disclaimer: string;
}

export interface SecondaryStructureEstimate {
  estimatedMfe: number;
  palindromicRegions: number;
  palindromePositions: number[];
  gcContent: number;
  hairpinRisk: "High" | "Medium" | "Low";
}

export interface ImmuneMotifHit {
  motif: string;
  label: string;
  start: number;
  end: number;
}

export interface GcCurvePoint {
  position: number;
  gc: number;
}

export interface NucleotideComposition {
  A: number;
  C: number;
  G: number;
  T: number;
  U: number;
}

export type Modality = "aso" | "sirna" | "mrna" | "sgrna";

export interface ModalityAnalysis {
  recommendations?: string[];
  recommendedChemistry?: string;
  optimalLength?: string;
  targetRegion?: string;
  strand?: string;
  thermodynamicBias?: string;
  needsCodonOptimization?: boolean;
  needsPolyA?: boolean;
  needsUTR?: boolean;
  nucleosideModifications?: string[];
  casProtein?: string;
  offTargetMitigation?: string;
}

export interface RestrictionSite {
  enzyme: string;
  recognitionSite: string;
  cutPosition: number;
  strand: "+" | "-";
  overhang: "5'" | "3'" | "blunt";
}

export interface MiRNATarget {
  mirnaId: string;
  seedSequence: string;
  start: number;
  end: number;
  /** null — no TargetScan context++ scoring is wired (F7). */
  bindingScore: number | null;
  /** GC fraction of the matched seed motif. Real. */
  seedGcContent?: number;
  conservationNote: string;
}

export interface Hairpin {
  start: number;
  end: number;
  stemLength: number;
  loopSize: number;
  stabilityScore: number;
  type: "hairpin" | "bulge" | "internal_loop";
}

export interface KmerFrequency {
  k: number;
  totalKmers: number;
  uniqueKmers: number;
  repeats: { kmer: string; count: number; positions: number[] }[];
  shannonEntropy: number;
}

export interface ThermoProfile {
  avgEnthalpy: number;
  avgEntropy: number;
  freeEnergy37: number;
  gcEnrichment: number;
  atEnrichment: number;
  stabilityClass: "stable" | "moderate" | "unstable";
  notes: string[];
}

export interface DotPlotPoint {
  x: number;
  y: number;
  matchLen: number;
}

export interface ModificationLandscapePoint {
  position: number;
  accessibilityScore: number;
  recommendedModification: string;
  confidenceLevel: "high" | "medium" | "low";
}

export interface RiskScores {
  specificity: number;
  stability: number;
  immunogenicity: number;
  delivery: number;
  toxicity: number;
  overall: number;
}

export interface PhysicochemicalProfile {
  molecularWeight: number;
  netCharge: number;
  hydrophobicityIndex: number;
  hydrophobicityProfile: { position: number; value: number }[];
  chargeProfile: { position: number; value: number }[];
}

export interface StabilityIndexPoint {
  position: number;
  rnaseH: number;
  duplexStability: number;
  singleStrandStability: number;
}

export interface AnalyzeResponse {
  sequence: string;
  sequenceType: "dna" | "rna" | "protein" | "unknown";
  length: number;
  gcContent: number | null;
  offTarget: SpecificityHeuristic;
  secondaryStructure: SecondaryStructureEstimate;
  immuneScreen: ImmuneMotifHit[];
  modality: ModalityAnalysis;
  gcCurve: GcCurvePoint[];
  composition: Record<string, number>;
  orfs: OrfHit[];
  meltingTemp?: {
    tmNearestNeighbor: number;
    tmBasicGC: number;
    length: number;
    gcContent?: number;
    method: string;
    note: string;
  } | null;
  complexity?: {
    dinucRepeats: { pattern: string; start: number; end: number; repeats: number }[];
    trinucRepeats: { pattern: string; count: number; positions: number[] }[];
    gcRichRegions: { start: number; end: number; length: number }[];
    atRichRegions: { start: number; end: number; length: number }[];
    selfComplementarity: { sequence: string; position: number; size: number }[];
    complexityScore: number;
  } | null;
  codonUsage?: {
    codons: { codon: string; position: number; adaptiveness: number; isRare: boolean }[];
    cai: number;
    rareCodons: { codon: string; position: number; adaptiveness: number }[];
    totalCodons: number;
    note: string;
  } | null;
  modificationScores?: {
    modality: string;
    scores: Record<string, { score: number; rationale: string }>;
    overallScore: number;
  } | null;
  energyProfile?: { position: number; energy: number }[];
  restrictionSites?: RestrictionSite[];
  mirnaTargets?: MiRNATarget[];
  hairpins?: Hairpin[];
  kmerFrequency?: KmerFrequency | null;
  thermoProfile?: ThermoProfile | null;
  dotPlot?: DotPlotPoint[];
  modificationLandscape?: ModificationLandscapePoint[];
  riskScores?: RiskScores;
  physicochemical?: PhysicochemicalProfile;
  stabilityIndex?: StabilityIndexPoint[];
  grnaCandidates?: GrnaCandidate[];
  proteinAnalysis?: {
    aminoAcidComposition: Record<string, number>;
    molecularWeight: number;
    length: number;
    hydrophobicFraction: number;
    hydrophilicFraction: number;
    chargedFraction: number;
    aromaticFraction: number;
  };
}

export interface GrnaCandidate {
  id: string;
  position: number;
  sequence: string;
  pam: string;
  strand: "+" | "-";
  score: number;
  gc: number;
  selfComplementarity: number;
  /** Fraction of 6-mers that repeat. NOT a genomic off-target count —
   *  no alignment is performed anywhere in this service. */
  internalRepetitiveness: number;
  polyT: boolean;
  color: string;
  specificityScore: number;
  efficiencyScore: number;
  mismatchDistribution: number[];
}
