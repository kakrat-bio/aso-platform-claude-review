export interface ValidationReport {
  valid: boolean;
  error?: string;
  sequence: string;
  sequenceType: "dna" | "rna" | "protein" | "unknown";
  length: number;
  gcContent: number | null;
  invalidChars: string[];
  features: string[];
  orfs: OrfInfo[];
  filename?: string;
  hasPolyA: boolean;
  hasPolyG: boolean;
}

export interface OrfInfo {
  strand: string;
  frame: number;
  start: number;
  end: number;
  length: number;
  proteinLength: number;
}

export interface OffTargetResult {
  lengthBasedRiskEstimate: "High" | "Medium" | "Low";
  note: string;
  internalRepetitiveness: number;
  recommendedMinLength: number;
  disclaimer: string;
}

export interface SecondaryStructureResult {
  estimatedMfe: number;
  palindromicRegions: number;
  palindromePositions: number[];
  gcContent: number;
  hairpinRisk: "High" | "Medium" | "Low";
}

export interface ImmuneMotif {
  motif: string;
  label: string;
  start: number;
  end: number;
}

export interface ModalityResult {
  recommendedChemistry?: string;
  recommendations: string[];
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

export interface GcWindow {
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

export interface MeltingTemp {
  tmNearestNeighbor: number;
  tmBasicGC: number;
  length: number;
  gcContent?: number;
  method: string;
  note: string;
}

export interface ComplexityRegion {
  start: number;
  end: number;
  length: number;
  pattern?: string;
  repeats?: number;
  count?: number;
  positions?: number[];
}

export interface SequenceComplexity {
  dinucRepeats: (ComplexityRegion & { pattern: string; repeats: number })[];
  trinucRepeats: (ComplexityRegion & { pattern: string; count: number; positions: number[] })[];
  gcRichRegions: ComplexityRegion[];
  atRichRegions: ComplexityRegion[];
  selfComplementarity: { sequence: string; position: number; size: number }[];
  complexityScore: number;
}

export interface CodonInfo {
  codon: string;
  position: number;
  adaptiveness: number;
  isRare: boolean;
}

export interface CodonUsage {
  codons: CodonInfo[];
  cai: number;
  rareCodons: { codon: string; position: number; adaptiveness: number }[];
  totalCodons: number;
  note: string;
}

export interface ModScore {
  score: number;
  rationale: string;
}

export interface ModificationScores {
  modality: string;
  scores: Record<string, ModScore>;
  overallScore: number;
}

export interface EnergyPoint {
  position: number;
  energy: number;
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

export interface AnalysisReport {
  sequence: string;
  sequenceType: string;
  length: number;
  gcContent: number | null;
  offTarget: OffTargetResult;
  secondaryStructure: SecondaryStructureResult;
  immuneScreen: ImmuneMotif[];
  modality: ModalityResult;
  gcCurve: GcWindow[];
  composition: Record<string, number>;
  orfs: OrfInfo[];
  meltingTemp?: MeltingTemp | null;
  complexity?: SequenceComplexity | null;
  codonUsage?: CodonUsage | null;
  modificationScores?: ModificationScores | null;
  energyProfile?: EnergyPoint[];
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
  offTargets: number;
  polyT: boolean;
  color: string;
  specificityScore: number;
  efficiencyScore: number;
  mismatchDistribution: number[];
}

export type Modality = "aso" | "sirna" | "mrna" | "sgrna";
