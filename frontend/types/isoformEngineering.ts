export interface IsoformEngineeringInputs {
  targetSymbol: string;
  isoformGoal: string;
  targetExonLocus: string;
  spliceElementTarget: string;
  stericChemistry: string;
  enforceInFrame: boolean;
  asoLength?: number;
  maxCandidates?: number;
}

// TG07 designs a steric-blocking oligonucleotide, not an mRNA construct.
// The previous shape carried CAI, uridine content, TLR scores and a predicted
// half-life — properties of a delivered mRNA (TG08) that an ASO does not have,
// and which the backend was inventing from a loop index. They are gone.
export interface IsoformOverview {
  targetGene: string;
  geneId: string | null;
  refSeq: string | null;
  transcriptLength: number;
  exonCount: number;
  targetExon: number;
  exonLength: number;
  targetWindow: string;
  windowStart: number;
  windowEnd: number;
  isoformGoal: string;
  spliceElementTarget: string;
  primaryMechanism: string;
  inFrameStatus: string;
  frameNote: string | null;
  spliceSiteStrength: number | null;
}

export interface IsoformCandidate {
  rank: number;
  constructId: string;
  modality: string;
  mechanismChemistry: string;
  sequence: string;
  targetSequence: string;
  transcriptStart: number;
  transcriptEnd: number;
  length: number;
  gcContent: number;
  meltingTempC: number | null;
  selfMfe: number | null;
  targetDuplexDg: number | null;
  targetWindow: string;
  exonNumber: number;
  exonLength: number;
  inFrameStatus: string;
  spliceSiteStrength: number | null;
  /** Null with a stated reason rather than a fabricated value. */
  predictedIsoformYield: null;
  tlrRisk: null;
  notComputed: Record<string, string>;
}

export interface IsoformRanking {
  orderedBy: string;
  caveat: string;
}

export interface IsoformEngineeringResponse {
  status: "OK" | "UNAVAILABLE";
  message?: string;
  overview: IsoformOverview;
  ranking?: IsoformRanking;
  dataProvenance?: Record<string, string>;
  candidates: IsoformCandidate[];
}

export interface DesignOptions {
  isoformGoals: { id: string; label: string; description: string }[];
  targetExonLoci: { id: string; label: string }[];
  spliceElementTargets: { id: string; label: string }[];
  stericChemistries: { id: string; label: string }[];
}
