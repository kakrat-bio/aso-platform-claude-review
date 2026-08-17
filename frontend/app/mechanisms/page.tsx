"use client";

import { useEffect, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import { AlertCircle, ArrowLeft, Loader2, Dna, ChevronDown, Beaker, Download, FlaskConical, FileText, X, Copy, Check } from "lucide-react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import MechanismCard from "@/components/MechanismCard";
import ProteinReplacementDashboard from "@/components/ProteinReplacementDashboard";
import { Card, SectionHeader, FieldLabel, InfoField, Pill } from "@/components/ui";
import { GeneTargetObject } from "@/types/gene";
import {
  MechanismOptions,
  MechanismRankingResponse,
  GeneFeaturesResponse,
  TherapeuticGoalId,
  THERAPEUTIC_GOALS,
} from "@/types/mechanism";
import {
  fetchMechanismOptions,
  fetchGeneFeatures,
  rankGeneSilencingMechanisms,
  rankGeneUpregulationMechanisms,
  rankRnaProcessingMechanisms,
  rankRnaEditingMechanisms,
  rankRnaNeutralizationMechanisms,
  rankTranslationalRegulationMechanisms,
  rankRnaEngineeringMechanisms,
  rankIsoformEngineeringMechanisms,
  getGoalLabel,
} from "@/lib/mechanismApi";
import { fetchRnaEditingClinVarVariants } from "@/lib/rnaEditingApi";
import {
  resolveOrganism,
  organismCapabilities,
  canOptIntoMechanisms,
  organismBannerMessage,
} from "@/lib/organisms";
import OrganismCapabilityBanner from "@/components/OrganismCapabilityBanner";
import { ClinVarVariant } from "@/types/geneSilencing";
import RnaEditingVariantSelector from "@/components/RnaEditingVariantSelector";
import { saveReport } from "@/lib/auth";
import { MechanismFeature } from "@/types/mechanism";
import { parseHgvsC } from "@/lib/hgvsParser";
import {
  ProteinReplacementInputs,
  ProteinReplacementResponse,
  DesignOptions,
  RnaCandidate,
} from "@/types/proteinReplacement";
import { fetchDesignOptions as fetchPrDesignOptions, generateConstructs } from "@/lib/proteinReplacementApi";
import {
  IsoformEngineeringInputs,
  IsoformEngineeringResponse,
  DesignOptions as IsoformDesignOptions,
  IsoformCandidate,
} from "@/types/isoformEngineering";
import { fetchDesignOptions as fetchIeDesignOptions, generateConstructs as generateIeConstructs } from "@/lib/isoformEngineeringApi";
import { useKeyboardShortcut } from "@/hooks/useKeyboardShortcut";

const CONFIRMED_TARGET_KEY = "aso:confirmedTarget";
const SELECTED_MECHANISM_KEY = "aso:selectedMechanism";
const SELECTED_GOAL_KEY = "aso:therapeuticGoal";
const PROJECT_TARGET_TISSUE_KEY = "aso:projectTargetTissue";
const GENE_FEATURES_BACKUP_PREFIX = "aso:geneFeaturesBackup:";

interface StoredGeneFeatures extends GeneFeaturesResponse {
  savedAt?: number;
}

function geneFeaturesBackupKey(organism: string, geneSymbol: string): string {
  return `${GENE_FEATURES_BACKUP_PREFIX}${organism}:${geneSymbol.toLowerCase()}`;
}

// Keep the last-known-good gene feature analysis locally so the Gene Function
// section still renders when the backend / Ensembl site is unreachable.
function saveGeneFeaturesBackup(organism: string, geneSymbol: string, data: GeneFeaturesResponse) {
  try {
    localStorage.setItem(
      geneFeaturesBackupKey(organism, geneSymbol),
      JSON.stringify({ ...data, savedAt: Date.now() }),
    );
  } catch {
    /* storage full or unavailable — ignore */
  }
}

function loadGeneFeaturesBackup(organism: string, geneSymbol: string): StoredGeneFeatures | null {
  try {
    const raw = localStorage.getItem(geneFeaturesBackupKey(organism, geneSymbol));
    if (!raw) return null;
    const data = JSON.parse(raw) as StoredGeneFeatures;
    return data && data.features ? data : null;
  } catch {
    return null;
  }
}

// Map free-text target tissue to deliveryContext dropdown values
function mapTargetTissueToDeliveryContext(tissue: string): string {
  const t = tissue.toLowerCase().trim();
  if (t.includes("liver") || t.includes("hepatic")) return "liver";
  if (t.includes("kidney") || t.includes("renal")) return "kidney";
  if (t.includes("brain") || t.includes("cns") || t.includes("central nervous")) return "cns";
  if (t.includes("muscle") || t.includes("skeletal") || t.includes("myocyte")) return "muscle";
  if (t.includes("heart") || t.includes("cardiac") || t.includes("myocard")) return "heart";
  if (t.includes("lung") || t.includes("pulmonary") || t.includes("respiratory")) return "lung";
  if (t.includes("eye") || t.includes("retina") || t.includes("ocular") || t.includes("vitreous")) return "eye";
  if (t.includes("tumor") || t.includes("cancer") || t.includes("neoplasm") || t.includes("malignan")) return "tumor";
  if (t.includes("blood") || t.includes("bone marrow") || t.includes("hematopoietic") || t.includes("leukemia") || t.includes("lymphoma")) return "blood";
  if (t.includes("skin") || t.includes("dermal") || t.includes("epidermal") || t.includes("cutaneous")) return "skin";
  if (t.includes("pancreas") || t.includes("pancreatic")) return "pancreas";
  if (t.includes("gut") || t.includes("intestine") || t.includes("intestinal") || t.includes("colon") || t.includes("bowel")) return "gut";
  if (t.includes("spinal") || t.includes("cord")) return "spinal cord";
  return "";
}

// Defect type → compatible mechanism IDs (mirrors UPREGULATION_DEFECT_COMPATIBILITY in backend)
const DEFECT_TO_MECHANISMS: Record<string, string[]> = {
  haploinsufficiency: ["A3", "A4", "A6", "A23"],
  poison_exon_inclusion: ["A3"],
  nat_mediated_repression: ["A4"],
  uorf_mediated_repression: ["A5"],
  mirna_mediated_repression: ["A6"],
  rbp_mediated_repression: ["A28"],
  epigenetic_promoter_silencing: ["A23"],
};

// Mechanism ID → feature key in gene features response
const MECHANISM_TO_FEATURE: Record<string, string> = {
  A3: "TANGO",
  A4: "NAT",
  A5: "uORF",
  A6: "miRNA_block",
  A23: "saRNA",
  A28: "RBP_block",
};

interface MechanismCard {
  key: string;
  label: string;
  mechanism: string;
  description: string;
}

const MECHANISM_CARDS: MechanismCard[] = [
  { key: "saRNA", label: "saRNA (Promoter)", mechanism: "A23", description: "Recruits RNA Pol II & AGO2 to boost transcription" },
  { key: "uORF", label: "uORF Blocking", mechanism: "A5", description: "Blocks inhibitory upstream ORFs to enhance translation" },
  { key: "TANGO", label: "Poison Exon (TANGO)", mechanism: "A3", description: "Prevents poison exon inclusion to restore functional mRNA" },
  { key: "NAT", label: "NAT / lncRNA Silencing", mechanism: "A4", description: "Degrades antisense lncRNAs that repress the gene" },
  { key: "miRNA_block", label: "miRNA Site Blocking", mechanism: "A6", description: "Blocks miRNA binding sites on the target mRNA" },
  { key: "RBP_block", label: "RBP Site Blocking", mechanism: "A28", description: "Blocks RNA-binding protein sites to relieve translational repression" },
];

function MechanismAvailabilityCards({
  features,
  selectedDefectType,
  geneSymbol,
}: {
  features: Record<string, MechanismFeature>;
  selectedDefectType: string;
  geneSymbol: string;
}) {
  // Which mechanisms are compatible with the selected defect type?
  const compatibleMechanisms = selectedDefectType
    ? new Set(DEFECT_TO_MECHANISMS[selectedDefectType] || [])
    : null;

  return (
    <div className="rounded-lg border border-[#E5E7EB] bg-slate-50 p-4">
      <div className="flex items-center justify-between mb-3">
        <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500">
          Mechanism Availability for {geneSymbol}
        </p>
        {compatibleMechanisms && (
          <p className="text-[10.5px] text-brand font-medium">
            {compatibleMechanisms.size} mechanism{compatibleMechanisms.size !== 1 ? "s" : ""} match this defect
          </p>
        )}
      </div>
      <div className="grid grid-cols-1 gap-2 md:grid-cols-3">
        {MECHANISM_CARDS.map(({ key, label, mechanism, description }) => {
          const feat = features[key];
          const structurallyAvailable = feat?.available ?? true;
          const isCompatible = compatibleMechanisms ? compatibleMechanisms.has(mechanism) : null;

          // Card state: eligible (green), compatible but structurally unavailable (yellow),
          // incompatible with selected defect (grey), or no defect selected (default)
          let borderColor = "border-[#E5E7EB] bg-white";
          let dotColor = "bg-slate-300";
          let labelColor = "text-slate-700";
          let badge = null;

          if (!selectedDefectType) {
            // No defect selected — show structural availability only
            if (structurallyAvailable) {
              borderColor = "border-emerald-200 bg-emerald-50/50";
              dotColor = "bg-emerald-500";
              labelColor = "text-emerald-800";
            } else {
              borderColor = "border-[#E5E7EB] bg-white opacity-60";
            }
          } else if (isCompatible && structurallyAvailable) {
            // Compatible AND structurally available — highlight as eligible
            borderColor = "border-brand bg-brand/5 ring-1 ring-brand/30";
            dotColor = "bg-brand";
            labelColor = "text-brand";
            badge = (
              <span className="ml-1.5 inline-flex items-center rounded-full bg-brand/10 px-1.5 py-0.5 text-[9px] font-bold text-brand">
                ELIGIBLE
              </span>
            );
          } else if (isCompatible && !structurallyAvailable) {
            // Compatible but gene lacks the feature — show warning
            borderColor = "border-amber-200 bg-amber-50/50";
            dotColor = "bg-amber-400";
            labelColor = "text-amber-700";
            badge = (
              <span className="ml-1.5 inline-flex items-center rounded-full bg-amber-100 px-1.5 py-0.5 text-[9px] font-bold text-amber-600">
                GENE LACKS FEATURE
              </span>
            );
          } else {
            // Not compatible with this defect type
            borderColor = "border-[#E5E7EB] bg-white opacity-50";
            dotColor = "bg-slate-300";
            labelColor = "text-slate-500";
          }

          return (
            <div
              key={key}
              className={`flex items-start gap-2 rounded-md border p-2.5 text-[11.5px] transition-all ${borderColor}`}
            >
              <span className={`mt-0.5 h-2 w-2 shrink-0 rounded-full ${dotColor}`} />
              <div className="flex-1 min-w-0">
                <div className="flex items-center flex-wrap gap-0.5">
                  <p className={`font-semibold ${labelColor}`}>
                    {label}
                    <span className="ml-1 font-mono text-[10px] opacity-60">({mechanism})</span>
                  </p>
                  {badge}
                </div>
                <p className="text-[10.5px] text-slate-500 mt-0.5 leading-snug">
                  {feat?.reason || description}
                </p>
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}

export default function MechanismSelectionPage() {
  const router = useRouter();

  const [gene, setGene] = useState<GeneTargetObject | null>(null);
  const [options, setOptions] = useState<MechanismOptions | null>(null);
  const [optionsError, setOptionsError] = useState<string | null>(null);

  const [selectedGoal, setSelectedGoal] = useState<TherapeuticGoalId | null>(null);

  // TG01 fields
  const [defectType, setDefectType] = useState("");
  const [silencingScope, setSilencingScope] = useState("");

  // TG02 fields
  const [upregDefectType, setUpregDefectType] = useState("");
  const [knownRegulatoryElement, setKnownRegulatoryElement] = useState("");
  const [geneFeatures, setGeneFeatures] = useState<GeneFeaturesResponse | null>(null);
  const [geneFeaturesLoading, setGeneFeaturesLoading] = useState(false);
  const [geneFeaturesNote, setGeneFeaturesNote] = useState<string | null>(null);

  // TG04 fields
  const [spliceDefectType, setSpliceDefectType] = useState("");
  const [targetExon, setTargetExon] = useState("");

  // TG03 fields
  const [editType, setEditType] = useState("");
  const [variantHgvs, setVariantHgvs] = useState("");
  const [enzymeRecruitment, setEnzymeRecruitment] = useState("");
  const [guideLength, setGuideLength] = useState<number>(71);
  const [mismatchPocket, setMismatchPocket] = useState("");
  const [maxBystanderEdits, setMaxBystanderEdits] = useState<number>(0);
  const [splicingDirection, setSplicingDirection] = useState("");
  const [intronSite, setIntronSite] = useState("");
  const [abdLength, setAbdLength] = useState<number>(150);

  const [clinvarVariants, setClinvarVariants] = useState<ClinVarVariant[]>([]);
  const [clinvarLoading, setClinvarLoading] = useState(false);
  const [clinvarError, setClinvarError] = useState<string | null>(null);
  const [useCustomVariant, setUseCustomVariant] = useState(false);
  const [selectedClinvarVariant, setSelectedClinvarVariant] = useState<ClinVarVariant | null>(null);

  // TG05 fields
  const [molecularDefect, setMolecularDefect] = useState("");
  const [neutralizationMode, setNeutralizationMode] = useState("");
  const [repeatUnit, setRepeatUnit] = useState("");
  const [estimatedRepeatCount, setEstimatedRepeatCount] = useState("");
  const [stericChemistry, setStericChemistry] = useState("");
  const [targetRbp, setTargetRbp] = useState("");
  const [oligoLength, setOligoLength] = useState<number | null>(null);
  const [targetGeneType, setTargetGeneType] = useState("");

  // TG06 fields
  const [translationalGoal, setTranslationalGoal] = useState("");
  const [targetElement, setTargetElement] = useState("");
  const [translationStericChemistry, setTranslationStericChemistry] = useState("");
  const [translationTargetRbp, setTranslationTargetRbp] = useState("");
  const [translationOligoLength, setTranslationOligoLength] = useState<number | null>(null);

  // TG07 fields
  const [ieOptions, setIeOptions] = useState<IsoformDesignOptions | null>(null);
  const [ieOptionsError, setIeOptionsError] = useState<string | null>(null);
  const [ieTargetSymbol, setIeTargetSymbol] = useState("");
  const [ieIsoformGoal, setIeIsoformGoal] = useState("");
  const [ieTargetExonLocus, setIeTargetExonLocus] = useState("");
  const [ieSpliceElementTarget, setIeSpliceElementTarget] = useState("");
  const [ieStericChemistry, setIeStericChemistry] = useState("");
  const [ieEnforceInFrame, setIeEnforceInFrame] = useState(true);
  const [ieResults, setIeResults] = useState<IsoformEngineeringResponse | null>(null);
  const [ieLoading, setIeLoading] = useState(false);
  const [ieError, setIeError] = useState<string | null>(null);
  const [selectedIeCandidate, setSelectedIeCandidate] = useState<IsoformCandidate | null>(null);
  const [ieCopied, setIeCopied] = useState(false);

  // TG08 fields
  const [prOptions, setPrOptions] = useState<DesignOptions | null>(null);
  const [prOptionsError, setPrOptionsError] = useState<string | null>(null);
  const [targetSymbol, setTargetSymbol] = useState("");
  const [rnaModality, setRnaModality] = useState("");
  const [codonStrategy, setCodonStrategy] = useState("");
  const [utrPair, setUtrPair] = useState("");
  const [iresSelection, setIresSelection] = useState("");
  const [nucleotideModification, setNucleotideModification] = useState("");
  const [prResults, setPrResults] = useState<ProteinReplacementResponse | null>(null);
  const [prLoading, setPrLoading] = useState(false);
  const [prError, setPrError] = useState<string | null>(null);
  const [selectedPrCandidate, setSelectedPrCandidate] = useState<RnaCandidate | null>(null);
  const [prCopied, setPrCopied] = useState(false);

  // TG09 fields
  const [tg09StructuralClass, setTg09StructuralClass] = useState("");
  const [tg09TargetType, setTg09TargetType] = useState("");
  const [tg09Scaffold, setTg09Scaffold] = useState("");
  const [tg09ChemStabilization, setTg09ChemStabilization] = useState("");
  const [tg09KdGoal, setTg09KdGoal] = useState("");

  // Shared fields
  const [deliveryContext, setDeliveryContext] = useState("");
  const [knownVariant, setKnownVariant] = useState("");

  const [ranking, setRanking] = useState<MechanismRankingResponse | null>(null);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const rankingResultsRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (ranking) {
      rankingResultsRef.current?.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  }, [ranking]);

  useEffect(() => {
    const stored = sessionStorage.getItem(CONFIRMED_TARGET_KEY);
    if (stored) {
      try {
        setGene(JSON.parse(stored));
      } catch {
        setGene(null);
      }
    }

    const savedGoal = sessionStorage.getItem(SELECTED_GOAL_KEY) as TherapeuticGoalId | null;
    if (savedGoal && THERAPEUTIC_GOALS.some((g) => g.id === savedGoal)) {
      setSelectedGoal(savedGoal);
    }

    // Load project target tissue and map to deliveryContext
    const projectTissue = sessionStorage.getItem(PROJECT_TARGET_TISSUE_KEY);
    if (projectTissue) {
      const mapped = mapTargetTissueToDeliveryContext(projectTissue);
      if (mapped) {
        setDeliveryContext(mapped);
      }
    }

    fetchMechanismOptions()
      .then(setOptions)
      .catch((e) => setOptionsError(e instanceof Error ? e.message : "Failed to load options."));
  }, []);

  // Fetch gene structural features when TG02 is selected
  useEffect(() => {
    if (selectedGoal !== "TG02" || !gene) {
      setGeneFeatures(null);
      setGeneFeaturesNote(null);
      return;
    }

    const organism = gene.organism || "homo_sapiens";
    setGeneFeaturesLoading(true);
    fetchGeneFeatures({
      geneSymbol: gene.geneSymbol,
      organism,
      ensemblId: gene.geneId,
      tissueTpm: gene.tissueTpm,
      exonCount: gene.exonCount,
      totalTranscripts: gene.totalTranscripts,
      geneType: gene.geneType || undefined,
    })
      .then((features) => {
        if (features.source !== "backup") {
          saveGeneFeaturesBackup(organism, gene.geneSymbol, features);
        }
        setGeneFeatures(features);
        setGeneFeaturesNote(
          features.source === "backup"
            ? `The Ensembl site is unreachable right now — showing the last saved analysis for ${gene.geneSymbol}.`
            : features.source === "fallback"
              ? `Could not verify ${gene.geneSymbol} structure from Ensembl — showing a conservative estimate that requires experimental validation.`
              : null,
        );
      })
      .catch(() => {
        // Backend / Ensembl entirely unreachable — replay the local backup.
        const backup = loadGeneFeaturesBackup(organism, gene.geneSymbol);
        if (backup) {
          setGeneFeatures({
            ...backup,
            source: "backup",
            backupTimestamp: backup.savedAt ?? Date.now(),
          });
          setGeneFeaturesNote(
            `The Ensembl site is unreachable right now — showing the last saved analysis for ${gene.geneSymbol}.`,
          );
        } else {
          setGeneFeatures(null);
          setGeneFeaturesNote(
            "Could not load gene features — structure-dependent mechanisms will be treated as potentially applicable.",
          );
        }
      })
      .finally(() => setGeneFeaturesLoading(false));
  }, [selectedGoal, gene]);

  // Fetch protein replacement design options when TG08 is selected
  useEffect(() => {
    if (selectedGoal !== "TG08") {
      setPrOptions(null);
      setPrOptionsError(null);
      return;
    }
    fetchPrDesignOptions()
      .then(setPrOptions)
      .catch((e) => setPrOptionsError(e instanceof Error ? e.message : "Failed to load options."));
  }, [selectedGoal]);

  // Fetch ClinVar variants for TG03 RNA editing
  useEffect(() => {
    if (selectedGoal !== "TG03" || !gene?.geneId) {
      setClinvarVariants([]);
      setClinvarLoading(false);
      setClinvarError(null);
      return;
    }
    setClinvarLoading(true);
    setClinvarError(null);
    fetchRnaEditingClinVarVariants(gene.geneId)
      .then(setClinvarVariants)
      .catch((e) => {
        setClinvarError(e instanceof Error ? e.message : "Failed to load ClinVar variants.");
        setClinvarVariants([]);
      })
      .finally(() => setClinvarLoading(false));
  }, [selectedGoal, gene?.geneId]);

  // Fetch isoform engineering design options when TG07 is selected
  useEffect(() => {
    if (selectedGoal !== "TG07") {
      setIeOptions(null);
      setIeOptionsError(null);
      return;
    }
    fetchIeDesignOptions()
      .then(setIeOptions)
      .catch((e) => setIeOptionsError(e instanceof Error ? e.message : "Failed to load options."));
  }, [selectedGoal]);

  // Pre-fill TG08 target protein replacement symbol from confirmed target when empty
  useEffect(() => {
    if (selectedGoal === "TG08" && gene?.geneSymbol && !targetSymbol) {
      setTargetSymbol(gene.geneSymbol);
    }
  }, [selectedGoal, gene, targetSymbol]);

  // Pre-fill TG07 target gene symbol from confirmed target when empty
  useEffect(() => {
    if (selectedGoal === "TG07" && gene?.geneSymbol && !ieTargetSymbol) {
      setIeTargetSymbol(gene.geneSymbol);
    }
  }, [selectedGoal, gene, ieTargetSymbol]);

  useKeyboardShortcut("enter", handleRank, { ctrl: true, meta: true });

  const goalIds = THERAPEUTIC_GOALS.map((g) => g.id);
  const inInput = () => {
    const el = document.activeElement;
    return el && (el.tagName === "INPUT" || el.tagName === "TEXTAREA" || el.tagName === "SELECT");
  };
  useKeyboardShortcut("1", () => { if (!inInput()) goalIds[0] && handleSelectGoal(goalIds[0]); });
  useKeyboardShortcut("2", () => { if (!inInput()) goalIds[1] && handleSelectGoal(goalIds[1]); });
  useKeyboardShortcut("3", () => { if (!inInput()) goalIds[2] && handleSelectGoal(goalIds[2]); });
  useKeyboardShortcut("4", () => { if (!inInput()) goalIds[3] && handleSelectGoal(goalIds[3]); });
  useKeyboardShortcut("5", () => { if (!inInput()) goalIds[4] && handleSelectGoal(goalIds[4]); });
  useKeyboardShortcut("6", () => { if (!inInput()) goalIds[5] && handleSelectGoal(goalIds[5]); });
  useKeyboardShortcut("7", () => { if (!inInput()) goalIds[6] && handleSelectGoal(goalIds[6]); });
  useKeyboardShortcut("8", () => { if (!inInput()) goalIds[7] && handleSelectGoal(goalIds[7]); });
  useKeyboardShortcut("9", () => { if (!inInput()) goalIds[8] && handleSelectGoal(goalIds[8]); });

  function handleSelectGoal(goalId: TherapeuticGoalId) {
    setSelectedGoal(goalId);
    sessionStorage.setItem(SELECTED_GOAL_KEY, goalId);
    setRanking(null);
    setSelectedId(null);
    setDefectType("");
    setSilencingScope("");
    setUpregDefectType("");
    setKnownRegulatoryElement("");
    setSpliceDefectType("");
    setTargetExon("");
    setEditType("");
    setVariantHgvs("");
    setEnzymeRecruitment("");
    setGuideLength(71);
    setMismatchPocket("");
    setMaxBystanderEdits(0);
    setSplicingDirection("");
    setIntronSite("");
    setAbdLength(150);
    setClinvarVariants([]);
    setClinvarLoading(false);
    setClinvarError(null);
    setUseCustomVariant(false);
    setSelectedClinvarVariant(null);
    setMolecularDefect("");
    setNeutralizationMode("");
    setRepeatUnit("");
    setEstimatedRepeatCount("");
    setStericChemistry("");
    setTargetRbp("");
    setOligoLength(null);
    setTargetGeneType("");
    setTranslationalGoal("");
    setTargetElement("");
    setTranslationStericChemistry("");
    setTranslationTargetRbp("");
    setTranslationOligoLength(null);
    setIeIsoformGoal("");
    setIeTargetExonLocus("");
    setIeSpliceElementTarget("");
    setIeStericChemistry("");
    setIeEnforceInFrame(true);
    setIeResults(null);
    setSelectedIeCandidate(null);
    setIeError(null);
    setTargetSymbol("");
    setRnaModality("");
    setCodonStrategy("");
    setUtrPair("");
    setIresSelection("");
    setNucleotideModification("");
    setPrResults(null);
    setSelectedPrCandidate(null);
    setPrError(null);
    setTg09StructuralClass("");
    setTg09TargetType("");
    setTg09Scaffold("");
    setTg09ChemStabilization("");
    setTg09KdGoal("");
    setDeliveryContext("");
    setKnownVariant("");
  }

  function clearRanking() {
    setRanking(null);
    setSelectedId(null);
  }

  async function handleRank() {
    if (!gene || !selectedGoal) return;
    setLoading(true);
    setError(null);
    setSelectedId(null);

    try {
      let result: MechanismRankingResponse;

      if (selectedGoal === "TG01") {
        if (!defectType || !silencingScope) return;
        result = await rankGeneSilencingMechanisms({
          geneSymbol: gene.geneSymbol,
          defectType,
          silencingScope,
          deliveryContext,
          knownVariant,
        });
      } else if (selectedGoal === "TG02") {
        if (!upregDefectType) return;
        result = await rankGeneUpregulationMechanisms({
          geneSymbol: gene.geneSymbol,
          defectType: upregDefectType,
          deliveryContext,
          knownRegulatoryElement,
          geneFeatures: geneFeatures as unknown as Record<string, unknown> | null,
        });
      } else if (selectedGoal === "TG04") {
        if (!spliceDefectType) return;
        result = await rankRnaProcessingMechanisms({
          geneSymbol: gene.geneSymbol,
          spliceDefectType,
          targetExon,
          deliveryContext,
          knownVariant,
        });
      } else if (selectedGoal === "TG03") {
        if (!editType || !variantHgvs) return;
        result = await rankRnaEditingMechanisms({
          geneSymbol: gene.geneSymbol,
          editType,
          variantHgvs,
          enzymeRecruitment,
          deliveryContext,
          guideLength,
          mismatchPocket,
          maxBystanderEdits,
          splicingDirection,
          intronSite,
          abdLength,
          exonCount: gene.exonCount,
          intronCount: gene.intronCount,
          totalTranscripts: gene.totalTranscripts,
        });
  } else if (selectedGoal === "TG05") {
        if (!molecularDefect || !neutralizationMode) return;
        result = await rankRnaNeutralizationMechanisms({
          geneSymbol: gene.geneSymbol,
          molecularDefect,
          neutralizationMode,
          repeatUnit: repeatUnit || undefined,
          estimatedRepeatCount: estimatedRepeatCount || undefined,
          stericChemistry: stericChemistry || undefined,
          targetRbp: targetRbp || undefined,
          oligoLength: oligoLength ?? undefined,
          deliveryContext,
          targetGeneType: targetGeneType || undefined,
        });
      } else if (selectedGoal === "TG06") {
        if (!translationalGoal || !targetElement) return;
        result = await rankTranslationalRegulationMechanisms({
          geneSymbol: gene.geneSymbol,
          translationalGoal,
          targetElement,
          stericChemistry: translationStericChemistry || undefined,
          targetRbp: translationTargetRbp || undefined,
          oligoLength: translationOligoLength ?? undefined,
          deliveryContext,
        });
      } else if (selectedGoal === "TG09") {
        if (!tg09StructuralClass || !tg09TargetType || !tg09Scaffold || !tg09ChemStabilization || !tg09KdGoal) return;
        result = await rankRnaEngineeringMechanisms({
          geneSymbol: gene.geneSymbol,
          structuralClass: tg09StructuralClass,
          targetType: tg09TargetType,
          scaffold: tg09Scaffold,
          chemStabilization: tg09ChemStabilization,
          kdGoal: tg09KdGoal,
          deliveryContext: deliveryContext || undefined,
        });
      } else if (selectedGoal === "TG07") {
        if (!ieIsoformGoal) return;
        result = await rankIsoformEngineeringMechanisms({
          geneSymbol: gene.geneSymbol,
          isoformGoal: ieIsoformGoal,
          targetExonLocus: ieTargetExonLocus || undefined,
          spliceElementTarget: ieSpliceElementTarget || undefined,
          stericChemistry: ieStericChemistry || undefined,
          deliveryContext: deliveryContext || undefined,
        });
      } else {
        return;
      }

      setRanking(result);
      saveReport({
        step: "mechanism",
        title: `Mechanism Analysis: ${gene.geneSymbol} (${selectedGoal || "N/A"})`,
        geneSymbol: gene.geneSymbol,
        disease: gene.disease || "",
        summary: `Ranked ${result.results.length} mechanisms for ${gene.geneSymbol}. Top: ${result.results[0]?.name || "N/A"}.`,
        data: { goal: selectedGoal, topMechanisms: result.results.slice(0, 5).map((m: any) => ({ id: m.id, name: m.name, score: m.score })) },
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Something went wrong.");
    } finally {
      setLoading(false);
    }
  }

  function handleSelectMechanism(id: string) {
    setSelectedId(id);
    const mechanism = ranking?.results.find((r) => r.id === id);
    if (mechanism) {
      const parsedVariant =
        selectedGoal === "TG01" && knownVariant.trim()
          ? parseHgvsC(knownVariant)
          : null;
      sessionStorage.setItem(
        SELECTED_MECHANISM_KEY,
        JSON.stringify({
          geneSymbol: gene?.geneSymbol,
          mechanism,
          silencingScope: selectedGoal === "TG01" ? silencingScope : null,
          defectType: selectedGoal === "TG01" ? defectType : null,
          therapeuticGoal: selectedGoal,
          knownVariant,
          ...(selectedGoal === "TG01" ? { parsedVariant } : {}),
          ...(selectedGoal === "TG02"
            ? {
                upregDefectType,
                knownRegulatoryElement,
                deliveryContext,
              }
            : {}),
          ...(selectedGoal === "TG04"
            ? {
                spliceDefectType,
                targetExon,
                deliveryContext,
              }
            : {}),
          ...(selectedGoal === "TG03"
            ? {
                editType,
                variantHgvs,
                enzymeRecruitment,
                deliveryContext,
                guideLength,
                mismatchPocket,
                maxBystanderEdits,
                splicingDirection,
                intronSite,
                abdLength,
              }
            : {}),
          ...(selectedGoal === "TG05"
            ? {
                molecularDefect,
                neutralizationMode,
                repeatUnit,
                estimatedRepeatCount,
                stericChemistry,
                targetRbp,
                oligoLength,
           deliveryContext,
           targetGeneType,
         }
             : {}),
          ...(selectedGoal === "TG06"
            ? {
                translationalGoal,
                targetElement,
                stericChemistry: translationStericChemistry,
                targetRbp: translationTargetRbp,
                oligoLength: translationOligoLength ?? null,
                deliveryContext,
              }
            : {}),
           ...(selectedGoal === "TG09"
             ? {
                 tg09StructuralClass,
                 tg09TargetType,
                 tg09Scaffold,
                 tg09ChemStabilization,
                 tg09KdGoal,
                 deliveryContext,
               }
             : {}),
           ...(selectedGoal === "TG07"
             ? {
                 ieTargetSymbol,
                 ieIsoformGoal,
                 ieTargetExonLocus,
                 ieSpliceElementTarget,
                 ieStericChemistry,
                 ieEnforceInFrame,
                 deliveryContext,
               }
             : {}),
         })
      );
    }
  }

  function isRankDisabled(): boolean {
    if (!gene || !selectedGoal || loading) return true;
    if (selectedGoal === "TG01") return !defectType || !silencingScope;
    if (selectedGoal === "TG02") return !upregDefectType;
    if (selectedGoal === "TG04") return !spliceDefectType;
    if (selectedGoal === "TG03") return !editType || !variantHgvs.trim();
    if (selectedGoal === "TG05") return !molecularDefect || !neutralizationMode;
    if (selectedGoal === "TG06") return !translationalGoal || !targetElement;
    if (selectedGoal === "TG07") return !ieIsoformGoal || !ieTargetExonLocus || !ieSpliceElementTarget || !ieStericChemistry;
    if (selectedGoal === "TG09") return !tg09StructuralClass || !tg09TargetType || !tg09Scaffold || !tg09ChemStabilization || !tg09KdGoal;
    if (selectedGoal === "TG08") return true;
    return true;
  }

  const mechanismOrganism = resolveOrganism(gene?.organism);
  const mechanismsEnabledHere = mechanismOrganism
    ? organismCapabilities(mechanismOrganism).mechanisms
    : true;

  // Reaching here for a Tier 4-6 organism means a direct URL or stale state:
  // Confirm & Proceed is not offered without the opt-in. Say so plainly and
  // offer the way back rather than redirecting, which would hide what
  // happened.
  if (gene && mechanismOrganism && !mechanismsEnabledHere) {
    return (
      <div className="flex min-h-screen bg-[#F8FAFC]">
        <Sidebar />
        <div className="flex min-h-screen flex-1 flex-col">
          <Topbar />
          <main className="flex flex-1 items-center justify-center px-6">
            <Card className="max-w-lg p-8">
              <AlertCircle className="mx-auto h-8 w-8 text-amber-400" />
              <p className="mt-3 text-center text-[14px] font-medium text-slate-700">
                Mechanism analysis is not enabled for {mechanismOrganism.commonName}
              </p>
              <p className="mt-2 text-[13px] leading-relaxed text-slate-600">
                {organismBannerMessage(mechanismOrganism.id)}
              </p>
              {canOptIntoMechanisms(mechanismOrganism) && (
                <p className="mt-3 rounded border border-amber-200 bg-amber-50 px-3 py-2 text-[12px] leading-snug text-amber-900">
                  To proceed anyway, go back to Basic Information and tick
                  &ldquo;Enable mechanism analysis anyway&rdquo; before
                  confirming the target.
                </p>
              )}
              <button
                onClick={() => router.push("/")}
                className="mt-4 inline-flex items-center gap-1.5 rounded-lg bg-brand px-4 py-2 text-[13px] font-medium text-white"
              >
                <ArrowLeft className="h-3.5 w-3.5" />
                Back to Basic Information
              </button>
            </Card>
          </main>
        </div>
      </div>
    );
  }

  if (!gene) {
    return (
      <div className="flex min-h-screen bg-[#F8FAFC]">
        <Sidebar />
        <div className="flex min-h-screen flex-1 flex-col">
          <Topbar />
          <main className="flex flex-1 items-center justify-center px-6">
            <Card className="max-w-md p-8 text-center">
              <AlertCircle className="mx-auto h-8 w-8 text-slate-300" />
              <p className="mt-3 text-[14px] font-medium text-slate-700">No confirmed target found</p>
              <p className="mt-1 text-[13px] text-slate-500">
                Go back to Basic Information, load a gene, and hit Confirm &amp; Proceed first.
              </p>
              <button
                onClick={() => router.push("/")}
                className="mt-4 inline-flex items-center gap-1.5 rounded-lg bg-brand px-4 py-2 text-[13px] font-medium text-white"
              >
                <ArrowLeft className="h-3.5 w-3.5" />
                Back to Basic Information
              </button>
            </Card>
          </main>
        </div>
      </div>
    );
  }

  return (
    <div className="flex min-h-screen bg-[#F8FAFC]">
      <Sidebar />
      <div className="flex min-h-screen flex-1 flex-col">
        <Topbar />
        <main className="flex-1 space-y-5 px-6 py-6">
          {/* Confirmed target banner */}
          <Card className="flex items-center gap-3 px-5 py-3">
            <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-gradient-to-br from-indigo-50 to-blue-50">
              <Dna className="h-4.5 w-4.5 text-indigo-500" />
            </span>
            <div>
              <p className="text-[13px] font-semibold text-slate-800">
                {gene.geneSymbol} <span className="font-normal text-slate-400">· {gene.organism}</span>
              </p>
              <p className="text-[12px] text-slate-500">
                {gene.geneName ?? "—"}
                {gene.diseaseName ? ` · ${gene.diseaseName}` : ""}
              </p>
            </div>
            <button
              onClick={() => router.push("/")}
              className="ml-auto text-[12.5px] font-medium text-brand hover:underline"
            >
              Change target
            </button>
          </Card>

          {/* Step 2: Therapeutic Goal Selection */}
          <Card>
            <SectionHeader step="2" title="Select Therapeutic Goal" />
            <div className="grid grid-cols-1 gap-3 px-6 pb-5 sm:grid-cols-2 lg:grid-cols-3">
              {THERAPEUTIC_GOALS.map((goal) => (
                <button
                  key={goal.id}
                  onClick={() => handleSelectGoal(goal.id)}
                  className={`rounded-lg border p-4 text-left transition-all duration-200 ${
                    selectedGoal === goal.id
                      ? "border-brand bg-brand/5 ring-1 ring-brand shadow-sm"
                      : "border-[#E5E7EB] bg-white hover:border-slate-300 hover:bg-slate-50 hover:shadow-sm hover:-translate-y-0.5"
                  }`}
                >
                  <div className="flex items-center gap-2">
                    <span className="text-[11px] font-mono text-slate-400">{goal.id}</span>
                    <h3 className="text-[14px] font-semibold text-slate-800">{goal.name}</h3>
                  </div>
                  <p className="mt-1 text-[12.5px] text-slate-500 leading-snug">{goal.description}</p>
                </button>
              ))}
            </div>
          </Card>

          {/* Step 3: Mechanism inputs */}
          {selectedGoal && (
            <Card>
              <SectionHeader step="3" title="Mechanism Selection" />
              {optionsError && (
                <p className="px-6 pb-2 text-[12.5px] text-red-600">{optionsError}</p>
              )}

              {selectedGoal === "TG01" && (
                <div className="grid grid-cols-1 gap-4 px-6 pb-4 md:grid-cols-3">
                  <div>
                    <FieldLabel hint="What kind of molecular defect are you trying to counteract?">
                      Molecular Defect Type <span className="text-red-500">*</span>
                    </FieldLabel>
                    <select
                      value={defectType}
                      onChange={(e) => {
                        setDefectType(e.target.value);
                        clearRanking();
                      }}
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    >
                      <option value="">Select defect type</option>
                      {options?.geneSilencing.defectTypes.map((o) => (
                        <option key={o.id} value={o.id}>
                          {o.label}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div>
                    <FieldLabel hint="Do you need to silence the whole transcript, or spare the wild-type allele?">
                      Silencing Scope <span className="text-red-500">*</span>
                    </FieldLabel>
                    <select
                      value={silencingScope}
                      onChange={(e) => {
                        setSilencingScope(e.target.value);
                        clearRanking();
                      }}
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    >
                      <option value="">Select scope</option>
                      {options?.geneSilencing.silencingScopes.map((o) => (
                        <option key={o.id} value={o.id}>
                          {o.label}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div>
                    <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                      Delivery / Tissue Context
                    </FieldLabel>
                    <select
                      value={deliveryContext}
                      onChange={(e) => {
                        setDeliveryContext(e.target.value);
                        clearRanking();
                      }}
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    >
                      <option value="">Not specified</option>
                      {options?.deliveryContexts.map((o) => (
                        <option key={o.id} value={o.id}>
                          {o.label}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div className="md:col-span-3">
                    <FieldLabel hint="Optional — a ClinVar ID, HGVS notation, or free-text description">
                      Known Variant (optional)
                    </FieldLabel>
                    <input
                      value={knownVariant}
                      onChange={(e) => {
                        setKnownVariant(e.target.value);
                        clearRanking();
                      }}
                      placeholder="e.g. c.1521_1523delCTT / p.Phe508del"
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    />
                    {knownVariant.trim() && (() => {
                      const parsed = parseHgvsC(knownVariant);
                      if (parsed.parsed) {
                        return (
                          <p className="mt-1.5 rounded-md bg-emerald-50 px-3 py-1.5 text-[11.5px] text-emerald-700">
                            Parsed position: CDS index {parsed.cdsStart}
                            {parsed.cdsStart !== parsed.cdsEnd ? `–${parsed.cdsEnd}` : ""} ({parsed.type})
                            {parsed.length != null && parsed.length > 1 ? `, ${parsed.length} bp` : ""}
                          </p>
                        );
                      }
                      return (
                        <p className="mt-1.5 rounded-md bg-amber-50 px-3 py-1.5 text-[11.5px] text-amber-700">
                          {parsed.reason}
                        </p>
                      );
                    })()}
                    {silencingScope === "allele_specific" && !knownVariant.trim() && (
                      <p className="mt-1.5 rounded-md bg-amber-50 px-3 py-1.5 text-[11.5px] text-amber-700">
                        Allele-specific ranking without a variant position can only confirm a mechanism
                        supports this approach in principle — it can't verify a specific candidate will
                        discriminate mutant from wild-type.
                      </p>
                    )}
                  </div>
                </div>
              )}

              {selectedGoal === "TG02" && (
                <div className="space-y-4 px-6 pb-4">
                  {/* Gene feature analysis loading */}
                  {geneFeaturesLoading && (
                    <div className="flex items-center gap-2 rounded-lg border border-blue-200 bg-blue-50 px-4 py-3 text-[12.5px] text-blue-700">
                      <Loader2 className="h-4 w-4 animate-spin" />
                      Analyzing gene structural features...
                    </div>
                  )}

                  {/* Defect type selection — this drives which mechanisms are eligible */}
                  <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                    <div>
                      <FieldLabel hint="Select the molecular defect to see which upregulation mechanisms apply">
                        Molecular Defect Type <span className="text-red-500">*</span>
                      </FieldLabel>
                      <select
                        value={upregDefectType}
                        onChange={(e) => {
                          setUpregDefectType(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Select defect type</option>
                        {options?.geneUpregulation.defectTypes.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>

                    <div>
                      <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                        Delivery / Tissue Context
                      </FieldLabel>
                      <select
                        value={deliveryContext}
                        onChange={(e) => {
                          setDeliveryContext(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Not specified</option>
                        {options?.deliveryContexts.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>

                    <div className="md:col-span-3">
                      <FieldLabel hint="Optional — known poison exon, NAT, uORF, or miRNA binding site for this gene">
                        Known Regulatory Element (optional)
                      </FieldLabel>
                      <input
                        value={knownRegulatoryElement}
                        onChange={(e) => {
                          setKnownRegulatoryElement(e.target.value);
                          clearRanking();
                        }}
                        placeholder="e.g. BDNF-AS antisense transcript / chr11:53886643 poison exon"
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      />
                    </div>
                  </div>

                  {/* Mechanism availability — filtered by selected defect type */}
                  {geneFeaturesNote && !geneFeaturesLoading && (
                    <div className="flex items-start gap-2 rounded-lg border border-amber-200 bg-amber-50 px-4 py-3 text-[12.5px] text-amber-800">
                      <AlertCircle className="mt-0.5 h-4 w-4 shrink-0" />
                      <span>{geneFeaturesNote}</span>
                    </div>
                  )}

                  {geneFeatures && !geneFeaturesLoading && (
                    <MechanismAvailabilityCards
                      features={geneFeatures.features}
                      selectedDefectType={upregDefectType}
                      geneSymbol={gene?.geneSymbol || ""}
                    />
                  )}

                  {/* Overexpression warnings */}
                  {geneFeatures?.warnings.map((w, i) => (
                    <div
                      key={i}
                      className={`rounded-lg border px-4 py-3 text-[12.5px] ${
                        w.severity === "high"
                          ? "border-amber-300 bg-amber-50 text-amber-800"
                          : "border-yellow-200 bg-yellow-50 text-yellow-700"
                      }`}
                    >
                      {w.message}
                    </div>
                  ))}

                  {/* Hint when no defect type selected */}
                  {!upregDefectType && (
                    <p className="text-[12px] text-slate-500 italic">
                      Select a Molecular Defect Type above to see which mechanisms are eligible for this gene.
                    </p>
                  )}
                </div>
              )}

              {selectedGoal === "TG04" && (
                <div className="grid grid-cols-1 gap-4 px-6 pb-4 md:grid-cols-3">
                  <div>
                    <FieldLabel hint="What type of RNA processing defect are you targeting?">
                      Splice Defect Type <span className="text-red-500">*</span>
                    </FieldLabel>
                    <select
                      value={spliceDefectType}
                      onChange={(e) => {
                        setSpliceDefectType(e.target.value);
                        clearRanking();
                      }}
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    >
                      <option value="">Select splice defect type</option>
                      {options?.rnaProcessing.spliceDefectTypes.map((o) => (
                        <option key={o.id} value={o.id}>
                          {o.label}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div>
                    <FieldLabel hint="Optional — the exon number affected by the mutation (e.g. exon 51 for DMD)">
                      Target Exon (optional)
                    </FieldLabel>
                    <input
                      value={targetExon}
                      onChange={(e) => {
                        setTargetExon(e.target.value);
                        clearRanking();
                      }}
                      placeholder="e.g. exon 51"
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    />
                  </div>

                  <div>
                    <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                      Delivery / Tissue Context
                    </FieldLabel>
                    <select
                      value={deliveryContext}
                      onChange={(e) => {
                        setDeliveryContext(e.target.value);
                        clearRanking();
                      }}
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    >
                      <option value="">Not specified</option>
                      {options?.deliveryContexts.map((o) => (
                        <option key={o.id} value={o.id}>
                          {o.label}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div className="md:col-span-3">
                    <FieldLabel hint="Optional — a ClinVar ID, HGVS notation, or free-text description of the known variant">
                      Known Variant (optional)
                    </FieldLabel>
                    <input
                      value={knownVariant}
                      onChange={(e) => {
                        setKnownVariant(e.target.value);
                        clearRanking();
                      }}
                      placeholder="e.g. c.1521_1523delCTT / p.Phe508del"
                      className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                    />
                  </div>
                </div>
              )}

              {selectedGoal === "TG03" && (
                <div className="space-y-4 px-6 pb-4">
                  <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                    <div>
                      <FieldLabel hint="Which RNA editing modality do you want to use to repair the transcript?">
                        Editing Modality <span className="text-red-500">*</span>
                      </FieldLabel>
                      <select
                        value={editType}
                        onChange={(e) => {
                          setEditType(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Select editing modality</option>
                        {options?.rnaEditing.editTypes.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>

                    <div className="md:col-span-2">
                      <FieldLabel hint="General delivery/tissue context — a soft tie-breaker and used to gauge endogenous enzyme expression">
                        Delivery / Tissue Context
                      </FieldLabel>
                      <select
                        value={deliveryContext}
                        onChange={(e) => {
                          setDeliveryContext(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Not specified</option>
                        {options?.deliveryContexts.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>
                  </div>

                  <div>
                    <FieldLabel hint="The exact variant to correct, in HGVS notation">
                      Target Variant (HGVS) <span className="text-red-500">*</span>
                    </FieldLabel>
                    {clinvarLoading && (
                      <p className="mt-1.5 text-[11.5px] text-slate-400">Loading ClinVar variants…</p>
                    )}
                    {clinvarError && (
                      <p className="mt-1.5 text-[11.5px] text-red-600">{clinvarError}</p>
                    )}
                    <RnaEditingVariantSelector
                      variants={clinvarVariants}
                      editType={editType}
                      selectedVariant={selectedClinvarVariant}
                      customVariant={variantHgvs}
                      useCustom={useCustomVariant}
                      onSelectVariant={(v) => {
                        setSelectedClinvarVariant(v);
                        setVariantHgvs(v?.hgvsc || "");
                        clearRanking();
                      }}
                      onCustomVariantChange={(v) => {
                        setVariantHgvs(v);
                        clearRanking();
                      }}
                      onToggleCustom={(useCustom) => {
                        setUseCustomVariant(useCustom);
                        if (!useCustom) {
                          setSelectedClinvarVariant(null);
                        }
                        clearRanking();
                      }}
                    />
                  </div>

                  {editType !== "trans_splicing" && (
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-4">
                      <div>
                        <FieldLabel hint="Which enzyme machinery should execute the edit? ADAR2 is CNS-enriched; APOBEC1 is gut/liver-enriched">
                          Recruiting Enzyme Machinery
                        </FieldLabel>
                        <select
                          value={enzymeRecruitment}
                          onChange={(e) => {
                            setEnzymeRecruitment(e.target.value);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.rnaEditing.enzymeRecruitment.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Length of the guide RNA / editase oligo (30–120 nt)">
                          Guide RNA Length (nt)
                        </FieldLabel>
                        <input
                          type="number"
                          min={30}
                          max={120}
                          value={guideLength}
                          onChange={(e) => {
                            setGuideLength(parseInt(e.target.value, 10));
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      {editType === "a_to_i" && (
                        <div>
                          <FieldLabel hint="The orphan base opposite the target adenosine on the guide RNA — A-C mismatch reports the highest editing efficiency">
                            Opposing Base Mismatch
                          </FieldLabel>
                          <select
                            value={mismatchPocket}
                            onChange={(e) => {
                              setMismatchPocket(e.target.value);
                              clearRanking();
                            }}
                            className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                          >
                            <option value="">Not specified</option>
                            {options?.rnaEditing.mismatchPocket.map((o) => (
                              <option key={o.id} value={o.id}>
                                {o.label}
                              </option>
                            ))}
                          </select>
                  </div>
                )}

                       <div>
                        <FieldLabel hint="Maximum allowable bystander A/C edits within ±20 bp of the target site">
                          Max Bystander Edits
                        </FieldLabel>
                        <input
                          type="number"
                          min={0}
                          value={maxBystanderEdits}
                          onChange={(e) => {
                            setMaxBystanderEdits(parseInt(e.target.value, 10) || 0);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>
                    </div>
                  )}

                  {editType === "trans_splicing" && (
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-4">
                      <div>
                        <FieldLabel hint="Which end of the transcript are you replacing?">
                          Splicing Direction
                        </FieldLabel>
                        <select
                          value={splicingDirection}
                          onChange={(e) => {
                            setSplicingDirection(e.target.value);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.rnaEditing.splicingDirections.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="The intron junction adjacent to the mutated region (e.g. Intron 12 Acceptor Junction)">
                          Intron Acceptor / Donor Site
                        </FieldLabel>
                        <select
                          value={intronSite}
                          onChange={(e) => {
                            setIntronSite(e.target.value);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.rnaEditing.intronSites.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Length of the antisense binding domain complementary to the target intron (100–300 bp)">
                          ABD Length (bp)
                        </FieldLabel>
                        <input
                          type="number"
                          min={100}
                          max={300}
                          value={abdLength}
                          onChange={(e) => {
                            setAbdLength(parseInt(e.target.value, 10));
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>
                    </div>
                  )}

                  {editType === "trans_splicing" && gene.exonCount !== null && gene.exonCount <= 1 && (
                    <div className="rounded-lg border border-amber-200 bg-amber-50 px-4 py-3 text-[12.5px] text-amber-800">
                      <strong>Warning:</strong> {gene.geneSymbol} appears to be a single-exon / intronless
                      gene ({gene.exonCount} exon). Trans-splicing relies on spliceosomal intron junctions and
                      is not applicable — the ranked mechanism will be marked ineligible.
                    </div>
                  )}
                </div>
              )}

              {selectedGoal === "TG05" && (
                <div className="space-y-4 px-6 pb-4">
                  <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                    <div>
                      <FieldLabel hint="What kind of molecular defect is the RNA itself driving?">
                        Molecular Defect Type <span className="text-red-500">*</span>
                      </FieldLabel>
                      <select
                        value={molecularDefect}
                        onChange={(e) => {
                          setMolecularDefect(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Select defect type</option>
                        {options?.rnaNeutralization.molecularDefects.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>

                    <div>
                      <FieldLabel hint="Which neutralization strategy do you want to use?">
                        Neutralization Mode <span className="text-red-500">*</span>
                      </FieldLabel>
                      <select
                        value={neutralizationMode}
                        onChange={(e) => {
                          setNeutralizationMode(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Select mode</option>
                        {options?.rnaNeutralization.neutralizationModes.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>

                    <div>
                      <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                        Delivery / Tissue Context
                      </FieldLabel>
                      <select
                        value={deliveryContext}
                        onChange={(e) => {
                          setDeliveryContext(e.target.value);
                          clearRanking();
                        }}
                        className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                      >
                        <option value="">Not specified</option>
                        {options?.deliveryContexts.map((o) => (
                          <option key={o.id} value={o.id}>
                            {o.label}
                          </option>
                        ))}
                      </select>
                    </div>
                  </div>

                  {molecularDefect === "loss_of_function" && (
                    <div className="rounded-lg border border-amber-200 bg-amber-50 px-4 py-3 text-[12.5px] text-amber-800">
                      <strong>Note:</strong> A pure loss-of-function defect cannot be addressed by RNA
                      neutralization — neutralizing the RNA cannot restore the missing protein. Ranked
                      mechanisms will be marked ineligible; consider TG02 (Gene Activation) or TG08
                      (Protein Replacement) instead.
                    </div>
                  )}

                  {neutralizationMode === "steric_repeat_masking" && (
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                      <div>
                        <FieldLabel hint="The pathogenic repeat motif in the transcript (e.g. CUG, CAG, GGGGCC)">
                          Repeat Unit (optional)
                        </FieldLabel>
                        <input
                          value={repeatUnit}
                          onChange={(e) => {
                            setRepeatUnit(e.target.value);
                            clearRanking();
                          }}
                          placeholder="e.g. CUG"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="Approximate number of repeat copies, e.g. '>50' or '55–200'">
                          Estimated Repeat Count (optional)
                        </FieldLabel>
                        <input
                          value={estimatedRepeatCount}
                          onChange={(e) => {
                            setEstimatedRepeatCount(e.target.value);
                            clearRanking();
                          }}
                          placeholder="e.g. &gt;50 copies"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="RNase H-independent chemistry for occupancy-only repeat masking">
                          Steric Chemistry
                        </FieldLabel>
                        <select
                          value={stericChemistry}
                          onChange={(e) => {
                            setStericChemistry(e.target.value);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.rnaNeutralization.stericChemistries.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="The RNA-binding protein sequestered by the toxic foci (e.g. MBNL1)">
                          Target RBP (optional)
                        </FieldLabel>
                        <input
                          value={targetRbp}
                          onChange={(e) => {
                            setTargetRbp(e.target.value);
                            clearRanking();
                          }}
                          placeholder="e.g. MBNL1"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="Oligo length in nucleotides (15–30 nt typical for repeat masking)">
                          Oligo Length (nt)
                        </FieldLabel>
                        <input
                          type="number"
                          min={15}
                          max={30}
                          value={oligoLength ?? ""}
                          onChange={(e) => {
                            setOligoLength(e.target.value === "" ? null : parseInt(e.target.value, 10));
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>
                    </div>
                  )}

                  {neutralizationMode === "microrna_antagomir" && (
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                      <div>
                        <FieldLabel hint="Chemistries suitable for anti-miR oligos (2'-O-MOE full-PS is the standard fully modified backbone)">
                          Oligo Chemistry
                        </FieldLabel>
                        <select
                          value={stericChemistry}
                          onChange={(e) => {
                            setStericChemistry(e.target.value);
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.rnaNeutralization.stericChemistries.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Antagomirs act on small non-coding RNAs — protein-coding mRNAs gate this mechanism ineligible">
                          Target Gene Type (optional)
                        </FieldLabel>
                        <input
                          value={targetGeneType}
                          onChange={(e) => {
                            setTargetGeneType(e.target.value);
                            clearRanking();
                          }}
                          placeholder="e.g. microRNA (non-coding RNA)"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="Length complementary to the mature miRNA (15–23 nt typical)">
                          Oligo Length (nt)
                        </FieldLabel>
                        <input
                          type="number"
                          min={15}
                          max={23}
                          value={oligoLength ?? ""}
                          onChange={(e) => {
                            setOligoLength(e.target.value === "" ? null : parseInt(e.target.value, 10));
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>
                    </div>
                  )}

                  {neutralizationMode === "aptamer_decoy" && (
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
                      <div>
                        <FieldLabel hint="The RNA-binding protein or toxic RNA the aptamer should sequester">
                          Target RBP / RNA (optional)
                        </FieldLabel>
                        <input
                          value={targetRbp}
                          onChange={(e) => {
                            setTargetRbp(e.target.value);
                            clearRanking();
                          }}
                          placeholder="e.g. MBNL1"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="Aptamer length in nucleotides (20–100 nt typical)">
                          Aptamer Length (nt)
                        </FieldLabel>
                        <input
                          type="number"
                          min={20}
                          max={100}
                          value={oligoLength ?? ""}
                          onChange={(e) => {
                            setOligoLength(e.target.value === "" ? null : parseInt(e.target.value, 10));
                            clearRanking();
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>
                    </div>
                  )}

                  {!molecularDefect && (
                    <p className="text-[12px] text-slate-500 italic">
                      Select a Molecular Defect Type above to see which neutralization mechanisms apply.
                    </p>
                  )}
                </div>
               )}

               {selectedGoal === "TG06" && (
                 <div className="space-y-4 px-6 pb-4">
                   <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                     <div>
                       <FieldLabel hint="Do you want to increase or decrease protein synthesis from the target transcript?">
                         Translational Goal <span className="text-red-500">*</span>
                       </FieldLabel>
                       <select
                         value={translationalGoal}
                         onChange={(e) => {
                           setTranslationalGoal(e.target.value);
                           clearRanking();
                         }}
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       >
                         <option value="">Select goal</option>
                         {options?.translationalRegulation.translationalGoals.map((o) => (
                           <option key={o.id} value={o.id}>
                             {o.label}
                           </option>
                         ))}
                       </select>
                     </div>

                     <div>
                       <FieldLabel hint="Which regulatory element in the transcript should the ASO target to modulate translation?">
                         Translational Target Element <span className="text-red-500">*</span>
                       </FieldLabel>
                       <select
                         value={targetElement}
                         onChange={(e) => {
                           setTargetElement(e.target.value);
                           clearRanking();
                         }}
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       >
                         <option value="">Select target element</option>
                         {options?.translationalRegulation.targetElements.map((o) => (
                           <option key={o.id} value={o.id}>
                             {o.label}
                           </option>
                         ))}
                       </select>
                     </div>

                     <div>
                       <FieldLabel hint="Steric-blocking chemistries only — DNA gapmers and siRNA are disabled">
                         Steric Chemistry Selection (Non-Cleaving Only)
                       </FieldLabel>
                       <select
                         value={translationStericChemistry}
                         onChange={(e) => {
                           setTranslationStericChemistry(e.target.value);
                           clearRanking();
                         }}
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       >
                         <option value="">Not specified</option>
                         {options?.translationalRegulation.stericChemistries.map((o) => (
                           <option key={o.id} value={o.id}>
                             {o.label}
                           </option>
                         ))}
                       </select>
                     </div>

                     <div>
                       <FieldLabel hint="The RNA-binding protein or miRNA whose interaction you want to block (e.g. MBNL1, eIF4E, miR-122)">
                         Target RBP / miRNA Symbol (optional)
                       </FieldLabel>
                       <input
                         value={translationTargetRbp}
                         onChange={(e) => {
                           setTranslationTargetRbp(e.target.value);
                           clearRanking();
                         }}
                         placeholder="e.g. miR-122"
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       />
                     </div>

                     <div>
                       <FieldLabel hint="Oligo length in nucleotides (18–25 nt typical for translational steric blocking)">
                         Oligo Length (nt)
                       </FieldLabel>
                       <input
                         type="number"
                         min={18}
                         max={25}
                         value={translationOligoLength ?? ""}
                         onChange={(e) => {
                           setTranslationOligoLength(e.target.value === "" ? null : parseInt(e.target.value, 10));
                           clearRanking();
                         }}
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       />
                     </div>

                     <div>
                       <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                         Delivery / Tissue Context
                       </FieldLabel>
                       <select
                         value={deliveryContext}
                         onChange={(e) => {
                           setDeliveryContext(e.target.value);
                           clearRanking();
                         }}
                         className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                       >
                         <option value="">Not specified</option>
                         {options?.deliveryContexts.map((o) => (
                           <option key={o.id} value={o.id}>
                             {o.label}
                           </option>
                         ))}
                       </select>
                     </div>
                   </div>

                   {!translationalGoal && !targetElement && (
                     <p className="text-[12px] text-slate-500 italic">
                       Select a Translational Goal and Target Element above to see which mechanisms apply.
                     </p>
                   )}
                  </div>
                )}

                {selectedGoal === "TG07" && (
                  <div className="space-y-4 px-6 pb-4">
                    {ieOptionsError && (
                      <p className="text-[12.5px] text-red-600">{ieOptionsError}</p>
                    )}

                    <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-3">
                      <div>
                        <FieldLabel hint="The gene symbol for the protein you want to engineer (e.g., DMD, CFTR, SMN2)">
                          Target Gene Symbol <span className="text-red-500">*</span>
                        </FieldLabel>
                        <input
                          value={ieTargetSymbol}
                          onChange={(e) => setIeTargetSymbol(e.target.value)}
                          placeholder="e.g. DMD, CFTR, SMN2"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="What type of isoform switching do you want to achieve?">
                          Isoform Goal <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={ieIsoformGoal}
                          onChange={(e) => setIeIsoformGoal(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select isoform goal</option>
                          {ieOptions?.isoformGoals.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                        {ieIsoformGoal && (
                          <p className="mt-1 text-[10.5px] text-slate-400">
                            {ieOptions?.isoformGoals.find((o) => o.id === ieIsoformGoal)?.description}
                          </p>
                        )}
                      </div>

                      <div>
                        <FieldLabel hint="The exon locus you want to target for isoform modulation">
                          Target Exon Locus <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={ieTargetExonLocus}
                          onChange={(e) => setIeTargetExonLocus(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select exon locus</option>
                          {ieOptions?.targetExonLoci.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Which splice regulatory element should the ASO target?">
                          Splice Element Target <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={ieSpliceElementTarget}
                          onChange={(e) => setIeSpliceElementTarget(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select splice element</option>
                          {ieOptions?.spliceElementTargets.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Steric-blocking chemistry for splice modulation">
                          Steric Chemistry <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={ieStericChemistry}
                          onChange={(e) => setIeStericChemistry(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select chemistry</option>
                          {ieOptions?.stericChemistries.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div className="flex items-end">
                        <label className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-3 cursor-pointer">
                          <input
                            type="checkbox"
                            checked={ieEnforceInFrame}
                            onChange={(e) => setIeEnforceInFrame(e.target.checked)}
                            className="h-4 w-4 rounded border-slate-300 text-brand focus:ring-brand"
                          />
                          <span className="text-[13px] text-slate-700">Enforce In-Frame splicing</span>
                        </label>
                      </div>
                    </div>

                    <div className="flex justify-end">
                      <button
                        onClick={async () => {
                          if (!ieTargetSymbol.trim() || !ieIsoformGoal || !ieTargetExonLocus || !ieSpliceElementTarget || !ieStericChemistry) return;
                          setIeLoading(true);
                          setIeError(null);
                          setSelectedIeCandidate(null);
                          try {
                            const res = await generateIeConstructs({
                              targetSymbol: ieTargetSymbol.trim(),
                              isoformGoal: ieIsoformGoal,
                              targetExonLocus: ieTargetExonLocus,
                              spliceElementTarget: ieSpliceElementTarget,
                              stericChemistry: ieStericChemistry,
                              enforceInFrame: ieEnforceInFrame,
                            });
                            setIeResults(res);
                          } catch (err) {
                            setIeError(err instanceof Error ? err.message : "Generation failed.");
                          } finally {
                            setIeLoading(false);
                          }
                        }}
                        disabled={!ieTargetSymbol.trim() || !ieIsoformGoal || !ieTargetExonLocus || !ieSpliceElementTarget || !ieStericChemistry || ieLoading}
                        className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                      >
                        {ieLoading && <Loader2 className="h-4 w-4 animate-spin" />}
                        {ieLoading ? "Generating Constructs..." : "Generate Constructs / Candidates"}
                      </button>
                    </div>
                  </div>
                )}

                {selectedGoal === "TG08" && (
                  <div className="space-y-4 px-6 pb-4">
                    {prOptionsError && (
                      <p className="text-[12.5px] text-red-600">{prOptionsError}</p>
                    )}

                    <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-3">
                      <div>
                        <FieldLabel hint="The gene symbol for the protein you want to replace (e.g., CFTR, PAH, GBA, F9)">
                          Target Protein Replacement Symbol <span className="text-red-500">*</span>
                        </FieldLabel>
                        <input
                          value={targetSymbol}
                          onChange={(e) => setTargetSymbol(e.target.value)}
                          placeholder="e.g. CFTR, PAH, GBA, F9"
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        />
                      </div>

                      <div>
                        <FieldLabel hint="Choose the RNA architecture best suited for your therapeutic goal">
                          RNA Modality <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={rnaModality}
                          onChange={(e) => {
                            setRnaModality(e.target.value);
                            if (e.target.value !== "circrna") setIresSelection("");
                          }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select modality</option>
                          {prOptions?.rnaModalities.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                        {rnaModality && (
                          <p className="mt-1 text-[10.5px] text-slate-400">
                            {prOptions?.rnaModalities.find((o) => o.id === rnaModality)?.description}
                          </p>
                        )}
                      </div>

                      <div>
                        <FieldLabel hint="Optimization strategy for the open reading frame">
                          Codon Optimization Strategy <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={codonStrategy}
                          onChange={(e) => setCodonStrategy(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select strategy</option>
                          {prOptions?.codonStrategies.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="5' and 3' UTR pair for optimal translation">
                          5' / 3' UTR Pair <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={utrPair}
                          onChange={(e) => setUtrPair(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select UTR pair</option>
                          {prOptions?.utrPairs.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="IRES element for cap-independent translation initiation (circRNA only)">
                          IRES Selection (circRNA only)
                        </FieldLabel>
                        <select
                          value={iresSelection}
                          onChange={(e) => setIresSelection(e.target.value)}
                          disabled={rnaModality !== "circrna"}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20 disabled:opacity-50 disabled:cursor-not-allowed"
                        >
                          <option value="">Select IRES</option>
                          {prOptions?.iresSelections.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                        {rnaModality !== "circrna" && (
                          <p className="mt-1 text-[10.5px] text-slate-400">
                            Select circRNA modality to enable IRES options
                          </p>
                        )}
                      </div>

                      <div>
                        <FieldLabel hint="Nucleotide chemical modification to reduce immunogenicity and enhance stability">
                          Nucleotide Chemical Modification <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={nucleotideModification}
                          onChange={(e) => setNucleotideModification(e.target.value)}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select modification</option>
                          {prOptions?.nucleotideModifications.map((o) => (
                            <option key={o.id} value={o.id}>
                              {o.label}
                            </option>
                          ))}
                        </select>
                      </div>
                    </div>

                    <div className="flex justify-end">
                      <button
                        onClick={async () => {
                          if (!targetSymbol.trim() || !rnaModality || !codonStrategy || !utrPair || !nucleotideModification) return;
                          if (rnaModality === "circrna" && !iresSelection) return;
                          setPrLoading(true);
                          setPrError(null);
                          setSelectedPrCandidate(null);
                          try {
                            const res = await generateConstructs({
                              targetSymbol: targetSymbol.trim(),
                              rnaModality,
                              codonStrategy,
                              utrPair,
                              iresSelection,
                              nucleotideModification,
                            });
                            setPrResults(res);
                          } catch (err) {
                            setPrError(err instanceof Error ? err.message : "Generation failed.");
                          } finally {
                            setPrLoading(false);
                          }
                        }}
                        disabled={!targetSymbol.trim() || !rnaModality || !codonStrategy || !utrPair || !nucleotideModification || (rnaModality === "circrna" && !iresSelection) || prLoading}
                        className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                      >
                        {prLoading && <Loader2 className="h-4 w-4 animate-spin" />}
                        {prLoading ? "Generating Constructs..." : "Generate Constructs / Candidates"}
                      </button>
                    </div>
                  </div>
                )}

                {selectedGoal === "TG09" && (
                  <div className="space-y-4 px-6 pb-4">
                    <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-3">
                      <div>
                        <FieldLabel hint="What class of structured RNA molecule do you want to design?">
                          Structural Class <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={tg09StructuralClass}
                          onChange={(e) => { setTg09StructuralClass(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select structural class</option>
                          {options?.rnaEngineering.structuralClasses.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="What type of molecule is the primary target?">
                          Target Molecule Type <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={tg09TargetType}
                          onChange={(e) => { setTg09TargetType(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select target type</option>
                          {options?.rnaEngineering.targetTypes.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Which structural scaffold architecture should the design use?">
                          Structural Scaffold Selection <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={tg09Scaffold}
                          onChange={(e) => { setTg09Scaffold(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select scaffold</option>
                          {options?.rnaEngineering.scaffolds.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="Which chemical base stabilization strategy should be applied?">
                          Chemical Base Stabilizations <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={tg09ChemStabilization}
                          onChange={(e) => { setTg09ChemStabilization(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select stabilization</option>
                          {options?.rnaEngineering.chemStabilizations.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="What is the target binding affinity threshold?">
                          Binding Threshold Parameter (Kd Goal) <span className="text-red-500">*</span>
                        </FieldLabel>
                        <select
                          value={tg09KdGoal}
                          onChange={(e) => { setTg09KdGoal(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Select Kd goal</option>
                          {options?.rnaEngineering.kdGoals.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <FieldLabel hint="General delivery/chemistry precedent — a soft tie-breaker, not a hard filter">
                          Delivery / Tissue Context
                        </FieldLabel>
                        <select
                          value={deliveryContext}
                          onChange={(e) => { setDeliveryContext(e.target.value); clearRanking(); }}
                          className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                        >
                          <option value="">Not specified</option>
                          {options?.deliveryContexts.map((o) => (
                            <option key={o.id} value={o.id}>{o.label}</option>
                          ))}
                        </select>
                      </div>
                    </div>
                  </div>
                )}

                {selectedGoal === "TG08" && prResults && (
                  <>
                    {/* Header & Construct Overview */}
                    <Card>
                      <SectionHeader step="2" title="Header & Construct Overview" />
                      <div className="px-6 pb-5">
                        <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                          <div className="space-y-3">
                            <InfoField label="Target Gene & Wild-Type RefSeq" value={`${prResults.overview.targetGene} · ${prResults.overview.refSeq}`} />
                            <InfoField label="Selected RNA Modality" value={prResults.overview.vectorTopology} />
                            <InfoField label="Native Protein Length" value={prResults.overview.nativeLength} />
                          </div>
                          <div className="space-y-3">
                            <InfoField label="Codon Adaptation Index (CAI)" value={prResults.overview.cai.toFixed(2)} valueClassName={prResults.overview.cai >= 0.92 ? "text-emerald-600" : "text-amber-600"} />
                            <InfoField label="Uridine Percentage (U%)" value={`${prResults.overview.uContent.toFixed(1)}%`} valueClassName={prResults.overview.uContent < 20 ? "text-emerald-600" : "text-amber-600"} />
                            <InfoField label="Predicted Intracellular Half-Life" value={prResults.overview.predictedHalfLife} />
                          </div>
                          <div className="space-y-3">
                            <InfoField label="Primary Mechanism Assigned" value={prResults.overview.primaryMechanism} />
                            <InfoField label="Expression Feasibility Score" value={`${prResults.overview.feasibilityScore}/100`} valueClassName={prResults.overview.feasibilityScore >= 80 ? "text-emerald-600" : "text-amber-600"} />
                          </div>
                        </div>
                      </div>
                    </Card>

                    {/* Candidate Construct Table */}
                    <Card>
                      <SectionHeader step="3" title="Candidate Construct Table" />
                      <div className="px-5 pb-5 overflow-x-auto">
                        <table className="w-full text-left">
                          <thead>
                            <tr className="border-b-2 border-slate-200">
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Rank</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Construct ID</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Modality</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">CAI</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">U%</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">MFE</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Initiation %</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Protein Yield</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">TLR Risk</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Signal Peptide</th>
                              <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Structure</th>
                              <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Action</th>
                            </tr>
                          </thead>
                          <tbody>
                            {prResults.candidates.map((c) => (
                              <tr
                                key={c.constructId}
                                className={`border-b border-slate-100 last:border-0 hover:bg-slate-50/50 transition-colors cursor-pointer ${selectedPrCandidate?.constructId === c.constructId ? "bg-brand/5" : ""}`}
                                onClick={() => setSelectedPrCandidate(c)}
                              >
                                <td className="py-3 pr-3">
                                  <span className="flex h-6 w-6 items-center justify-center rounded-full bg-brand/10 text-[11px] font-bold text-brand">
                                    {c.rank}
                                  </span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className="text-[11.5px] font-semibold text-slate-800 font-mono">{c.constructId}</span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className="text-[11px] text-slate-600">{c.modality}</span>
                                </td>
                                <td className="py-3 pr-3 text-right">
                                  <span className={`text-[11.5px] font-semibold ${c.cai >= 0.92 ? "text-emerald-600" : "text-amber-600"}`}>
                                    {c.cai.toFixed(2)}
                                  </span>
                                </td>
                                <td className="py-3 pr-3 text-right">
                                  <span className={`text-[11.5px] font-semibold ${c.uContent < 20 ? "text-emerald-600" : "text-amber-600"}`}>
                                    {c.uContent.toFixed(1)}%
                                  </span>
                                </td>
                                <td className="py-3 pr-3 text-right">
                                  <span className="text-[11.5px] font-mono text-slate-700">{c.mfe.toFixed(1)}</span>
                                </td>
                                <td className="py-3 pr-3 text-right">
                                  <span className="text-[11.5px] font-semibold text-slate-700">{c.initiationEfficiency}%</span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${c.predictedProteinYield.includes("Ultra") ? "border-emerald-200 bg-emerald-50 text-emerald-700" : c.predictedProteinYield.includes("High") ? "border-blue-200 bg-blue-50 text-blue-700" : "border-amber-200 bg-amber-50 text-amber-700"}`}>
                                    {c.predictedProteinYield}
                                  </span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${c.tlrRisk.includes("Very Low") ? "border-emerald-200 bg-emerald-50 text-emerald-700" : c.tlrRisk.includes("Low") ? "border-blue-200 bg-blue-50 text-blue-700" : c.tlrRisk.includes("Moderate") ? "border-amber-200 bg-amber-50 text-amber-700" : "border-red-200 bg-red-50 text-red-700"}`}>
                                    {c.tlrRisk}
                                  </span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className="text-[10.5px] text-slate-600">{c.signalPeptideStatus}</span>
                                </td>
                                <td className="py-3 pr-3">
                                  <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${c.secondaryStructureFlag === "PASSED" ? "border-emerald-200 bg-emerald-50 text-emerald-700" : "border-red-200 bg-red-50 text-red-700"}`}>
                                    {c.secondaryStructureFlag}
                                  </span>
                                </td>
                                <td className="py-3">
                                  <button
                                    onClick={(e) => { e.stopPropagation(); setSelectedPrCandidate(c); }}
                                    className="flex items-center gap-1 rounded-md border border-slate-200 px-2.5 py-1.5 text-[11px] font-medium text-slate-600 hover:bg-white transition-colors"
                                  >
                                    Select
                                  </button>
                                </td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </Card>

                    {/* Candidate Analysis & Visualizations */}
                    <Card>
                      <SectionHeader step="3a" title="Candidate Analysis & Visualizations" />
                      <div className="px-6 pb-5">
                        <ProteinReplacementDashboard candidates={prResults.candidates} />
                      </div>
                    </Card>

                    {/* Inspection Drawer */}
                    {selectedPrCandidate && (
                      <Card className="overflow-hidden">
                        <div className="flex items-center justify-between px-6 pt-4 pb-3 border-b border-slate-100">
                           <SectionHeader step="3b" title={`Inspection: ${selectedPrCandidate.constructId}`} />
                          <button onClick={() => setSelectedPrCandidate(null)} className="flex h-6 w-6 items-center justify-center rounded-full border border-slate-200 text-slate-400 hover:text-slate-600">
                            <X className="h-3.5 w-3.5" />
                          </button>
                        </div>
                        <div className="px-6 pb-5 space-y-5">
                          <div>
                            <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">Full Transcript Feature Map</p>
                            <div className="rounded-lg border border-slate-200 bg-slate-50 p-4 overflow-x-auto">
                              <div className="flex items-center gap-1 min-w-[600px]">
                                {selectedPrCandidate.features.map((f, i) => (
                                  <div key={i} className="flex items-center gap-1">
                                    <div className={`flex flex-col items-center rounded-md px-2 py-1.5 text-[10px] font-medium text-center min-w-[70px] ${f.type === "cap" ? "bg-indigo-100 text-indigo-700" : f.type === "ires" ? "bg-purple-100 text-purple-700" : f.type === "utr" ? "bg-blue-50 text-blue-600" : f.type === "kozak" ? "bg-emerald-50 text-emerald-600" : f.type === "orf" ? "bg-slate-200 text-slate-700" : f.type === "signal" ? "bg-amber-50 text-amber-700" : f.type === "utr3" ? "bg-blue-50 text-blue-600" : f.type === "polyA" ? "bg-rose-50 text-rose-600" : f.type === "scarsplice" ? "bg-teal-50 text-teal-700" : "bg-slate-100 text-slate-600"}`}>
                                      <span className="font-mono text-[9px] opacity-70">{f.start}-{f.end}</span>
                                      <span className="leading-tight">{f.name}</span>
                                    </div>
                                    {i < selectedPrCandidate.features.length - 1 && <span className="text-slate-300 text-[10px]">→</span>}
                                  </div>
                                ))}
                              </div>
                            </div>
                          </div>

                          <div>
                            <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">Construct Transcript Sequence (5' → 3')</p>
                            <div className="flex items-center gap-2">
                              <div className="flex-1 rounded-lg bg-slate-900 px-4 py-2.5 font-mono text-[11px] text-emerald-400 break-all leading-relaxed select-all">
                                {selectedPrCandidate.sequence}
                              </div>
                              <button onClick={() => { navigator.clipboard.writeText(selectedPrCandidate.sequence); setPrCopied(true); setTimeout(() => setPrCopied(false), 1500); }} className="flex shrink-0 items-center gap-1 rounded-md border border-slate-200 px-2.5 py-1.5 text-[11px] text-slate-500 hover:bg-white transition-colors">
                                {prCopied ? <><Check className="h-3 w-3 text-emerald-500" /> Copied</> : <><Copy className="h-3 w-3" /> Copy</>}
                              </button>
                            </div>
                          </div>

                          <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                            <div className="rounded-lg border border-slate-200 bg-white p-3">
                              <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">Amino Acid Sequence Identity</p>
                              <p className="text-[14px] font-bold text-emerald-600">{selectedPrCandidate.diagnostics.aminoAcidIdentity}%</p>
                              <p className="text-[10px] text-slate-500 mt-0.5">100% match with wild-type functional protein</p>
                            </div>
                            <div className="rounded-lg border border-slate-200 bg-white p-3">
                              <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">TLR3 / dsRNA Risk</p>
                              <p className="text-[14px] font-bold text-slate-700">{selectedPrCandidate.diagnostics.tlr3Score}</p>
                              <p className="text-[10px] text-slate-500 mt-0.5">Double-stranded RNA contaminant risk</p>
                            </div>
                            <div className="rounded-lg border border-slate-200 bg-white p-3">
                              <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">TLR7 / TLR8 Risk</p>
                              <p className="text-[14px] font-bold text-slate-700">{selectedPrCandidate.diagnostics.tlr7Score} / {selectedPrCandidate.diagnostics.tlr8Score}</p>
                              <p className="text-[10px] text-slate-500 mt-0.5">Uridine motif cluster activation risk</p>
                            </div>
                          </div>

                          <div>
                            <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">5' UTR / Ribosome Entry Secondary Structure (ViennaRNA MFE)</p>
                            <div className="rounded-lg border border-slate-200 bg-slate-900 px-4 py-3 font-mono text-[11px] text-slate-300 break-all leading-relaxed">
                              {selectedPrCandidate.diagnostics.mfePlot}
                            </div>
                            <div className="mt-2">
                              <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${selectedPrCandidate.diagnostics.fiveUtrHairpin ? "border-red-200 bg-red-50 text-red-700" : "border-emerald-200 bg-emerald-50 text-emerald-700"}`}>
                                {selectedPrCandidate.diagnostics.fiveUtrHairpin ? "FAILED — 5' UTR Hairpin Obstacle Detected" : "PASSED — No 5' UTR Hairpin Obstacles"}
                              </span>
                            </div>
                          </div>
                        </div>
                      </Card>
                    )}

                    {/* Downstream Actions */}
                    <Card>
                      <SectionHeader step="4" title="Downstream Action & Export" />
                      <div className="px-6 pb-5 flex flex-wrap items-center gap-3">
                        <button onClick={() => { const blob = new Blob([prResults.candidates.map((c) => `>${c.constructId}\n${c.sequence}`).join("\n\n")], { type: "text/plain" }); const url = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = url; a.download = `${targetSymbol}_constructs.fasta`; a.click(); URL.revokeObjectURL(url); }} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                          <Download className="h-4 w-4" /> Export Candidate Sequences to FASTA
                        </button>
                        <button onClick={() => alert("LNP Formulation & Payload Optimization module is under development.")} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                          <FlaskConical className="h-4 w-4" /> Proceed to LNP Formulation & Payload Optimization
                        </button>
                        <button onClick={() => alert("IVT Plasmid Template & Synthesis Protocol generation is under development.")} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                          <FileText className="h-4 w-4" /> Generate IVT Plasmid Template & Synthesis Protocol
                        </button>
                      </div>
                    </Card>
                  </>
                )}

                {selectedGoal !== "TG07" && (
                  <div className="px-6 pb-4">
                    
                  </div>
                )}

              <div className="flex items-center justify-between gap-3 px-6 pb-5">
                <p className="flex items-center gap-1 text-[11.5px] text-slate-400">
                  Results will appear below <ChevronDown className="h-3.5 w-3.5" />
                </p>
                <button
                  onClick={handleRank}
                  disabled={isRankDisabled()}
                  className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                >
                  {loading && <Loader2 className="h-4 w-4 animate-spin" />}
                  {loading ? "Ranking mechanisms..." : "Rank Mechanisms"}
                </button>
              </div>
            </Card>
          )}

          {error && (
            <div className="flex items-center gap-2 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-[13px] text-red-600">
              <AlertCircle className="h-4 w-4 shrink-0" />
              {error}
            </div>
          )}

          {ranking && (
            <div ref={rankingResultsRef} className="scroll-mt-5 space-y-3">
              <p className="text-[13px] text-slate-500">
                Ranked {ranking.results.length} {ranking.therapeuticGoal} mechanisms for {ranking.geneSymbol}.
                Eligible mechanisms match your selected defect type and scope; poor-fit mechanisms are shown below for reference.
              </p>
              {(() => {
                const eligible = ranking.results.filter((m) => m.eligible);
                const poorFit = ranking.results.filter((m) => !m.eligible);
                return (
                  <>
                    {eligible.length > 0 && (
                      <div className="space-y-3">
                        <p className="text-[11px] font-semibold uppercase tracking-wider text-emerald-600">
                          Eligible Mechanisms ({eligible.length})
                        </p>
                        {eligible.map((m) => (
                          <MechanismCard
                            key={m.id}
                            mechanism={m}
                            selected={selectedId === m.id}
                            onSelect={() => handleSelectMechanism(m.id)}
                          />
                        ))}
                      </div>
                    )}
                    {poorFit.length > 0 && (
                      <div className="space-y-3 rounded-xl border border-dashed border-[#E5E7EB] bg-slate-50/50 p-4">
                        <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-400">
                          Poor Fit ({poorFit.length})
                        </p>
                        {poorFit.map((m) => (
                          <MechanismCard
                            key={m.id}
                            mechanism={m}
                            selected={selectedId === m.id}
                            onSelect={() => handleSelectMechanism(m.id)}
                          />
                        ))}
                      </div>
                    )}
                  </>
                );
              })()}
            </div>
          )}

          {(selectedId || (selectedGoal === "TG08" && prResults)) && (
            <div className="flex justify-end pt-2">
              {selectedGoal === "TG03" ? (
                <div className="flex items-center gap-3">
                  <p className="text-[12px] text-slate-500">
                    The ASO design pipeline for RNA editing mechanisms (guide RNA / PTM design) is under development.
                  </p>
                  <span className="flex items-center gap-2 rounded-lg border border-[#E5E7EB] bg-slate-50 px-5 py-3 text-[13px] font-medium text-slate-400">
                    Design pipeline coming soon
                  </span>
                </div>
              ) : selectedGoal === "TG02" ? (
                <button
                  onClick={() => router.push("/gene-upregulation")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  <Beaker className="h-4 w-4" />
                  Proceed to Gene Upregulation Design
                </button>
              ) : selectedGoal === "TG04" ? (
                <button
                  onClick={() => router.push("/gene-silencing")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  <Beaker className="h-4 w-4" />
                  Proceed to RNA Processing Design
                </button>
              ) : selectedGoal === "TG05" ? (
                <button
                  onClick={() => router.push("/rna-neutralization")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  Proceed to RNA Neutralization
                </button>
               ) : selectedGoal === "TG06" ? (
                 <button
                   onClick={() => router.push("/translational-regulation")}
                   className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                 >
                   Proceed to Translational Regulation
                 </button>
               ) : selectedGoal === "TG07" ? (
                 <button
                   onClick={() => router.push("/isoform-engineering")}
                   className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                 >
                   <Beaker className="h-4 w-4" />
                   Proceed to Isoform Engineering Design
                 </button>
               ) : selectedGoal === "TG08" ? (
                <button
                  onClick={() => router.push("/protein-replacement")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  <Beaker className="h-4 w-4" />
                  Proceed to Protein Replacement Design
                </button>
              ) : selectedGoal === "TG09" ? (
                <button
                  onClick={() => router.push("/rna-engineering")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  <Beaker className="h-4 w-4" />
                  Proceed to RNA Engineering
                </button>
              ) : (
                <button
                  onClick={() => router.push("/gene-silencing")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-6 py-3 text-[14px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                >
                  Proceed to Gene Silencing
                </button>
              )}
            </div>
          )}
        </main>
      </div>
    </div>
  );
}
