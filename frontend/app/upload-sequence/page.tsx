"use client";

import { useCallback, useRef, useState } from "react";
import {
  UploadCloud,
  FileText,
  AlertCircle,
  Loader2,
  CheckCircle2,
  ArrowRight,
  ArrowLeft,
  Dna,
  Scissors,
  BookOpen,
  Syringe,
  Zap,
  Shield,
  Thermometer,
  Target,
  BarChart3,
  Star,
  ChevronDown,
  ChevronUp,
  Info,
  LayoutDashboard,
  AlignLeft,
  List,
  GitBranch,
  FlaskConical,
} from "lucide-react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import { Card, SectionHeader } from "@/components/ui";
import { ValidationReport, AnalysisReport, GrnaCandidate } from "@/types/upload";
import { Modality } from "@/types/upload";
import type { AnalyzeResponse } from "@/types/upload-types";
import { validateSequence, analyzeSequence } from "@/lib/uploadApi";
import { saveReport } from "@/lib/auth";
import SequenceTrackViewer from "@/components/SequenceTrackViewer";
import GcContentChart from "@/components/GcContentChart";
import NucleotideCompositionChart from "@/components/NucleotideCompositionChart";
import ExportMenu from "@/components/ExportMenu";
import MeltingTemperatureCard from "@/components/MeltingTemperatureCard";
import CodonUsageCard from "@/components/CodonUsageCard";
import SequenceComplexityCard from "@/components/SequenceComplexityCard";
import ModificationScorecard from "@/components/ModificationScorecard";
import StackingEnergyChart from "@/components/StackingEnergyChart";
import RestrictionSiteMap from "@/components/RestrictionSiteMap";
import MiRNATargetingCard from "@/components/MiRNATargetingCard";
import HairpinDiagram from "@/components/HairpinDiagram";
import KmerFrequencyChart from "@/components/KmerFrequencyChart";
import ThermodynamicProfile from "@/components/ThermodynamicProfile";
import SequenceAnnotationBar from "@/components/SequenceAnnotationBar";
import BasePairDotPlot from "@/components/BasePairDotPlot";
import ModificationLandscapeCard from "@/components/ModificationLandscapeCard";
import RiskScoreDashboard from "@/components/RiskScoreDashboard";
import PhysicochemicalCard from "@/components/PhysicochemicalCard";
import StabilityIndexChart from "@/components/StabilityIndexChart";
import AnalysisTabs from "@/components/AnalysisTabs";
import PairwiseAlignmentViewer from "@/components/PairwiseAlignmentViewer";
import FeatureTable from "@/components/FeatureTable";
import AnnotatedSequenceViewer from "@/components/AnnotatedSequenceViewer";
import CrisprScanTrack from "@/components/CrisprScanTrack";
import CrisprCandidateTable from "@/components/CrisprCandidateTable";
import CrisprPrimerDesignCard from "@/components/CrisprPrimerDesignCard";
import InfoTooltip from "@/components/InfoTooltip";

type Step = "upload" | "validate" | "modality" | "analysis";

const MODALITIES = [
  {
    id: "aso",
    name: "ASO / Splice Switching",
    icon: Scissors,
    color: "indigo",
    description: "Exon skipping/inclusion or RNase H-mediated degradation.",
  },
  {
    id: "sirna",
    name: "siRNA / shRNA Design",
    icon: Zap,
    color: "emerald",
    description: "Target knock-down via RISC pathway.",
  },
  {
    id: "mrna",
    name: "mRNA Design / Optimization",
    icon: BookOpen,
    color: "blue",
    description: "Codon optimization, UTR selection, and nucleoside modification.",
  },
  {
    id: "sgrna",
    name: "sgRNA / CRISPR Targeting",
    icon: Target,
    color: "purple",
    description: "Off-target analysis and guide selection.",
  },
];

export default function UploadSequencePage() {
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [step, setStep] = useState<Step>("upload");
  const [rawInput, setRawInput] = useState("");
  const [filename, setFilename] = useState<string | undefined>();
  const [validation, setValidation] = useState<ValidationReport | null>(null);
  const [selectedModality, setSelectedModality] = useState<string | null>(null);
  const [analysis, setAnalysis] = useState<AnalysisReport | null>(null);
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [expandedRecs, setExpandedRecs] = useState<Record<string, boolean>>({});
  const [activeTab, setActiveTab] = useState("overview");

  const handleFileUpload = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setFilename(file.name);
    const reader = new FileReader();
    reader.onload = (ev) => {
      const text = ev.target?.result as string;
      setRawInput(text);
    };
    reader.readAsText(file);
  }, []);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    const file = e.dataTransfer.files?.[0];
    if (!file) return;
    setFilename(file.name);
    const reader = new FileReader();
    reader.onload = (ev) => {
      const text = ev.target?.result as string;
      setRawInput(text);
    };
    reader.readAsText(file);
  }, []);

  function clientSideValidate(raw: string, fname?: string) {
    const seq = raw.replace(/^>.*$/gm, "").replace(/[^A-Za-z]/g, "").toUpperCase();
    if (!seq) return null;

    const hasT = seq.includes("T");
    const hasU = seq.includes("U");
    let seqType: "dna" | "rna" | "protein" = "dna";
    if (hasT && !hasU) seqType = "dna";
    else if (hasU && !hasT) seqType = "rna";
    else if (!hasT && !hasU) {
      const proteinLetters = new Set("DEFHIKLMNPQRSTVWY");
      if (seq.split("").some((b) => proteinLetters.has(b))) seqType = "protein";
      else seqType = "dna";
    }

    const invalid = [...new Set(seq.split("").filter((b) => !"ACGTRYSWKMBDHVNUEFHIKLMNPQRSTVWY".includes(b)))].sort();
    const features: string[] = [];

    if (seqType === "protein") {
      features.push("Protein sequence detected");
      const aaCounts: Record<string, number> = {};
      for (const aa of "ACDEFGHIKLMNPQRSTVWY") aaCounts[aa] = 0;
      for (const b of seq) if (aaCounts[b] !== undefined) aaCounts[b]++;
      const top = Object.entries(aaCounts).sort((a, b) => b[1] - a[1])[0];
      if (top && top[1] > 0) features.push(`Most common residue: ${top[0]} (${top[1]}x)`);
    } else {
      const gc = seq.length > 0 ? Math.round((seq.split("").filter((b) => "GC".includes(b)).length / seq.length) * 1000) / 10 : 0;
      features.push(`GC content: ${gc}%`);
      if (/A{6,}$/.test(seq)) features.push("Poly-A tail detected");
      if (/G{4,}/.test(seq)) features.push("Poly-G tract detected");

      const isRna = hasU && !hasT;
      const startCodons = isRna ? ["AUG"] : ["ATG"];
      const stopCodons = isRna ? ["UAA", "UAG", "UGA"] : ["TAA", "TAG", "TGA"];
      const orfs: { strand: string; frame: number; start: number; end: number; length: number; proteinLength: number }[] = [];
      function revComp(s: string): string {
        const comp: Record<string, string> = isRna ? { A: "U", U: "A", G: "C", C: "G" } : { A: "T", T: "A", G: "C", C: "G" };
        return s.split("").reverse().map((b) => comp[b] ?? "N").join("");
      }
      function scanOrfs(s: string, strand: string) {
        for (let frame = 0; frame < 3; frame++) {
          let i = frame;
          while (i < s.length - 2) {
            const codon = s.slice(i, i + 3);
            if (startCodons.includes(codon)) {
              const start = i;
              let j = i + 3;
              while (j < s.length - 2) {
                if (stopCodons.includes(s.slice(j, j + 3))) {
                  orfs.push({ strand, frame: frame + 1, start: start + 1, end: j + 3, length: j + 3 - start, proteinLength: (j - start) / 3 });
                  break;
                }
                j += 3;
              }
              i = j + 3 < s.length ? j + 3 : s.length;
            } else { i += 3; }
          }
        }
      }
      scanOrfs(seq, "+");
      scanOrfs(revComp(seq), "-");
      if (orfs.length > 0) {
        const best = orfs.reduce((a, b) => (b.proteinLength > a.proteinLength ? b : a));
        features.push(`Longest ORF: ${best.proteinLength} aa (${best.strand} strand, frame ${best.frame})`);
      }
    }

    return {
      valid: invalid.length === 0,
      sequence: seq,
      sequenceType: seqType,
      length: seq.length,
      gcContent: seqType === "protein" ? null : (seq.length > 0 ? Math.round((seq.split("").filter((b) => "GC".includes(b)).length / seq.length) * 1000) / 10 : 0),
      invalidChars: invalid,
      features,
      orfs: seqType === "protein" ? [] : [],
      filename: fname,
      hasPolyA: seqType !== "protein" && /A{6,}$/.test(seq),
      hasPolyG: seqType !== "protein" && /G{4,}/.test(seq),
    };
  }

  function generateGrnaCandidates(seq: string): GrnaCandidate[] {
    const upper = seq.toUpperCase();
    const length = upper.length;
    const candidates: GrnaCandidate[] = [];

    const scanForward = (s: string, strand: "+" | "-") => {
      const searchSeq = strand === "-" ? revCompLocal(s) : s;
      for (let i = 0; i <= searchSeq.length - 23; i++) {
        const pam = searchSeq.slice(i + 20, i + 23);
        if (!/^NGG$/.test(pam)) continue;

        const spacer = searchSeq.slice(i, i + 20);
        const gc = spacer.split("").filter((b) => "GC".includes(b)).length;
        const gcPct = gc / 20;

        let score = 0;
        if (gcPct >= 0.4 && gcPct <= 0.8) score += 40;
        else if (gcPct >= 0.3 && gcPct <= 0.85) score += 25;
        else score += 10;

        if (spacer.startsWith("G") || spacer.startsWith("C")) score += 10;
        if (!/TTTT/.test(spacer)) score += 10;

        const selfComp = selfComplementarityScore(spacer);
        if (selfComp < 0.2) score += 15;
        else if (selfComp < 0.4) score += 8;

        const offTargets = estimateOffTargets(spacer);
        if (offTargets <= 2) score += 15;
        else if (offTargets <= 5) score += 8;
        else score += 3;

        score = Math.min(100, Math.max(0, score));

        const color = score >= 70 ? "emerald" : score >= 40 ? "amber" : "rose";

        // Specificity score (MIT-style: 100 / (100 + sum of off-target penalties))
        const specificityScore = Math.min(100, Math.max(0, 100 - offTargets * 8 - selfComp * 15));
        // Efficiency score (Doench 2016-inspired: favors GC 40-80%, G-start, low self-complementarity)
        let effScore = 0;
        if (gcPct >= 0.4 && gcPct <= 0.8) effScore += 35;
        else if (gcPct >= 0.3 && gcPct <= 0.9) effScore += 20;
        if (spacer.startsWith("G")) effScore += 20;
        else if (spacer.startsWith("A") || spacer.startsWith("C")) effScore += 10;
        if (!/TTTT/.test(spacer)) effScore += 15;
        effScore = Math.max(0, Math.min(100, Math.round(effScore + (30 - selfComp * 30))));
        // Mismatch distribution: counts at 0,1,2,3,4 mismatches (CRISPOR-style)
        const mismatchDistribution = computeMismatchDistribution(spacer, searchSeq);

        // Map position to the forward strand for consistent display
        const fwdPos = strand === "+" ? i + 21 : length - i - 22;

        candidates.push({
          id: `gRNA-${candidates.length + 1}`,
          position: fwdPos,
          sequence: spacer,
          pam,
          strand,
          score: Math.round(score * 10) / 10,
          gc: Math.round(gcPct * 1000) / 10,
          selfComplementarity: Math.round(selfComp * 1000) / 1000,
          offTargets,
          polyT: /TTTT/.test(spacer),
          color,
          specificityScore: Math.round(specificityScore),
          efficiencyScore: Math.round(effScore),
          mismatchDistribution,
        });
      }
    };

    function revCompLocal(s: string): string {
      const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };
      return s.split("").reverse().map((b) => comp[b] ?? "N").join("");
    }

    function computeMismatchDistribution(spacer: string, genome: string): number[] {
      const counts = [0, 0, 0, 0, 0];
      for (let i = 0; i <= genome.length - 20; i++) {
        const candidate = genome.slice(i, i + 20);
        let mm = 0;
        for (let j = 0; j < 20 && mm <= 4; j++) {
          if (candidate[j] !== spacer[j]) mm++;
        }
        if (mm <= 4) {
          counts[mm]++;
        }
      }
      return counts;
    }

    scanForward(upper, "+");
    scanForward(upper, "-");

    return candidates;
  }

  function selfComplementarityScore(seq: string): number {
    const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", U: "A" };
    const revComp = seq.split("").reverse().map((b) => comp[b] ?? "N").join("");
    let matches = 0;
    const len = Math.min(seq.length, revComp.length);
    for (let i = 0; i < len; i++) {
      if (seq[i] === revComp[i]) matches++;
    }
    return len > 0 ? matches / len : 0;
  }

  function estimateOffTargets(seq: string): number {
    const k = 6;
    const kmers = new Set<string>();
    for (let i = 0; i <= seq.length - k; i++) {
      kmers.add(seq.slice(i, i + k));
    }
    const unique = kmers.size;
    const total = seq.length - k + 1;
    const repetitiveness = total > 0 ? 1 - unique / total : 0;
    return Math.round(repetitiveness * 12);
  }

  function clientSideAnalyze(seq: string, modality: string) {
    const hasT = seq.includes("T");
    const hasU = seq.includes("U");
    let seqType: "dna" | "rna" | "protein" = "dna";
    if (hasT && !hasU) seqType = "dna";
    else if (hasU && !hasT) seqType = "rna";
    else if (!hasT && !hasU) {
      const proteinLetters = new Set("DEFHIKLMNPQRSTVWY");
      if (seq.split("").some((b) => proteinLetters.has(b))) seqType = "protein";
      else seqType = "dna";
    }

    if (seqType === "protein") {
      const aaCounts: Record<string, number> = {};
      for (const aa of "ACDEFGHIKLMNPQRSTVWY") aaCounts[aa] = 0;
      for (const b of seq) if (aaCounts[b] !== undefined) aaCounts[b]++;
      const length = seq.length;
      const hydrophobic = "AILMFWYV".split("").reduce((s, a) => s + (aaCounts[a] || 0), 0);
      const hydrophilic = "RNDQEKHP".split("").reduce((s, a) => s + (aaCounts[a] || 0), 0);
      const charged = "RDEKHP".split("").reduce((s, a) => s + (aaCounts[a] || 0), 0);
      const aromatic = "FWY".split("").reduce((s, a) => s + (aaCounts[a] || 0), 0);

      return {
        sequence: seq,
        sequenceType: "protein",
        length,
        gcContent: null,
        offTarget: { lengthBasedRiskEstimate: "N/A" as const, note: "Protein sequence — nucleotide-level off-target screening not applicable", internalRepetitiveness: 0, recommendedMinLength: 0, disclaimer: "Not applicable for protein sequences." },
        secondaryStructure: { estimatedMfe: null, palindromicRegions: 0, palindromePositions: [], gcContent: null, hairpinRisk: "N/A" as const },
        immuneScreen: [],
        modality: { recommendations: ["Protein sequence uploaded — select a nucleic-acid modality for ASO/siRNA/mRNA/sgRNA design against the coding transcript."], recommendedChemistry: undefined, optimalLength: undefined, targetRegion: undefined, strand: undefined, needsCodonOptimization: undefined, needsPolyA: undefined, needsUTR: undefined, nucleosideModifications: undefined, casProtein: undefined, offTargetMitigation: undefined },
        gcCurve: [],
        composition: aaCounts,
        orfs: [],
        meltingTemp: undefined,
        complexity: undefined,
        codonUsage: undefined,
        modificationScores: undefined,
        energyProfile: [],
        restrictionSites: [],
        mirnaTargets: [],
        hairpins: [],
        kmerFrequency: undefined,
        thermoProfile: undefined,
        dotPlot: [],
        modificationLandscape: [],
        riskScores: { specificity: 0, stability: 0, immunogenicity: 0, delivery: 0, toxicity: 0, overall: 0 },
        physicochemical: { molecularWeight: 0, netCharge: 0, hydrophobicityIndex: 0, hydrophobicityProfile: [], chargeProfile: [] },
        stabilityIndex: [],
        grnaCandidates: [],
        proteinAnalysis: {
          aminoAcidComposition: aaCounts,
          molecularWeight: 0,
          length,
          hydrophobicFraction: +(hydrophobic / length).toFixed(3),
          hydrophilicFraction: +(hydrophilic / length).toFixed(3),
          chargedFraction: +(charged / length).toFixed(3),
          aromaticFraction: +(aromatic / length).toFixed(3),
        },
      };
    }

    const _gc = (s: string) => s.length > 0 ? Math.round((s.split("").filter((b) => "GC".includes(b)).length / s.length) * 1000) / 10 : 0;
    const gc = _gc(seq);
    const length = seq.length;
    const offRisk = length < 18 ? "High" : length < 20 ? "Medium" : "Low";
    const k = 6;
    const kmers = length >= k ? Array.from({ length: length - k + 1 }, (_, i) => seq.slice(i, i + k)) : [];
    const repetitiveness = kmers.length > 0 ? 1 - new Set(kmers).size / kmers.length : 0;

    // Palindrome positions
    const palindromePositions: number[] = [];
    for (let i = 0; i < length - 5; i++) {
      const chunk = seq.slice(i, i + 6);
      if (chunk === chunk.split("").reverse().join("")) palindromePositions.push(i + 1);
    }
    const palindromes = palindromePositions.length;
    const mfe = Math.round((gc / 100 * -1.5 + (100 - gc) / 100 * -0.9) * length / 2 * 10) / 10;

    // Immune hits with positions
    const immune: { motif: string; label: string; start: number; end: number }[] = [];
    const guRichRe = /[GU]{2,}U[GU]{2,}/gi;
    let m;
    while ((m = guRichRe.exec(seq)) !== null && immune.length < 40) {
      immune.push({ motif: m[0], label: "GU-rich stretch (literature-associated with TLR7/8 sensing; not a confirmed motif)", start: m.index + 1, end: m.index + m[0].length });
    }
    const homoRe = /(.)\1{3,}/g;
    while ((m = homoRe.exec(seq)) !== null && immune.length < 40) {
      immune.push({ motif: m[0].slice(0, 6), label: "Homopolymer run (4+ repeats; general repetitive-element flag)", start: m.index + 1, end: m.index + m[0].length });
    }
    if (seqType === "dna") {
      const cpgRe = /[AG][AG]CG[CT][CT]/gi;
      while ((m = cpgRe.exec(seq)) !== null && immune.length < 40) {
        immune.push({ motif: m[0], label: "Unmethylated CpG in a purine-purine-CG-pyrimidine-pyrimidine context (literature TLR9 motif pattern; not a confirmed assay)", start: m.index + 1, end: m.index + m[0].length });
      }
    }

    // ORFs (both strands)
    const orfs: { strand: string; frame: number; start: number; end: number; length: number; proteinLength: number }[] = [];
    const isRna = hasU && !hasT;
    const startCodons = isRna ? ["AUG"] : ["ATG"];
    const stopCodons = isRna ? ["UAA", "UAG", "UGA"] : ["TAA", "TAG", "TGA"];
    function revComp(s: string): string {
      const comp: Record<string, string> = isRna ? { A: "U", U: "A", G: "C", C: "G" } : { A: "T", T: "A", G: "C", C: "G" };
      return s.split("").reverse().map((b) => comp[b] ?? "N").join("");
    }
    function scanOrfs(s: string, strand: string) {
      for (let frame = 0; frame < 3; frame++) {
        let i = frame;
        while (i < s.length - 2) {
          if (startCodons.includes(s.slice(i, i + 3))) {
            const start = i;
            let j = i + 3;
            while (j < s.length - 2) {
              if (stopCodons.includes(s.slice(j, j + 3))) {
                orfs.push({ strand, frame: frame + 1, start: start + 1, end: j + 3, length: j + 3 - start, proteinLength: (j - start) / 3 });
                break;
              }
              j += 3;
            }
            i = j + 3 < s.length ? j + 3 : s.length;
          } else { i += 3; }
        }
      }
    }
    scanOrfs(seq, "+");
    scanOrfs(revComp(seq), "-");

    // GC sliding window
    const windowSize = 10;
    const step = 2;
    const gcCurve: { position: number; gc: number }[] = [];
    for (let i = 0; i <= length - windowSize; i += step) {
      const chunk = seq.slice(i, i + windowSize);
      const gcCount = chunk.split("").filter((b) => "GC".includes(b)).length;
      gcCurve.push({ position: i + 1, gc: Math.round((gcCount / windowSize) * 1000) / 10 });
    }

    // Nucleotide composition
    const composition = { A: (seq.match(/A/g) || []).length, C: (seq.match(/C/g) || []).length, G: (seq.match(/G/g) || []).length, T: (seq.match(/T/g) || []).length, U: (seq.match(/U/g) || []).length };

    const offNote = offRisk === "Low" ? "Adequate length for specificity in general, not verified against any genome" : offRisk === "Medium" ? "Moderate length — not verified against any genome" : "Short sequence — generally correlates with higher off-target probability, not verified against any genome";
    const modRecs: string[] = [];
    let modDetails: Record<string, unknown> = { recommendations: modRecs };
    if (modality === "aso") {
      if (gc < 30) modRecs.push("Low GC% — consider LNA or 2'-OMe modifications to boost Tm");
      else if (gc > 70) modRecs.push("High GC% — risk of G-quadruplexes; consider shorter ASO");
      else modRecs.push("GC content in optimal range for RNase H recruitment");
      if (length < 15) modRecs.push("Very short — high off-target risk; minimum 18 nt recommended");
      else if (length > 25) modRecs.push("Long ASO — may have reduced cellular uptake; consider gapmer design");
      modDetails = { ...modDetails, recommendedChemistry: gc >= 35 ? "gapmer" : "pmo", optimalLength: "18-22 nt", targetRegion: "Exon junction or mutated region recommended" };
    } else if (modality === "sirna") {
      if (length < 19 || length > 25) modRecs.push("Optimal siRNA length is 19-25 nt");
      if (gc < 30 || gc > 52) modRecs.push("Optimal GC content for siRNA is 30-52%");
      modRecs.push("Guide strand + Passenger strand design");
      modDetails = { ...modDetails, strand: "Guide strand (antisense) + Passenger strand", optimalLength: "21 nt with 2-nt 3' overhangs" };
    } else if (modality === "mrna") {
      modRecs.push("Consider 5' Cap analog (Anti-Reverse Cap ARCA)");
      modRecs.push("Evaluate codon optimization for human expression");
      if (!/A{6,}$/.test(seq)) modRecs.push("Add 100-150 nt poly(A) for stability");
      modDetails = { ...modDetails, needsCodonOptimization: true, needsPolyA: !/A{6,}$/.test(seq), needsUTR: true, nucleosideModifications: ["N1-methylpseudouridine (m1Ψ)", "5-methylcytidine (m5C)"] };
    } else if (modality === "sgrna") {
      if (length < 17 || length > 21) modRecs.push("Optimal sgRNA spacer length is 17-21 nt");
      modRecs.push("Requires NGG PAM adjacent to target site (SpCas9)");
      if (/TTTT/.test(seq)) modRecs.push("Poly-T tract detected — may cause premature transcription termination");
      modDetails = { ...modDetails, casProtein: "SpCas9 (NGG PAM)", optimalLength: "20 nt spacer + PAM" };
    }

    // --- Melting temperature (nearest-neighbor simplified) ---
    const DNA_NN: Record<string, [number, number]> = {
      "AA": [-7.9, -22.2], "TT": [-7.9, -22.2], "AT": [-7.2, -20.4], "TA": [-7.2, -21.3],
      "CA": [-8.5, -22.7], "TG": [-8.5, -22.7], "GT": [-8.4, -22.4], "AC": [-8.4, -22.4],
      "CT": [-7.8, -21.0], "AG": [-7.8, -21.0], "GA": [-8.2, -22.2], "TC": [-8.2, -22.2],
      "CG": [-10.6, -27.2], "GC": [-9.8, -24.4], "GG": [-8.0, -19.9], "CC": [-8.0, -19.9],
    };
    const nnSeq = seq.replace(/U/g, "T");
    let dH = 0, dS = 0, nnCount = 0;
    for (let i = 0; i < nnSeq.length - 1; i++) {
      const dinuc = nnSeq.slice(i, i + 2);
      if (DNA_NN[dinuc]) { dH += DNA_NN[dinuc][0]; dS += DNA_NN[dinuc][1]; nnCount++; }
    }
    let tmNN = 0;
    if (nnCount > 0 && dS !== 0) {
      const R = 1.987, C = 250e-6;
      tmNN = Math.round((dH * 1000 / (dS + R * Math.log(C / 4)) - 273.15) * 10) / 10;
    }
    const tmBasicGC = length > 0 ? Math.round((64.9 + 41 * (gc - 16.4) / length) * 10) / 10 : 0;

    // --- Sequence complexity ---
    const dinucRepeats: { pattern: string; start: number; end: number; repeats: number }[] = [];
    for (let i = 0; i < length - 3; i++) {
      const dinuc = seq.slice(i, i + 2);
      let run = 1, j = i + 2;
      while (j <= length - 2 && seq.slice(j, j + 2) === dinuc) { run++; j += 2; }
      if (run >= 3) { dinucRepeats.push({ pattern: dinuc, start: i + 1, end: i + run * 2, repeats: run }); if (dinucRepeats.length >= 20) break; }
    }
    const gcRichRegions: { start: number; end: number; length: number }[] = [];
    const atRichRegions: { start: number; end: number; length: number }[] = [];
    let m2: RegExpExecArray | null;
    const gcRichRe = /[GC]{5,}/g;
    while ((m2 = gcRichRe.exec(seq)) !== null && gcRichRegions.length < 10) gcRichRegions.push({ start: m2.index + 1, end: m2.index + m2[0].length, length: m2[0].length });
    const atRichRe = /[AT]{5,}/g;
    while ((m2 = atRichRe.exec(seq)) !== null && atRichRegions.length < 10) atRichRegions.push({ start: m2.index + 1, end: m2.index + m2[0].length, length: m2[0].length });
    const complexityScore = Math.round((1 - dinucRepeats.length / Math.max(length, 1)) * 1000) / 1000;

    // --- Codon usage ---
    const HUMAN_CODON_ADAPT: Record<string, number> = {
      "UUU": 0.52, "UUC": 0.48, "UUA": 0.07, "UUG": 0.13, "CUU": 0.13, "CUC": 0.20, "CUA": 0.07, "CUG": 0.40,
      "AUU": 0.36, "AUC": 0.47, "AUA": 0.18, "AUG": 1.00, "GUU": 0.18, "GUC": 0.24, "GUA": 0.12, "GUG": 0.46,
      "UCU": 0.19, "UCC": 0.22, "UCA": 0.15, "UCG": 0.06, "CCU": 0.19, "CCC": 0.20, "CCA": 0.20, "CCG": 0.06,
      "ACU": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11, "GCU": 0.21, "GCC": 0.27, "GCA": 0.23, "GCG": 0.09,
      "UAU": 0.44, "UAC": 0.56, "UAA": 0.30, "UAG": 0.24, "CAU": 0.42, "CAC": 0.58, "CAA": 0.27, "CAG": 0.73,
      "AAU": 0.47, "AAC": 0.53, "AAA": 0.43, "AAG": 0.57, "GAU": 0.46, "GAC": 0.54, "GAA": 0.42, "GAG": 0.58,
      "UGU": 0.45, "UGC": 0.55, "UGA": 0.26, "UGG": 1.00, "CGU": 0.08, "CGC": 0.19, "CGA": 0.06, "CGG": 0.21,
      "AGU": 0.15, "AGC": 0.22, "AGA": 0.21, "AGG": 0.20, "GGU": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
    };
    const codonRna = seq.replace(/T/g, "U");
    const codons: { codon: string; position: number; adaptiveness: number; isRare: boolean }[] = [];
    const rareCodons: { codon: string; position: number; adaptiveness: number }[] = [];
    let caiSum = 0, caiCount = 0;
    let codingStart: number | null = null;
    for (let i = 0; i < length - 2; i++) { if (codonRna.slice(i, i + 3) === "AUG") { codingStart = i; break; } }
    if (codingStart !== null) {
      let ci = codingStart;
      while (ci < length - 2) {
        const codon = codonRna.slice(ci, ci + 3);
        if (codon === "UAA" || codon === "UAG" || codon === "UGA") break;
        const adapt = HUMAN_CODON_ADAPT[codon] ?? 0.5;
        const isRare = adapt < 0.2;
        codons.push({ codon, position: ci + 1, adaptiveness: adapt, isRare });
        if (isRare) rareCodons.push({ codon, position: ci + 1, adaptiveness: adapt });
        caiSum += adapt; caiCount++;
        ci += 3;
      }
    }

    // --- Modification scores ---
    const modScores: Record<string, { score: number; rationale: string }> = {};
    if (modality === "aso") {
      modScores.lnaBoosting = { score: Math.max(0, Math.min(100, Math.round((70 - gc) * 1.5))), rationale: gc < 40 ? "Low GC benefits most from LNA-mediated Tm boost" : "GC already adequate — fewer LNA substitutions needed" };
      modScores.gapmerSuitability = { score: Math.max(0, Math.min(100, length >= 18 ? Math.round(80 + (gc - 40) * 0.5) : 30)), rationale: length >= 18 ? "Central gap of DNA flanked by modified wings" : "Too short for typical gapmer design" };
      modScores.psBackbone = { score: Math.max(0, Math.min(100, Math.round(length >= 12 ? 50 + length * 1.5 : 20))), rationale: "PS bonds increase nuclease resistance and protein binding" };
      modScores.cellUptake = { score: Math.max(0, Math.min(100, Math.round(90 - Math.abs(length - 20) * 3))), rationale: "18-22 nt optimal for cellular uptake of ASOs" };
    } else if (modality === "sirna") {
      const fgc = _gc(seq.slice(0, Math.floor(length / 2)));
      const sgc = _gc(seq.slice(Math.floor(length / 2)));
      modScores.thermodynamicBias = { score: Math.max(0, Math.min(100, Math.round(50 + Math.abs(fgc - sgc) * 2))), rationale: fgc < sgc ? "5'-end less stable — favors guide strand loading" : "Consider strand polarity" };
      modScores.riscLoading = { score: 30 <= gc && gc <= 52 ? 90 : 50, rationale: "GC 30-52% optimal for RISC loading efficiency" };
      modScores.specificity = { score: Math.max(0, Math.min(100, Math.round(length * 4))), rationale: "19-25 nt length balances potency and specificity" };
    } else if (modality === "mrna") {
      modScores.capEfficiency = { score: 70, rationale: "5' cap required for ribosome recruitment" };
      modScores.polyAStability = { score: /A{10,}$/.test(seq) ? 85 : 30, rationale: /A{10,}$/.test(seq) ? "Poly-A tail present" : "No poly-A tail detected" };
      modScores.mrnaStability = { score: 40 <= gc && gc <= 60 ? 80 : 50, rationale: "GC 40-60% optimal for mRNA half-life" };
      modScores.nucleosideMod = { score: 90, rationale: "m1Ψ and m5C modifications reduce innate immune activation" };
    } else if (modality === "sgrna") {
      modScores.gcOptimal = { score: 40 <= gc && gc <= 80 ? 90 : 40, rationale: "GC 40-80% optimal for sgRNA on-target activity" };
      modScores.pamProximity = { score: seq.slice(-5).includes("GG") ? 80 : 50, rationale: "NGG PAM motif detected near 3' end" };
      modScores.offTargetScore = { score: Math.max(0, Math.min(100, Math.round(length * 5))), rationale: "20 nt spacer optimal for specificity" };
    }
    const overallScore = Math.round(Object.values(modScores).reduce((s, v) => s + v.score, 0) / Math.max(Object.keys(modScores).length, 1));

    // --- Stacking energy profile ---
    const STACK_ENERGY: Record<string, number> = {
      "AA": -1.0, "TT": -1.0, "AT": -0.88, "TA": -0.58, "CA": -1.45, "TG": -1.45,
      "GT": -1.44, "AC": -1.44, "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
      "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
    };
    const nnSeq2 = seq.replace(/U/g, "T");
    const energyProfile: { position: number; energy: number }[] = [];
    const eWindow = 10, eStep = 2;
    for (let i = 0; i <= length - eWindow; i += eStep) {
      const chunk = nnSeq2.slice(i, i + eWindow);
      let sum = 0;
      for (let j = 0; j < chunk.length - 1; j++) sum += STACK_ENERGY[chunk.slice(j, j + 2)] ?? -1.0;
      energyProfile.push({ position: i + 1, energy: Math.round(sum / (chunk.length - 1) * 1000) / 1000 });
    }

    // === NEW FIELD CALCULATIONS ===

    // --- Restriction enzyme sites ---
    const RESTRICTION_ENZYMES: { enzyme: string; site: string; overhang: "5'" | "3'" | "blunt" }[] = [
      { enzyme: "EcoRI", site: "GAATTC", overhang: "5'" },
      { enzyme: "BamHI", site: "GGATCC", overhang: "5'" },
      { enzyme: "HindIII", site: "AAGCTT", overhang: "5'" },
      { enzyme: "NotI", site: "GCGGCCGC", overhang: "5'" },
      { enzyme: "XhoI", site: "CTCGAG", overhang: "5'" },
      { enzyme: "SacI", site: "GAGCTC", overhang: "3'" },
      { enzyme: "KpnI", site: "GGTACC", overhang: "3'" },
      { enzyme: "SpeI", site: "ACTAGT", overhang: "5'" },
      { enzyme: "NdeI", site: "CATATG", overhang: "5'" },
      { enzyme: "SmaI", site: "CCCGGG", overhang: "blunt" },
      { enzyme: "XbaI", site: "TCTAGA", overhang: "5'" },
      { enzyme: "PstI", site: "CTGCAG", overhang: "3'" },
      { enzyme: "SalI", site: "GTCGAC", overhang: "5'" },
      { enzyme: "Apal", site: "GGGCCC", overhang: "3'" },
      { enzyme: "NcoI", site: "CCATGG", overhang: "5'" },
    ];
    const restrictionSites: { enzyme: string; recognitionSite: string; cutPosition: number; strand: "+" | "-"; overhang: "5'" | "3'" | "blunt" }[] = [];
    for (const enz of RESTRICTION_ENZYMES) {
      const site = enz.site;
      for (let i = 0; i <= length - site.length; i++) {
        if (seq.slice(i, i + site.length) === site) {
          restrictionSites.push({ enzyme: enz.enzyme, recognitionSite: site, cutPosition: i + Math.floor(site.length / 2) + 1, strand: "+", overhang: enz.overhang });
        }
      }
    }

    // --- miRNA seed-motif matches -------------------------------------
    // These are exact string matches for common seed hexamers. They are NOT
    // miRNA target predictions:
    //   * the motifs are a fixed generic list, not a miRBase query
    //   * the previous "miR-seed-N" labels were invented identifiers that
    //     read as specific miRNAs
    //   * the previous bindingScore was
    //       0.5 + gc/200 + Math.random() * 0.15
    //     i.e. mostly noise, rendered as a percentage bar and sorted on
    //
    // Real target prediction needs TargetScan context++ scoring (feature F7),
    // which is not wired. What is reported here is exactly what was measured:
    // a motif, where it occurs, and its GC content.
    const SEED_HEXA_MERS = [
      "AACCCU", "AGCACCA", "GGAGCUA", "UAAGGCA", "CUCCAGA",
      "GAGGUUG", "UGCACUU", "AACAGUC", "GGCUGCA", "UCUACAG",
      "AAUGCCC", "UUCCGGA", "CCAGUGA", "GGCUGAU", "AUUGCCU",
    ];
    const mirnaTargets: {
      mirnaId: string;
      seedSequence: string;
      start: number;
      end: number;
      bindingScore: number | null;
      seedGcContent: number;
      conservationNote: string;
    }[] = [];
    for (let si = 0; si < SEED_HEXA_MERS.length && mirnaTargets.length < 15; si++) {
      const seed = SEED_HEXA_MERS[si];
      for (let i = 0; i <= length - seed.length; i++) {
        if (seq.slice(i, i + seed.length) === seed) {
          mirnaTargets.push({
            mirnaId: `seed motif ${si + 1}`,
            seedSequence: seed,
            start: i + 1,
            end: i + seed.length,
            bindingScore: null,
            seedGcContent: Math.round(_gc(seed) * 100) / 100,
            conservationNote:
              "Exact seed-motif match only. Not a miRNA target prediction — " +
              "no TargetScan context++ scoring is wired (feature F7).",
          });
          break;
        }
      }
    }

    // --- Hairpin structures ---
    const hairpins: { start: number; end: number; stemLength: number; loopSize: number; stabilityScore: number; type: "hairpin" | "bulge" | "internal_loop" }[] = [];
    for (let i = 0; i < length - 12; i += 3) {
      for (let j = i + 12; j < Math.min(i + 40, length); j += 2) {
        const seg = seq.slice(i, j + 1);
        const half = Math.floor(seg.length / 2);
        const stem1 = seg.slice(0, half);
        const stem2 = seg.slice(seg.length - half);
        const rc = revComp(stem2);
        let matches = 0;
        for (let k = 0; k < Math.min(stem1.length, rc.length); k++) {
          if (stem1[k] === rc[k]) matches++;
          else break;
        }
        if (matches >= 4) {
          const loopStart = i + matches;
          const loopEnd = j - matches;
          const loopSize = Math.max(1, loopEnd - loopStart + 1);
          const score = Math.round((matches / half) * (1 - loopSize / 20) * 100) / 100;
          if (score > 0.3 && hairpins.length < 10) {
            hairpins.push({
              start: i + 1,
              end: j + 1,
              stemLength: matches,
              loopSize,
              stabilityScore: Math.max(0, Math.min(1, score)),
              type: loopSize <= 4 ? "hairpin" : loopSize <= 8 ? "bulge" : "internal_loop",
            });
          }
        }
      }
    }
    // Deduplicate overlapping hairpins
    const uniqueHairpins = hairpins.filter((h, i) => {
      return !hairpins.slice(0, i).some((prev) => Math.abs(prev.start - h.start) < 5 && Math.abs(prev.end - h.end) < 5);
    });

    // --- Kmer frequency (k=6) ---
    const km = 6;
    const kmerMap: Record<string, number[]> = {};
    for (let i = 0; i <= length - km; i++) {
      const kmer = seq.slice(i, i + km);
      if (!kmerMap[kmer]) kmerMap[kmer] = [];
      kmerMap[kmer].push(i + 1);
    }
    const totalKmers = Math.max(length - km + 1, 1);
    const uniqueKmers = Object.keys(kmerMap).length;
    const shannonH = Object.values(kmerMap).reduce((sum, positions) => {
      const p = positions.length / totalKmers;
      return sum - p * Math.log2(p);
    }, 0);
    const kmerRepeats = Object.entries(kmerMap)
      .filter(([, positions]) => positions.length > 1)
      .sort((a, b) => b[1].length - a[1].length)
      .slice(0, 20)
      .map(([kmer, positions]) => ({ kmer, count: positions.length, positions }));

    // --- Thermodynamic profile ---
    const avgEnthalpy = energyProfile.length > 0
      ? Math.round(energyProfile.reduce((s, e) => s + e.energy, 0) / energyProfile.length * 1000) / 1000
      : -1.0;
    const avgEntropy = Math.round(avgEnthalpy * 30 * 1000) / 1000;
    const freeEnergy37 = Math.round((avgEnthalpy - 310.15 * avgEntropy / 1000) * 100) / 100;
    const stabilityClass = freeEnergy37 < -10 ? "stable" as const : freeEnergy37 < -3 ? "moderate" as const : "unstable" as const;
    const thermoProfile = {
      avgEnthalpy,
      avgEntropy,
      freeEnergy37,
      gcEnrichment: gc,
      atEnrichment: 100 - gc,
      stabilityClass,
      notes: [
        stabilityClass === "stable" ? "Sequence has favorable free energy for duplex formation." : stabilityClass === "unstable" ? "Low duplex stability — may require chemical modifications." : "Moderate thermodynamic stability.",
      ],
    };

    // --- Dot plot (self-complementarity k=6) ---
    const dotPlotPoints: { x: number; y: number; matchLen: number }[] = [];
    for (let i = 0; i <= length - k; i += 2) {
      for (let j = i + k; j <= length - k; j += 2) {
        const fwd = seq.slice(i, i + k);
        const rev = seq.slice(j, j + k);
        if (fwd === rev || fwd === revComp(rev)) {
          let matchLen = 0;
          for (let m = 0; m < Math.min(k, length - Math.max(i, j)); m++) {
            if (seq[i + m] === seq[j + m] || seq[i + m] === ({ A: "T", T: "A", G: "C", C: "G", U: "A" }[seq[j + m]] ?? "N")) {
              matchLen++;
            } else break;
          }
          if (matchLen >= 4 && Math.abs(i - j) > 5) {
            dotPlotPoints.push({ x: i + 1, y: j + 1, matchLen });
          }
        }
      }
    }
    const uniqueDotPlot = dotPlotPoints.filter((p, idx) => {
      return !dotPlotPoints.slice(0, idx).some((prev) => Math.abs(prev.x - p.x) < 3 && Math.abs(prev.y - p.y) < 3);
    }).slice(0, 200);

    // --- Modification landscape ---
    const modLandscape: { position: number; accessibilityScore: number; recommendedModification: string; confidenceLevel: "high" | "medium" | "low" }[] = [];
    const windowMod = 6;
    for (let i = 0; i <= length - windowMod; i += 3) {
      const chunk = seq.slice(i, i + windowMod);
      const localGc = _gc(chunk);
      const localG = (chunk.match(/G/g) || []).length;
      const accessibility = Math.round((localGc / 100 * 0.6 + (1 - Math.abs(localGc - 50) / 50) * 0.4) * 100) / 100;
      const recMod = modality === "aso"
        ? (localGc < 35 ? "LNA" : localGc < 55 ? "2'-MOE" : "PS")
        : modality === "sirna"
        ? (localGc < 40 ? "2'-OMe" : "PS")
        : modality === "mrna"
        ? (localG > 2 ? "N1-methylpseudouridine" : "5-methylcytidine")
        : (localGc < 40 ? "2'-OMe" : "LNA");
      const conf = Math.abs(localGc - 50) < 15 ? "high" as const : Math.abs(localGc - 50) < 25 ? "medium" as const : "low" as const;
      modLandscape.push({ position: i + 1, accessibilityScore: accessibility, recommendedModification: recMod, confidenceLevel: conf });
    }

    // --- Risk scores ---
    const specificityScore = Math.max(0, Math.min(100, length >= 18 ? 70 + (length - 18) * 2 : length * 3 + 10));
    const stabilityScore = Math.max(0, Math.min(100, Math.round(50 + (gc - 50) * 0.8 - Math.abs(length - 20) * 1.5)));
    const immunogenicityScore = Math.max(0, Math.min(100, 100 - immune.length * 8));
    const deliveryScore = Math.max(0, Math.min(100, length >= 15 && length <= 25 ? 80 : length < 15 ? 40 : 50));
    const toxicityScore = Math.max(0, Math.min(100, 80 - (palindromes > 3 ? 20 : 0) - (immune.length > 5 ? 15 : 0)));
    const riskScores = {
      specificity: Math.round(specificityScore),
      stability: Math.round(stabilityScore),
      immunogenicity: Math.round(immunogenicityScore),
      delivery: Math.round(deliveryScore),
      toxicity: Math.round(toxicityScore),
      overall: Math.round((specificityScore + stabilityScore + immunogenicityScore + deliveryScore + toxicityScore) / 5),
    };

    // --- Physicochemical profile ---
    const BASE_MW: Record<string, number> = { A: 331.2, C: 307.2, G: 347.2, T: 322.2, U: 324.2 };
    const BASE_CHARGE: Record<string, number> = { A: -1, C: -1, G: -1, T: -1, U: -1 };
    const BASE_HYDRO: Record<string, number> = { A: 0.5, C: -0.2, G: 0.1, T: 0.8, U: 0.3 };
    let mw = 0;
    let charge = 0;
    let hydroSum = 0;
    const hydroProfile: { position: number; value: number }[] = [];
    const chargeProfileArr: { position: number; value: number }[] = [];
    const mwWindow = 6;
    for (let i = 0; i <= length - mwWindow; i += 3) {
      const chunk = seq.slice(i, i + mwWindow);
      let chunkMw = 0;
      let chunkHydro = 0;
      for (const b of chunk) {
        chunkMw += BASE_MW[b] ?? 330;
        chunkHydro += BASE_HYDRO[b] ?? 0;
      }
      mw += chunkMw;
      charge += chunk.length * -1;
      hydroSum += chunkHydro;
      hydroProfile.push({ position: i + 1, value: Math.round(chunkHydro / chunk.length * 100) / 100 });
      chargeProfileArr.push({ position: i + 1, value: Math.round(chunk.length * -1 * 100) / 100 });
    }
    const physicochemical = {
      molecularWeight: Math.round(mw * 100) / 100,
      netCharge: charge,
      hydrophobicityIndex: Math.round(hydroSum / Math.max(length, 1) * 100) / 100,
      hydrophobicityProfile: hydroProfile,
      chargeProfile: chargeProfileArr,
    };

    // --- Stability index ---
    const stabilityIndex: { position: number; rnaseH: number; duplexStability: number; singleStrandStability: number }[] = [];
    const siWindow = 8;
    for (let i = 0; i <= length - siWindow; i += 3) {
      const chunk = seq.slice(i, i + siWindow);
      const localGc = _gc(chunk);
      const rnaseH = Math.round((localGc / 100 * 0.6 + (chunk.includes("T") || chunk.includes("U") ? 0.3 : 0) + 0.1) * 1000) / 1000;
      const duplex = Math.round((localGc * 0.01 + 0.3) * 1000) / 1000;
      const ss = Math.round((1 - localGc / 100 * 0.4 + 0.2) * 1000) / 1000;
      stabilityIndex.push({ position: i + 1, rnaseH: Math.min(1, Math.max(0, rnaseH)), duplexStability: Math.min(1, Math.max(0, duplex)), singleStrandStability: Math.min(1, Math.max(0, ss)) });
    }

    const grnaCandidates = modality === "sgrna" ? generateGrnaCandidates(seq) : [];

    return {
      sequence: seq,
      sequenceType: seqType,
      length,
      gcContent: gc,
      offTarget: {
        lengthBasedRiskEstimate: offRisk,
        note: offNote,
        internalRepetitiveness: Math.round(repetitiveness * 1000) / 1000,
        recommendedMinLength: 18,
        disclaimer: "This is a length/repetitiveness heuristic only — it does not check the sequence against any real genome or transcriptome. Use a real alignment tool (e.g. BLAST) for actual off-target screening.",
      },
      secondaryStructure: { estimatedMfe: mfe, palindromicRegions: palindromes, palindromePositions: palindromePositions.slice(0, 50), gcContent: gc, hairpinRisk: palindromes > 3 ? "High" : palindromes > 1 ? "Medium" : "Low" },
      immuneScreen: immune,
      modality: modDetails,
      gcCurve,
      composition,
      orfs: orfs.slice(0, 20),
      meltingTemp: { tmNearestNeighbor: tmNN, tmBasicGC, length, gcContent: gc, method: "Nearest-neighbor (SantaLucia 1998) at 50 mM Na+, 250 µM oligo", note: "Estimates only — actual Tm depends on salt, DMSO, and oligo concentration. Validate experimentally." },
      complexity: { dinucRepeats, trinucRepeats: [], gcRichRegions, atRichRegions, selfComplementarity: [], complexityScore },
      codonUsage: { codons: codons.slice(0, 50), cai: caiCount > 0 ? Math.round(caiSum / caiCount * 1000) / 1000 : 0, rareCodons: rareCodons.slice(0, 20), totalCodons: caiCount, note: "CAI ranges 0-1; higher = more optimized for human expression." },
      modificationScores: { modality, scores: modScores, overallScore },
      energyProfile,
      restrictionSites,
      mirnaTargets,
      hairpins: uniqueHairpins,
      kmerFrequency: { k: km, totalKmers, uniqueKmers, repeats: kmerRepeats, shannonEntropy: Math.round(shannonH * 1000) / 1000 },
      thermoProfile,
      dotPlot: uniqueDotPlot,
      modificationLandscape: modLandscape,
      riskScores,
      physicochemical,
      stabilityIndex,
      grnaCandidates,
    };
  }

  async function handleValidate() {
    if (!rawInput.trim()) return;
    setLoading(true);
    setError(null);
    try {
      let result: ValidationReport | null = null;
      try {
        result = await validateSequence(rawInput, filename) as unknown as ValidationReport;
      } catch {
        result = clientSideValidate(rawInput, filename) as unknown as ValidationReport;
      }
      if (!result) {
        setError("Could not parse sequence.");
        setLoading(false);
        return;
      }
      setValidation(result);
      if (!result.valid) {
        setError(`Invalid characters found: ${result.invalidChars.join(", ")}`);
      }
      setStep("validate");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Validation failed.");
    } finally {
      setLoading(false);
    }
  }

  async function handleAnalyze() {
    if (!validation?.sequence || !selectedModality) return;
    setLoading(true);
    setProgress(0);
    setError(null);

    // Simulated progress steps: stalls at 90% until real work finishes
    const steps = [
      { pct: 12, delay: 200 },
      { pct: 28, delay: 350 },
      { pct: 45, delay: 400 },
      { pct: 62, delay: 350 },
      { pct: 78, delay: 300 },
      { pct: 88, delay: 250 },
      { pct: 90, delay: 200 },
    ];
    let cancelled = false;
    let timerIdx = 0;
    const advance = () => {
      if (cancelled || timerIdx >= steps.length) return;
      const step = steps[timerIdx++];
      setProgress(step.pct);
      if (timerIdx < steps.length) {
        setTimeout(advance, step.delay);
      }
    };
    setTimeout(advance, 150);

    try {
      let result: AnalysisReport | null = null;
      try {
        result = await analyzeSequence(validation.sequence, selectedModality as Modality) as unknown as AnalysisReport;
      } catch {
        result = clientSideAnalyze(validation.sequence, selectedModality) as unknown as AnalysisReport;
      }
      cancelled = true;
      if (!result) {
        setError("Analysis failed.");
        setLoading(false);
        setProgress(0);
        return;
      }
      setProgress(100);
      // Brief pause at 100% so user sees the bar complete
      await new Promise((r) => setTimeout(r, 350));
      setAnalysis(result);
      setStep("analysis");
      saveReport({
        step: "sequence_analysis",
        title: `Sequence Analysis: ${filename || "Uploaded Sequence"}`,
        geneSymbol: "",
        disease: "",
        summary: `Analyzed sequence (${rawInput.length} bp). GC: ${((result.gcContent || 0) * 100).toFixed(1)}%. Tm: ${result.meltingTemp?.tmNearestNeighbor?.toFixed(1) || "N/A"}°C.`,
        data: { gcContent: result.gcContent, tm: result.meltingTemp?.tmNearestNeighbor ?? null, length: rawInput.length },
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Analysis failed.");
    } finally {
      setLoading(false);
      setProgress(0);
    }
  }

  function handleReset() {
    setStep("upload");
    setRawInput("");
    setFilename(undefined);
    setValidation(null);
    setSelectedModality(null);
    setAnalysis(null);
    setError(null);
  }

  function toggleRec(key: string) {
    setExpandedRecs((prev) => ({ ...prev, [key]: !prev[key] }));
  }

  return (
    <div className="flex min-h-screen bg-[#F8FAFC]">
      <Sidebar />
      <div className="flex min-h-screen flex-1 flex-col">
        <Topbar />
        <main className="flex-1 space-y-5 px-6 py-6">
          {/* Header */}
          <Card className="flex items-center gap-3 px-5 py-3">
            <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-gradient-to-br from-violet-50 to-purple-50">
              <UploadCloud className="h-4.5 w-4.5 text-violet-500" />
            </span>
            <div>
              <p className="text-[13px] font-semibold text-slate-800">Upload Sequence</p>
              <p className="text-[12px] text-slate-500">
                Upload FASTA files or paste raw sequences for custom analysis
              </p>
            </div>
            {step !== "upload" && (
              <button
                onClick={handleReset}
                className="ml-auto text-[12.5px] font-medium text-brand hover:underline"
              >
                Start over
              </button>
            )}
          </Card>

          {/* Step indicators */}
          <div className="flex items-center gap-2 text-[11.5px] font-medium text-slate-400">
            {(["upload", "validate", "modality", "analysis"] as Step[]).map((s, i) => (
              <span key={s} className="flex items-center gap-1.5">
                {i > 0 && <span className="text-slate-300">→</span>}
                <span
                  className={`rounded-md px-2 py-0.5 ${
                    step === s
                      ? "bg-brand text-white"
                      : ["upload", "validate", "modality", "analysis"].indexOf(step) > i
                      ? "bg-emerald-100 text-emerald-700"
                      : "bg-slate-100 text-slate-400"
                  }`}
                >
                  {i + 1}. {s.charAt(0).toUpperCase() + s.slice(1)}
                </span>
              </span>
            ))}
          </div>

          {error && (
            <div className="flex items-center gap-2 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-[13px] text-red-600">
              <AlertCircle className="h-4 w-4 shrink-0" /> {error}
            </div>
          )}

          {/* ===== STEP 1: UPLOAD ===== */}
          {step === "upload" && (
            <Card>
              <SectionHeader step="1" title="Input Sequence" />
              <div className="px-6 pb-6 space-y-4">
                {/* Paste area */}
                <div>
                  <label className="mb-1.5 block text-[12.5px] font-medium text-slate-600">
                    Paste Sequence (FASTA or raw)
                  </label>
                  <textarea
                    value={rawInput}
                    onChange={(e) => setRawInput(e.target.value)}
                    placeholder={`>header optional\ndna or rna sequence here...\n\nOr paste raw:\nATGCGATCGATCGATCG...`}
                    rows={8}
                    className="w-full rounded-lg border border-slate-300 bg-white px-4 py-3 text-[13px] font-mono text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20 resize-y"
                  />
                </div>

                {/* File drop zone */}
                <div
                  onDragOver={(e) => e.preventDefault()}
                  onDrop={handleDrop}
                  onClick={() => fileInputRef.current?.click()}
                  className="flex cursor-pointer flex-col items-center gap-2 rounded-xl border-2 border-dashed border-slate-300 bg-slate-50 p-8 transition-colors hover:border-brand/40 hover:bg-brand/[0.02]"
                >
                  <UploadCloud className="h-8 w-8 text-slate-300" />
                  <p className="text-[13px] font-medium text-slate-500">
                    Drop a FASTA/GenBank file here, or click to browse
                  </p>
                  <p className="text-[11px] text-slate-400">
                    Supports .fasta, .fa, .txt, .genbank, .gb
                  </p>
                  <input
                    ref={fileInputRef}
                    type="file"
                    accept=".fasta,.fa,.txt,.genbank,.gb,.fastq,.fq"
                    onChange={handleFileUpload}
                    className="hidden"
                  />
                </div>

                {filename && rawInput && (
                  <div className="flex items-center gap-2 rounded-lg bg-emerald-50 px-3 py-2 text-[12px] text-emerald-600">
                    <FileText className="h-3.5 w-3.5" />
                    Loaded: {filename} ({rawInput.length.toLocaleString()} chars)
                  </div>
                )}

                <div className="flex justify-end">
                  <button
                    onClick={handleValidate}
                    disabled={!rawInput.trim() || loading}
                    className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                  >
                    {loading && <Loader2 className="h-4 w-4 animate-spin" />}
                    {loading ? "Validating..." : "Validate Sequence"}
                    <ArrowRight className="h-3.5 w-3.5" />
                  </button>
                </div>
              </div>
            </Card>
          )}

          {/* ===== STEP 2: VALIDATION ===== */}
          {step === "validate" && validation && (
            <Card>
              <SectionHeader
                step="2"
                title="Sequence Validation"
                right={
                  <button
                    onClick={() => setStep("upload")}
                    className="flex items-center gap-1 text-[12px] font-medium text-brand hover:underline"
                  >
                    <ArrowLeft className="h-3 w-3" /> Edit
                  </button>
                }
              />
              <div className="px-6 pb-6 space-y-4">
                {/* Status banner */}
                <div
                  className={`flex items-center gap-3 rounded-lg px-4 py-3 ${
                    validation.valid
                      ? "border border-emerald-200 bg-emerald-50"
                      : "border border-amber-200 bg-amber-50"
                  }`}
                >
                  {validation.valid ? (
                    <CheckCircle2 className="h-5 w-5 text-emerald-500" />
                  ) : (
                    <AlertCircle className="h-5 w-5 text-amber-500" />
                  )}
                  <div>
                    <p className={`text-[13px] font-medium ${validation.valid ? "text-emerald-700" : "text-amber-700"}`}>
                      {validation.valid ? "Sequence is valid" : "Sequence has issues"}
                    </p>
                    <p className={`text-[11.5px] ${validation.valid ? "text-emerald-600" : "text-amber-600"}`}>
                      {validation.valid
                        ? "All characters are valid IUPAC nucleotides"
                        : `Invalid characters: ${validation.invalidChars.join(", ")}`}
                    </p>
                  </div>
                </div>

                {/* Stats grid */}
                <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
                  <div className="rounded-lg border border-slate-100 bg-slate-50 p-3 text-center">
                    <p className="text-[11px] text-slate-400 uppercase tracking-wider">Length</p>
                    <p className="text-[18px] font-bold text-slate-800 mt-0.5">
                      {validation.length.toLocaleString()}
                    </p>
                    <p className="text-[10px] text-slate-400">nucleotides</p>
                  </div>
                  <div className="rounded-lg border border-slate-100 bg-slate-50 p-3 text-center">
                    <p className="text-[11px] text-slate-400 uppercase tracking-wider">Type</p>
                    <p className="text-[18px] font-bold text-slate-800 mt-0.5 uppercase">
                      {validation.sequenceType}
                    </p>
                    <p className="text-[10px] text-slate-400">detected</p>
                  </div>
                  <div className="rounded-lg border border-slate-100 bg-slate-50 p-3 text-center">
                    <p className="text-[11px] text-slate-400 uppercase tracking-wider">GC Content</p>
                    <p className="text-[18px] font-bold text-slate-800 mt-0.5">
                      {validation.gcContent != null ? `${validation.gcContent}%` : "N/A"}
                    </p>
                    <p className="text-[10px] text-slate-400">
                      {validation.gcContent != null && validation.gcContent >= 40 && validation.gcContent <= 60 ? "optimal" : validation.gcContent != null ? "suboptimal" : "—"}
                    </p>
                  </div>
                  <div className="rounded-lg border border-slate-100 bg-slate-50 p-3 text-center">
                    <p className="text-[11px] text-slate-400 uppercase tracking-wider">Features</p>
                    <p className="text-[18px] font-bold text-slate-800 mt-0.5">
                      {validation.features.length}
                    </p>
                    <p className="text-[10px] text-slate-400">detected</p>
                  </div>
                </div>

                {/* Features */}
                {validation.features.length > 0 && (
                  <div>
                    <p className="mb-2 text-[12.5px] font-medium text-slate-600">Detected Features</p>
                    <div className="flex flex-wrap gap-1.5">
                      {validation.features.map((f) => (
                        <span key={f} className="inline-flex items-center rounded-full bg-indigo-50 px-2.5 py-1 text-[11px] font-medium text-indigo-600">
                          {f}
                        </span>
                      ))}
                    </div>
                  </div>
                )}

                {/* ORFs */}
                {validation.orfs.length > 0 && (
                  <div>
                    <p className="mb-2 text-[12.5px] font-medium text-slate-600">Open Reading Frames ({validation.orfs.length} found, both strands)</p>
                    <div className="space-y-1.5">
                      {validation.orfs.map((orf, i) => (
                        <div key={`${orf.frame}-${orf.start}-${i}`} className="flex items-center gap-3 rounded-lg bg-slate-50 px-3 py-2 text-[12px]">
                          <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-semibold ${
                            orf.strand === "-" ? "bg-amber-100 text-amber-700" : "bg-blue-100 text-blue-700"
                          }`}>
                            {orf.strand} strand
                          </span>
                          <span className="font-mono text-slate-500">Frame {orf.frame}</span>
                          <span className="text-slate-400">|</span>
                          <span className="text-slate-600">
                            {orf.start}–{orf.end} ({orf.proteinLength} aa)
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {/* Preview */}
                <div>
                  <p className="mb-1.5 text-[12.5px] font-medium text-slate-600">Sequence Preview</p>
                  <p className="max-h-20 overflow-y-auto rounded-lg bg-slate-50 p-3 font-mono text-[11px] text-slate-500 break-all">
                    {validation.sequence.slice(0, 200)}
                    {validation.sequence.length > 200 && "..."}
                  </p>
                </div>

                <div className="flex justify-end">
                  <button
                    onClick={() => setStep("modality")}
                    className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm hover:bg-brand-dark"
                  >
                    Choose Modality <ArrowRight className="h-3.5 w-3.5" />
                  </button>
                </div>
              </div>
            </Card>
          )}

          {/* ===== STEP 3: MODALITY ===== */}
          {step === "modality" && (
            <Card>
              <SectionHeader
                step="3"
                title="Select Therapeutic Modality"
                right={
                  <button
                    onClick={() => setStep("validate")}
                    className="flex items-center gap-1 text-[12px] font-medium text-brand hover:underline"
                  >
                    <ArrowLeft className="h-3 w-3" /> Back
                  </button>
                }
              />
              <div className="px-6 pb-6 space-y-4">
                <p className="text-[12.5px] text-slate-500">
                  Choose the therapeutic modality for your uploaded sequence.
                </p>

                <div className="grid grid-cols-1 gap-3 sm:grid-cols-2">
                  {MODALITIES.map((m) => {
                    const Icon = m.icon;
                    const isSelected = selectedModality === m.id;
                    return (
                      <button
                        key={m.id}
                        onClick={() => setSelectedModality(m.id)}
                        className={`flex items-start gap-3 rounded-xl border p-4 text-left transition-colors ${
                          isSelected
                            ? "border-brand bg-brand/5 ring-1 ring-brand"
                            : "border-[#E5E7EB] bg-white hover:border-slate-300 hover:bg-slate-50"
                        }`}
                      >
                        <span
                          className={`flex h-9 w-9 shrink-0 items-center justify-center rounded-lg ${
                            isSelected
                              ? "bg-brand text-white"
                              : "bg-slate-100 text-slate-400"
                          }`}
                        >
                          <Icon className="h-4 w-4" />
                        </span>
                        <div>
                          <p className="text-[13px] font-semibold text-slate-800">{m.name}</p>
                          <p className="text-[11.5px] text-slate-500 mt-0.5">{m.description}</p>
                        </div>
                      </button>
                    );
                  })}
                </div>

                <div className="flex justify-end">
                  <button
                    onClick={handleAnalyze}
                    disabled={!selectedModality || loading}
                    className="relative flex items-center gap-2 overflow-hidden rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                  >
                    {loading && (
                      <div
                        className="absolute inset-0 bg-white/20 transition-all duration-300 ease-out"
                        style={{ width: `${progress}%` }}
                      />
                    )}
                    <span className="relative flex items-center gap-2">
                      {loading ? (
                        <>
                          <Loader2 className="h-4 w-4 animate-spin" />
                          Analyzing… {progress}%
                        </>
                      ) : (
                        <>
                          Run Analysis
                          <ArrowRight className="h-3.5 w-3.5" />
                        </>
                      )}
                    </span>
                  </button>
                </div>
                {loading && (
                  <div className="mt-3">
                    <div className="h-1.5 w-full overflow-hidden rounded-full bg-slate-200">
                      <div
                        className="h-full rounded-full bg-brand transition-all duration-300 ease-out"
                        style={{ width: `${progress}%` }}
                      />
                    </div>
                    <p className="mt-1.5 text-[11px] text-slate-400">
                      {progress < 15
                        ? "Preparing sequence data…"
                        : progress < 50
                        ? "Computing thermodynamic properties…"
                        : progress < 75
                        ? "Analyzing modality-specific features…"
                        : progress < 95
                        ? "Running motif and complexity scans…"
                        : "Finalizing results…"}
                    </p>
                  </div>
                )}
              </div>
            </Card>
          )}

          {/* ===== STEP 4: ANALYSIS ===== */}
          {step === "analysis" && analysis && (
            <div className="space-y-5">
              <div className="flex items-center justify-between">
                <SectionHeader step="4" title="Analysis Results" />
                <div className="flex items-center gap-3">
                  <ExportMenu
                    sequence={analysis.sequence}
                    validation={{
                      valid: true,
                      sequenceType: (validation?.sequenceType ?? analysis.sequenceType) as "dna" | "rna" | "unknown",
                      length: validation?.length ?? analysis.length,
                      gcContent: validation?.gcContent ?? analysis.gcContent,
                      features: validation?.features ?? [],
                      orfs: (validation?.orfs ?? analysis.orfs).map((orf) => ({
                        ...orf,
                        strand: orf.strand as "+" | "-",
                      })),
                      invalidChars: [],
                      hasPolyA: validation?.hasPolyA ?? false,
                      hasPolyG: validation?.hasPolyG ?? false,
                    }}
                    analysis={analysis as unknown as AnalyzeResponse}
                    modalityName={MODALITIES.find((m) => m.id === selectedModality)?.name ?? selectedModality ?? ""}
                  />
                  <button
                    onClick={() => setStep("modality")}
                    className="flex items-center gap-1 text-[12px] font-medium text-brand hover:underline"
                  >
                    <ArrowLeft className="h-3 w-3" /> Change modality
                  </button>
                </div>
              </div>

              {/* Summary bar */}
              <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
                <div className="rounded-xl border border-[#E5E7EB] bg-white p-3 shadow-card">
                  <p className="text-[10px] uppercase tracking-wider text-slate-400">Sequence Type</p>
                  <p className="text-[15px] font-bold text-slate-800 mt-0.5 uppercase">{analysis.sequenceType}</p>
                </div>
                <div className="rounded-xl border border-[#E5E7EB] bg-white p-3 shadow-card">
                  <p className="text-[10px] uppercase tracking-wider text-slate-400">Length</p>
                  <p className="text-[15px] font-bold text-slate-800 mt-0.5">{analysis.length.toLocaleString()} nt</p>
                </div>
                <div className="rounded-xl border border-[#E5E7EB] bg-white p-3 shadow-card">
                  <p className="text-[10px] uppercase tracking-wider text-slate-400">GC Content</p>
                  <p className="text-[15px] font-bold text-slate-800 mt-0.5">{analysis.gcContent}%</p>
                </div>
                <div className="rounded-xl border border-[#E5E7EB] bg-white p-3 shadow-card">
                  <p className="text-[10px] uppercase tracking-wider text-slate-400">Modality</p>
                  <p className="text-[15px] font-bold text-slate-800 mt-0.5 uppercase">{MODALITIES.find((m) => m.id === selectedModality)?.name ?? selectedModality}</p>
                </div>
              </div>

              {/* Tabbed analysis interface */}
              <AnalysisTabs
                activeTab={activeTab}
                onTabChange={setActiveTab}
                tabs={[
                  { id: "overview", label: "Overview", icon: LayoutDashboard },
                  { id: "alignments", label: "Alignments", icon: AlignLeft },
                  {
                    id: "features",
                    label: "Features",
                    icon: List,
                    badge: (analysis.restrictionSites?.length ?? 0) + (analysis.mirnaTargets?.length ?? 0) + analysis.immuneScreen.length + (analysis.orfs?.length ?? 0),
                  },
                  { id: "structure", label: "Structure", icon: GitBranch },
                  { id: "properties", label: "Properties", icon: FlaskConical },
                ]}
              >
                {/* ===== OVERVIEW TAB ===== */}
                {activeTab === "overview" && (
                  <div className="p-5 space-y-5">
                     {/* Sequence Map */}
                     <div>
                       <p className="text-[14px] font-semibold text-slate-800 mb-3 flex items-center gap-1.5">
                         Sequence Map
                         <InfoTooltip content="Genomic coordinates of predicted features including ORFs, immune motifs, restriction sites, miRNA targets, and PAM sites for CRISPR." />
                       </p>
                      <SequenceTrackViewer
                        seqLength={analysis.length}
                        orfs={analysis.orfs}
                        immuneHits={analysis.immuneScreen}
                        palindromePositions={analysis.secondaryStructure.palindromePositions ?? []}
                        restrictionSites={analysis.restrictionSites}
                        mirnaTargets={analysis.mirnaTargets}
                        grnaCandidates={selectedModality === "sgrna" ? analysis.grnaCandidates : []}
                      />
                    </div>

                     {/* gRNA Candidates — CRISPRscan-style track + table */}
                    {selectedModality === "sgrna" && analysis.grnaCandidates && analysis.grnaCandidates.length > 0 && (
                      <div className="space-y-5">
                        <Card className="p-5">
                          <p className="text-[14px] font-semibold text-slate-800 mb-1 flex items-center gap-1.5">
                            CRISPRscan gRNA Score Track
                            <InfoTooltip content="Distribution of gRNA on-target scores (0–100) along the sequence for NGG-PAM adjacent 20nt spacers." />
                          </p>
                          <p className="text-[11px] text-slate-400 mb-3">Score distribution across the sequence for NGG-PAM adjacent spacers</p>
                          <CrisprScanTrack candidates={analysis.grnaCandidates} seqLength={analysis.length} />
                        </Card>

                        <Card className="p-5">
                          <p className="text-[14px] font-semibold text-slate-800 mb-3 flex items-center gap-1.5">
                            gRNA Candidate Table
                            <InfoTooltip content="Sortable table of guide RNAs with specificity (green >50, yellow >30, red), efficiency scores, and off-target mismatch distribution." />
                          </p>
                          <CrisprCandidateTable candidates={analysis.grnaCandidates} />
                        </Card>

                        {/* CRISPR Primer Design Card */}
                        <Card className="p-5">
                          <CrisprPrimerDesignCard
                            candidates={analysis.grnaCandidates ?? []}
                            sequence={analysis.sequence ?? ""}
                          />
                        </Card>
                      </div>
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                      {/* GC Content Chart */}
                       <Card className="p-5 lg:col-span-2">
                         <p className="text-[14px] font-semibold text-slate-800 mb-1 flex items-center gap-1.5">
                           GC Content Distribution
                           <InfoTooltip content="10 nt sliding window GC content. Green band = optimal 40–60% range for oligonucleotide stability." />
                         </p>
                        <p className="text-[11px] text-slate-400 mb-3">10 nt sliding window — green band marks the 40–60% optimal range</p>
                        <GcContentChart data={analysis.gcCurve} seqLength={analysis.length} />
                      </Card>

                       {/* Nucleotide Composition */}
                       {analysis.sequenceType !== "protein" && (
                         <Card className="p-5">
                           <p className="text-[14px] font-semibold text-slate-800 mb-3 flex items-center gap-1.5">
                             Nucleotide Composition
                             <InfoTooltip content="Absolute counts of A, C, G, T/U bases in the input sequence." />
                           </p>
                          <NucleotideCompositionChart composition={analysis.composition} />
                        </Card>
                      )}

                       {/* Specificity Heuristic */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Shield className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Specificity Heuristic
                             <InfoTooltip content="Heuristic estimate of off-target risk based on sequence length and internal repetitiveness. Not a genome-wide alignment." />
                           </p>
                         </div>
                        <div className="flex items-center gap-3 mb-2">
                          <span
                            className={`inline-flex items-center rounded-full px-2.5 py-1 text-[11px] font-semibold ${
                              analysis.offTarget.lengthBasedRiskEstimate === "Low"
                                ? "bg-emerald-100 text-emerald-700"
                                : analysis.offTarget.lengthBasedRiskEstimate === "Medium"
                                ? "bg-amber-100 text-amber-700"
                                : "bg-red-100 text-red-700"
                            }`}
                          >
                            {analysis.offTarget.lengthBasedRiskEstimate} Risk
                          </span>
                        </div>
                        <p className="text-[12px] text-slate-500 mb-2">{analysis.offTarget.note}</p>
                        <div className="rounded-lg bg-slate-50 p-3 mb-2">
                          <div className="flex justify-between text-[11px]">
                            <span className="text-slate-400">Internal repetitiveness</span>
                            <span className="font-medium text-slate-600">{(analysis.offTarget.internalRepetitiveness * 100).toFixed(1)}%</span>
                          </div>
                          <div className="mt-1 h-1.5 rounded-full bg-slate-200">
                            <div
                              className={`h-full rounded-full ${analysis.offTarget.internalRepetitiveness > 0.3 ? "bg-red-400" : "bg-emerald-400"}`}
                              style={{ width: `${Math.min(analysis.offTarget.internalRepetitiveness * 100, 100)}%` }}
                            />
                          </div>
                        </div>
                        <div className="flex items-start gap-1.5 rounded-lg bg-amber-50 px-3 py-2 text-[11.5px] text-amber-700">
                          <AlertCircle className="h-3.5 w-3.5 shrink-0 mt-0.5" />
                          {analysis.offTarget.disclaimer}
                        </div>
                      </Card>

                       {/* Secondary Structure */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Thermometer className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Secondary Structure
                             <InfoTooltip content="Estimated minimum free energy (MFE) and palindrome density. Hairpin risk rises with more palindromic regions." />
                           </p>
                         </div>
                        <div className="grid grid-cols-2 gap-3 mb-3">
                          <div className="rounded-lg bg-slate-50 p-3 text-center">
                            <p className="text-[10px] uppercase text-slate-400">Est. ΔG</p>
                            <p className="text-[16px] font-bold text-slate-800 mt-0.5">{analysis.secondaryStructure.estimatedMfe} kcal/mol</p>
                          </div>
                          <div className="rounded-lg bg-slate-50 p-3 text-center">
                            <p className="text-[10px] uppercase text-slate-400">Palindromes</p>
                            <p className="text-[16px] font-bold text-slate-800 mt-0.5">{analysis.secondaryStructure.palindromicRegions}</p>
                          </div>
                        </div>
                        <div className="flex items-center gap-2 mb-2">
                          <span className="text-[11px] text-slate-400">Hairpin risk:</span>
                          <span
                            className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-semibold ${
                              analysis.secondaryStructure.hairpinRisk === "Low"
                                ? "bg-emerald-100 text-emerald-700"
                                : analysis.secondaryStructure.hairpinRisk === "Medium"
                                ? "bg-amber-100 text-amber-700"
                                : "bg-red-100 text-red-700"
                            }`}
                          >
                            {analysis.secondaryStructure.hairpinRisk}
                          </span>
                        </div>
                        <div className="flex items-start gap-1.5 rounded-lg bg-slate-50 px-3 py-2 text-[11.5px] text-slate-500">
                          <Info className="h-3.5 w-3.5 shrink-0 mt-0.5" />
                          Composition-based MFE estimate, not a real folding prediction (e.g. RNAfold). Treat as a rough proxy only.
                        </div>
                      </Card>

                       {/* Immune Screen */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Syringe className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Immune Sensing Patterns
                             <InfoTooltip content="Sequence motifs associated with innate immune activation (TLR7/8, TLR9). A negative result does not rule out immunogenicity." />
                           </p>
                        </div>
                        {analysis.immuneScreen.length === 0 ? (
                          <div className="flex items-center gap-2 rounded-lg bg-emerald-50 px-3 py-2.5">
                            <CheckCircle2 className="h-4 w-4 text-emerald-500" />
                            <p className="text-[12px] text-emerald-600">No immunostimulatory motifs detected</p>
                          </div>
                        ) : (
                          <div className="space-y-1.5 max-h-40 overflow-y-auto">
                            {analysis.immuneScreen.map((m, i) => (
                              <div key={i} className="flex items-center gap-2 rounded-lg bg-amber-50 border border-amber-200 px-3 py-1.5 text-[11px]">
                                <span className="font-mono font-semibold text-amber-700 shrink-0">{m.motif}</span>
                                <span className="text-amber-600 shrink-0">@ {m.start}–{m.end}</span>
                                <span className="text-slate-500 truncate">{m.label}</span>
                              </div>
                            ))}
                          </div>
                        )}
                        <div className="mt-3 flex items-start gap-1.5 rounded-lg bg-amber-50 px-3 py-2 text-[11.5px] text-amber-700">
                          <AlertCircle className="h-3.5 w-3.5 shrink-0 mt-0.5" />
                          Pattern-matching against a short literature-informed list, not a validated immunogenicity assay. Positions are exact regex matches; biological labels are literature-informed guesses.
                        </div>
                      </Card>

                       {/* Secondary Structure Summary */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Thermometer className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Secondary Structure Summary
                             <InfoTooltip content="Summary of estimated thermodynamic stability metrics." />
                           </p>
                         </div>
                        <div className="grid grid-cols-3 gap-3 mb-3">
                          <div className="rounded-lg bg-slate-50 p-3 text-center">
                            <p className="text-[10px] uppercase text-slate-400">ΔG</p>
                            <p className="text-[14px] font-bold text-slate-800 mt-0.5">{analysis.secondaryStructure.estimatedMfe}</p>
                            <p className="text-[9px] text-slate-400">kcal/mol</p>
                          </div>
                          <div className="rounded-lg bg-slate-50 p-3 text-center">
                            <p className="text-[10px] uppercase text-slate-400">Hairpins</p>
                            <p className="text-[14px] font-bold text-slate-800 mt-0.5">{analysis.hairpins?.length ?? 0}</p>
                            <p className="text-[9px] text-slate-400">detected</p>
                          </div>
                          <div className="rounded-lg bg-slate-50 p-3 text-center">
                            <p className="text-[10px] uppercase text-slate-400">Palindrome</p>
                            <p className="text-[14px] font-bold text-slate-800 mt-0.5">{analysis.secondaryStructure.palindromicRegions}</p>
                            <p className="text-[9px] text-slate-400">regions</p>
                          </div>
                        </div>
                        {analysis.thermoProfile && (
                          <div className="flex items-center gap-3 text-[11px] text-slate-500">
                            <span className="font-medium">Stability:</span>
                            <span
                              className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-semibold ${
                                analysis.thermoProfile.stabilityClass === "stable"
                                  ? "bg-emerald-100 text-emerald-700"
                                  : analysis.thermoProfile.stabilityClass === "moderate"
                                  ? "bg-amber-100 text-amber-700"
                                  : "bg-red-100 text-red-700"
                              }`}
                            >
                              {analysis.thermoProfile.stabilityClass}
                            </span>
                            <span>ΔG free energy: {analysis.thermoProfile.freeEnergy37} kcal/mol</span>
                          </div>
                        )}
                       </Card>

                       {/* Sequence Risk Summary */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Shield className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Sequence Risk Summary
                             <InfoTooltip content="Composite risk scores across specificity, stability, immunogenicity, delivery, and toxicity dimensions." />
                           </p>
                         </div>
                         <div className="grid grid-cols-2 gap-2">
                           <div className="rounded-lg bg-slate-50 p-2.5">
                             <p className="text-[9px] uppercase text-slate-400">Specificity</p>
                             <div className="flex items-center gap-1.5 mt-1">
                               <div className="h-1.5 flex-1 rounded-full bg-slate-200">
                                 <div className="h-full rounded-full bg-blue-400" style={{ width: `${analysis.riskScores?.specificity ?? 50}%` }} />
                               </div>
                               <span className="text-[10px] font-bold text-slate-700">{analysis.riskScores?.specificity ?? "—"}</span>
                             </div>
                           </div>
                           <div className="rounded-lg bg-slate-50 p-2.5">
                             <p className="text-[9px] uppercase text-slate-400">Stability</p>
                             <div className="flex items-center gap-1.5 mt-1">
                               <div className="h-1.5 flex-1 rounded-full bg-slate-200">
                                 <div className="h-full rounded-full bg-emerald-400" style={{ width: `${analysis.riskScores?.stability ?? 50}%` }} />
                               </div>
                               <span className="text-[10px] font-bold text-slate-700">{analysis.riskScores?.stability ?? "—"}</span>
                             </div>
                           </div>
                           <div className="rounded-lg bg-slate-50 p-2.5">
                             <p className="text-[9px] uppercase text-slate-400">Immunogenicity</p>
                             <div className="flex items-center gap-1.5 mt-1">
                               <div className="h-1.5 flex-1 rounded-full bg-slate-200">
                                 <div className="h-full rounded-full bg-amber-400" style={{ width: `${analysis.riskScores?.immunogenicity ?? 50}%` }} />
                               </div>
                               <span className="text-[10px] font-bold text-slate-700">{analysis.riskScores?.immunogenicity ?? "—"}</span>
                             </div>
                           </div>
                           <div className="rounded-lg bg-slate-50 p-2.5">
                             <p className="text-[9px] uppercase text-slate-400">Delivery</p>
                             <div className="flex items-center gap-1.5 mt-1">
                               <div className="h-1.5 flex-1 rounded-full bg-slate-200">
                                 <div className="h-full rounded-full bg-purple-400" style={{ width: `${analysis.riskScores?.delivery ?? 50}%` }} />
                               </div>
                               <span className="text-[10px] font-bold text-slate-700">{analysis.riskScores?.delivery ?? "—"}</span>
                             </div>
                           </div>
                        </div>
                        {analysis.riskScores && (
                          <div className="mt-3 flex items-center justify-between rounded-lg bg-slate-50 px-3 py-2">
                            <span className="text-[11px] font-medium text-slate-500">Overall Risk Score</span>
                            <span className={`text-[14px] font-bold ${
                              analysis.riskScores.overall >= 70 ? "text-emerald-600" : analysis.riskScores.overall >= 40 ? "text-amber-600" : "text-red-600"
                            }`}>
                              {analysis.riskScores.overall}/100
                            </span>
                          </div>
                        )}
                        <div className="mt-2 flex items-start gap-1.5 rounded-lg bg-blue-50 px-3 py-2 text-[11px] text-blue-700">
                          <Info className="h-3 w-3 shrink-0 mt-0.5" />
                          Composite heuristic scores based on sequence composition, length, and motif analysis. Not a clinical risk assessment.
                        </div>
                      </Card>

                       {/* Modification Recommendations Summary */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <Star className="h-4 w-4 text-slate-500" />
                           <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                             Modification Recommendations
                             <InfoTooltip content="Per-position modification suggestions (LNA, 2′-MOE, PS, etc.) based on local GC content and modality." />
                           </p>
                        </div>
                        {analysis.modificationLandscape && analysis.modificationLandscape.length > 0 ? (
                          <div className="space-y-2">
                            {(() => {
                              const modCounts: Record<string, number> = {};
                              analysis.modificationLandscape.forEach((m) => {
                                modCounts[m.recommendedModification] = (modCounts[m.recommendedModification] ?? 0) + 1;
                              });
                              const sorted = Object.entries(modCounts).sort((a, b) => b[1] - a[1]).slice(0, 5);
                              const total = analysis.modificationLandscape.length;
                              return sorted.map(([mod, count]) => (
                                <div key={mod} className="flex items-center gap-2">
                                  <span className="text-[11px] font-mono text-slate-600 w-32 truncate">{mod}</span>
                                  <div className="flex-1 h-1.5 rounded-full bg-slate-200">
                                    <div className="h-full rounded-full bg-brand" style={{ width: `${(count / total) * 100}%` }} />
                                  </div>
                                  <span className="text-[10px] text-slate-400 w-10 text-right">{count}×</span>
                                </div>
                              ));
                            })()}
                          </div>
                        ) : (
                          <p className="text-[12px] text-slate-400">No modification landscape data available.</p>
                        )}
                        <div className="mt-3 flex items-start gap-1.5 rounded-lg bg-slate-50 px-3 py-2 text-[11px] text-slate-500">
                          <Info className="h-3.5 w-3.5 shrink-0 mt-0.5" />
                          Recommendations based on local GC content and modality-specific guidelines. See Properties tab for full landscape.
                        </div>
                      </Card>

                       {/* Modality Recommendations */}
                       <Card className="p-5">
                         <div className="flex items-center gap-2 mb-3">
                           <BarChart3 className="h-4 w-4 text-slate-500" />
                            <p className="text-[14px] font-semibold text-slate-800 flex items-center gap-1.5">
                              {MODALITIES.find((m) => m.id === selectedModality)?.name ?? "Modality"} Recommendations
                              <InfoTooltip content="Modality-specific design recommendations for chemistry, length, and target region." />
                            </p>
                          </div>
                          <div className="space-y-2">
                          {(analysis.modality.recommendations ?? []).map((rec, i) => (
                            <div key={i} className="flex items-start gap-2 rounded-lg bg-slate-50 px-3 py-2">
                              <Star className="h-3 w-3 mt-0.5 shrink-0 text-brand" />
                              <p className="text-[12px] text-slate-600">{rec}</p>
                            </div>
                          ))}
                        </div>
                        {analysis.modality.optimalLength && (
                          <div className="mt-3 flex items-center gap-2 text-[12px] text-slate-500">
                            <span className="font-medium">Optimal length:</span>
                            <span className="rounded-md bg-indigo-50 px-2 py-0.5 text-[11px] font-medium text-indigo-600">{analysis.modality.optimalLength}</span>
                          </div>
                        )}
                        {analysis.modality.recommendedChemistry && (
                          <div className="flex items-center gap-2 text-[12px] text-slate-500 mt-1">
                            <span className="font-medium">Chemistry:</span>
                            <span className="rounded-md bg-indigo-50 px-2 py-0.5 text-[11px] font-medium text-indigo-600">{analysis.modality.recommendedChemistry}</span>
                          </div>
                        )}
                        {analysis.modality.casProtein && (
                          <div className="flex items-center gap-2 text-[12px] text-slate-500 mt-1">
                            <span className="font-medium">Cas protein:</span>
                            <span className="rounded-md bg-purple-50 px-2 py-0.5 text-[11px] font-medium text-purple-600">{analysis.modality.casProtein}</span>
                          </div>
                        )}
                        {analysis.modality.nucleosideModifications && (
                          <div className="mt-3">
                            <p className="text-[11px] font-medium text-slate-500 mb-1">Suggested modifications:</p>
                            <div className="flex flex-wrap gap-1.5">
                              {analysis.modality.nucleosideModifications.map((mod) => (
                                <span key={mod} className="rounded-full bg-blue-50 px-2.5 py-0.5 text-[10px] font-medium text-blue-600">{mod}</span>
                              ))}
                            </div>
                          </div>
                        )}
                      </Card>
                    </div>

                    {/* Melting Temperature */}
                    {analysis.meltingTemp && (
                      <Card className="p-5">
                        <MeltingTemperatureCard tm={analysis.meltingTemp} />
                      </Card>
                    )}

                    {/* Modification Scores */}
                    {analysis.modificationScores && (
                      <Card className="p-5">
                        <ModificationScorecard scores={analysis.modificationScores} />
                      </Card>
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                      {/* Stacking Energy Profile */}
                      {analysis.energyProfile && analysis.energyProfile.length > 0 && (
                           <Card className="p-5">
                           <p className="text-[14px] font-semibold text-slate-800 mb-1 flex items-center gap-1.5">
                             Base-Stacking Energy Profile
                             <InfoTooltip content="Nearest-neighbor stacking ΔG (kcal/mol) averaged over a 10 nt sliding window." />
                           </p>
                           <p className="text-[11px] text-slate-400 mb-3">10 nt sliding window average nearest-neighbor ΔG (kcal/mol)</p>
                           <StackingEnergyChart data={analysis.energyProfile} seqLength={analysis.length} />
                         </Card>
                       )}

                       {/* Sequence Complexity */}
                       {analysis.complexity && (
                         <Card className="p-5">
                           <SequenceComplexityCard complexity={analysis.complexity} />
                         </Card>
                       )}
                    </div>

                    {/* Codon Usage (relevant for mRNA, show for all) */}
                    {analysis.codonUsage && analysis.codonUsage.totalCodons > 0 && (
                      <Card className="p-5">
                        <CodonUsageCard codonUsage={analysis.codonUsage} />
                      </Card>
                    )}

                    {/* Risk Score Dashboard */}
                    {analysis.riskScores && (
                      <Card className="p-5">
                        <RiskScoreDashboard riskScores={analysis.riskScores} />
                      </Card>
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                      {/* Restriction Enzyme Map */}
                      {analysis.restrictionSites && analysis.restrictionSites.length > 0 && (
                        <Card className="p-5">
                          <RestrictionSiteMap
                            sites={analysis.restrictionSites}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}

                      {/* miRNA Targeting */}
                      {analysis.mirnaTargets && analysis.mirnaTargets.length > 0 && (
                        <Card className="p-5">
                          <MiRNATargetingCard
                            targets={analysis.mirnaTargets}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}

                      {/* Thermodynamic Profile */}
                      {analysis.thermoProfile && (
                        <Card className="p-5">
                          <ThermodynamicProfile profile={analysis.thermoProfile} />
                        </Card>
                      )}

                      {/* Kmer Frequency */}
                      {analysis.kmerFrequency && (
                        <Card className="p-5">
                          <KmerFrequencyChart kmerData={analysis.kmerFrequency} />
                        </Card>
                      )}
                    </div>

                    {/* Hairpin Diagram */}
                    {analysis.hairpins && analysis.hairpins.length > 0 && (
                      <Card className="p-5">
                        <HairpinDiagram
                          hairpins={analysis.hairpins}
                          seqLength={analysis.length}
                        />
                      </Card>
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                      {/* Modification Landscape */}
                      {analysis.modificationLandscape && analysis.modificationLandscape.length > 0 && (
                        <Card className="p-5">
                          <ModificationLandscapeCard
                            landscape={analysis.modificationLandscape}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}

                      {/* Base Pair Dot Plot */}
                      {analysis.dotPlot && analysis.dotPlot.length > 0 && (
                        <Card className="p-5">
                          <BasePairDotPlot
                            points={analysis.dotPlot}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}

                      {/* Physicochemical Properties */}
                      {analysis.physicochemical && (
                        <Card className="p-5">
                          <PhysicochemicalCard profile={analysis.physicochemical} />
                        </Card>
                      )}

                      {/* Stability Index */}
                      {analysis.stabilityIndex && analysis.stabilityIndex.length > 0 && (
                        <Card className="p-5">
                          <StabilityIndexChart
                            data={analysis.stabilityIndex}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}
                    </div>

                    {/* Sequence Annotation Bar (full width) */}
                    <Card className="p-5">
                      <SequenceAnnotationBar
                        annotations={[
                          ...(analysis.restrictionSites ?? []).map((s, i) => ({
                            id: `restr-${i}`,
                            label: s.enzyme,
                            start: s.cutPosition,
                            end: s.cutPosition,
                            type: "restriction" as const,
                          })),
                          ...(analysis.mirnaTargets ?? []).map((t, i) => ({
                            id: `mirna-${i}`,
                            label: t.mirnaId,
                            start: t.start,
                            end: t.end,
                            type: "mirna" as const,
                          })),
                          ...(analysis.immuneScreen ?? []).slice(0, 15).map((m, i) => ({
                            id: `imm-${i}`,
                            label: m.motif,
                            start: m.start,
                            end: m.end,
                            type: "immune" as const,
                          })),
                          ...(analysis.orfs ?? []).map((o, i) => ({
                            id: `orf-${i}`,
                            label: `ORF ${o.strand} f${o.frame}`,
                            start: o.start,
                            end: o.end,
                            type: "orfs" as const,
                          })),
                          ...(analysis.complexity?.gcRichRegions ?? []).map((r, i) => ({
                            id: `gc-${i}`,
                            label: `GC-rich`,
                            start: r.start,
                            end: r.end,
                            type: "complexity" as const,
                          })),
                          ...(analysis.hairpins ?? []).map((h, i) => ({
                            id: `hp-${i}`,
                            label: h.type.replace("_", " "),
                            start: h.start,
                            end: h.end,
                            type: "structure" as const,
                          })),
                        ]}
                        seqLength={analysis.length}
                      />
                    </Card>
                  </div>
                )}

                {/* ===== ALIGNMENTS TAB ===== */}
                {activeTab === "alignments" && (
                  <div className="p-5 space-y-5">
                    <Card className="p-5">
                      <PairwiseAlignmentViewer
                        sequence={analysis.sequence}
                        features={[
                          ...(analysis.restrictionSites ?? []).map((s) => ({
                            start: s.cutPosition,
                            end: s.cutPosition + s.recognitionSite.length - 1,
                            label: s.enzyme,
                            type: "restriction",
                            color: "#ec4899",
                          })),
                          ...(analysis.mirnaTargets ?? []).map((t) => ({
                            start: t.start,
                            end: t.end,
                            label: t.mirnaId,
                            type: "mirna",
                            color: "#8b5cf6",
                          })),
                          ...(analysis.immuneScreen ?? []).slice(0, 10).map((m) => ({
                            start: m.start,
                            end: m.end,
                            label: m.motif,
                            type: "immune",
                            color: "#f59e0b",
                          })),
                          ...(analysis.orfs ?? []).slice(0, 5).map((o) => ({
                            start: o.start,
                            end: o.end,
                            label: `ORF ${o.strand} f${o.frame}`,
                            type: "orf",
                            color: "#3b82f6",
                          })),
                        ]}
                      />
                    </Card>

                    {/* Restriction Enzyme Map */}
                    {analysis.restrictionSites && analysis.restrictionSites.length > 0 && (
                      <Card className="p-5">
                        <RestrictionSiteMap
                          sites={analysis.restrictionSites}
                          seqLength={analysis.length}
                        />
                      </Card>
                    )}

                    {/* miRNA Targeting */}
                    {analysis.mirnaTargets && analysis.mirnaTargets.length > 0 && (
                      <Card className="p-5">
                        <MiRNATargetingCard
                          targets={analysis.mirnaTargets}
                          seqLength={analysis.length}
                        />
                      </Card>
                    )}
                  </div>
                )}

                {/* ===== FEATURES TAB ===== */}
                {activeTab === "features" && (
                  <div className="p-5 space-y-5">
                    <FeatureTable
                      restrictionSites={analysis.restrictionSites}
                      mirnaTargets={analysis.mirnaTargets}
                      immuneHits={analysis.immuneScreen}
                      orfs={analysis.orfs}
                      hairpins={analysis.hairpins}
                    />

                        <AnnotatedSequenceViewer
                          sequence={analysis.sequence}
                          features={[
                            ...(analysis.restrictionSites ?? []).map((s) => ({
                              start: s.cutPosition,
                              end: s.cutPosition + s.recognitionSite.length - 1,
                              label: s.enzyme,
                              type: "restriction",
                              color: "#ec4899",
                            })),
                            ...(analysis.mirnaTargets ?? []).map((t) => ({
                              start: t.start,
                              end: t.end,
                              label: t.mirnaId,
                              type: "mirna",
                              color: "#8b5cf6",
                            })),
                            ...(analysis.immuneScreen ?? []).slice(0, 10).map((m) => ({
                              start: m.start,
                              end: m.end,
                              label: m.motif,
                              type: "immune",
                              color: "#f59e0b",
                            })),
                            ...(analysis.orfs ?? []).slice(0, 5).map((o) => ({
                              start: o.start,
                              end: o.end,
                              label: `ORF ${o.strand} f${o.frame}`,
                              type: "orfs",
                              color: "#3b82f6",
                            })),
                            ...(analysis.hairpins ?? []).map((h) => ({
                              start: h.start,
                              end: h.end,
                              label: h.type.replace("_", " "),
                              type: "structure",
                              color: "#10b981",
                            })),
                            ...(analysis.grnaCandidates ?? []).map((g) => ({
                              start: g.position,
                              end: g.position + 22,
                              label: `PAM ${g.pam}`,
                              type: "pam",
                              color: "#f97316",
                            })),
                          ]}
                        />

                    {/* Modification Scores */}
                    {analysis.modificationScores && (
                      <Card className="p-5">
                        <ModificationScorecard scores={analysis.modificationScores} />
                      </Card>
                    )}

                    {/* Codon Usage */}
                    {analysis.codonUsage && analysis.codonUsage.totalCodons > 0 && (
                      <Card className="p-5">
                        <CodonUsageCard codonUsage={analysis.codonUsage} />
                      </Card>
                    )}
                  </div>
                )}

                {/* ===== STRUCTURE TAB ===== */}
                {activeTab === "structure" && (
                  <div className="p-5 space-y-5">
                    {/* Thermodynamic Profile */}
                    {analysis.thermoProfile && (
                      <Card className="p-5">
                        <ThermodynamicProfile profile={analysis.thermoProfile} />
                      </Card>
                    )}

                    {/* Hairpin Diagram */}
                    {analysis.hairpins && analysis.hairpins.length > 0 && (
                      <Card className="p-5">
                        <HairpinDiagram
                          hairpins={analysis.hairpins}
                          seqLength={analysis.length}
                        />
                      </Card>
                    )}

                    {/* Base Pair Dot Plot */}
                    {analysis.dotPlot && analysis.dotPlot.length > 0 && (
                      <Card className="p-5">
                        <BasePairDotPlot
                          points={analysis.dotPlot}
                          seqLength={analysis.length}
                        />
                      </Card>
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                       {/* Stacking Energy Profile (Structure tab) */}
                       {analysis.energyProfile && analysis.energyProfile.length > 0 && (
                         <Card className="p-5">
                           <p className="text-[14px] font-semibold text-slate-800 mb-1 flex items-center gap-1.5">
                             Base-Stacking Energy Profile
                             <InfoTooltip content="Nearest-neighbor stacking ΔG (kcal/mol) averaged over a 10 nt sliding window." />
                           </p>
                           <p className="text-[11px] text-slate-400 mb-3">10 nt sliding window average nearest-neighbor ΔG (kcal/mol)</p>
                           <StackingEnergyChart data={analysis.energyProfile} seqLength={analysis.length} />
                         </Card>
                       )}

                       {/* Stability Index */}
                       {analysis.stabilityIndex && analysis.stabilityIndex.length > 0 && (
                         <Card className="p-5">
                           <p className="text-[14px] font-semibold text-slate-800 mb-1 flex items-center gap-1.5">
                             Stability Index
                             <InfoTooltip content="Per-window RNase H cleavage potential and duplex/single-strand stability metrics." />
                           </p>
                           <StabilityIndexChart
                             data={analysis.stabilityIndex}
                             seqLength={analysis.length}
                           />
                         </Card>
                       )}
                    </div>

                    {/* Sequence Complexity */}
                    {analysis.complexity && (
                      <Card className="p-5">
                        <SequenceComplexityCard complexity={analysis.complexity} />
                      </Card>
                    )}

                    {/* Sequence Annotation Bar */}
                    <Card className="p-5">
                      <SequenceAnnotationBar
                        annotations={[
                          ...(analysis.restrictionSites ?? []).map((s, i) => ({
                            id: `restr-${i}`,
                            label: s.enzyme,
                            start: s.cutPosition,
                            end: s.cutPosition,
                            type: "restriction" as const,
                          })),
                          ...(analysis.mirnaTargets ?? []).map((t, i) => ({
                            id: `mirna-${i}`,
                            label: t.mirnaId,
                            start: t.start,
                            end: t.end,
                            type: "mirna" as const,
                          })),
                          ...(analysis.immuneScreen ?? []).slice(0, 15).map((m, i) => ({
                            id: `imm-${i}`,
                            label: m.motif,
                            start: m.start,
                            end: m.end,
                            type: "immune" as const,
                          })),
                          ...(analysis.orfs ?? []).map((o, i) => ({
                            id: `orf-${i}`,
                            label: `ORF ${o.strand} f${o.frame}`,
                            start: o.start,
                            end: o.end,
                            type: "orfs" as const,
                          })),
                          ...(analysis.complexity?.gcRichRegions ?? []).map((r, i) => ({
                            id: `gc-${i}`,
                            label: `GC-rich`,
                            start: r.start,
                            end: r.end,
                            type: "complexity" as const,
                          })),
                          ...(analysis.hairpins ?? []).map((h, i) => ({
                            id: `hp-${i}`,
                            label: h.type.replace("_", " "),
                            start: h.start,
                            end: h.end,
                            type: "structure" as const,
                          })),
                        ]}
                        seqLength={analysis.length}
                      />
                    </Card>
                  </div>
                )}

                {/* ===== PROPERTIES TAB ===== */}
                {activeTab === "properties" && (
                  <div className="p-5 space-y-5">
                    {/* Risk Score Dashboard */}
                    {analysis.riskScores && (
                      <RiskScoreDashboard riskScores={analysis.riskScores} />
                    )}

                    <div className="grid grid-cols-1 gap-5 lg:grid-cols-2">
                      {/* Modification Landscape */}
                      {analysis.modificationLandscape && analysis.modificationLandscape.length > 0 && (
                        <Card className="p-5">
                          <ModificationLandscapeCard
                            landscape={analysis.modificationLandscape}
                            seqLength={analysis.length}
                          />
                        </Card>
                      )}

                      {/* Kmer Frequency */}
                      {analysis.kmerFrequency && (
                        <Card className="p-5">
                          <KmerFrequencyChart kmerData={analysis.kmerFrequency} />
                        </Card>
                      )}

                      {/* Physicochemical Properties */}
                      {analysis.physicochemical && (
                        <Card className="p-5">
                          <PhysicochemicalCard profile={analysis.physicochemical} />
                        </Card>
                      )}
                    </div>
                  </div>
                )}
              </AnalysisTabs>

              {/* Action buttons */}
              <div className="flex justify-end gap-3">
                <button
                  onClick={handleReset}
                  className="rounded-lg border border-slate-300 px-4 py-2 text-[13px] font-medium text-slate-600 hover:bg-slate-50"
                >
                  Upload New Sequence
                </button>
                <button
                  onClick={() => setStep("modality")}
                  className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm hover:bg-brand-dark"
                >
                  Try Different Modality <ArrowRight className="h-3.5 w-3.5" />
                </button>
              </div>
            </div>
          )}
        </main>
      </div>
    </div>
  );
}
