"use client";

import { useState, useEffect } from "react";
import { RnaNeutralizationCandidate } from "@/types/rnaNeutralization";
import {
  CheckCircle,
  ChevronRight,
  ArrowLeft,
  Beaker,
  Shield,
  Wind,
} from "lucide-react";
import { useKeyboardShortcut } from "@/hooks/useKeyboardShortcut";

interface AsoDesignPipelineProps {
  candidate: RnaNeutralizationCandidate;
  deliveryContext: string;
  onComplete: (design: {
    candidate: RnaNeutralizationCandidate;
    modifications: string[];
    conjugation: string;
    secondaryStructurePassed: boolean;
    selfDimerPassed: boolean;
  }) => void;
  onBack: () => void;
}

export default function AsoDesignPipeline({
  candidate,
  deliveryContext,
  onComplete,
  onBack,
}: AsoDesignPipelineProps) {
  const [step, setStep] = useState(1);
  const [selectedModPattern, setSelectedModPattern] = useState<string>("");
  const [selectedConjugation, setSelectedConjugation] = useState<string>("");
  const [secondaryChecked, setSecondaryChecked] = useState(false);

  useKeyboardShortcut("escape", onBack);

  const chemistryModPatterns: Record<string, { id: string; label: string; description: string }[]> = {
    "2-o-moe-full-ps": [
      { id: "uniform-2moess", label: "Uniform 2′-O-MOE Full PS", description: "Standard pattern: every nucleotide carries 2′-O-MOE + phosphorothioate backbone. Maximizes nuclease resistance and RNase H independence." },
      { id: "alternating-2moess-2omeps", label: "Alternating 2′-MOE / 2′-OMe + PS", description: "Positions 1,3,5... 2′-O-MOE; positions 2,4,6... 2′-O-Me. Reduces synthesis cost while maintaining >90% nuclease resistance." },
    ],
    "pmo": [
      { id: "full-pmo", label: "Uniform Phosphorodiamidate Morpholino", description: "Full PMO backbone with no phosphate groups. Completely non-ionic; requires delivery conjugation for cellular uptake." },
    ],
    "lna-dna-mixmer": [
      { id: "lnana-flank-dna-core", label: "LNA Flanks / Modified Core", description: "3-nt LNA wings for Tm elevation, modified 2′-O-MOE core for RNase H independence. Locked sugar pucker enhances binding affinity." },
      { id: "uniform-lna-ps", label: "Uniform LNA/PS Mixmer", description: "Alternating LNA and DNA nucleotides with PS backbone. Each LNA adds ~2–8°C to Tm per modification." },
    ],
  };

  const deliveryConjugations: Record<string, { id: string; label: string; description: string }[]> = {
    muscle: [
      { id: "p-pmo", label: "P-PMO (Cell-Penetrating Peptide)", description: "Conjugated arginine-rich peptide enables efficient muscle uptake. Proven in DMD exon-skipping trials." },
      { id: "palmitic-acid", label: "Palmitic Acid Conjugation", description: "Lipid conjugation enhances membrane penetration in cardiac and skeletal muscle tissue." },
    ],
    heart: [
      { id: "p-pmo", label: "P-PMO (Cell-Penetrating Peptide)", description: "Optimal for cardiomyocyte uptake. Enhanced by tissue ischemia-mediated uptake." },
      { id: "palmitic-acid", label: "Palmitic Acid Conjugation", description: "Lipid-mediated uptake in cardiac myocytes. Supports systemic IV delivery." },
    ],
    cns: [
      { id: "unconjugated", label: "Unconjugated 2′-MOE ASO", description: "IT (intrathecal) delivery bypasses BBB. No conjugate needed for CNS parenchymal uptake." },
      { id: "cholesterol", label: "Cholesterol Conjugation", description: "Enhances BBB penetration and broad CNS distribution after systemic administration." },
    ],
    default: [
      { id: "unconjugated", label: "Unconjugated", description: "Relies on naked ASO uptake via gymnosis or endogenous uptake pathways." },
      { id: "cholesterol", label: "Cholesterol Conjugation", description: "Broad tissue distribution enhancement; supports systemic delivery." },
    ],
  };

  const currentMods = chemistryModPatterns[candidate.chemistry] ?? chemistryModPatterns["2-o-moe-full-ps"];
  const currentConjugations = deliveryConjugations[deliveryContext] ?? deliveryConjugations["default"];

  function handleComplete() {
    onComplete({
      candidate,
      modifications: selectedModPattern ? [selectedModPattern] : [],
      conjugation: selectedConjugation,
      secondaryStructurePassed: secondaryChecked,
      selfDimerPassed: candidate.realMetrics.selfStructureMfe > -3,
    });
  }

  return (
    <div className="space-y-5">
      <div className="flex items-center gap-3">
        <button
          onClick={onBack}
          className="flex items-center gap-1 rounded-lg border border-slate-200 px-3 py-1.5 text-[12px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
        >
          <ArrowLeft className="h-3.5 w-3.5" /> Back to Results
        </button>
        <div className="flex items-center gap-2 text-[12px] text-slate-500">
          <span className="font-medium text-slate-700">Candidate #{candidate.rank}</span>
          <span>·</span>
          <span className="font-mono text-[11px]">{candidate.sequence}</span>
        </div>
      </div>

      {/* Step indicators */}
      <div className="flex items-center gap-2">
        {[1, 2, 3].map((s) => (
          <div key={s} className="flex items-center gap-2">
            <button
              onClick={() => (s < step ? setStep(s) : null)}
              className={`flex h-7 w-7 items-center justify-center rounded-full text-[11px] font-bold transition-colors ${
                s === step
                  ? "bg-brand text-white"
                  : s < step
                    ? "bg-emerald-100 text-emerald-700"
                    : "bg-slate-100 text-slate-400"
              }`}
            >
              {s < step ? <CheckCircle className="h-3.5 w-3.5" /> : s}
            </button>
            <span className={`text-[12px] font-medium ${s === step ? "text-slate-700" : "text-slate-400"}`}>
              {s === 1 ? "Chem Mods" : s === 2 ? "Delivery" : "Structure Filter"}
            </span>
            {s < 3 && <ChevronRight className="h-3.5 w-3.5 text-slate-300" />}
          </div>
        ))}
      </div>

      {/* Step 1: Chemical Modification Optimization */}
      {step === 1 && (
        <div className="space-y-4">
          <div className="flex items-center gap-2">
            <Beaker className="h-4 w-4 text-brand" />
            <h3 className="text-[14px] font-semibold text-slate-800">
              Step 1: Chemical Modification Optimization
            </h3>
          </div>
          <p className="text-[12px] text-slate-500">
            Apply full-length chemical modifications to prevent nuclease degradation while maintaining non-ionic/steric properties.
            The selected chemistry ({candidate.chemistry}) supports RNase H-independent steric blockade.
          </p>

          <div className="space-y-2">
            {currentMods.map((mod) => (
              <button
                key={mod.id}
                onClick={() => setSelectedModPattern(mod.id)}
                className={`w-full rounded-xl border p-4 text-left transition-colors ${
                  selectedModPattern === mod.id
                    ? "border-brand bg-brand/5 ring-2 ring-brand/20"
                    : "border-slate-200 hover:border-slate-300 hover:bg-slate-50"
                }`}
              >
                <p className="text-[13px] font-semibold text-slate-700">{mod.label}</p>
                <p className="mt-1 text-[11px] text-slate-500 leading-relaxed">{mod.description}</p>
              </button>
            ))}
          </div>

          <div className="flex justify-end">
            <button
              onClick={() => setStep(2)}
              disabled={!selectedModPattern}
              className="flex items-center gap-1.5 rounded-lg bg-brand px-4 py-2 text-[12px] font-medium text-white shadow-sm hover:bg-brand-dark transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
            >
              Continue <ChevronRight className="h-3.5 w-3.5" />
            </button>
          </div>
        </div>
      )}

      {/* Step 2: Tissue Delivery Conjugation */}
      {step === 2 && (
        <div className="space-y-4">
          <div className="flex items-center gap-2">
            <Shield className="h-4 w-4 text-brand" />
            <h3 className="text-[14px] font-semibold text-slate-800">
              Step 2: Tissue Delivery Conjugation Selection
            </h3>
          </div>
          <p className="text-[12px] text-slate-500">
            Based on the <strong>{deliveryContext || "systemic"}</strong> delivery context, select the optimal conjugation
            strategy for cellular uptake in the target tissue.
          </p>

          <div className="space-y-2">
            {currentConjugations.map((conj) => (
              <button
                key={conj.id}
                onClick={() => setSelectedConjugation(conj.id)}
                className={`w-full rounded-xl border p-4 text-left transition-colors ${
                  selectedConjugation === conj.id
                    ? "border-brand bg-brand/5 ring-2 ring-brand/20"
                    : "border-slate-200 hover:border-slate-300 hover:bg-slate-50"
                }`}
              >
                <p className="text-[13px] font-semibold text-slate-700">{conj.label}</p>
                <p className="mt-1 text-[11px] text-slate-500 leading-relaxed">{conj.description}</p>
              </button>
            ))}
          </div>

          <div className="flex justify-between">
            <button
              onClick={() => setStep(1)}
              className="flex items-center gap-1.5 rounded-lg border border-slate-200 px-4 py-2 text-[12px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
            >
              <ArrowLeft className="h-3.5 w-3.5" /> Back
            </button>
            <button
              onClick={() => setStep(3)}
              disabled={!selectedConjugation}
              className="flex items-center gap-1.5 rounded-lg bg-brand px-4 py-2 text-[12px] font-medium text-white shadow-sm hover:bg-brand-dark transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
            >
              Continue <ChevronRight className="h-3.5 w-3.5" />
            </button>
          </div>
        </div>
      )}

      {/* Step 3: Secondary Structure & Self-Dimer Filter */}
      {step === 3 && (
        <div className="space-y-4">
          <div className="flex items-center gap-2">
            <Wind className="h-4 w-4 text-brand" />
            <h3 className="text-[14px] font-semibold text-slate-800">
              Step 3: Secondary Structure &amp; Self-Dimer Filter
            </h3>
          </div>
          <p className="text-[12px] text-slate-500">
            Run the sequence through RNAfold/ViennaRNA to ensure high-concentration formulations
            won&apos;t aggregate or form hairpin structures before binding the toxic repeat targets.
          </p>

          <div className="rounded-xl border border-slate-200 p-4 space-y-3">
            <div className="flex items-center justify-between">
              <span className="text-[12px] font-medium text-slate-700">Hairpin Structure Prediction</span>
              <span className={`rounded-full px-2 py-0.5 text-[10px] font-semibold ${
                false
                  ? "bg-emerald-50 text-emerald-600 border border-emerald-200"
                  : false
                    ? "bg-amber-50 text-amber-600 border border-amber-200"
                    : "bg-red-50 text-red-600 border border-red-200"
              }`}>
                {"not computed"} Risk
              </span>
            </div>

            <div className="rounded-lg bg-slate-900 p-3 font-mono text-[11px] text-emerald-400 leading-relaxed">
              <span className="text-slate-500">$ RNAfold |</span> {candidate.sequence}
              <br />
              <span className="text-slate-400">{candidate.sequence}</span>
              <br />
              <span className="text-amber-400">(((...)))</span>
              <span className="text-slate-500">  ({candidate.realMetrics.selfStructureMfe.toFixed(2)} kcal/mol)</span>
            </div>

            <div className="flex items-center justify-between">
              <span className="text-[12px] font-medium text-slate-700">Self-Dimer MFE Threshold</span>
              <span className={`text-[12px] font-semibold ${candidate.realMetrics.selfStructureMfe > -3 ? "text-emerald-600" : "text-amber-600"}`}>
                {candidate.realMetrics.selfStructureMfe.toFixed(1)} kcal/mol {candidate.realMetrics.selfStructureMfe > -3 ? "✓ (> -3 threshold)" : "⚠ (at risk)"}
              </span>
            </div>

            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={secondaryChecked}
                onChange={(e) => setSecondaryChecked(e.target.checked)}
                className="h-4 w-4 rounded border-slate-300 text-brand focus:ring-brand"
              />
              <span className="text-[12px] text-slate-600">
                Confirm secondary structure passes aggregation filter for target concentration
              </span>
            </label>
          </div>

          <div className="flex justify-between">
            <button
              onClick={() => setStep(2)}
              className="flex items-center gap-1.5 rounded-lg border border-slate-200 px-4 py-2 text-[12px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
            >
              <ArrowLeft className="h-3.5 w-3.5" /> Back
            </button>
            <button
              onClick={handleComplete}
              disabled={!secondaryChecked}
              className="flex items-center gap-1.5 rounded-lg bg-emerald-600 px-5 py-2 text-[12px] font-medium text-white shadow-sm hover:bg-emerald-700 transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
            >
              <CheckCircle className="h-3.5 w-3.5" /> Finalize Design
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
