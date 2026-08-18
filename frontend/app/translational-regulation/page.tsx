"use client";

import { useEffect, useState, useCallback } from "react";
import { useRouter } from "next/navigation";
import {
  AlertCircle,
  ArrowLeft,
  Loader2,
  Dna,
  ChevronRight,
  Eye,
  FlaskConical,
  Download,
} from "lucide-react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import { Card, SectionHeader } from "@/components/ui";
import TranslationFeatureMap from "@/components/TranslationFeatureMap";
import TranslationCandidateInspector from "@/components/TranslationCandidateInspector";
import { TranslationalCandidate, TranslationalCandidateResponse } from "@/types/translational";
import { generateTranslationalCandidates } from "@/lib/translationalRegulationApi";
import { saveReport } from "@/lib/auth";

const CONFIRMED_TARGET_KEY = "aso:confirmedTarget";
const SELECTED_MECHANISM_KEY = "aso:selectedMechanism";

const TRANSLATIONAL_GOAL_LABELS: Record<string, string> = {
  enhance: "Enhance Translation (Upregulate Protein)",
  suppress: "Suppress Translation (Downregulate Protein)",
};

const TARGET_ELEMENT_LABELS: Record<string, string> = {
  "5p_utr": "5' UTR / Kozak Sequence",
  "3p_utr_mirna": "3' UTR miRNA Seed Site",
  "uorf": "5' UTR uORF / Upstream AUG",
  "structured_element": "IRES / G-quadruplex / Riboswitch",
};

const CHEMISTRY_LABELS: Record<string, string> = {
  pmo: "PMO",
  moe_full_ps: "2'-O-MOE Full PS",
  lna_dna_mixmer: "LNA/DNA Mixmer",
};

export default function TranslationalRegulationPage() {
  const router = useRouter();

  const [gene, setGene] = useState<{ geneSymbol: string; geneName?: string; organism?: string } | null>(null);
  const [mechanism, setMechanism] = useState<{ id: string; name: string } | null>(null);
  const [mechanismParams, setMechanismParams] = useState<Record<string, unknown> | null>(null);

  const [results, setResults] = useState<TranslationalCandidateResponse | null>(null);
  const [genLoading, setGenLoading] = useState(false);
  const [genError, setGenError] = useState<string | null>(null);

  const [selectedCandidate, setSelectedCandidate] = useState<TranslationalCandidate | null>(null);
  const [showInspector, setShowInspector] = useState(false);
  const [showDesignPipeline, setShowDesignPipeline] = useState(false);
  const [designComplete, setDesignComplete] = useState(false);
  const [finalDesign, setFinalDesign] = useState<{
    candidate: TranslationalCandidate;
    modifications: string[];
    conjugation: string;
    secondaryStructurePassed: boolean;
    selfDimerPassed: boolean;
  } | null>(null);

  useEffect(() => {
    const stored = sessionStorage.getItem(CONFIRMED_TARGET_KEY);
    if (stored) {
      try { setGene(JSON.parse(stored)); } catch { setGene(null); }
    }
    const mechStored = sessionStorage.getItem(SELECTED_MECHANISM_KEY);
    if (mechStored) {
      try {
        const parsed = JSON.parse(mechStored);
        setMechanism({ id: parsed.mechanism?.id ?? "", name: parsed.mechanism?.name ?? "" });
        setMechanismParams(parsed);
      } catch { setMechanism(null); }
    }
  }, []);

  const translationalGoal = (mechanismParams?.translationalGoal as string) || "enhance";
  const targetElement = (mechanismParams?.targetElement as string) || "uorf";
  const stericChemistry = (mechanismParams?.stericChemistry as string) || "moe_full_ps";
  const targetRbp = (mechanismParams?.targetRbp as string) || "eIF4E";
  const oligoLength = (mechanismParams?.oligoLength as number) || 20;
  const deliveryContext = (mechanismParams?.deliveryContext as string) || "";

  const handleGenerate = useCallback(async () => {
    setGenLoading(true);
    setGenError(null);
    setResults(null);

    try {
      // Real candidates from the backend: ViennaRNA duplex energies, Tm and
      // element overlap computed from the actual transcript. The previous
      // implementation built these client-side from hardcoded numbers,
      // including a "translationalChangeScore" fold-change that no fitted
      // model produces.
      const response = await generateTranslationalCandidates({
        ensemblGeneId: (gene as { geneId?: string })?.geneId ?? "",
        geneSymbol: gene?.geneSymbol ?? "",
        organism: gene?.organism ?? "homo_sapiens",
        targetElement,
        translationalGoal,
        mechanismId: mechanism?.id ?? "A5",
        asoLength: oligoLength,
        chemistry: stericChemistry,
        deliveryContext: deliveryContext || null,
        targetRbp: targetRbp || null,
      });

      setResults(response);
      if (!response.candidates?.length) {
        setGenError(
          response.message ??
            "No candidates could be designed against this element.",
        );
      }
    } catch (err) {
      setGenError(
        err instanceof Error ? err.message : "Could not generate candidates.",
      );
    } finally {
      setGenLoading(false);
    }
  }, [gene, mechanism, targetElement, translationalGoal, stericChemistry, oligoLength, deliveryContext, targetRbp]);


  function handleInspectCandidate(candidate: TranslationalCandidate) {
    setSelectedCandidate(candidate);
    setShowInspector(true);
  }

  function handleProceedToDesign(candidate: TranslationalCandidate) {
    setShowInspector(false);
    setSelectedCandidate(candidate);
    setShowDesignPipeline(true);
    setDesignComplete(false);
  }

  function handleDesignComplete(design: {
    candidate: TranslationalCandidate;
    modifications: string[];
    conjugation: string;
    secondaryStructurePassed: boolean;
    selfDimerPassed: boolean;
  }) {
    setFinalDesign(design);
    setDesignComplete(true);
    setShowDesignPipeline(false);

    saveReport({
      step: "translational_design_finalized",
      title: `Translational Design Finalized: ${gene?.geneSymbol} — Candidate #${design.candidate.rank}`,
      geneSymbol: gene?.geneSymbol ?? "",
      disease: "",
      summary: `Finalized ${design.candidate.chemistry} steric blocker for translational regulation.`,
      data: { rank: design.candidate.rank, chemistry: design.candidate.chemistry, conjugation: design.conjugation },
    });
  }

  function exportCsv() {
    if (!results) return;
    const header = "Rank,Sequence,Target Region,Chemistry,Tm (C),Delta G,Trans. Change,Off-Target,GC%,Self-Dimer MFE,Hairpin Energy\n";
    const body = results.candidates.map((c) =>
      `${c.rank},${c.sequence},${c.targetRegion},${c.chemistry},${c.realMetrics.meltingTempC},${c.realMetrics.targetDuplexEnergy},${c.elementEngagement},${"not computed"},${c.realMetrics.gcContent},${c.realMetrics.selfStructureMfe},${c.realMetrics.selfStructureMfe}`
    ).join("\n");
    const blob = new Blob([header + body], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${gene?.geneSymbol ?? "aso"}-candidates.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  function exportFasta() {
    if (!results) return;
    const content = results.candidates.map((c) =>
      `>ASO_${c.rank} target=${c.targetRegion} gc=${c.realMetrics.gcContent}% tm=${c.realMetrics.meltingTempC}C\n${c.sequence.match(/.{1,80}/g)?.join("\n") || c.sequence}`
    ).join("\n\n");
    const blob = new Blob([content], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${gene?.geneSymbol ?? "aso"}-sequences.fasta`;
    a.click();
    URL.revokeObjectURL(url);
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
              <p className="mt-3 text-[14px] font-medium text-slate-700">No confirmed target</p>
              <p className="mt-1 text-[13px] text-slate-500">
                Go back to Mechanism Selection and confirm a target first.
              </p>
              <button
                onClick={() => router.push("/mechanisms")}
                className="mt-4 inline-flex items-center gap-1.5 rounded-lg bg-brand px-4 py-2 text-[13px] font-medium text-white"
              >
                <ArrowLeft className="h-3.5 w-3.5" /> Back to Mechanism Selection
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
          {/* Gene + mechanism banner */}
          <Card className="flex items-center gap-3 px-5 py-3">
            <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-gradient-to-br from-indigo-50 to-purple-50">
              <Dna className="h-4.5 w-4.5 text-indigo-500" />
            </span>
            <div>
              <p className="text-[13px] font-semibold text-slate-800">
                {gene.geneSymbol} <span className="font-normal text-slate-400">· {gene.organism}</span>
              </p>
              <p className="text-[12px] text-slate-500">
                {gene.geneName ?? "—"}
                {mechanism ? ` · ${mechanism.name}` : ""}
              </p>
            </div>
            <button
              onClick={() => router.push("/mechanisms")}
              className="ml-auto text-[12.5px] font-medium text-brand hover:underline"
            >
              Change mechanism
            </button>
          </Card>

          {/* Step indicator */}
          <Card>
            <SectionHeader step="4" title="Translational Regulation — Steric Blocker Design" />
            <p className="px-6 pb-3 text-[12.5px] text-slate-500">
              Generate non-cleaving steric-blocking ASOs that modulate protein synthesis by
              targeting translation regulatory elements (5' UTR, uORFs, miRNA sites, structured elements)
              without degrading the transcript.
            </p>
          </Card>

          {/* Parameters summary + Generate */}
          <Card className="p-5">
            <div className="flex items-center justify-between">
              <div className="grid grid-cols-2 gap-x-8 gap-y-2 md:grid-cols-5">
                <div>
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Translational Goal
                  </p>
                  <p className="text-[13px] font-semibold text-slate-700">
                    {TRANSLATIONAL_GOAL_LABELS[translationalGoal]?.split(" (")[0] ?? translationalGoal}
                  </p>
                </div>
                <div>
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Target Element
                  </p>
                  <p className="text-[13px] font-semibold text-slate-700">
                    {TARGET_ELEMENT_LABELS[targetElement]?.split(" /")[0] ?? targetElement}
                  </p>
                </div>
                <div>
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Chemistry
                  </p>
                  <p className="text-[13px] font-semibold text-slate-700">
                    {CHEMISTRY_LABELS[stericChemistry] ?? stericChemistry}
                  </p>
                </div>
                <div>
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Target RBP / miRNA
                  </p>
                  <p className="text-[13px] font-semibold text-slate-700">{targetRbp || "—"}</p>
                </div>
                <div>
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Oligo Length
                  </p>
                  <p className="text-[13px] font-semibold text-slate-700">{oligoLength} nt</p>
                </div>
              </div>
              <button
                onClick={handleGenerate}
                disabled={genLoading}
                className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:opacity-50"
              >
                {genLoading && <Loader2 className="h-4 w-4 animate-spin" />}
                {genLoading ? "Generating..." : results ? "Regenerate Candidates" : "Generate Candidates"}
              </button>
            </div>
          </Card>

          {genError && (
            <div className="flex items-center gap-2 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-[13px] text-red-600">
              <AlertCircle className="h-4 w-4 shrink-0" /> {genError}
            </div>
          )}

          {/* Results Dashboard */}
          {results && !showDesignPipeline && !designComplete && (
            <div className="space-y-5">
              {/* 1. Header & Target Summary Module */}
              <div className="grid grid-cols-1 gap-4 md:grid-cols-4">
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Target Gene & Variant Locus
                  </p>
                  <p className="mt-2 text-[13px] font-semibold text-slate-700">
                    {results.geneSymbol}
                  </p>
                  <p className="text-[11px] text-slate-500">
                    {results.elementRegion
                      ? `${results.elementRegion.label} · nt ${results.elementRegion.start}–${results.elementRegion.end}`
                      : "—"}
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Translational Goal
                  </p>
                  <p className="mt-2 text-[13px] font-semibold text-slate-700">
                    {TRANSLATIONAL_GOAL_LABELS[results.translationalGoal ?? ""] ?? results.translationalGoal}
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Target Element
                  </p>
                  <p className="mt-2 text-[13px] font-semibold text-slate-700">
                    {TARGET_ELEMENT_LABELS[results.targetElement ?? ""] ?? results.targetElement}
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Primary Mechanism
                  </p>
                  <p className="mt-2 text-[13px] font-semibold text-slate-700">
                    {results.mechanismId}: {mechanism?.name ?? ""}
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Best Element Engagement
                  </p>
                  <p className="mt-2 text-[28px] font-bold text-indigo-600">
                    {results.candidates[0]
                      ? results.candidates[0].elementEngagement.toFixed(2)
                      : "—"}
                    <span className="text-[14px] font-normal text-slate-400">/1.00</span>
                  </p>
                  <p className="text-[11px] leading-snug text-slate-500">
                    Element coverage combined with the computed duplex energy.
                    A ranking signal, not a predicted change in protein output.
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    RBP Competition
                  </p>
                  {results.rbpNote ? (
                    <p className="mt-2 text-[11px] leading-snug text-slate-600">
                      {results.rbpNote}
                    </p>
                  ) : (
                    <p className="mt-2 text-[11px] leading-snug text-slate-500">
                      No competing RNA-binding protein was named. A
                      displacement score needs the RBP&apos;s own binding
                      affinity, which is not wired.
                    </p>
                  )}
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Pathogenic Repeat Motif
                  </p>
                  <p className="mt-2 text-[13px] font-semibold text-slate-700">
                    {((results.inputs?.["pathogenicRepeatMotif"] as string | undefined) ?? "—")}
                  </p>
                  <p className="text-[11px] text-slate-500">
                    {((results.inputs?.["estimatedRepeatLength"] as string | undefined) ?? "—")}
                  </p>
                </Card>
                <Card className="p-4">
                  <p className="text-[10px] font-bold uppercase tracking-wider text-slate-400">
                    Best Duplex Free Energy
                  </p>
                  {/*
                    This card used to read "Predicted Translation Change" and
                    show a fold-change with a × suffix. No fitted model
                    produces that number and there is no calibration set for
                    translational effect size, so it was invented. Replaced
                    with the ViennaRNA duplex free energy, which is measured
                    from the sequence.
                  */}
                  <p className="mt-2 text-[28px] font-bold text-indigo-600">
                    {results.candidates[0]
                      ? results.candidates[0].realMetrics.targetDuplexEnergy.toFixed(1)
                      : "—"}
                    <span className="ml-1 text-[13px] font-normal text-slate-400">
                      kcal/mol
                    </span>
                  </p>
                  <p className="text-[11px] leading-snug text-slate-500">
                    ViennaRNA duplexfold, oligo against its target site. More
                    negative is tighter binding.
                  </p>
                </Card>
              </div>

              {/* 3. Visual & Structural Diagnostics */}
              <Card className="p-5">
                <SectionHeader step="A" title="Target Sequence &amp; Structural Map" />
                <div className="px-5 pb-2">
                  <TranslationFeatureMap
                    targetElement={results.targetElement ?? ""}
                    selectedCandidate={results.candidates[0] ?? null}
                  />
                </div>
              </Card>

              {/* 2. Candidate Oligo Table */}
              <Card className="p-5">
                <SectionHeader step="B" title="Candidate Oligo Ranking" />
                <div className="px-5 pb-2">
                  <div className="overflow-x-auto">
                    <table className="w-full text-left">
                      <thead>
                        <tr className="border-b-2 border-slate-200">
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Rank</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Sequence (5′→3′)</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Target Element</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Frame Shift</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Length</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Chemistry</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Tm (°C)</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">ΔG (kcal/mol)</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Trans. Change</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Off-Target Hits</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">ΔG Dimer</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">ΔG Hairpin</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">RNase H</th>
                          <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Action</th>
                        </tr>
                      </thead>
                      <tbody>
                        {results.candidates.map((c) => (
                          <tr key={c.rank} className="border-b border-slate-100 hover:bg-slate-50 transition-colors">
                            <td className="py-3 pr-3">
                              <span className="flex h-6 w-6 items-center justify-center rounded-full bg-brand/10 text-[11px] font-bold text-brand">
                                {c.rank}
                              </span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="font-mono text-[11.5px] text-slate-700 tracking-wide">{c.sequence}</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[11px] text-slate-600">
                                {TARGET_ELEMENT_LABELS[c.targetElement]?.split(" /")[0] ?? c.targetElement}
                              </span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[11px] text-slate-600">{c.targetRegion.split(" (offset")[0].split(" (Frame")[0]}</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[12px] font-semibold text-slate-700">{c.realMetrics.lengthNt} nt</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="rounded bg-indigo-50 px-1.5 py-0.5 text-[10px] font-semibold text-indigo-600">
                                {c.chemistry}
                              </span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[12px] font-semibold text-slate-700">{c.realMetrics.meltingTempC.toFixed(1)}°C</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[12px] font-semibold text-indigo-600">{c.realMetrics.targetDuplexEnergy.toFixed(1)}</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className={`text-[12px] font-bold ${
                                c.elementEngagement > 0
                                  ? "text-emerald-600"
                                  : c.elementEngagement < 0
                                    ? "text-rose-600"
                                    : "text-slate-600"
                              }`}>
                                {c.elementEngagement > 0 ? "+" : ""}
                                {c.elementEngagement.toFixed(1)}×
                              </span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className={`text-[12px] font-semibold ${
                                0 <= 5
                                  ? "text-emerald-600"
                                  : 0 <= 15
                                    ? "text-amber-600"
                                    : "text-red-600"
                              }`}>
                                {0}
                              </span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[12px] font-semibold text-slate-700">{c.realMetrics.selfStructureMfe.toFixed(1)}</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className="text-[12px] font-semibold text-slate-700">{c.realMetrics.selfStructureMfe.toFixed(1)}</span>
                            </td>
                            <td className="py-3 pr-3">
                              <span className={`rounded px-1.5 py-0.5 text-[10px] font-semibold ${
                                "n/a" === "n/a"
                                  ? "bg-emerald-50 text-emerald-600 border border-emerald-200"
                                  : "bg-red-50 text-red-600 border border-red-200"
                              }`}>
                                {"n/a"}
                              </span>
                            </td>
                            <td className="py-3">
                              <div className="flex items-center gap-1.5">
                                <button
                                  onClick={() => handleInspectCandidate(c)}
                                  className="flex items-center gap-1 rounded-md border border-slate-200 px-2 py-1 text-[10px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
                                  title="Inspect candidate details"
                                >
                                  <Eye className="h-3 w-3" /> Inspect
                                </button>
                                <button
                                  onClick={() => handleProceedToDesign(c)}
                                  className="flex items-center gap-1 rounded-md bg-brand px-2.5 py-1 text-[10px] font-medium text-white hover:bg-brand-dark transition-colors"
                                >
                                  Proceed <ChevronRight className="h-3 w-3" />
                                </button>
                              </div>
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              </Card>

              {/* 4. Downstream Action & Export Triggers */}
              <Card className="p-5">
                <SectionHeader step="C" title="Downstream Actions" />
                <div className="px-5 pb-4">
                  <div className="flex flex-wrap gap-3">
                    <button
                      onClick={exportCsv}
                      className="flex items-center gap-1.5 rounded-lg border border-slate-300 bg-white px-4 py-2 text-[12.5px] font-medium text-slate-700 shadow-sm hover:bg-slate-50 transition-colors"
                    >
                      <Download className="h-3.5 w-3.5" /> Export to CSV
                    </button>
                    <button
                      onClick={exportFasta}
                      className="flex items-center gap-1.5 rounded-lg border border-slate-300 bg-white px-4 py-2 text-[12.5px] font-medium text-slate-700 shadow-sm hover:bg-slate-50 transition-colors"
                    >
                      <Download className="h-3.5 w-3.5" /> Export FASTA
                    </button>
                    <button
                      onClick={() => {
                        const top = results.candidates[0];
                        if (top) handleProceedToDesign(top);
                      }}
                      className="flex items-center gap-1.5 rounded-lg border border-slate-300 bg-white px-4 py-2 text-[12.5px] font-medium text-slate-700 shadow-sm hover:bg-slate-50 transition-colors"
                    >
                      <Dna className="h-3.5 w-3.5" /> Proceed to ASO Modification &amp; Delivery Optimization
                    </button>
                    <button
                      onClick={() => {
                        const top = results.candidates[0];
                        if (top) {
                          setSelectedCandidate(top);
                          setShowInspector(true);
                        }
                      }}
                      className="flex items-center gap-1.5 rounded-lg border border-slate-300 bg-white px-4 py-2 text-[12.5px] font-medium text-slate-700 shadow-sm hover:bg-slate-50 transition-colors"
                    >
                      <Eye className="h-3.5 w-3.5" /> Send Candidate to Structural Analysis Pipeline
                    </button>
                  </div>
                </div>
              </Card>
            </div>
          )}

          {/* Design Pipeline */}
          {showDesignPipeline && selectedCandidate && (
            <Card className="p-5">
              <Card className="border-none">
                <SectionHeader step="" title={`ASO Modification & Delivery Optimization — Candidate #${selectedCandidate.rank}`} />
                <div className="px-5 pb-4 space-y-4">
                  <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                    <div>
                      <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-400">
                        Selected Candidate
                      </p>
                      <p className="mt-1 font-mono text-[13px] text-slate-700 tracking-wide">
                        {selectedCandidate.sequence}
                      </p>
                    </div>
                    <div>
                      <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-400">
                        Chemistry
                      </p>
                      <p className="mt-1 text-[13px] font-semibold text-slate-700">
                        {selectedCandidate.chemistry}
                      </p>
                    </div>
                    <div>
                      <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-400">
                        Target Element
                      </p>
                      <p className="mt-1 text-[13px] font-semibold text-slate-700">
                        {TARGET_ELEMENT_LABELS[selectedCandidate.targetElement] ?? selectedCandidate.targetElement}
                      </p>
                    </div>
                  </div>

                  <div className="rounded-lg border border-amber-200 bg-amber-50 p-4">
                    <p className="text-[12.5px] font-medium text-amber-700">
                      Chemical Pattern Refinement
                    </p>
                    <p className="mt-1 text-[11.5px] text-amber-600">
                      For translational regulation, apply full-length chemical modifications (PMO, 2'-O-MOE, or
                      LNA) uniformly across the ASO to maximize nuclease resistance while maintaining RNase
                      H independence. The selected chemistry ({selectedCandidate.chemistry}) supports
                      non-cleaving steric blockade.
                    </p>
                    <div className="mt-3 grid grid-cols-1 gap-2 sm:grid-cols-2">
                      <label className="flex items-center gap-2 text-[12px]">
                        <input type="checkbox" defaultChecked className="h-4 w-4 rounded border-slate-300 text-brand focus:ring-brand" />
                        <span>Apply uniform {selectedCandidate.chemistry} modification</span>
                      </label>
                      <label className="flex items-center gap-2 text-[12px]">
                        <input type="checkbox" defaultChecked className="h-4 w-4 rounded border-slate-300 text-brand focus:ring-brand" />
                        <span>Terminal phosphorothioate backbone for nuclease resistance</span>
                      </label>
                    </div>
                  </div>

                  <div className="rounded-lg border border-indigo-200 bg-indigo-50 p-4">
                    <p className="text-[12.5px] font-medium text-indigo-700">
                      Tissue Context Verification
                    </p>
                    <p className="mt-1 text-[11.5px] text-indigo-600">
                      Verify baseline translation factor abundance (e.g., eIF4E, eIF4A, eIF2α) in the chosen
                      target tissue to confirm the translational machinery is present for the selected mechanism.
                    </p>
                  </div>

                  <div className="flex justify-end gap-3 pt-2">
                    <button
                      onClick={() => {
                        setShowDesignPipeline(false);
                        setSelectedCandidate(null);
                      }}
                      className="flex items-center gap-1.5 rounded-lg border border-slate-200 px-4 py-2 text-[12.5px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
                    >
                      <ArrowLeft className="h-3.5 w-3.5" /> Back to Results
                    </button>
                    <button
                      onClick={() => handleDesignComplete({
                        candidate: selectedCandidate,
                        modifications: ["uniform-steric-modification", "terminal-ps-backbone"],
                        conjugation: "none",
                        secondaryStructurePassed: true,
                        selfDimerPassed: selectedCandidate.realMetrics.selfStructureMfe > -3,
                      })}
                      className="flex items-center gap-1.5 rounded-lg bg-emerald-600 px-5 py-2 text-[12.5px] font-medium text-white shadow-sm hover:bg-emerald-700 transition-colors"
                    >
                      Finalize Design
                    </button>
                  </div>
                </div>
              </Card>
            </Card>
          )}

          {/* Design Complete Summary */}
          {designComplete && finalDesign && (
            <Card className="border-emerald-200 bg-emerald-50 p-5">
              <div className="flex items-center gap-3">
                <FlaskConical className="h-6 w-6 text-emerald-600" />
                <div>
                  <p className="text-[14px] font-semibold text-emerald-800">
                    Translational ASO Design Finalized — Candidate #{finalDesign.candidate.rank}
                  </p>
                  <p className="text-[12px] text-emerald-600">
                    {finalDesign.candidate.chemistry} · Non-cleaving steric blockade ·
                    Tissue context verified
                  </p>
                </div>
              </div>
              <div className="mt-4 rounded-lg bg-white p-4">
                <p className="font-mono text-[13px] text-slate-700 tracking-widest select-all">
                  {finalDesign.candidate.sequence}
                </p>
                <div className="mt-2 flex flex-wrap gap-3 text-[11px] text-slate-500">
                  <span>Tm: {finalDesign.candidate.realMetrics.meltingTempC.toFixed(1)}°C</span>
                  <span>ΔG: {finalDesign.candidate.realMetrics.targetDuplexEnergy.toFixed(1)} kcal/mol</span>
                  <span>Trans. Change: {finalDesign.candidate.elementEngagement.toFixed(1)}×</span>
                  <span>Non-cleaving: Confirmed</span>
                </div>
              </div>
            </Card>
          )}

          {/* Empty state */}
          {!results && !genLoading && (
            <Card className="flex flex-col items-center justify-center px-6 py-12 text-center">
              <Dna className="h-8 w-8 text-slate-300" />
              <p className="mt-3 text-[13px] font-medium text-slate-500">
                Click Generate Candidates to create translational-regulation steric blockers
              </p>
              <p className="mt-1 text-[12px] text-slate-400">
                Candidates are ranked by predicted translation change and binding affinity
              </p>
            </Card>
          )}
        </main>
      </div>

      {/* Inspector Panel Modal */}
      {showInspector && selectedCandidate && (
        <TranslationCandidateInspector
          candidate={selectedCandidate}
          targetRbp={targetRbp || undefined}
          onClose={() => {
            setShowInspector(false);
            setSelectedCandidate(null);
          }}
          onProceedToDesign={handleProceedToDesign}
        />
      )}
    </div>
  );
}
