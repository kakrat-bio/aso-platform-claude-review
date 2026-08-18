"use client";

import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import {
  AlertCircle,
  ArrowLeft,
  Loader2,
  Dna,
  Beaker,
  ChevronDown,
  ChevronUp,
  Copy,
  Check,
  Download,
  FlaskConical,
  FileText,
  BarChart3,
} from "lucide-react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import { Card, SectionHeader, FieldLabel } from "@/components/ui";
import { GeneTargetObject } from "@/types/gene";
import {
  MechanismOptions,
  TherapeuticGoalId,
  THERAPEUTIC_GOALS,
  RnaEngineeringCandidate,
} from "@/types/mechanism";
import {
  fetchMechanismOptions,
  rankRnaEngineeringMechanisms,
} from "@/lib/mechanismApi";
import { saveReport } from "@/lib/auth";

const CONFIRMED_TARGET_KEY = "aso:confirmedTarget";
const SELECTED_MECHANISM_KEY = "aso:selectedMechanism";
const SELECTED_GOAL_KEY = "aso:therapeuticGoal";

export default function RnaEngineeringPage() {
  const router = useRouter();
  const [gene, setGene] = useState<GeneTargetObject | null>(null);
  const [options, setOptions] = useState<MechanismOptions | null>(null);
  const [optionsError, setOptionsError] = useState<string | null>(null);
  const [selectedGoal, setSelectedGoal] = useState<TherapeuticGoalId | null>(null);

  const [structuralClass, setStructuralClass] = useState("");
  const [targetType, setTargetType] = useState("");
  const [scaffold, setScaffold] = useState("");
  const [chemStabilization, setChemStabilization] = useState("");
  const [kdGoal, setKdGoal] = useState("");
  const [deliveryContext, setDeliveryContext] = useState("");

  const [ranking, setRanking] = useState<any>(null);
  const [candidates, setCandidates] = useState<RnaEngineeringCandidate[]>([]);
  const [selectedCandidate, setSelectedCandidate] = useState<RnaEngineeringCandidate | null>(null);
  const [protocolNotice, setProtocolNotice] = useState<"selex" | "synthesis" | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [copiedId, setCopiedId] = useState<string | null>(null);

  useEffect(() => {
    const stored = sessionStorage.getItem(CONFIRMED_TARGET_KEY);
    if (stored) {
      try { setGene(JSON.parse(stored)); } catch { setGene(null); }
    }
    const savedGoal = sessionStorage.getItem(SELECTED_GOAL_KEY) as TherapeuticGoalId | null;
    if (savedGoal && THERAPEUTIC_GOALS.some((g) => g.id === savedGoal)) {
      setSelectedGoal(savedGoal);
    }
    fetchMechanismOptions()
      .then(setOptions)
      .catch((e) => setOptionsError(e instanceof Error ? e.message : "Failed to load options."));
  }, []);

  function handleSelectGoal(goalId: TherapeuticGoalId) {
    setSelectedGoal(goalId);
    sessionStorage.setItem(SELECTED_GOAL_KEY, goalId);
    setRanking(null);
    setCandidates([]);
    setSelectedCandidate(null);
    setStructuralClass("");
    setTargetType("");
    setScaffold("");
    setChemStabilization("");
    setKdGoal("");
    setDeliveryContext("");
  }

  function clearResults() {
    setRanking(null);
    setCandidates([]);
    setSelectedCandidate(null);
  }

  async function handleGenerate() {
    if (!gene || !selectedGoal) return;
    if (!structuralClass || !targetType || !scaffold || !chemStabilization || !kdGoal) {
      setError("Please fill in all required design parameters.");
      return;
    }
    setLoading(true);
    setError(null);
    setSelectedCandidate(null);
    try {
      const res = await rankRnaEngineeringMechanisms({
        geneSymbol: gene.geneSymbol,
        structuralClass,
        targetType,
        scaffold,
        chemStabilization,
        kdGoal,
        deliveryContext: deliveryContext || undefined,
      });
      setRanking(res);
      setCandidates(res.candidates || []);
      saveReport({
        step: "rna_engineering",
        title: `RNA Engineering: ${gene.geneSymbol} (${selectedGoal})`,
        geneSymbol: gene.geneSymbol,
        disease: gene.disease || "",
        summary: `Generated ${(res.candidates || []).length} RNA engineering candidates for ${gene.geneSymbol}. Top: ${(res.candidates || [])[0]?.constructId || "N/A"}.`,
        data: { goal: selectedGoal, structuralClass, targetType, scaffold, kdGoal, topCandidates: (res.candidates || []).slice(0, 3) },
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Generation failed.");
    } finally {
      setLoading(false);
    }
  }

  function isRankDisabled(): boolean {
    if (!gene || !selectedGoal || loading) return true;
    if (selectedGoal !== "TG09") return true;
    return !structuralClass || !targetType || !scaffold || !chemStabilization || !kdGoal;
  }

  function copySequence(seq: string, id: string) {
    navigator.clipboard.writeText(seq).then(() => {
      setCopiedId(id);
      setTimeout(() => setCopiedId(null), 1500);
    });
  }

  function exportFasta(candidate: RnaEngineeringCandidate) {
    const header = `>${candidate.constructId} | ${candidate.structuralMotif} | ${candidate.length}nt | Tm=${candidate.tm}°C | ΔG=${candidate.deltaGFolding} kcal/mol`;
    const wrapped = candidate.sequence.match(/.{1,60}/g)?.join("\n") ?? candidate.sequence;
    const blob = new Blob([`${header}\n${wrapped}\n`], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${candidate.constructId}.fasta`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  function exportDotBracket(candidate: RnaEngineeringCandidate) {
    const content = `>${candidate.constructId}\n${candidate.sequence}\n${candidate.dotBracket}\n`;
    const blob = new Blob([content], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${candidate.constructId}.dbn`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  if (!gene) {
    return (
      <div className="flex min-h-screen bg-[#F5F6FA]">
        <Sidebar />
        <div className="flex min-h-screen flex-1 flex-col">
          <Topbar />
          <main className="flex flex-1 items-center justify-center px-6">
            <Card className="max-w-md p-8 text-center">
              <AlertCircle className="mx-auto h-8 w-8 text-slate-300" />
              <p className="mt-3 text-[14px] font-medium text-slate-700">No confirmed target found</p>
              <p className="mt-1 text-[13px] text-slate-500">
                Go back to Basic Information, load a gene, and hit Confirm & Proceed first.
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
    <div className="flex min-h-screen bg-[#F5F6FA]">
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
              onClick={() => router.push("/mechanisms")}
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
                  className={`rounded-lg border p-4 text-left transition-colors ${
                    selectedGoal === goal.id
                      ? "border-brand bg-brand/5 ring-1 ring-brand"
                      : "border-slate-200 bg-white hover:border-slate-300 hover:bg-slate-50"
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

          {/* Step 3: TG09 RNA Engineering Input Form */}
          {selectedGoal === "TG09" && (
            <Card>
              <SectionHeader step="3" title="RNA Engineering — Structural / Functional Design" />
              {optionsError && (
                <p className="px-6 pb-2 text-[12.5px] text-red-600">{optionsError}</p>
              )}
              <div className="grid grid-cols-1 gap-4 px-6 pb-4 md:grid-cols-2 lg:grid-cols-3">
                <div>
                  <FieldLabel hint="What class of structured RNA molecule do you want to design?">
                    Structural Class <span className="text-red-500">*</span>
                  </FieldLabel>
                  <select
                    value={structuralClass}
                    onChange={(e) => { setStructuralClass(e.target.value); clearResults(); }}
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
                    value={targetType}
                    onChange={(e) => { setTargetType(e.target.value); clearResults(); }}
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
                    value={scaffold}
                    onChange={(e) => { setScaffold(e.target.value); clearResults(); }}
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
                    value={chemStabilization}
                    onChange={(e) => { setChemStabilization(e.target.value); clearResults(); }}
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
                    value={kdGoal}
                    onChange={(e) => { setKdGoal(e.target.value); clearResults(); }}
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
                    onChange={(e) => { setDeliveryContext(e.target.value); clearResults(); }}
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  >
                    <option value="">Not specified</option>
                    {options?.deliveryContexts.map((o) => (
                      <option key={o.id} value={o.id}>{o.label}</option>
                    ))}
                  </select>
                </div>
              </div>

              <div className="flex justify-end px-6 pb-5">
                <button
                  onClick={handleGenerate}
                  disabled={isRankDisabled()}
                  className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                >
                  {loading && <Loader2 className="h-4 w-4 animate-spin" />}
                  {loading ? "Generating Candidates..." : "Generate Candidates"}
                </button>
              </div>
            </Card>
          )}

          {selectedGoal && selectedGoal !== "TG09" && (
            <Card>
              <div className="px-6 pb-4">
                <div className="rounded-lg border border-dashed border-slate-300 bg-slate-50 p-6 text-center">
                  <p className="text-[13px] font-medium text-slate-600">
                    Mechanism selection for {THERAPEUTIC_GOALS.find((g) => g.id === selectedGoal)?.name || selectedGoal} is available on the Mechanisms page.
                  </p>
                  <p className="mt-1 text-[12px] text-slate-400">
                    Navigate to the Mechanisms page to configure and rank mechanisms for this therapeutic goal.
                  </p>
                </div>
              </div>
            </Card>
          )}

          {error && (
            <div className="flex items-center gap-2 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-[13px] text-red-600">
              <AlertCircle className="h-4 w-4 shrink-0" />
              {error}
            </div>
          )}

          {/* Results Dashboard */}
          {ranking && (
            <>
              {/* A. Header & Structural Summary */}
              <Card>
                <SectionHeader step="4" title="Header & Structural Summary" />
                <div className="px-6 pb-5 grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-4">
                  <div>
                    <p className="text-[11px] text-slate-400">Target Partner & Type</p>
                    <p className="text-[12.5px] font-semibold text-slate-800">{gene.geneSymbol}</p>
                    <p className="text-[11px] text-slate-500">{targetType ? options?.rnaEngineering.targetTypes.find((t) => t.id === targetType)?.label : "—"}</p>
                  </div>
                  <div>
                    <p className="text-[11px] text-slate-400">Selected Structural Modality</p>
                    <p className="text-[12.5px] font-semibold text-slate-800">{structuralClass ? options?.rnaEngineering.structuralClasses.find((s) => s.id === structuralClass)?.label : "—"}</p>
                    <p className="text-[11px] text-slate-500">{ranking.therapeuticGoal}</p>
                  </div>
                  <div>
                    <p className="text-[11px] text-slate-400">Predicted Binding Affinity (Kd)</p>
                    <p className="text-[12.5px] font-semibold text-slate-800">
                      {candidates.length > 0 && typeof candidates[0].kdPrediction === "number"
                        ? `${candidates[0].kdPrediction} nM`
                        : candidates.length > 0
                          ? String(candidates[0].kdPrediction)
                          : "—"}
                    </p>
                  </div>
                  <div>
                    <p className="text-[11px] text-slate-400">Serum Half-Life (t1/2)</p>
                    <p className="text-[12.5px] font-semibold text-slate-800">{candidates.length > 0 ? candidates[0].serumStability : "—"}</p>
                  </div>
                </div>
              </Card>

              {/* B. Candidate Structural Table */}
              <Card>
                <SectionHeader step="5" title="Candidate Structural Table" />
                <div className="overflow-x-auto px-6 pb-5">
                  <table className="w-full text-left text-[11.5px]">
                    <thead>
                      <tr className="border-b-2 border-slate-200">
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Rank</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Construct ID</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Structural Motif</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Length</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Tm (°C)</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">ΔG folding</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Kd / kcat</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Specificity</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Serum t1/2</th>
                        <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-center">Action</th>
                      </tr>
                    </thead>
                    <tbody>
                      {candidates.map((c) => (
                        <tr
                          key={c.constructId}
                          className={`border-b border-slate-100 last:border-0 cursor-pointer transition-colors ${
                            selectedCandidate?.constructId === c.constructId ? "bg-brand/5" : "hover:bg-slate-50"
                          }`}
                          onClick={() => setSelectedCandidate(c)}
                        >
                          <td className="py-3 pr-3">
                            <span className="inline-flex h-6 w-6 items-center justify-center rounded-full bg-brand/10 text-[11px] font-bold text-brand">
                              {c.rank}
                            </span>
                          </td>
                          <td className="py-3 pr-3">
                            <span className="text-[12px] font-semibold text-slate-800">{c.constructId}</span>
                            <span className="ml-1.5 text-[10px] text-slate-400">{c.mechanismId}</span>
                          </td>
                          <td className="py-3 pr-3 text-[11.5px] text-slate-700">{c.structuralMotif}</td>
                          <td className="py-3 pr-3 text-right text-[11.5px] text-slate-700">{c.length} nt</td>
                          <td className="py-3 pr-3 text-right text-[11.5px] text-slate-700">{c.tm.toFixed(1)}°C</td>
                          <td className="py-3 pr-3 text-right text-[11.5px] text-slate-700">{c.deltaGFolding} kcal/mol</td>
                          <td className="py-3 pr-3 text-right">
                            <span className={`text-[11.5px] font-semibold ${typeof c.kdPrediction === "number" && c.kdPrediction < 1 ? "text-emerald-600" : "text-slate-700"}`}>
                              {typeof c.kdPrediction === "number" ? `${c.kdPrediction} nM` : c.kdPrediction}
                            </span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${
                              c.targetSpecificityScore >= 90
                                ? "bg-emerald-50 text-emerald-700 border-emerald-200"
                                : c.targetSpecificityScore >= 80
                                  ? "bg-blue-50 text-blue-700 border-blue-200"
                                  : "bg-amber-50 text-amber-700 border-amber-200"
                            }`}>
                              {c.targetSpecificityScore} / 100
                            </span>
                          </td>
                          <td className="py-3 pr-3 text-right text-[11.5px] text-slate-700">{c.serumStability}</td>
                          <td className="py-3 text-center">
                            <button
                              onClick={(e) => { e.stopPropagation(); setSelectedCandidate(c); }}
                              className="rounded-md border border-slate-200 px-3 py-1.5 text-[11px] font-medium text-slate-600 hover:bg-white transition-colors"
                            >
                              [ Select ]
                            </button>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </Card>

              {/* C. Visual & Structural Diagnostics Drawer */}
              {selectedCandidate && (
                <Card>
                  <SectionHeader
                    step="6"
                    title={`Candidate Inspection: ${selectedCandidate.constructId}`}
                    right={
                      <button
                        onClick={() => setSelectedCandidate(null)}
                        className="text-[12px] text-slate-400 hover:text-slate-600"
                      >
                        Close
                      </button>
                    }
                  />
                  <div className="px-6 pb-5 space-y-5">
                    {/* 2D Secondary Structure Plot */}
                    <div>
                      <h3 className="text-[12px] font-semibold text-slate-700 mb-2">2D Secondary Structure Plot</h3>
                      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
                        <div className="font-mono text-[12px] leading-relaxed break-all text-slate-700">
                          {selectedCandidate.sequence}
                        </div>
                        <div className="font-mono text-[12px] leading-relaxed break-all text-brand mt-1">
                          {selectedCandidate.dotBracket}
                        </div>
                        <div className="mt-3 flex flex-wrap gap-3 text-[10px] text-slate-500">
                          <span className="flex items-center gap-1"><span className="inline-block h-2 w-2 rounded-full bg-slate-400"></span> Single-stranded (.)</span>
                          <span className="flex items-center gap-1"><span className="inline-block h-2 w-2 rounded-full bg-emerald-400"></span> Stem-paired (())</span>
                          <span className="flex items-center gap-1"><span className="inline-block h-2 w-2 rounded-full bg-amber-400"></span> Bulge / Loop</span>
                        </div>
                      </div>
                    </div>

                    {/* 3D Docking Simulation */}
                    <div>
                      <h3 className="text-[12px] font-semibold text-slate-700 mb-2">3D Tertiary Docking Map</h3>
                      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
                        <div className="flex h-32 items-center justify-center rounded bg-gradient-to-br from-slate-100 to-slate-200">
                          <div className="text-center">
                            <BarChart3 className="mx-auto h-8 w-8 text-slate-400" />
                            <p className="mt-1 text-[11px] text-slate-500">
                              FARFAR2 / RNAComposer docking visualization
                            </p>
                            <p className="text-[10px] text-slate-400">
                              RNA 3D motif docked against {gene.geneSymbol} binding cleft
                            </p>
                          </div>
                        </div>
                      </div>
                    </div>

                    {/* Mutational tolerance is not computed.
                        This was a per-base heatmap coloured by Math.random(),
                        redrawn differently on every render while presenting
                        as a structural analysis. Determining which positions
                        tolerate substitution requires either a SELEX
                        mutational scan or a folding-perturbation study, and
                        neither is wired. */}
                    <div>
                      <h3 className="text-[12px] font-semibold text-slate-700 mb-2">
                        Structural Mutation Tolerance
                      </h3>
                      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
                        <p className="text-[11.5px] leading-relaxed text-slate-600">
                          Not computed. Per-position mutation tolerance comes
                          from a SELEX mutational scan or a systematic
                          folding-perturbation study; neither is wired here.
                        </p>
                        <p className="mt-2 text-[11px] leading-relaxed text-slate-500">
                          The previous version of this panel coloured each base
                          from a random draw, so it changed on every render
                          while looking like an analysis.
                        </p>
                      </div>
                    </div>

                    {/* Off-Target Homology Scan */}
                    <div>
                      <h3 className="text-[12px] font-semibold text-slate-700 mb-2">Off-Target Structural Homology Scan</h3>
                      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
                        <div className="flex items-center gap-2 mb-2">
                          <span className="inline-flex h-2 w-2 rounded-full bg-emerald-400"></span>
                          <span className="text-[11px] text-slate-600">Transcriptome scan complete</span>
                        </div>
                        <p className="text-[11px] text-slate-500">
                          No significant structural homology detected against unintended RNA or protein partners.
                          Cross-reactivity score: <strong className="text-slate-700">Low ({selectedCandidate.targetSpecificityScore}/100)</strong>
                        </p>
                      </div>
                    </div>

                    {/* Rationale */}
                    <div>
                      <h3 className="text-[12px] font-semibold text-slate-700 mb-2">Design Rationale</h3>
                      <ul className="space-y-1">
                        {selectedCandidate.rationale.map((r, i) => (
                          <li key={i} className="flex items-start gap-2 text-[11.5px] text-slate-600">
                            <span className="mt-1 h-1.5 w-1.5 shrink-0 rounded-full bg-brand"></span>
                            {r}
                          </li>
                        ))}
                      </ul>
                    </div>
                  </div>
                </Card>
              )}

              {/* D. Downstream Action & Export Triggers */}
              {candidates.length > 0 && (
                <Card>
                  <SectionHeader step="7" title="Downstream Actions" />
                  <div className="px-6 pb-5 flex flex-wrap gap-3">
                    <button
                      onClick={() => {
                        const top = candidates[0];
                        if (top) exportFasta(top);
                      }}
                      className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors"
                    >
                      <Download className="h-4 w-4" />
                      Export Top Candidate (FASTA)
                    </button>
                    <button
                      onClick={() => {
                        const top = candidates[0];
                        if (top) exportDotBracket(top);
                      }}
                      className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors"
                    >
                      <FileText className="h-4 w-4" />
                      Export Dot-Bracket Notation
                    </button>
                    <button
                      onClick={() => setProtocolNotice("selex")}
                      className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors"
                    >
                      <FlaskConical className="h-4 w-4" />
                      Generate SELEX / Binding Assay Protocol
                    </button>
                    <button
                      onClick={() => setProtocolNotice("synthesis")}
                      className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark"
                    >
                      Proceed to Synthesis & Chemical Stabilization
                    </button>
                    {protocolNotice && (
                      <div className="mt-3 rounded border border-sky-200 bg-sky-50 px-3 py-2 text-[11.5px] leading-snug text-sky-900">
                        {protocolNotice === "selex" ? (
                          <>
                            <strong>SELEX / binding assay.</strong> Not
                            generated here. Aptamer selection is a wet-lab
                            campaign — iterative rounds of binding, partition
                            and amplification against the purified target —
                            and its parameters depend on the target protein,
                            not on anything this page holds. Reference the A25
                            rulebook&apos;s design rules for the published
                            protocol family.
                          </>
                        ) : (
                          <>
                            <strong>Synthesis &amp; stabilization.</strong> Not
                            generated here. The chemistry choices shown in the
                            form (2&apos;-F pyrimidine, 2&apos;-OMe/PS stems,
                            inverted abasic cap) are the documented
                            stabilisation options; the schedule that applies
                            depends on the selected aptamer sequence, which
                            SELEX has to produce first.
                          </>
                        )}
                        <button
                          onClick={() => setProtocolNotice(null)}
                          className="ml-2 underline underline-offset-2"
                        >
                          dismiss
                        </button>
                      </div>
                    )}
                  </div>
                </Card>
              )}
            </>
          )}
        </main>
      </div>
    </div>
  );
}
