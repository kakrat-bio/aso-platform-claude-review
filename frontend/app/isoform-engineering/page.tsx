"use client";

import { useEffect, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import {
  AlertCircle,
  ArrowLeft,
  Loader2,
  Dna,
  Download,
  FlaskConical,
  FileText,
  X,
  Copy,
  Check,
} from "lucide-react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import { Card, SectionHeader, FieldLabel, InfoField, Pill } from "@/components/ui";
import IsoformAnalysisDashboard from "@/components/IsoformAnalysisDashboard";
import { GeneTargetObject } from "@/types/gene";
import { saveReport } from "@/lib/auth";
import {
  IsoformEngineeringInputs,
  IsoformEngineeringResponse,
  DesignOptions,
  IsoformCandidate,
} from "@/types/isoformEngineering";
import {
  fetchDesignOptions,
  generateConstructs,
} from "@/lib/isoformEngineeringApi";

const CONFIRMED_TARGET_KEY = "aso:confirmedTarget";

export default function IsoformEngineeringPage() {
  const router = useRouter();
  const [gene, setGene] = useState<GeneTargetObject | null>(null);
  const [options, setOptions] = useState<DesignOptions | null>(null);
  const [optionsError, setOptionsError] = useState<string | null>(null);

  const [targetSymbol, setTargetSymbol] = useState("");
  const [isoformGoal, setIsoformGoal] = useState("");
  const [targetExonLocus, setTargetExonLocus] = useState("");
  const [spliceElementTarget, setSpliceElementTarget] = useState("");
  const [stericChemistry, setStericChemistry] = useState("");
  const [enforceInFrame, setEnforceInFrame] = useState(true);

  const [results, setResults] = useState<IsoformEngineeringResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const [selectedCandidate, setSelectedCandidate] = useState<IsoformCandidate | null>(null);
  const [copied, setCopied] = useState(false);
  const selectedCandidateRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (selectedCandidate && selectedCandidateRef.current) {
      selectedCandidateRef.current.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  }, [selectedCandidate]);

  useEffect(() => {
    const stored = sessionStorage.getItem(CONFIRMED_TARGET_KEY);
    if (stored) {
      try {
        const g = JSON.parse(stored) as GeneTargetObject;
        setGene(g);
        setTargetSymbol(g.geneSymbol || "");
      } catch {
        setGene(null);
      }
    }
    fetchDesignOptions()
      .then(setOptions)
      .catch((e) => setOptionsError(e instanceof Error ? e.message : "Failed to load options."));
  }, []);

  function handleReset() {
    setResults(null);
    setSelectedCandidate(null);
    setError(null);
  }

  async function handleGenerate() {
    if (!targetSymbol.trim()) return;
    setLoading(true);
    setError(null);
    setSelectedCandidate(null);
    try {
      const res = await generateConstructs({
        targetSymbol: targetSymbol.trim(),
        isoformGoal,
        targetExonLocus,
        spliceElementTarget,
        stericChemistry,
        enforceInFrame,
      });
      setResults(res);
      saveReport({
        step: "isoform_engineering",
        title: `Isoform Engineering: ${targetSymbol.toUpperCase()}`,
        geneSymbol: targetSymbol.toUpperCase(),
        disease: gene?.disease || "",
        summary: `Generated ${res.candidates.length} isoform engineering constructs for ${targetSymbol.toUpperCase()}. Top: ${res.candidates[0]?.constructId || "N/A"}.`,
        data: {
          isoformGoal,
          topCandidates: res.candidates.slice(0, 3).map((c) => ({
            constructId: c.constructId,
            sequence: c.sequence,
            targetDuplexDg: c.targetDuplexDg,
            meltingTempC: c.meltingTempC,
            inFrameStatus: c.inFrameStatus,
          })),
        },
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : "Generation failed.");
    } finally {
      setLoading(false);
    }
  }

  function handleSelectCandidate(candidate: IsoformCandidate) {
    setSelectedCandidate(candidate);
  }

  function copySequence(seq: string) {
    navigator.clipboard.writeText(seq).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    });
  }

  function isFormValid(): boolean {
    return (
      !!targetSymbol.trim() &&
      !!isoformGoal &&
      !!targetExonLocus &&
      !!spliceElementTarget &&
      !!stericChemistry
    );
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
    <div className="flex min-h-screen bg-[#F5F6FA]">
      <Sidebar />
      <div className="flex min-h-screen flex-1 flex-col">
        <Topbar />
        <main className="flex-1 space-y-5 px-6 py-6">
          <Card className="flex items-center gap-3 px-5 py-3">
            <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-gradient-to-br from-indigo-50 to-blue-50">
              <Dna className="h-4.5 w-4.5 text-indigo-500" />
            </span>
            <div>
              <p className="text-[13px] font-semibold text-slate-800">
                {gene.geneSymbol}{" "}
                <span className="font-normal text-slate-400">· {gene.organism}</span>
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

          {!results && (
            <Card>
              <SectionHeader step="1" title="Isoform Engineering Design" />
              <p className="px-6 pb-3 text-[12.5px] text-slate-500">
                Configure RNA payload parameters for isoform engineering. Select your isoform goal, target exon locus, splice element, and steric chemistry to generate optimized constructs.
              </p>
              {optionsError && (
                <p className="px-6 pb-2 text-[12.5px] text-red-600">{optionsError}</p>
              )}
              <div className="grid grid-cols-1 gap-4 px-6 pb-4 md:grid-cols-2 lg:grid-cols-3">
                <div>
                  <FieldLabel hint="The gene symbol for the protein you want to engineer (e.g., DMD, CFTR, SMN2)">
                    Target Gene Symbol <span className="text-red-500">*</span>
                  </FieldLabel>
                  <input
                    value={targetSymbol}
                    onChange={(e) => setTargetSymbol(e.target.value)}
                    placeholder="e.g. DMD, CFTR, SMN2"
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  />
                </div>

                <div>
                  <FieldLabel hint="What type of isoform switching do you want to achieve?">
                    Isoform Goal <span className="text-red-500">*</span>
                  </FieldLabel>
                  <select
                    value={isoformGoal}
                    onChange={(e) => setIsoformGoal(e.target.value)}
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  >
                    <option value="">Select isoform goal</option>
                    {options?.isoformGoals.map((o) => (
                      <option key={o.id} value={o.id}>
                        {o.label}
                      </option>
                    ))}
                  </select>
                  {isoformGoal && (
                    <p className="mt-1 text-[10.5px] text-slate-400">
                      {options?.isoformGoals.find((o) => o.id === isoformGoal)?.description}
                    </p>
                  )}
                </div>

                <div>
                  <FieldLabel hint="The exon locus you want to target for isoform modulation">
                    Target Exon Locus <span className="text-red-500">*</span>
                  </FieldLabel>
                  <select
                    value={targetExonLocus}
                    onChange={(e) => setTargetExonLocus(e.target.value)}
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  >
                    <option value="">Select exon locus</option>
                    {options?.targetExonLoci.map((o) => (
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
                    value={spliceElementTarget}
                    onChange={(e) => setSpliceElementTarget(e.target.value)}
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  >
                    <option value="">Select splice element</option>
                    {options?.spliceElementTargets.map((o) => (
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
                    value={stericChemistry}
                    onChange={(e) => setStericChemistry(e.target.value)}
                    className="w-full rounded-lg border border-slate-300 bg-white py-2.5 px-3 text-[13.5px] text-slate-700 focus:border-brand focus:outline-none focus:ring-2 focus:ring-brand/20"
                  >
                    <option value="">Select chemistry</option>
                    {options?.stericChemistries.map((o) => (
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
                      checked={enforceInFrame}
                      onChange={(e) => setEnforceInFrame(e.target.checked)}
                      className="h-4 w-4 rounded border-slate-300 text-brand focus:ring-brand"
                    />
                    <span className="text-[13px] text-slate-700">Enforce In-Frame splicing</span>
                  </label>
                </div>
              </div>

              <div className="flex justify-end px-6 pb-5">
                <button
                  onClick={handleGenerate}
                  disabled={!isFormValid() || loading}
                  className="flex items-center gap-2 rounded-lg bg-brand px-5 py-2.5 text-[13.5px] font-medium text-white shadow-sm transition-colors hover:bg-brand-dark disabled:cursor-not-allowed disabled:opacity-50"
                >
                  {loading && <Loader2 className="h-4 w-4 animate-spin" />}
                  {loading ? "Generating Constructs..." : "Generate Constructs / Candidates"}
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

          {results && (
            <>
              <Card>
                <SectionHeader
                  step="2"
                  title="Header & Construct Overview"
                  right={
                    <button
                      onClick={handleReset}
                      className="text-[12px] font-medium text-slate-500 hover:text-slate-700"
                    >
                      Modify Parameters
                    </button>
                  }
                />
                <div className="px-6 pb-5">
                  <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                    <div className="space-y-3">
                      <InfoField label="Target Gene & Canonical Transcript" value={`${results.overview.targetGene} · ${results.overview.refSeq ?? "—"}`} />
                      <InfoField label="Selected Isoform Goal" value={isoformGoal.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase())} />
                      <InfoField label="Transcript Length" value={`${results.overview.transcriptLength.toLocaleString()} nt · ${results.overview.exonCount} exons`} />
                    </div>
                    <div className="space-y-3">
                      <InfoField label="Target Exon" value={`Exon ${results.overview.targetExon} · ${results.overview.exonLength} nt`} />
                      <InfoField
                        label="Reading Frame If Skipped"
                        value={results.overview.inFrameStatus}
                        valueClassName={results.overview.inFrameStatus === "In-Frame" ? "text-emerald-600" : "text-amber-600"}
                      />
                      <InfoField
                        label="Splice-Site Strength"
                        value={results.overview.spliceSiteStrength !== null ? results.overview.spliceSiteStrength.toFixed(3) : "Not fetched for this exon"}
                        valueClassName={results.overview.spliceSiteStrength === null ? "text-slate-400" : undefined}
                      />
                    </div>
                    <div className="space-y-3">
                      <InfoField label="Primary Mechanism Assigned" value={results.overview.primaryMechanism} />
                      <InfoField label="Design Window" value={`${results.overview.targetWindow} (${results.overview.windowStart}–${results.overview.windowEnd})`} />
                    </div>
                  </div>
                </div>
              </Card>

              <Card>
                <SectionHeader step="3" title="Candidate Construct Table" />
                <div className="px-5 pb-5 overflow-x-auto">
                  <table className="w-full text-left">
                    <thead>
                      <tr className="border-b-2 border-slate-200">
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Rank</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Construct ID</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Modality</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">ASO Sequence (5′→3′)</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Len</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">GC</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Tm</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Duplex ΔG</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400 text-right">Self MFE</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Transcript Pos</th>
                        <th className="pb-2 pr-3 text-[10px] font-bold uppercase tracking-wider text-slate-400">Frame If Skipped</th>
                        <th className="pb-2 text-[10px] font-bold uppercase tracking-wider text-slate-400">Action</th>
                      </tr>
                    </thead>
                    <tbody>
                      {results.candidates.map((c) => (
                        <tr
                          key={c.constructId}
                          className={`border-b border-slate-100 last:border-0 hover:bg-slate-50/50 transition-colors cursor-pointer ${selectedCandidate?.constructId === c.constructId ? "bg-brand/5" : ""}`}
                          onClick={() => handleSelectCandidate(c)}
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
                          <td className="py-3 pr-3">
                            <span className="text-[10.5px] font-mono text-slate-700">{c.sequence}</span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className="text-[11.5px] font-mono text-slate-700">{c.length}</span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className="text-[11.5px] font-semibold text-slate-700">{c.gcContent.toFixed(1)}%</span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className="text-[11.5px] font-mono text-slate-700">
                              {c.meltingTempC !== null ? `${c.meltingTempC.toFixed(1)}°C` : "—"}
                            </span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className="text-[11.5px] font-mono text-slate-700">
                              {c.targetDuplexDg !== null ? c.targetDuplexDg.toFixed(1) : "—"}
                            </span>
                          </td>
                          <td className="py-3 pr-3 text-right">
                            <span className="text-[11.5px] font-mono text-slate-700">
                              {c.selfMfe !== null ? c.selfMfe.toFixed(1) : "—"}
                            </span>
                          </td>
                          <td className="py-3 pr-3">
                            <span className="text-[11px] text-slate-600">{c.transcriptStart}–{c.transcriptEnd}</span>
                          </td>
                          <td className="py-3 pr-3">
                            <span className={`inline-flex items-center rounded-full border px-2 py-0.5 text-[10px] font-semibold ${c.inFrameStatus === "In-Frame" ? "border-emerald-200 bg-emerald-50 text-emerald-700" : "border-amber-200 bg-amber-50 text-amber-700"}`}>
                              {c.inFrameStatus}
                            </span>
                          </td>
                          <td className="py-3">
                            <button
                              onClick={(e) => { e.stopPropagation(); handleSelectCandidate(c); }}
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

              <Card className="p-5">
                <SectionHeader step="3a" title="Candidate Analysis & Visualizations" />
                <div className="px-5 pb-5">
                  <IsoformAnalysisDashboard candidates={results.candidates} />
                </div>
              </Card>

              {selectedCandidate && (
                <div ref={selectedCandidateRef}>
                <Card className="overflow-hidden">
                  <div className="flex items-center justify-between px-6 pt-4 pb-3 border-b border-slate-100">
                    <SectionHeader step="3b" title={`Inspection: ${selectedCandidate.constructId}`} />
                    <button onClick={() => setSelectedCandidate(null)} className="flex h-6 w-6 items-center justify-center rounded-full border border-slate-200 text-slate-400 hover:text-slate-600">
                      <X className="h-3.5 w-3.5" />
                    </button>
                  </div>
                  <div className="px-6 pb-5 space-y-5">
                    <div>
                      <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">Target Location</p>
                      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4 text-[11px] text-slate-700 space-y-1">
                        <p>
                          Exon <span className="font-semibold">{selectedCandidate.exonNumber}</span> ({selectedCandidate.exonLength} nt),{" "}
                          {selectedCandidate.targetWindow}
                        </p>
                        <p className="font-mono text-[10.5px] text-slate-500">
                          transcript positions {selectedCandidate.transcriptStart}–{selectedCandidate.transcriptEnd}
                        </p>
                      </div>
                    </div>

                    <div>
                      <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">ASO Sequence (5&apos; &rarr; 3&apos;)</p>
                      <div className="flex items-center gap-2">
                        <div className="flex-1 rounded-lg bg-slate-900 px-4 py-2.5 font-mono text-[11px] text-emerald-400 break-all leading-relaxed select-all">
                          {selectedCandidate.sequence}
                        </div>
                        <button onClick={() => { navigator.clipboard.writeText(selectedCandidate.sequence); setCopied(true); setTimeout(() => setCopied(false), 1500); }} className="flex shrink-0 items-center gap-1 rounded-md border border-slate-200 px-2.5 py-1.5 text-[11px] text-slate-500 hover:bg-white transition-colors">
                          {copied ? <><Check className="h-3 w-3 text-emerald-500" /> Copied</> : <><Copy className="h-3 w-3" /> Copy</>}
                        </button>
                      </div>
                      <p className="mt-1.5 text-[10px] text-slate-500">
                        Reverse complement of transcript window{" "}
                        <span className="font-mono">{selectedCandidate.targetSequence}</span>
                      </p>
                    </div>

                    <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                      <div className="rounded-lg border border-slate-200 bg-white p-3">
                        <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">Melting Temperature</p>
                        <p className="text-[14px] font-bold text-slate-700">
                          {selectedCandidate.meltingTempC !== null ? `${selectedCandidate.meltingTempC.toFixed(1)} °C` : "—"}
                        </p>
                        <p className="text-[10px] text-slate-500 mt-0.5">primer3, SantaLucia nearest-neighbour (unmodified backbone)</p>
                      </div>
                      <div className="rounded-lg border border-slate-200 bg-white p-3">
                        <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">Target Duplex &Delta;G</p>
                        <p className="text-[14px] font-bold text-slate-700">
                          {selectedCandidate.targetDuplexDg !== null ? `${selectedCandidate.targetDuplexDg.toFixed(2)} kcal/mol` : "—"}
                        </p>
                        <p className="text-[10px] text-slate-500 mt-0.5">ViennaRNA duplexfold against the real transcript window</p>
                      </div>
                      <div className="rounded-lg border border-slate-200 bg-white p-3">
                        <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400 mb-1">Splice-Site Strength</p>
                        <p className={`text-[14px] font-bold ${selectedCandidate.spliceSiteStrength !== null ? "text-slate-700" : "text-slate-400"}`}>
                          {selectedCandidate.spliceSiteStrength !== null ? selectedCandidate.spliceSiteStrength.toFixed(3) : "Not fetched"}
                        </p>
                        <p className="text-[10px] text-slate-500 mt-0.5">From Ensembl flanking genomic sequence, when available</p>
                      </div>
                    </div>

                    {selectedCandidate.notComputed && Object.keys(selectedCandidate.notComputed).length > 0 && (
                      <div>
                        <p className="text-[11px] font-semibold uppercase tracking-wider text-slate-500 mb-2">Not Computed</p>
                        <div className="rounded-lg border border-amber-200 bg-amber-50 p-3 space-y-1.5">
                          {Object.entries(selectedCandidate.notComputed).map(([field, reason]) => (
                            <p key={field} className="text-[10.5px] text-amber-900">
                              <span className="font-semibold font-mono">{field}</span>: {reason}
                            </p>
                          ))}
                        </div>
                      </div>
                    )}

                    {results.ranking && (
                      <p className="text-[10px] text-slate-500 leading-relaxed">
                        Ranked by <span className="font-mono">{results.ranking.orderedBy}</span>. {results.ranking.caveat}
                      </p>
                    )}
                  </div>
                </Card>
                </div>
              )}

              <Card>
                <SectionHeader step="4" title="Downstream Action & Export" />
                <div className="px-6 pb-5 flex flex-wrap items-center gap-3">
                  <button onClick={() => { const blob = new Blob([results.candidates.map((c) => `>${c.constructId}\n${c.sequence}`).join("\n\n")], { type: "text/plain" }); const url = URL.createObjectURL(blob); const a = document.createElement("a"); a.href = url; a.download = `${targetSymbol}_isoform_constructs.fasta`; a.click(); URL.revokeObjectURL(url); }} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                    <Download className="h-4 w-4" /> Export Candidate Sequences to FASTA
                  </button>
                  <button onClick={() => alert("Splice Modulation & ASO Synthesis module is under development.")} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                    <FlaskConical className="h-4 w-4" /> Proceed to Splice Modulation & ASO Synthesis
                  </button>
                  <button onClick={() => alert("IVT Plasmid Template & Synthesis Protocol generation is under development.")} className="flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-4 py-2.5 text-[13px] font-medium text-slate-700 hover:bg-slate-50 transition-colors">
                    <FileText className="h-4 w-4" /> Generate IVT Plasmid Template & Synthesis Protocol
                  </button>
                </div>
              </Card>
            </>
          )}
        </main>
      </div>
    </div>
  );
}
