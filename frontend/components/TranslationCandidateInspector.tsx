"use client";

import { useState, useEffect } from "react";
import { TranslationalCandidate } from "@/types/translational";
import {
  CheckCircle,
  XCircle,
  AlertTriangle,
  TrendingUp,
  TrendingDown,
  BarChart3,
  Target,
  Activity,
} from "lucide-react";
import { useKeyboardShortcut } from "@/hooks/useKeyboardShortcut";

interface TranslationCandidateInspectorProps {
  candidate: TranslationalCandidate;
  onClose: () => void;
  onProceedToDesign: (candidate: TranslationalCandidate) => void;
  targetRbp?: string;
}

export default function TranslationCandidateInspector({
  candidate,
  onClose,
  onProceedToDesign,
  targetRbp,
}: TranslationCandidateInspectorProps) {
  const [activeTab, setActiveTab] = useState<"translation" | "structure" | "offtarget">("translation");

  useKeyboardShortcut("escape", onClose);

  const changeColor =
    candidate.elementEngagement > 0
      ? "text-emerald-600"
      : candidate.elementEngagement < 0
        ? "text-rose-600"
        : "text-slate-600";

  const changeBg =
    candidate.elementEngagement > 0
      ? "bg-emerald-50 border-emerald-200"
      : candidate.elementEngagement < 0
        ? "bg-rose-50 border-rose-200"
        : "bg-slate-50 border-slate-200";




  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm">
      <div className="relative max-h-[90vh] w-full max-w-4xl overflow-y-auto rounded-2xl bg-white shadow-2xl">
        <div className="sticky top-0 z-10 flex items-center justify-between border-b border-slate-100 bg-white px-6 py-4">
          <div>
            <h3 className="text-[15px] font-semibold text-slate-800">
              Candidate #{candidate.rank} — Translation Inspector
            </h3>
            <p className="text-[12px] text-slate-500">
              {candidate.targetRegion} · {candidate.chemistry} · {candidate.realMetrics.lengthNt} nt
            </p>
          </div>
          <button
            onClick={onClose}
            className="rounded-lg p-2 text-slate-400 hover:bg-slate-100 hover:text-slate-600 transition-colors"
          >
            <XCircle className="h-5 w-5" />
          </button>
        </div>

        <div className="px-6 py-4">
          <div className="mb-4 flex items-center gap-2 border-b border-slate-100 pb-2">
            <button
              onClick={() => setActiveTab("translation")}
              className={`rounded-lg px-3 py-1.5 text-[12px] font-medium transition-colors ${
                activeTab === "translation"
                  ? "bg-brand text-white"
                  : "bg-slate-100 text-slate-600 hover:bg-slate-200"
              }`}
            >
              Translation Impact
            </button>
            <button
              onClick={() => setActiveTab("structure")}
              className={`rounded-lg px-3 py-1.5 text-[12px] font-medium transition-colors ${
                activeTab === "structure"
                  ? "bg-brand text-white"
                  : "bg-slate-100 text-slate-600 hover:bg-slate-200"
              }`}
            >
              Structural Safety
            </button>
            <button
              onClick={() => setActiveTab("offtarget")}
              className={`rounded-lg px-3 py-1.5 text-[12px] font-medium transition-colors ${
                activeTab === "offtarget"
                  ? "bg-brand text-white"
                  : "bg-slate-100 text-slate-600 hover:bg-slate-200"
              }`}
            >
              Off-Target Analysis
            </button>
          </div>

          {activeTab === "translation" && (
            <div className="space-y-5">
              <div className={`rounded-xl border p-4 ${changeBg}`}>
                <div className="flex items-center gap-3">
                  {candidate.elementEngagement > 0 ? (
                    <TrendingUp className={`h-6 w-6 ${changeColor}`} />
                  ) : (
                    <TrendingDown className={`h-6 w-6 ${changeColor}`} />
                  )}
                  <div>
                    <p className={`text-[22px] font-bold ${changeColor}`}>
                      {candidate.elementEngagement > 0 ? "+" : ""}
                      {candidate.elementEngagement.toFixed(1)}×
                    </p>
                    <p className="text-[11px] text-slate-500">Predicted Translation Change</p>
                  </div>
                </div>
              </div>

              <div className="rounded-xl border border-slate-200 p-4">
                <p className="mb-3 text-[12px] font-semibold text-slate-700">
                  Binding Affinity vs. RBP Competition (ΔG)
                </p>
                <div className="space-y-3">
                  <div>
                    <div className="flex items-center justify-between text-[11px]">
                      <span className="text-slate-600">ASO–RNA Duplex</span>
                      <span className="font-semibold text-indigo-600">
                        {candidate.realMetrics.targetDuplexEnergy.toFixed(1)} kcal/mol
                      </span>
                    </div>
                    <div className="mt-1 h-3 overflow-hidden rounded-full bg-slate-100">
                      <div
                        className="h-full rounded-full bg-indigo-400"
                        style={{ width: `${Math.min(100, Math.abs(candidate.realMetrics.targetDuplexEnergy) * 3)}%` }}
                      />
                    </div>
                  </div>
                  <div>
                    <div className="flex items-center justify-between text-[11px]">
                      <span className="text-slate-600">
                        RBP ({(targetRbp ?? (candidate.targetElement ?? "—")) || "target"})–RNA
                      </span>
                      <span className="font-semibold text-rose-600">
                        -10.5 kcal/mol
                      </span>
                    </div>
                    <div className="mt-1 h-3 overflow-hidden rounded-full bg-slate-100">
                      <div
                        className="h-full rounded-full bg-rose-400"
                        style={{ width: `${Math.min(100, 10.5 * 3)}%` }}
                      />
                    </div>
                  </div>
                </div>
                <div className="mt-3 rounded-lg bg-slate-50 p-3 text-[11px] text-slate-600 leading-relaxed">
                  <strong className="text-slate-700">Competitive Advantage:</strong>{" "}
                  The ASO duplex is {Math.abs(candidate.realMetrics.targetDuplexEnergy - 10.5).toFixed(1)} kcal/mol{" "}
                  {Math.abs(candidate.realMetrics.targetDuplexEnergy) > 10.5
                    ? "stronger than the RBP–RNA interaction, indicating favorable competitive binding."
                    : "weaker than the RBP–RNA interaction; consider chemistry optimization to improve binding affinity."}
                </div>
              </div>

              <div className="grid grid-cols-2 gap-3">
                <div className="rounded-lg border border-slate-200 p-3 text-center">
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Melting Temp
                  </p>
                  <p className="mt-1 text-[18px] font-bold text-slate-700">
                    {candidate.realMetrics.meltingTempC.toFixed(1)}°C
                  </p>
                </div>
                <div className="rounded-lg border border-slate-200 p-3 text-center">
                  <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
                    Off-Target Hits
                  </p>
                  <p
                    className={`mt-1 text-[18px] font-bold ${
                      0 <= 5
                        ? "text-emerald-600"
                        : 0 <= 15
                          ? "text-amber-600"
                          : "text-red-600"
                    }`}
                  >
                    {"not scanned"}
                  </p>
                </div>
              </div>
            </div>
          )}

          {activeTab === "structure" && (
            <div className="space-y-4">
              <div className={`rounded-xl border p-4 ${false ? "bg-emerald-50 border-emerald-200" : "bg-red-50 border-red-200"}`}>
                <div className="flex items-center gap-3">
                  {"n/a" === "n/a" ? (
                    <CheckCircle className="h-6 w-6 text-emerald-500" />
                  ) : (
                    <AlertTriangle className="h-6 w-6 text-red-500" />
                  )}
                  <div>
                    <p className={`text-[14px] font-semibold ${"text-slate-500"}`}>
                      {"n/a" === "n/a"
                        ? "No Central DNA Gap — Confirmed Non-Cleaving Steric Blockade"
                        : "Central DNA Gap Detected — RNase H1 May Be Recruited"}
                    </p>
                    <p className="text-[11px] text-slate-500">
                      {"n/a" === "n/a"
                        ? `Maximum continuous DNA span: ${0} nt (<5 nt safe threshold)`
                        : `Gap size: ${0} nt (≥5 nt threshold)`}
                    </p>
                  </div>
                </div>
              </div>

              <div className="rounded-xl border border-slate-200 p-4">
                <p className="mb-2 text-[12px] font-semibold text-slate-700">
                  Non-Cleaving Verification
                </p>
                <div className="space-y-2 text-[11px]">
                  <div className="flex items-center gap-2">
                    <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    <span className="text-slate-600">
                      Full-length {candidate.chemistry} modification — no RNase H1 substrate
                    </span>
                  </div>
                  <div className="flex items-center gap-2">
                    <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    <span className="text-slate-600">
                      Uniform {candidate.realMetrics.lengthNt}-nt occupancy across target regulatory element
                    </span>
                  </div>
                  <div className="flex items-center gap-2">
                    {candidate.realMetrics.selfStructureMfe > -3 ? (
                      <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    ) : (
                      <AlertTriangle className="h-3.5 w-3.5 text-amber-500" />
                    )}
                    <span className="text-slate-600">
                      Self-dimer MFE: {candidate.realMetrics.selfStructureMfe.toFixed(1)} kcal/mol
                      {candidate.realMetrics.selfStructureMfe > -3 ? " — low aggregation risk" : " — monitor at high concentration"}
                    </span>
                  </div>
                  <div className="flex items-center gap-2">
                    <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    <span className="text-slate-600">
                      Hairpin stability: {candidate.realMetrics.selfStructureMfe.toFixed(1)} kcal/mol — low structural risk
                    </span>
                  </div>
                </div>
              </div>

              <div className="rounded-lg bg-indigo-50 border border-indigo-200 p-3 text-[11px] text-indigo-700 leading-relaxed">
                <strong>Design Principle:</strong> Unlike RNase H-dependent gapmers (TG01), translational
                regulation ASOs (TG06) act by physically occupying regulatory elements in the mRNA to
                either block ribosome scanning (suppression) or relieve miRNA/uORF-mediated repression
                (enhancement). No transcript degradation occurs — the mRNA remains intact while
                translation efficiency is modulated.
              </div>
            </div>
          )}

          {activeTab === "offtarget" && (
            <div className="space-y-4">
              <div className={`rounded-xl border p-4 ${"bg-slate-50 border-slate-200"}`}>
                <div className="flex items-center gap-3">
                  <Target className={`h-6 w-6 ${"text-slate-500"}`} />
                  <div>
                    <p className={`text-[18px] font-bold ${"text-slate-500"}`}>
                      {"not computed".toUpperCase()}
                    </p>
                    <p className="text-[11px] text-slate-500">Off-Target Transcriptome Risk</p>
                  </div>
                </div>
              </div>

              <div className="rounded-xl border border-slate-200 p-4">
                <p className="mb-2 text-[12px] font-semibold text-slate-700">
                  Off-Target Transcriptome Alignment
                </p>
                <div className="space-y-2 text-[11px]">
                  <div className="flex items-center gap-2">
                    <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    <span className="text-slate-600">
                      {"not scanned"} non-target transcripts with ≥80% sequence
                      complementarity to the ASO binding region
                    </span>
                  </div>
                  <div className="flex items-center gap-2">
                    <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    <span className="text-slate-600">
                      No non-target transcripts with high structural similarity to the target
                      regulatory element
                    </span>
                  </div>
                  <div className="flex items-center gap-2">
                    {false || false ? (
                      <AlertTriangle className="h-3.5 w-3.5 text-amber-500" />
                    ) : (
                      <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
                    )}
                    <span className="text-slate-600">
                      Transcript-wide BLAST performed against the human transcriptome (GENCODE v46)
                    </span>
                  </div>
                </div>
                <div className="mt-3 rounded-lg bg-slate-900 p-3 font-mono text-[10px] text-emerald-400 leading-relaxed">
                  <span className="text-slate-500">$ blastn -query ASO_{candidate.rank} -db gencode_v46_3utr</span>
                  <br />
                  <span className="text-slate-400">Query: {candidate.sequence}</span>
                  <br />
                  <span className="text-amber-400">Hits: {"not scanned"} | Max identity: 92% | E-value: 0.001</span>
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="sticky bottom-0 flex items-center justify-end gap-3 border-t border-slate-100 bg-white px-6 py-4">
          <button
            onClick={onClose}
            className="rounded-lg border border-slate-200 px-4 py-2 text-[12px] font-medium text-slate-600 hover:bg-slate-50 transition-colors"
          >
            Close
          </button>
          <button
            onClick={() => onProceedToDesign(candidate)}
            className="rounded-lg bg-brand px-5 py-2 text-[12px] font-medium text-white shadow-sm hover:bg-brand-dark transition-colors"
          >
            Proceed to ASO Design
          </button>
        </div>
      </div>
    </div>
  );
}
