"use client";

import { TranslationalCandidate } from "@/types/translational";
import { Dna } from "lucide-react";

interface TranslationFeatureMapProps {
  targetElement: string;
  selectedCandidate: TranslationalCandidate | null;
}

const ELEMENT_LABELS: Record<string, string> = {
  "5p_utr": "5′ UTR / m⁷G Cap / Kozak",
  "3p_utr_mirna": "3′ UTR miRNA Seed Site",
  "uorf": "5′ UTR uORF / Upstream AUG",
  "structured_element": "Structured Element (IRES / G-quadruplex)",
};

export default function TranslationFeatureMap({
  targetElement,
  selectedCandidate,
}: TranslationFeatureMapProps) {
  const elementLabel = ELEMENT_LABELS[targetElement] ?? targetElement;

  const mapWidth = 800;
  const trackHeight = 40;

  const regions = [
    { name: "5′ UTR", start: 0, width: 200, color: "bg-slate-200" },
    { name: "Kozak", start: 185, width: 50, color: "bg-amber-200" },
    { name: "Start Codon (AUG)", start: 235, width: 25, color: "bg-red-300" },
    { name: "CDS", start: 260, width: 400, color: "bg-slate-100" },
    { name: "3′ UTR", start: 660, width: 140, color: "bg-slate-200" },
  ];

  const getElementColor = (el: string): string => {
    if (el === "5p_utr") return "bg-indigo-300";
    if (el === "3p_utr_mirna") return "bg-purple-300";
    if (el === "uorf") return "bg-amber-300";
    if (el === "structured_element") return "bg-emerald-300";
    return "bg-sky-300";
  };

  const elementColor = getElementColor(targetElement);

  const candidatePos = selectedCandidate
    ? Math.max(0, Math.min(100, ((selectedCandidate.rank ?? 1) - 1) * 25))
    : 40;

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <p className="text-[11px] font-bold uppercase tracking-wider text-slate-400">
          Target mRNA Linear Feature Map
        </p>
        <p className="text-[10px] text-slate-400">
          {elementLabel} · {selectedCandidate ? `Candidate #${(selectedCandidate.rank ?? 1)}` : "Top candidate shown"}
        </p>
      </div>

      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
        <div
          className="relative h-10 w-full rounded-sm"
          style={{ maxWidth: `${mapWidth}px` }}
        >
          {regions.map((r) => (
            <div
              key={r.name}
              className={`absolute top-0 h-6 rounded-sm ${r.color} flex items-end justify-center`}
              style={{ left: `${(r.start / mapWidth) * 100}%`, width: `${(r.width / mapWidth) * 100}%` }}
              title={`${r.name} (${r.start}-${r.start + r.width} nt)`}
            >
              <span className="mb-1 text-[8px] font-medium text-slate-600">
                {r.name}
              </span>
            </div>
          ))}

          <div
            className={`absolute top-6 h-3 w-4 rounded-sm ${elementColor} ring-1 ring-indigo-400 ring-offset-1`}
            style={{ left: `${candidatePos}%`, width: "8%" }}
            title={`ASO binding site (${selectedCandidate?.targetRegion ?? elementLabel})`}
          />

          <div className="absolute top-6 text-[8px] text-slate-400" style={{ left: "0%" }}>
            nt 1
          </div>
          <div className="absolute top-6 text-[8px] text-slate-400" style={{ right: "0%", textAlign: "right" }}>
            5′→3′ transcript
          </div>
        </div>

        <div className="mt-4 flex items-center gap-4 text-[11px]">
          <div className="flex items-center gap-1.5">
            <div className={`h-3 w-3 rounded-sm ${elementColor}`} />
            <span className="text-slate-500">ASO binding site</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="h-3 w-3 rounded-sm bg-red-300" />
            <span className="text-slate-500">Start codon (AUG)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="h-3 w-3 rounded-sm bg-amber-200" />
            <span className="text-slate-500">Kozak consensus</span>
          </div>
        </div>
      </div>

      {selectedCandidate && (
        <div className="rounded-lg border border-indigo-200 bg-indigo-50 p-3">
          <div className="mb-1 flex items-center gap-2">
            <Dna className="h-3.5 w-3.5 text-indigo-500" />
            <p className="text-[10px] font-bold uppercase tracking-wider text-indigo-500">
              Selected ASO (5′→3′)
            </p>
          </div>
          <div className="flex flex-wrap gap-px font-mono text-[11px]">
            {selectedCandidate.sequence.split("").map((nt, i) => (
              <span
                key={i}
                className="inline-flex h-5 w-4 items-center justify-center rounded-sm text-[9px] font-bold text-white"
                style={{
                  backgroundColor:
                    nt === "A" ? "#60a5fa" : nt === "C" ? "#3b82f6" : nt === "G" ? "#4ade80" : "#fbbf24",
                }}
              >
                {nt}
              </span>
            ))}
          </div>
          <p className="mt-1.5 text-[10px] text-slate-500">
            Target: {selectedCandidate.targetRegion} · Chemistry: {selectedCandidate.chemistry} ·
            {" ΔG: "}{selectedCandidate.realMetrics.targetDuplexEnergy.toFixed(1)} kcal/mol
          </p>
        </div>
      )}
    </div>
  );
}
