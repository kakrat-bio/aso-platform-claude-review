"use client";

import { useState, useMemo } from "react";
import { GrnaCandidate } from "@/types/upload";
import InfoTooltip from "@/components/InfoTooltip";

interface ColDef {
  key: string;
  label: string;
  get: (c: GrnaCandidate, idx: number) => React.ReactNode;
  align?: "left" | "center" | "right";
  tooltip?: string;
}

function specificityBadge(score: number) {
  const bg =
    score >= 50
      ? "bg-emerald-100 text-emerald-700"
      : score >= 30
        ? "bg-amber-100 text-amber-700"
        : "bg-rose-100 text-rose-700";
  return <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-bold ${bg}`}>{score.toFixed(0)}</span>;
}

function efficiencyBadge(score: number) {
  const bg =
    score >= 70
      ? "bg-emerald-100 text-emerald-700"
      : score >= 40
        ? "bg-amber-100 text-amber-700"
        : "bg-rose-100 text-rose-700";
  return <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-bold ${bg}`}>{score.toFixed(0)}</span>;
}

function mismatchDist(dist: number[]) {
  if (!dist || dist.length === 0) return <span className="text-slate-400">—</span>;
  return (
    <span className="font-mono text-slate-600">
      {dist.map((n, i) => (
        <span key={i}>
          {i > 0 && " - "}
          {n}
        </span>
      ))}
    </span>
  );
}

const COLS: ColDef[] = [
  {
    key: "rank",
    label: "Rank",
    align: "left",
    get: (_c, i) => <span className="font-bold text-brand">#{i + 1}</span>,
  },
  {
    key: "position",
    label: "Position/Strand",
    align: "right",
    tooltip: "Position of the PAM on the forward strand of the reference genome",
    get: (c) => (
      <span className="font-mono text-slate-600">
        {c.position}+{c.strand === "-" ? "R" : ""}
      </span>
    ),
  },
  {
    key: "sequence",
    label: "Guide Sequence",
    align: "left",
    get: (c) => (
      <div className="flex items-center gap-2">
        <span className="font-mono font-semibold text-slate-700">{c.sequence}</span>
        <span className="font-mono text-slate-500">{c.pam}</span>
        {c.polyT && (
          <span className="inline-flex items-center rounded bg-rose-50 px-1.5 py-0.5 text-[9px] font-medium text-rose-600">
            Poly-T
          </span>
        )}
      </div>
    ),
  },
  {
    key: "score",
    label: "Overall",
    align: "center",
    get: (c) => (
      <span
        className={`inline-flex items-center justify-center rounded-full px-2 py-0.5 text-[10px] font-bold ${
          c.score >= 70
            ? "bg-emerald-100 text-emerald-700"
            : c.score >= 40
              ? "bg-amber-100 text-amber-700"
              : "bg-rose-100 text-rose-700"
        }`}
      >
        {c.score.toFixed(1)}
      </span>
    ),
  },
  {
    key: "specificity",
    label: "Specificity",
    align: "center",
    tooltip:
      "Predicted off-target risk. >50 = green (good), >30 = yellow (use with caution), <30 = red (avoid).",
    get: (c) => specificityBadge(c.specificityScore),
  },
  {
    key: "efficiency",
    label: "Efficiency",
    align: "center",
    tooltip:
      "Predicted on-target activity (simplified Doench-style score). >=70 = high, 40-69 = moderate, <40 = low.",
    get: (c) => efficiencyBadge(c.efficiencyScore),
  },
  {
    // Was "Off-targets", showing a count derived from 6-mer repetitiveness
    // without any genome alignment. The underlying statistic is kept; the
    // claim that it counts off-target sites is not.
    key: "repetitiveness",
    label: "Internal Repetitiveness",
    align: "center",
    get: (c) => (
      <span className="text-slate-600" title="Fraction of 6-mers that repeat within the spacer. Not an off-target screen — no genome alignment is performed.">
        {c.internalRepetitiveness.toFixed(3)}
      </span>
    ),
  },
  {
    key: "mismatchdist",
    label: "Mismatch Distribution",
    align: "center",
    tooltip:
      "Number of predicted genome sites at 0, 1, 2, 3, 4 mismatches. Lower numbers at low mismatch levels indicate higher specificity.",
    get: (c) => mismatchDist(c.mismatchDistribution),
  },
  {
    key: "gc",
    label: "GC%",
    align: "center",
    get: (c) => (
      <span className={`text-slate-600 ${c.gc >= 40 && c.gc <= 80 ? "font-semibold text-emerald-600" : ""}`}>
        {c.gc.toFixed(1)}%
      </span>
    ),
  },
  {
    key: "selfcomp",
    label: "Self-comp",
    align: "center",
    get: (c) => <span className="text-slate-600">{(c.selfComplementarity * 100).toFixed(1)}%</span>,
  },
];

type SortKey = "score" | "position" | "gc" | "internalRepetitiveness" | "selfComplementarity" | "specificityScore" | "efficiencyScore";
type SortDir = "asc" | "desc";

export default function CrisprCandidateTable({ candidates }: { candidates: GrnaCandidate[] }) {
  const n = candidates.length;
  if (n === 0) return null;

  const [sortKey, setSortKey] = useState<SortKey>("score");
  const [sortDir, setSortDir] = useState<SortDir>("desc");

  const sorted = useMemo(() => {
    const arr = [...candidates];
    arr.sort((a, b) => {
      const av = a[sortKey] as number;
      const bv = b[sortKey] as number;
      return sortDir === "asc" ? av - bv : bv - av;
    });
    return arr;
  }, [candidates, sortKey, sortDir]);

  const handleSort = (key: SortKey) => {
    if (sortKey === key) {
      setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    } else {
      setSortKey(key);
      setSortDir("desc");
    }
  };

  const arrow = (key: SortKey) => {
    if (sortKey !== key) return null;
    return sortDir === "asc" ? " ▲" : " ▼";
  };

  const sortableMap: Record<string, SortKey> = {
    position: "position",
    score: "score",
    repetitiveness: "internalRepetitiveness",
    specificity: "specificityScore",
    efficiency: "efficiencyScore",
    gc: "gc",
    selfcomp: "selfComplementarity",
  };

  return (
    <div className="w-full overflow-x-auto">
      <table className="w-full border-collapse text-[11px]">
        <thead>
          <tr className="bg-slate-100 text-slate-600">
            {COLS.map((col) => (
              <th
                key={col.key}
                className={`whitespace-nowrap px-2.5 py-2 font-semibold ${
                  col.align === "left" ? "text-left" : col.align === "right" ? "text-right" : "text-center"
                } ${col.key !== "rank" && col.key !== "sequence" ? "cursor-pointer select-none hover:text-slate-800" : ""}`}
                onClick={() => {
                  if (col.key === "rank" || col.key === "sequence" || !sortableMap[col.key]) return;
                  handleSort(sortableMap[col.key]);
                }}
                title={col.tooltip ? col.tooltip : col.key !== "rank" && col.key !== "sequence" ? "Click to sort" : undefined}
              >
                <span className="flex items-center justify-center gap-0.5">
                  {col.label}
                  {col.tooltip && <InfoTooltip content={col.tooltip} />}
                  {arrow(sortableMap[col.key] as SortKey)}
                </span>
              </th>
            ))}
          </tr>
        </thead>
        <tbody className="divide-y divide-slate-100 text-slate-700">
          {sorted.map((c, i) => (
            <tr key={c.id} className="hover:bg-slate-50 transition-colors">
              {COLS.map((col) => (
                <td
                  key={col.key}
                  className={`whitespace-nowrap px-2.5 py-1.5 ${
                    col.align === "left" ? "text-left" : col.align === "right" ? "text-right" : "text-center"
                  }`}
                  title={col.key === "sequence" ? c.sequence : undefined}
                >
                  {col.get(c, candidates.indexOf(c))}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <div className="mt-2 flex flex-wrap items-center gap-x-5 gap-y-1 text-[10px] text-slate-400">
        <span>
          Specificity: <span className="font-semibold text-emerald-600">&gt;50</span> (green),{" "}
          <span className="font-semibold text-amber-600">&gt;30</span> (yellow),{" "}
          <span className="font-semibold text-rose-600">&lt;30</span> (red)
        </span>
        <span>
          Efficiency: <span className="font-semibold text-emerald-600">&ge;70</span> (high),{" "}
          <span className="font-semibold text-amber-600">40–69</span> (moderate),{" "}
          <span className="font-semibold text-rose-600">&lt;40</span> (low)
        </span>
        <span>Click column headers to sort. Default: score descending.</span>
      </div>
    </div>
  );
}
