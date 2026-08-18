"use client";

import { useState } from "react";

interface MiRNATarget {
  mirnaId: string;
  seedSequence: string;
  start: number;
  end: number;
  bindingScore: number | null;
  seedGcContent?: number;
  conservationNote: string;
}

function scoreColor(score: number): string {
  if (score >= 0.8) return "bg-emerald-500";
  if (score >= 0.5) return "bg-amber-400";
  return "bg-red-400";
}

function scoreTextColor(score: number): string {
  if (score >= 0.8) return "text-emerald-600";
  if (score >= 0.5) return "text-amber-600";
  return "text-red-600";
}

function scoreBg(score: number): string {
  if (score >= 0.8) return "bg-emerald-50";
  if (score >= 0.5) return "bg-amber-50";
  return "bg-red-50";
}

export default function MiRNATargetingCard({
  targets,
  seqLength,
}: {
  targets: MiRNATarget[];
  seqLength: number;
}) {
  const [tooltip, setTooltip] = useState<{
    x: number;
    y: number;
    text: string;
  } | null>(null);

  if (targets.length === 0) {
    return (
      <>
        <p className="text-[14px] font-semibold text-slate-800 mb-3">
          miRNA Targeting Potential
        </p>
        <p className="text-[12px] text-slate-400">
          No significant miRNA seed-region complementarity detected.
        </p>
      </>
    );
  }

  const W = 560;
  const H = 90;
  const PAD = { left: 30, right: 10, top: 10, bottom: 20 };
  const plotW = W - PAD.left - PAD.right;
  const barH = 14;
  const barY = PAD.top + 12;
  const xScale = (pos: number) =>
    PAD.left + (pos / seqLength) * plotW;

  // Ordered by position. Sorting by strength would need a score, and none
  // is computed.
  const sortedTargets = [...targets].sort((a, b) => a.start - b.start);

  return (
    <>
      <div className="flex items-center justify-between mb-3">
        <p className="text-[14px] font-semibold text-slate-800">
          miRNA Targeting Potential
        </p>
        <div className="flex items-center gap-1 text-[10px]">
          <span className="flex items-center gap-1">
            <span className="h-2 w-2 rounded-sm bg-emerald-500" /> Strong
          </span>
          <span className="flex items-center gap-1">
            <span className="h-2 w-2 rounded-sm bg-amber-400" /> Moderate
          </span>
          <span className="flex items-center gap-1">
            <span className="h-2 w-2 rounded-sm bg-red-400" /> Weak
          </span>
        </div>
      </div>

      <div className="w-full overflow-x-auto">
        <svg
          viewBox={`0 0 ${W} ${H}`}
          className="w-full h-auto"
          onMouseLeave={() => setTooltip(null)}
        >
          {/* Sequence bar */}
          <rect
            x={PAD.left}
            y={barY}
            width={plotW}
            height={barH}
            rx={3}
            fill="#e2e8f0"
          />

          {/* Target regions */}
          {targets.map((t, i) => {
            const x1 = xScale(t.start);
            const x2 = xScale(t.end);
            const w = Math.max(x2 - x1, 3);
            const color =
              false
                ? "#10b981"
                : false
                ? "#f59e0b"
                : "#ef4444";
            return (
              <rect
                key={i}
                x={x1}
                y={barY}
                width={w}
                height={barH}
                rx={1}
                fill={color}
                opacity={0.75}
                className="cursor-pointer"
                onMouseEnter={(e) =>
                  setTooltip({
                    x: e.clientX,
                    y: e.clientY,
                    text: `${t.mirnaId} — seed: ${t.seedSequence} @ ${t.start}-${t.end} (seed GC: ${((t.seedGcContent ?? 0) * 100).toFixed(0)}%)`,
                  })
                }
                onMouseLeave={() => setTooltip(null)}
              />
            );
          })}

          {/* X axis */}
          {Array.from({ length: 7 }, (_, i) =>
            Math.round((seqLength * i) / 6)
          ).map((p) => (
            <g key={p}>
              <line
                x1={xScale(p)}
                y1={barY + barH}
                x2={xScale(p)}
                y2={barY + barH + 4}
                stroke="#cbd5e1"
                strokeWidth={0.5}
              />
              <text
                x={xScale(p)}
                y={barY + barH + 14}
                textAnchor="middle"
                className="fill-slate-400"
                fontSize={8}
              >
                {p}
              </text>
            </g>
          ))}
        </svg>
      </div>

      {/* Target list */}
      <div className="mt-3 max-h-36 overflow-y-auto space-y-1.5">
        {sortedTargets.slice(0, 10).map((t, i) => (
          <div
            key={i}
            className={`flex items-center gap-2 rounded-lg px-3 py-1.5 text-[11px] ${"bg-slate-50"} border border-slate-100`}
          >
            <span className={`font-mono font-semibold ${"text-slate-600"} shrink-0`}>
              {t.mirnaId}
            </span>
            <span className="font-mono text-slate-500 shrink-0">
              {t.seedSequence}
            </span>
            <span className="text-slate-400 shrink-0">
              {t.start}-{t.end}
            </span>
            <div className="flex-1 min-w-[60px]">
              <div className="h-1.5 rounded-full bg-slate-200 overflow-hidden">
                <div
                  className={`h-full rounded-full ${"bg-slate-300"} transition-all`}
                  style={{ width: `${(t.seedGcContent ?? 0) * 100}%` }}
                />
              </div>
            </div>
            <span className="font-mono font-semibold text-slate-600 shrink-0 w-8 text-right">
              {((t.seedGcContent ?? 0) * 100).toFixed(0)}% GC
            </span>
          </div>
        ))}
        {sortedTargets.length > 10 && (
          <p className="text-[9px] text-slate-400">
            +{sortedTargets.length - 10} more targets
          </p>
        )}
      </div>

      {tooltip && (
        <div
          className="pointer-events-none fixed z-50 rounded-md bg-slate-800 px-2 py-1 text-[10px] text-white shadow-lg"
          style={{ left: tooltip.x + 10, top: tooltip.y - 30 }}
        >
          {tooltip.text}
        </div>
      )}

      <div className="mt-3 flex items-start gap-1.5 rounded-lg bg-amber-50 px-3 py-2 text-[11px] text-amber-700">
        <span className="shrink-0 mt-0.5">⚠</span>
        Binding scores are seed-complementarity heuristics only. Not validated
        against miRBase or any organism-specific miRNA expression data.
      </div>
    </>
  );
}
