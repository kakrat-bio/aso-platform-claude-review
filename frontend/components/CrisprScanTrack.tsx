"use client";

import { useState } from "react";
import { GrnaCandidate } from "@/types/upload";

interface CrisprScanTrackProps {
  candidates: GrnaCandidate[];
  seqLength: number;
}

const SCORE_COLORS: Record<string, string> = {
  emerald: "#10b981",
  amber: "#f59e0b",
  rose: "#f43f5e",
};

export default function CrisprScanTrack({ candidates, seqLength }: CrisprScanTrackProps) {
  const [tooltip, setTooltip] = useState<{ x: number; y: number; text: string } | null>(null);

  if (seqLength < 1) return null;

  const W = 640;
  const PAD_L = 50;
  const PAD_R = 20;
  const PAD_TOP = 24;
  const PAD_BOT = 40;
  const plotW = W - PAD_L - PAD_R;
  const plotH = 140;
  const maxScore = 100;

  const xScale = (pos: number) => PAD_L + (Math.max(1, Math.min(pos, seqLength)) / seqLength) * plotW;
  const yScale = (score: number) => PAD_TOP + plotH - (Math.max(0, Math.min(score, maxScore)) / maxScore) * plotH;

  const xTicks = 6;
  const tickPositions = Array.from({ length: xTicks + 1 }, (_, i) => Math.round((seqLength * i) / xTicks));

  const showTip = (x: number, y: number, text: string) => setTooltip({ x, y, text });
  const hideTip = () => setTooltip(null);

  return (
    <div className="w-full overflow-x-auto">
      <svg viewBox={`0 0 ${W} ${PAD_TOP + plotH + PAD_BOT}`} className="w-full h-auto min-w-[400px]">
        <rect x={0} y={0} width={W} height={PAD_TOP + plotH + PAD_BOT} fill="#f8fafc" rx={8} />
        <rect x={PAD_L} y={PAD_TOP} width={plotW} height={plotH} fill="#ffffff" stroke="#e2e8f0" strokeWidth={0.5} rx={4} />

        <text x={PAD_L} y={PAD_TOP - 6} className="fill-slate-500" fontSize={9} fontWeight={600}>
          Score
        </text>

        <text x={W / 2} y={PAD_TOP - 6} textAnchor="middle" className="fill-slate-500" fontSize={9} fontWeight={600}>
          gRNA Position
        </text>

        <line x1={PAD_L} y1={PAD_TOP} x2={PAD_L} y2={PAD_TOP + plotH} stroke="#cbd5e1" strokeWidth={0.5} />
        <line x1={PAD_L} y1={PAD_TOP + plotH} x2={W - PAD_R} y2={PAD_TOP + plotH} stroke="#cbd5e1" strokeWidth={0.5} />

        <text x={PAD_L - 6} y={PAD_TOP + plotH + 4} textAnchor="end" className="fill-slate-400" fontSize={8}>0</text>
        <text x={PAD_L - 6} y={PAD_TOP + 4} textAnchor="end" className="fill-slate-400" fontSize={8}>100</text>

        <line x1={PAD_L} y1={PAD_TOP} x2={W - PAD_R} y2={PAD_TOP} stroke="#e2e8f0" strokeWidth={0.3} strokeDasharray="2 2" />

        {candidates.map((c, i) => {
          const x = xScale(c.position);
          const y = yScale(c.score);
          const barH = PAD_TOP + plotH - y;
          const color = SCORE_COLORS[c.color] ?? "#6366f1";

          return (
            <g key={c.id}>
              <rect
                x={x - 3}
                y={y}
                width={6}
                height={Math.max(barH, 4)}
                rx={2}
                fill={color}
                opacity={0.85}
                className="cursor-pointer"
                onMouseEnter={(e) =>
                  showTip(
                    e.clientX,
                    e.clientY,
                    `#${i + 1} | pos ${c.position} | ${c.sequence} | ${c.pam} | score ${c.score.toFixed(1)} | repetitiveness ${c.internalRepetitiveness.toFixed(3)}`
                  )
                }
                onMouseLeave={hideTip}
              />
              <polygon
                points={`${x},${PAD_TOP + plotH + 2} ${x - 3},${PAD_TOP + plotH + 10} ${x + 3},${PAD_TOP + plotH + 10}`}
                fill={color}
                opacity={0.6}
                className="cursor-pointer"
                onMouseEnter={(e) =>
                  showTip(
                    e.clientX,
                    e.clientY,
                    `#${i + 1} | pos ${c.position} | ${c.sequence} | ${c.pam} | score ${c.score.toFixed(1)} | repetitiveness ${c.internalRepetitiveness.toFixed(3)}`
                  )
                }
                onMouseLeave={hideTip}
              />
            </g>
          );
        })}

        {tickPositions.map((p) => (
          <g key={p}>
            <line x1={xScale(p)} y1={PAD_TOP + plotH} x2={xScale(p)} y2={PAD_TOP + plotH + 4} stroke="#94a3b8" strokeWidth={0.5} />
            <text x={xScale(p)} y={PAD_TOP + plotH + 14} textAnchor="middle" className="fill-slate-400" fontSize={8}>
              {p}
            </text>
          </g>
        ))}
      </svg>

      {tooltip && (
        <div
          className="pointer-events-none fixed z-50 rounded-md bg-slate-800 px-2 py-1.5 text-[10px] text-white shadow-lg max-w-[280px]"
          style={{ left: tooltip.x + 10, top: tooltip.y - 30 }}
        >
          {tooltip.text}
        </div>
      )}

      <div className="mt-2 flex flex-wrap items-center gap-4 text-[10px] text-slate-500">
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-3 rounded-sm" style={{ backgroundColor: SCORE_COLORS.emerald }} />
          Score ≥ 70 (high)
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-3 rounded-sm" style={{ backgroundColor: SCORE_COLORS.amber }} />
          Score 40–69 (moderate)
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-3 rounded-sm" style={{ backgroundColor: SCORE_COLORS.rose }} />
          Score &lt; 40 (low)
        </span>
      </div>
    </div>
  );
}
