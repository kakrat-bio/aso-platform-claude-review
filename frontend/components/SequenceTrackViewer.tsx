"use client";

import { useState } from "react";

interface OrfInfo {
  strand: string;
  frame: number;
  start: number;
  end: number;
  length: number;
  proteinLength: number;
}

interface ImmuneHit {
  motif: string;
  label: string;
  start: number;
  end: number;
}

interface RestrictionSite {
  enzyme: string;
  recognitionSite: string;
  cutPosition: number;
  strand: "+" | "-";
  overhang: "5'" | "3'" | "blunt";
}

interface MiRNATarget {
  mirnaId: string;
  seedSequence: string;
  start: number;
  end: number;
  bindingScore: number | null;
  seedGcContent?: number;
  conservationNote: string;
}

interface TrackProps {
  seqLength: number;
  orfs: OrfInfo[];
  immuneHits: ImmuneHit[];
  palindromePositions: number[];
  restrictionSites?: RestrictionSite[];
  mirnaTargets?: MiRNATarget[];
  grnaCandidates?: GrnaCandidate[];
}

interface GrnaCandidate {
  position: number;
  sequence: string;
  pam: string;
  strand: "+" | "-";
}

const ORF_COLORS = ["#6366f1", "#8b5cf6", "#a78bfa", "#c084fc", "#e879f9", "#f472b6"];
const IMMUNE_COLOR = "#f59e0b";
const PALINDROME_COLOR = "#10b981";
const RESTRICTION_COLOR = "#ec4899";
const MIRNA_COLOR = "#3b82f6";
const PAM_COLOR_FWD = "#f97316";
const PAM_COLOR_REV = "#ec4899";

export default function SequenceTrackViewer({
  seqLength,
  orfs,
  immuneHits,
  palindromePositions,
  restrictionSites = [],
  mirnaTargets = [],
  grnaCandidates = [],
}: TrackProps) {
  const [tooltip, setTooltip] = useState<{ x: number; y: number; text: string } | null>(null);

  if (seqLength < 1) return null;

  const W = 600;
  const PAD_L = 40;
  const PAD_R = 10;
  const plotW = W - PAD_L - PAD_R;
  const TRACK_H = 18;
  const TRACK_GAP = 6;
  const TRACKS: { label: string; y: number; items: React.ReactNode }[] = [];

  const xScale = (pos: number) => PAD_L + (Math.max(1, Math.min(pos, seqLength)) / seqLength) * plotW;

  function showTip(x: number, y: number, text: string) {
    setTooltip({ x, y, text });
  }

  // ORF track
  const orfItems: React.ReactNode[] = [];
  if (orfs.length > 0) {
    orfs.forEach((orf, i) => {
      const x1 = xScale(orf.start);
      const x2 = xScale(orf.end);
      const w = Math.max(x2 - x1, 2);
      const yOff = orf.strand === "-" ? TRACK_H + 2 : 0;
      const color = ORF_COLORS[i % ORF_COLORS.length];
      orfItems.push(
        <g key={`orf-${i}`}>
          <rect
            x={x1}
            y={yOff}
            width={w}
            height={TRACK_H / 2 - 1}
            rx={2}
            fill={color}
            opacity={0.8}
            className="cursor-pointer"
            onMouseEnter={(e) =>
              showTip(
                e.clientX,
                e.clientY,
                `ORF ${orf.strand} frame ${orf.frame}: ${orf.start}–${orf.end} (${orf.proteinLength} aa)`
              )
            }
            onMouseLeave={() => setTooltip(null)}
          />
        </g>
      );
    });
  }

  const orfY = 10;
  if (orfs.length > 0) {
    TRACKS.push({
      label: "ORFs",
      y: orfY,
      items: <g transform={`translate(0,${orfY})`}>{orfItems}</g>,
    });
  }

  // Immune hits track
  const immuneY = orfY + TRACK_H + TRACK_GAP + 8;
  if (immuneHits.length > 0) {
    const immuneItems = immuneHits.slice(0, 30).map((hit, i) => {
      const x1 = xScale(hit.start);
      const x2 = xScale(hit.end);
      const w = Math.max(x2 - x1, 2);
      return (
        <rect
          key={`imm-${i}`}
          x={x1}
          y={0}
          width={w}
          height={TRACK_H}
          rx={2}
          fill={IMMUNE_COLOR}
          opacity={0.7}
          className="cursor-pointer"
          onMouseEnter={(e) =>
            showTip(e.clientX, e.clientY, `${hit.motif} @ ${hit.start}–${hit.end}`)
          }
          onMouseLeave={() => setTooltip(null)}
        />
      );
    });
    TRACKS.push({
      label: "Immune",
      y: immuneY,
      items: <g transform={`translate(0,${immuneY})`}>{immuneItems}</g>,
    });
  }

  // Palindrome track
  const palY = immuneY + TRACK_H + TRACK_GAP + 8;
  if (palindromePositions.length > 0) {
    const palItems = palindromePositions.slice(0, 50).map((pos, i) => (
      <rect
        key={`pal-${i}`}
        x={xScale(pos)}
        y={0}
        width={3}
        height={TRACK_H}
        rx={1}
        fill={PALINDROME_COLOR}
        opacity={0.6}
        className="cursor-pointer"
        onMouseEnter={(e) => showTip(e.clientX, e.clientY, `Palindrome @ pos ${pos}`)}
        onMouseLeave={() => setTooltip(null)}
      />
    ));
    TRACKS.push({
      label: "Palindromes",
      y: palY,
      items: <g transform={`translate(0,${palY})`}>{palItems}</g>,
    });
  }

  // Restriction sites track
  const restrY = palY + TRACK_H + TRACK_GAP + 8;
  if (restrictionSites.length > 0) {
    const restrItems = restrictionSites.slice(0, 40).map((site, i) => {
      const x = xScale(site.cutPosition);
      return (
        <g key={`restr-${i}`}>
          <line
            x1={x}
            y1={-2}
            x2={x}
            y2={TRACK_H + 2}
            stroke={RESTRICTION_COLOR}
            strokeWidth={2}
            opacity={0.8}
          />
          <circle
            cx={x}
            cy={TRACK_H / 2}
            r={3}
            fill={RESTRICTION_COLOR}
            stroke="#fff"
            strokeWidth={1}
            className="cursor-pointer"
            onMouseEnter={(e) =>
              showTip(
                e.clientX,
                e.clientY,
                `${site.enzyme} (${site.recognitionSite}) @ pos ${site.cutPosition} [${site.strand}, ${site.overhang}]`
              )
            }
            onMouseLeave={() => setTooltip(null)}
          />
        </g>
      );
    });
    TRACKS.push({
      label: "Restriction",
      y: restrY,
      items: <g transform={`translate(0,${restrY})`}>{restrItems}</g>,
    });
  }

  // miRNA targets track
  const mirnaY = restrictionSites.length > 0
    ? restrY + TRACK_H + TRACK_GAP + 8
    : palY + TRACK_H + TRACK_GAP + 8;
  if (mirnaTargets.length > 0) {
    const mirnaItems = mirnaTargets.slice(0, 20).map((target, i) => {
      const x1 = xScale(target.start);
      const x2 = xScale(target.end);
      const w = Math.max(x2 - x1, 3);
      return (
        <rect
          key={`mirna-${i}`}
          x={x1}
          y={0}
          width={w}
          height={TRACK_H}
          rx={2}
          fill={MIRNA_COLOR}
          opacity={0.7}
          className="cursor-pointer"
          onMouseEnter={(e) =>
            showTip(
              e.clientX,
              e.clientY,
              `${target.mirnaId} seed:${target.seedSequence} @ ${target.start}–${target.end} (seed GC: ${((target.seedGcContent ?? 0) * 100).toFixed(0)}%)`
            )
          }
          onMouseLeave={() => setTooltip(null)}
        />
      );
    });
    const adjustedMirnaY = restrictionSites.length > 0
      ? restrY + TRACK_H + TRACK_GAP + 8
      : palY + TRACK_H + TRACK_GAP + 8;
    TRACKS.push({
      label: "miRNA",
      y: adjustedMirnaY,
      items: <g transform={`translate(0,${adjustedMirnaY})`}>{mirnaItems}</g>,
    });
  }

  // PAM sites track (CRISPR)
  const lastTrackY = TRACKS.length > 0 ? TRACKS[TRACKS.length - 1].y + TRACK_H : 0;
  const pamY = lastTrackY + TRACK_GAP + 8;
  if (grnaCandidates.length > 0) {
    const pamItems = grnaCandidates.slice(0, 60).map((c, i) => {
      const x = xScale(c.position);
      const isFwd = c.strand === "+";
      const color = isFwd ? PAM_COLOR_FWD : PAM_COLOR_REV;
      const arrowDy = isFwd ? 0 : TRACK_H / 2 + 2;
      return (
        <g key={`pam-${i}`}>
          <line x1={x} y1={arrowDy} x2={x} y2={arrowDy + TRACK_H / 2 - 2} stroke={color} strokeWidth={2} opacity={0.85} />
          <polygon
            points={`${x - 2},${arrowDy + TRACK_H / 2 - 2} ${x + 2},${arrowDy + TRACK_H / 2 - 2} ${x},${arrowDy + TRACK_H / 2 + 4}`}
            fill={color}
            opacity={isFwd ? 1 : 0.7}
          />
        </g>
      );
    });
    TRACKS.push({
      label: "PAM Sites",
      y: pamY,
      items: <g transform={`translate(0,${pamY})`}>{pamItems}</g>,
    });
  }

  const totalH = TRACKS.length > 0 ? TRACKS[TRACKS.length - 1].y + TRACK_H + 28 : 60;

  // X axis ticks
  const xTicks = 6;
  const tickPositions = Array.from({ length: xTicks + 1 }, (_, i) => Math.round((seqLength * i) / xTicks));

  return (
    <div className="w-full overflow-x-auto">
      <svg viewBox={`0 0 ${W} ${totalH}`} className="w-full h-auto min-w-[400px]">
        {/* X axis */}
        <line x1={PAD_L} y1={totalH - 14} x2={W - PAD_R} y2={totalH - 14} stroke="#cbd5e1" strokeWidth={0.5} />
        {tickPositions.map((p) => (
          <g key={p}>
            <line x1={xScale(p)} y1={totalH - 18} x2={xScale(p)} y2={totalH - 12} stroke="#cbd5e1" strokeWidth={0.5} />
            <text x={xScale(p)} y={totalH - 2} textAnchor="middle" className="fill-slate-400" fontSize={8}>
              {p}
            </text>
          </g>
        ))}

        {/* Track labels + tracks */}
        {TRACKS.map((track) => (
          <g key={track.label}>
            <text x={PAD_L - 4} y={track.y + TRACK_H / 2 + 3} textAnchor="end" className="fill-slate-500" fontSize={8} fontWeight={600}>
              {track.label}
            </text>
            {track.items}
          </g>
        ))}

        {TRACKS.length === 0 && (
          <text x={W / 2} y={30} textAnchor="middle" className="fill-slate-400" fontSize={10}>
            No tracks to display
          </text>
        )}
      </svg>

      {/* Tooltip */}
      {tooltip && (
        <div
          className="pointer-events-none fixed z-50 rounded-md bg-slate-800 px-2 py-1 text-[10px] text-white shadow-lg"
          style={{ left: tooltip.x + 10, top: tooltip.y - 30 }}
        >
          {tooltip.text}
        </div>
      )}

      {/* Legend */}
      <div className="mt-2 flex flex-wrap items-center gap-4 text-[10px] text-slate-500">
        <span className="flex items-center gap-1">
          <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: "#6366f1" }} />
          ORF (+)
        </span>
        <span className="flex items-center gap-1">
          <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: "#e879f9" }} />
          ORF (−)
        </span>
        <span className="flex items-center gap-1">
          <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: IMMUNE_COLOR }} />
          Immune hit
        </span>
        <span className="flex items-center gap-1">
          <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: PALINDROME_COLOR }} />
          Palindrome
        </span>
        {restrictionSites.length > 0 && (
          <span className="flex items-center gap-1">
            <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: RESTRICTION_COLOR }} />
            Restriction site
          </span>
        )}
        {mirnaTargets.length > 0 && (
          <span className="flex items-center gap-1">
            <span className="inline-block h-2 w-4 rounded-sm" style={{ backgroundColor: MIRNA_COLOR }} />
            miRNA target
          </span>
        )}
        {grnaCandidates.length > 0 && (
          <span className="flex items-center gap-1">
            <span className="inline-block h-2 w-2 rotate-45" style={{ backgroundColor: PAM_COLOR_FWD }} />
            Forward PAM
          </span>
        )}
        {grnaCandidates.length > 0 && (
          <span className="flex items-center gap-1">
            <span className="inline-block h-2 w-2 rotate-45" style={{ backgroundColor: PAM_COLOR_REV }} />
            Reverse PAM
          </span>
        )}
      </div>
    </div>
  );
}
