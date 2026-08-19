"use client";

import { IsoformCandidate } from "@/types/isoformEngineering";

/**
 * Summary of a TG07 steric-blocking ASO panel.
 *
 * This used to render five charts — splice efficiency, predicted isoform
 * yield, CAI vs uridine content, MFE, and TLR innate-immune risk. Four of
 * those plotted numbers the backend produced from a loop index rather than
 * from the transcript, and CAI / uridine / TLR describe a delivered mRNA
 * construct, which is TG08's modality, not an antisense oligonucleotide.
 * They are gone along with the fields behind them.
 *
 * What is plotted now is what is measured: duplex free energy against the
 * real transcript window, melting temperature, GC and position along the
 * exon.
 */
export default function IsoformAnalysisDashboard({
  candidates,
}: {
  candidates: IsoformCandidate[];
}) {
  if (!candidates.length) return null;

  const dgs = candidates
    .map((c) => c.targetDuplexDg)
    .filter((v): v is number => v !== null);
  const tms = candidates
    .map((c) => c.meltingTempC)
    .filter((v): v is number => v !== null);

  const mean = (xs: number[]) =>
    xs.length ? xs.reduce((s, x) => s + x, 0) / xs.length : null;
  const fmt = (v: number | null, digits = 1, suffix = "") =>
    v === null ? "—" : `${v.toFixed(digits)}${suffix}`;

  const meanGc = mean(candidates.map((c) => c.gcContent));
  const meanTm = mean(tms);
  const meanDg = mean(dgs);
  const strongest = dgs.length ? Math.min(...dgs) : null;
  const allInFrame = candidates.every((c) => c.inFrameStatus === "In-Frame");

  // Bar length is scaled against the strongest binder in the panel, so the
  // chart is a comparison within this design set and nothing more.
  const dgFloor = dgs.length ? Math.min(...dgs) : -1;

  return (
    <div id="isoform-analysis-dashboard" className="space-y-4">
      <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 lg:grid-cols-5">
        <Stat label="Candidates" value={String(candidates.length)} />
        <Stat label="Mean GC" value={fmt(meanGc, 1, "%")} />
        <Stat label="Mean Tm" value={fmt(meanTm, 1, " °C")} />
        <Stat label="Mean Duplex ΔG" value={fmt(meanDg, 2)} />
        <Stat
          label="Strongest ΔG"
          value={fmt(strongest, 2)}
          tone="text-emerald-600"
        />
      </div>

      <div className="rounded-xl border border-blue-100 bg-blue-50 p-4">
        <p className="text-[12px] font-semibold text-blue-800 mb-1">Key Findings</p>
        <ul className="space-y-1 text-[11.5px] text-blue-700">
          <li>
            • Top candidate <strong>{candidates[0].constructId}</strong> binds its
            target window at{" "}
            <strong>{fmt(candidates[0].targetDuplexDg, 2, " kcal/mol")}</strong>{" "}
            with Tm <strong>{fmt(candidates[0].meltingTempC, 1, " °C")}</strong>{" "}
            and GC <strong>{candidates[0].gcContent.toFixed(1)}%</strong>.
          </li>
          <li>
            • Target: exon <strong>{candidates[0].exonNumber}</strong> (
            {candidates[0].exonLength} nt), {candidates[0].targetWindow}.
          </li>
          <li>
            • Skipping this exon is{" "}
            <strong>{allInFrame ? "in-frame" : "out-of-frame"}</strong>{" "}
            ({candidates[0].exonLength} nt
            {allInFrame ? " is" : " is not"} a multiple of 3).
          </li>
          <li>
            • Ordering is thermodynamic. It is not a validated activity model —
            see the caveat below the candidate table.
          </li>
        </ul>
      </div>

      <div className="rounded-xl border border-slate-200 bg-white p-5">
        <p className="text-[12px] font-semibold text-slate-700 mb-3">
          Target Duplex ΔG by Candidate (more negative binds more tightly)
        </p>
        <div className="space-y-1.5">
          {candidates.map((c) => {
            const v = c.targetDuplexDg;
            const pct = v === null || dgFloor === 0 ? 0 : (v / dgFloor) * 100;
            return (
              <div key={c.constructId} className="flex items-center gap-2">
                <span className="w-40 shrink-0 truncate font-mono text-[10px] text-slate-500">
                  {c.constructId}
                </span>
                <div className="h-3 flex-1 overflow-hidden rounded-full bg-slate-100">
                  <div
                    className="h-full rounded-full bg-brand/70"
                    style={{ width: `${Math.max(0, Math.min(100, pct))}%` }}
                  />
                </div>
                <span className="w-16 shrink-0 text-right font-mono text-[10px] text-slate-600">
                  {fmt(v, 1)}
                </span>
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
}

function Stat({
  label,
  value,
  tone = "text-slate-800",
}: {
  label: string;
  value: string;
  tone?: string;
}) {
  return (
    <div className="rounded-xl border border-slate-200 bg-white p-3 text-center">
      <p className="text-[10px] font-semibold uppercase tracking-wider text-slate-400">
        {label}
      </p>
      <p className={`text-[18px] font-bold ${tone}`}>{value}</p>
    </div>
  );
}
