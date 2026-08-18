"use client";

import { useState, useMemo } from "react";
import { ChevronUp, ChevronDown, Search, Download } from "lucide-react";

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

interface ImmuneHit {
  motif: string;
  label: string;
  start: number;
  end: number;
}

interface OrfHit {
  strand: string;
  frame: number;
  start: number;
  end: number;
  length: number;
  proteinLength: number;
}

interface Hairpin {
  start: number;
  end: number;
  stemLength: number;
  loopSize: number;
  stabilityScore: number;
  type: "hairpin" | "bulge" | "internal_loop";
}

interface TableRow {
  type: string;
  typeColor: string;
  name: string;
  start: number;
  end: number;
  length: number;
  score: number | null;
  strand: string;
  details: string;
}

const TYPE_COLORS: Record<string, string> = {
  Restriction: "#ec4899",
  "miRNA Target": "#8b5cf6",
  "Immune Motif": "#f59e0b",
  ORF: "#3b82f6",
  Structure: "#10b981",
};

function downloadCsv(rows: TableRow[]) {
  const header = "Type,Name,Start,End,Length,Score,Strand,Details\n";
  const body = rows
    .map(
      (r) =>
        `"${r.type}","${r.name}",${r.start},${r.end},${r.length},${r.score ?? ""},${r.strand},"${r.details}"`
    )
    .join("\n");
  const blob = new Blob([header + body], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "feature-table.csv";
  a.click();
  URL.revokeObjectURL(url);
}

export default function FeatureTable({
  restrictionSites = [],
  mirnaTargets = [],
  immuneHits = [],
  orfs = [],
  hairpins = [],
}: {
  restrictionSites?: RestrictionSite[];
  mirnaTargets?: MiRNATarget[];
  immuneHits?: ImmuneHit[];
  orfs?: OrfHit[];
  hairpins?: Hairpin[];
}) {
  const [sortCol, setSortCol] = useState<keyof TableRow>("start");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");
  const [filterType, setFilterType] = useState<string | null>(null);
  const [search, setSearch] = useState("");

  const allRows: TableRow[] = useMemo(() => {
    const rows: TableRow[] = [];

    restrictionSites.forEach((s) =>
      rows.push({
        type: "Restriction",
        typeColor: TYPE_COLORS["Restriction"],
        name: s.enzyme,
        start: s.cutPosition,
        end: s.cutPosition,
        length: s.recognitionSite.length,
        score: null,
        strand: s.strand,
        details: `${s.recognitionSite} · ${s.overhang} overhang`,
      })
    );

    mirnaTargets.forEach((t) =>
      rows.push({
        type: "miRNA Target",
        typeColor: TYPE_COLORS["miRNA Target"],
        name: t.mirnaId,
        start: t.start,
        end: t.end,
        length: t.end - t.start + 1,
        score: t.seedGcContent ?? null,
        strand: "+",
        details: `seed: ${t.seedSequence}`,
      })
    );

    immuneHits.forEach((h) =>
      rows.push({
        type: "Immune Motif",
        typeColor: TYPE_COLORS["Immune Motif"],
        name: h.motif,
        start: h.start,
        end: h.end,
        length: h.end - h.start + 1,
        score: null,
        strand: "+",
        details: h.label.slice(0, 60) + (h.label.length > 60 ? "..." : ""),
      })
    );

    orfs.forEach((o) =>
      rows.push({
        type: "ORF",
        typeColor: TYPE_COLORS["ORF"],
        name: `ORF ${o.strand} f${o.frame}`,
        start: o.start,
        end: o.end,
        length: o.end - o.start + 1,
        score: null,
        strand: o.strand,
        details: `Frame ${o.frame} · ${o.proteinLength} aa`,
      })
    );

    hairpins.forEach((h) =>
      rows.push({
        type: "Structure",
        typeColor: TYPE_COLORS["Structure"],
        name: h.type.replace("_", " "),
        start: h.start,
        end: h.end,
        length: h.end - h.start + 1,
        score: h.stabilityScore,
        strand: "+",
        details: `stem: ${h.stemLength} bp · loop: ${h.loopSize} nt`,
      })
    );

    return rows;
  }, [restrictionSites, mirnaTargets, immuneHits, orfs, hairpins]);

  const filteredRows = useMemo(() => {
    let rows = allRows;
    if (filterType) rows = rows.filter((r) => r.type === filterType);
    if (search) {
      const q = search.toLowerCase();
      rows = rows.filter(
        (r) =>
          r.name.toLowerCase().includes(q) ||
          r.details.toLowerCase().includes(q)
      );
    }
    rows.sort((a, b) => {
      let va = a[sortCol];
      let vb = b[sortCol];
      if (va === null || va === undefined) va = "";
      if (vb === null || vb === undefined) vb = "";
      if (typeof va === "number" && typeof vb === "number") {
        return sortDir === "asc" ? va - vb : vb - va;
      }
      return sortDir === "asc"
        ? String(va).localeCompare(String(vb))
        : String(vb).localeCompare(String(va));
    });
    return rows;
  }, [allRows, filterType, search, sortCol, sortDir]);

  const typeCounts = useMemo(() => {
    const counts: Record<string, number> = {};
    allRows.forEach((r) => {
      counts[r.type] = (counts[r.type] || 0) + 1;
    });
    return counts;
  }, [allRows]);

  function toggleSort(col: keyof TableRow) {
    if (sortCol === col) {
      setSortDir(sortDir === "asc" ? "desc" : "asc");
    } else {
      setSortCol(col);
      setSortDir("asc");
    }
  }

  function SortIcon({ col }: { col: keyof TableRow }) {
    if (sortCol !== col) return <span className="text-slate-300 ml-0.5">↕</span>;
    return sortDir === "asc" ? (
      <ChevronUp className="h-3 w-3 inline ml-0.5" />
    ) : (
      <ChevronDown className="h-3 w-3 inline ml-0.5" />
    );
  }

  return (
    <>
      <div className="flex items-center justify-between mb-3">
        <p className="text-[14px] font-semibold text-slate-800">
          Feature Table
        </p>
        <div className="flex items-center gap-2">
          <span className="text-[10px] text-slate-400">
            {filteredRows.length} of {allRows.length} features
          </span>
          <button
            onClick={() => downloadCsv(filteredRows)}
            className="flex items-center gap-1 rounded-md border border-[#E5E7EB] px-2 py-1 text-[10px] font-medium text-slate-500 hover:bg-slate-50"
          >
            <Download className="h-3 w-3" /> CSV
          </button>
        </div>
      </div>

      {/* Filter chips + search */}
      <div className="flex flex-wrap items-center gap-2 mb-3">
        <button
          onClick={() => setFilterType(null)}
          className={`rounded-full px-2.5 py-0.5 text-[10px] font-medium transition-colors ${
            filterType === null
              ? "bg-brand text-white"
              : "bg-slate-100 text-slate-500 hover:bg-slate-200"
          }`}
        >
          All ({allRows.length})
        </button>
        {Object.entries(typeCounts).map(([type, count]) => (
          <button
            key={type}
            onClick={() => setFilterType(filterType === type ? null : type)}
            className={`inline-flex items-center gap-1 rounded-full px-2.5 py-0.5 text-[10px] font-medium transition-colors ${
              filterType === type ? "text-white" : "bg-slate-100 text-slate-500 hover:bg-slate-200"
            }`}
            style={filterType === type ? { backgroundColor: TYPE_COLORS[type] } : undefined}
          >
            <span className="h-1.5 w-1.5 rounded-full" style={{ backgroundColor: TYPE_COLORS[type] }} />
            {type} ({count})
          </button>
        ))}
        <div className="relative ml-auto">
          <Search className="absolute left-2 top-1/2 -translate-y-1/2 h-3 w-3 text-slate-400" />
          <input
            type="text"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            placeholder="Search..."
            className="rounded-md border border-[#E5E7EB] bg-white pl-6 pr-2 py-1 text-[11px] text-slate-600 placeholder:text-slate-400 focus:border-brand focus:outline-none focus:ring-1 focus:ring-brand/20 w-40"
          />
        </div>
      </div>

      {/* Table */}
      <div className="overflow-x-auto rounded-lg border border-[#E5E7EB]">
        <table className="w-full text-[11px]">
          <thead>
            <tr className="bg-slate-50 text-left text-slate-500">
              {(
                [
                  ["type", "Type"],
                  ["name", "Name"],
                  ["start", "Start"],
                  ["end", "End"],
                  ["length", "Length"],
                  ["score", "Score"],
                  ["strand", "Strand"],
                  ["details", "Details"],
                ] as [keyof TableRow, string][]
              ).map(([col, label]) => (
                <th
                  key={col}
                  onClick={() => toggleSort(col)}
                  className="px-3 py-2 font-medium cursor-pointer select-none hover:bg-slate-100 whitespace-nowrap"
                >
                  {label}
                  <SortIcon col={col} />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {filteredRows.slice(0, 100).map((row, i) => (
              <tr
                key={i}
                className="border-t border-slate-100 hover:bg-slate-50/50"
              >
                <td className="px-3 py-1.5">
                  <span
                    className="inline-flex items-center rounded-full px-2 py-0.5 text-[9px] font-semibold text-white"
                    style={{ backgroundColor: row.typeColor }}
                  >
                    {row.type}
                  </span>
                </td>
                <td className="px-3 py-1.5 font-mono font-medium text-slate-700">
                  {row.name}
                </td>
                <td className="px-3 py-1.5 font-mono text-slate-600">{row.start}</td>
                <td className="px-3 py-1.5 font-mono text-slate-600">{row.end}</td>
                <td className="px-3 py-1.5 text-slate-600">{row.length}</td>
                <td className="px-3 py-1.5 font-mono text-slate-600">
                  {row.score !== null ? row.score.toFixed(2) : "—"}
                </td>
                <td className="px-3 py-1.5 text-slate-600">{row.strand}</td>
                <td className="px-3 py-1.5 text-slate-500 max-w-[200px] truncate">
                  {row.details}
                </td>
              </tr>
            ))}
            {filteredRows.length === 0 && (
              <tr>
                <td colSpan={8} className="px-3 py-6 text-center text-slate-400">
                  No features match your filter.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {filteredRows.length > 100 && (
        <p className="mt-2 text-[9px] text-slate-400 text-center">
          Showing first 100 of {filteredRows.length} results. Use filters or export CSV for the full list.
        </p>
      )}
    </>
  );
}
