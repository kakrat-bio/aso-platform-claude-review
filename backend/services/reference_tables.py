"""Static reference tables for features that cannot be computed from sequence.

Specified in `docs/planning/data_sources_halted_flagged.md`. See
`backend/data/reference/README.md` for why every table currently ships with
its header row only.

THE CONTRACT
------------
A table that is missing, header-only, or has no row for the gene returns
None. Callers turn that into an UNRESOLVED feature, which halts the mechanism
or withholds the flag. Absence is never zero and never a negative finding:
"this gene is not in a sixty-row catalogue of known repeat expansions" is a
real answer, but it is a different answer from "this gene has no repeat".

No external API is called at request time. Every source here is a periodic
download, because runtime calls would put latency, availability and silent
version drift into a system whose whole claim is reproducibility.
"""

from __future__ import annotations

import csv
import os
import threading

REFERENCE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "reference",
)

# Table name -> the column its rows are keyed by.
TABLES: dict[str, str] = {
    "repeat_expansion_loci": "gene_symbol",
    "tissue_expression": "gene_symbol",
    "dominant_negative_genes": "gene_symbol",
    "protein_localisation": "gene_symbol",
    "clingen_dosage": "gene_symbol",
    "rbp_repressor_sites": "gene_symbol",
    "polyadenylation_sites": "gene_symbol",
    "apa_therapeutic_benefit": "gene_symbol",
    "alt_promoters": "gene_symbol",
    "alt_promoter_benefit": "gene_symbol",
    "intron_retention_potential": "gene_symbol",
    "intron_retention_benefit": "gene_symbol",
}

_cache: dict[str, dict[str, list[dict]]] = {}
_lock = threading.Lock()


def _load(name: str) -> dict[str, list[dict]]:
    """Read one TSV into {key: [row, ...]}. Missing or empty yields {}."""
    path = os.path.join(REFERENCE_DIR, f"{name}.tsv")
    index: dict[str, list[dict]] = {}
    if not os.path.exists(path):
        return index
    key_col = TABLES[name]
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row.get(key_col) or "").strip().upper()
            if key:
                index.setdefault(key, []).append(row)
    return index


def rows_for(table: str, gene_symbol: str | None) -> list[dict]:
    """Every row for a gene, or [] when the table or the gene is absent."""
    if table not in TABLES or not gene_symbol:
        return []
    with _lock:
        if table not in _cache:
            _cache[table] = _load(table)
    return _cache[table].get(gene_symbol.strip().upper(), [])


def row_for(table: str, gene_symbol: str | None) -> dict | None:
    rows = rows_for(table, gene_symbol)
    return rows[0] if rows else None


def provenance_of(row: dict | None) -> str | None:
    """The `source vversion` string that belongs in the audit trail.

    A result produced against one release of a reference table and one
    produced against a later release are not the same result, so the version
    travels with the finding rather than sitting in a config file.
    """
    if not row:
        return None
    source = (row.get("source") or "").strip()
    version = (row.get("source_version") or "").strip()
    if source and version:
        return f"{source} {version}"
    return source or None


def status() -> dict[str, dict]:
    """Which tables are populated. For the /scope endpoint and diagnostics."""
    out: dict[str, dict] = {}
    for name in TABLES:
        path = os.path.join(REFERENCE_DIR, f"{name}.tsv")
        rows = rows_for(name, None)  # warms the cache without matching
        with _lock:
            index = _cache.get(name, {})
        out[name] = {
            "present": os.path.exists(path),
            "genes": len(index),
            "populated": bool(index),
        }
        del rows
    return out


def reset_cache() -> None:
    """Drop cached tables. Call after replacing a file on disk."""
    with _lock:
        _cache.clear()
