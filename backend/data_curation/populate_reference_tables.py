"""Populate the reference tables by FETCHING from the live sources.

WHY THIS EXISTS RATHER THAN A HAND-WRITTEN TABLE
------------------------------------------------
`data_sources_halted_flagged.md` is explicit that no accession, version or
row may be written from recall: "Do not accept a URL, accession or version
number from any AI assistant, including me, without checking it. Standing
project rule — in Session 5, eight of nine recalled PMIDs were wrong."

So every row this script writes is retrieved from the source at run time and
stamped with where it came from and when. Nothing is typed in from memory.
Re-running refreshes the tables and updates `source_version`.

WHAT IT CAN AND CANNOT FILL
---------------------------
Two sources answer directly over HTTP and are populated here:

  protein_localisation   UniProt REST, per gene           -> B1, aptamer flag
  clingen_dosage         ClinGen dosage-sensitivity file  -> P6, replacement flag

The rest need a bulk download, a licence decision, or curation that is not a
lookup at all, and the script says so per table rather than filling them with
something plausible:

  tissue_expression            GTEx / HPA bulk matrix (SO-DATA-04 licence)
  repeat_expansion_loci        STRipy / ExpansionHunter catalogue
  polyadenylation_sites        PolyASite / PolyA_DB
  alt_promoters                CAGE / FANTOM atlas
  rbp_repressor_sites          literature curation, PMID per row
  *_benefit tables             per-gene therapeutic judgement, not a lookup
  curated_transcript_parts     seeded automatically by real_data_cache

GENE SET
--------
Defaults to the genes this platform actually reasons about — the benchmark's
targets plus the curated repeat-expansion loci — rather than the whole
genome. Pass --genes to extend it.

Run:
    python -m backend.data_curation.populate_reference_tables --dry-run
    python -m backend.data_curation.populate_reference_tables
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import io
import sys
import time
from pathlib import Path

import requests

REFERENCE_DIR = Path(__file__).resolve().parents[1] / "data" / "reference"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
CLINGEN_GENE_LIST = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"

# Genes the platform actually reasons about: every target in the
# mechanism-recovery benchmark, plus the repeat-expansion loci the TG05
# designer needs. Not a genome-wide dump.
DEFAULT_GENES = [
    # mechanism-recovery benchmark targets
    "TTR", "SOD1", "APOB", "DMD", "SMN2", "MFSD8",
    "ALAS1", "HAO1", "PCSK9",
    # repeat-expansion loci
    "DMPK", "HTT", "CNBP", "FMR1", "FXN", "C9orf72", "ATXN1", "ATXN2",
    # commonly discussed ASO / aptamer targets
    "SCN1A", "MAPT", "GRN", "STMN2", "FN1", "VEGFA", "CFTR", "HBB",
]

def _today() -> str:
    return _dt.date.today().isoformat()


# UniProt localisation strings -> the B1 classes feature_service expects.
#
# "Membrane" here must mean the PLASMA membrane. An aptamer cannot reach the
# ER, Golgi, mitochondrial, nuclear, endosomal or lysosomal membrane any more
# than it can reach the cytosol, so those count as intracellular.
_SURFACE = ("cell membrane", "cell surface", "apical", "basolateral",
            "plasma membrane", "postsynaptic cell membrane",
            "presynaptic cell membrane")
_OUTSIDE = ("secreted", "extracellular")
_INSIDE = ("cytoplasm", "nucleus", "nucleolus", "mitochondri", "cytosol",
           "endoplasmic", "golgi", "lysosome", "peroxisome", "endosome",
           "autophagosome", "autolysosome", "chromosome", "perikaryon",
           "p-body", "stress granule", "ribonucleoprotein", "cytoskeleton",
           "vesicle", "endomembrane", "synaptosome")


def _classify_localisation(locations: list[str], has_signal: bool,
                           n_transmembrane: int = 0) -> tuple[str, str]:
    """Classify by weight of evidence, not by first match.

    Returns (class, evidence).

    A single mention must not decide this. UniProt lists every compartment a
    protein has been observed in, so "any secreted annotation wins" calls
    C9orf72 secreted off one entry among eighteen, and "any membrane wins"
    calls FMR1 — a cytoplasmic RNA-binding protein — a surface target off two
    synaptic annotations among twenty-seven. B1 gates whether to suggest an
    aptamer, and suggesting one for an intracellular target is exactly the
    error it exists to prevent.

    So: a signal peptide plus any extracellular annotation is decisive, since
    that is what a signal peptide means. Transmembrane segments plus any
    surface annotation are likewise decisive — a multi-pass channel or
    receptor sits at the surface however many internal compartments it is
    also seen in during trafficking, and pure counting gets CFTR (a canonical
    apical chloride channel) wrong because its ER, endosome and recycling
    annotations outnumber the one that matters.

    Otherwise the dominant class wins, and ties fall to intracellular — the
    conservative direction, because the cost of withholding a flag is a
    missing prompt and the cost of raising a wrong one is a suggested therapy
    that cannot reach its target.
    """
    lowered = [loc.lower() for loc in locations]
    n_out = sum(1 for loc in lowered if any(k in loc for k in _OUTSIDE))
    n_surf = sum(1 for loc in lowered if any(k in loc for k in _SURFACE))
    n_in = sum(1 for loc in lowered
               if any(k in loc for k in _INSIDE)
               and not any(k in loc for k in _OUTSIDE + _SURFACE))
    evidence = (f"secreted={n_out} surface={n_surf} intracellular={n_in} "
                f"tm_segments={n_transmembrane}")

    if not (n_out or n_surf or n_in):
        return "", evidence + " (no usable annotation)"
    # A signal peptide is the strongest single indicator of secretion.
    if has_signal and n_out and not n_transmembrane:
        return "secreted", evidence + " +signal_peptide"
    # Transmembrane segments plus a surface annotation: an integral membrane
    # protein at the cell surface, whatever else it is seen in.
    if n_transmembrane and n_surf:
        return "membrane", evidence + " +transmembrane"
    if has_signal and n_out:
        return "secreted", evidence + " +signal_peptide"
    if n_out > n_in and n_out >= n_surf:
        return "secreted", evidence
    if n_surf > n_in:
        return "membrane", evidence
    return "intracellular", evidence


def fetch_uniprot_localisation(genes: list[str], pause: float = 0.4) -> list[dict]:
    """Real subcellular localisation, one reviewed human entry per gene."""
    rows = []
    for gene in genes:
        try:
            resp = requests.get(
                UNIPROT_SEARCH,
                params={
                    "query": (f"gene_exact:{gene} AND organism_id:9606 "
                              f"AND reviewed:true"),
                    "fields": ("accession,gene_primary,cc_subcellular_location,"
                               "ft_signal,ft_transmem"),
                    "format": "json",
                    # More than one, because the top hit is not always the
                    # gene asked for — see the symbol-collision guard below.
                    "size": 10,
                },
                timeout=30,
            )
            resp.raise_for_status()
            results = resp.json().get("results", [])
        except Exception as exc:  # noqa: BLE001
            print(f"  {gene:10} FETCH FAILED: {type(exc).__name__}")
            continue
        if not results:
            print(f"  {gene:10} no reviewed human entry")
            continue

        # UniProt ranks by relevance, not by symbol identity, and legacy
        # symbols collide: `gene_exact:HTT` returns SLC6A4 (the serotonin
        # transporter, historically "HTT") ABOVE huntingtin. Taking
        # results[0] silently writes a 12-transmembrane surface protein into
        # the row for a cytoplasmic scaffold, which would then fire the
        # aptamer flag for a target no aptamer can reach.
        #
        # So require the entry's own primary gene name to be the symbol we
        # asked for, and skip the gene entirely rather than guess if none
        # matches.
        entry = None
        for cand in results:
            primary = ((cand.get("genes") or [{}])[0]
                       .get("geneName", {}).get("value", ""))
            if primary.upper() == gene.upper():
                entry = cand
                break
        if entry is None:
            got = ", ".join(
                ((c.get("genes") or [{}])[0].get("geneName", {}).get("value", "?"))
                for c in results[:3]
            )
            print(f"  {gene:10} SKIPPED — no reviewed entry whose primary gene "
                  f"symbol is {gene} (top hits: {got})")
            continue

        locations = [
            loc["location"]["value"]
            for c in entry.get("comments", [])
            if c.get("commentType") == "SUBCELLULAR LOCATION"
            for loc in c.get("subcellularLocations", [])
            if loc.get("location")
        ]
        feats = entry.get("features", [])
        has_signal = any(f.get("type") == "Signal" for f in feats)
        n_tm = sum(1 for f in feats if f.get("type") == "Transmembrane")
        cls, evidence = _classify_localisation(locations, has_signal, n_tm)
        if not cls:
            print(f"  {gene:10} localisation not classifiable: {locations[:2]}")
            continue
        rows.append({
            "gene_symbol": gene,
            "uniprot_id": entry.get("primaryAccession", ""),
            "has_signal_peptide": "true" if has_signal else "false",
            "localisation_class": cls,
            "evidence": evidence,
            "source": "UniProt REST (reviewed, organism 9606)",
            "source_version": f"retrieved {_today()}",
        })
        print(f"  {gene:10} {entry.get('primaryAccession',''):8} {cls:14} "
              f"signal={'y' if has_signal else 'n'}  {evidence}")
        time.sleep(pause)
    return rows


def fetch_clingen_dosage(genes: list[str] | None = None) -> list[dict]:
    """ClinGen haploinsufficiency / triplosensitivity scores.

    The score is kept as the STRING the file carries. It is a code, not an
    ordinal: 30 means "autosomal recessive" and 40 means "dosage sensitivity
    unlikely", so comparing the column numerically would read 40 as the
    strongest evidence available (SO-DATA-02).
    """
    try:
        resp = requests.get(CLINGEN_GENE_LIST, timeout=60)
        resp.raise_for_status()
    except Exception as exc:  # noqa: BLE001
        print(f"  ClinGen fetch FAILED: {type(exc).__name__}: {exc}")
        return []

    text = resp.text
    version = f"retrieved {_today()}"
    for line in text.splitlines()[:20]:
        if line.lower().startswith("#") and "date" in line.lower():
            version = line.lstrip("# ").strip()
            break

    # The file carries '#'-prefixed preamble lines; the last one is the header.
    lines = text.splitlines()
    header_idx = max(i for i, ln in enumerate(lines[:20]) if ln.startswith("#"))
    reader = csv.DictReader(
        io.StringIO("\n".join([lines[header_idx].lstrip("#")] + lines[header_idx + 1:])),
        delimiter="\t",
    )
    wanted = {g.upper() for g in genes} if genes else None
    rows = []
    for rec in reader:
        sym = (rec.get("Gene Symbol") or rec.get("gene_symbol") or "").strip()
        if not sym or (wanted and sym.upper() not in wanted):
            continue
        hi = (rec.get("Haploinsufficiency Score") or "").strip()
        ts = (rec.get("Triplosensitivity Score") or "").strip()
        if not hi and not ts:
            continue
        rows.append({
            "gene_symbol": sym,
            "haploinsufficiency_score": hi,
            "triplosensitivity_score": ts,
            "source": "ClinGen Dosage Sensitivity (GRCh38 gene curation list)",
            "source_version": version,
        })
    for r in rows:
        print(f"  {r['gene_symbol']:10} HI={r['haploinsufficiency_score'] or '-':4} "
              f"TS={r['triplosensitivity_score'] or '-'}")
    return rows


def write_table(name: str, rows: list[dict]) -> int:
    """Rewrite a table, preserving its existing header order."""
    path = REFERENCE_DIR / f"{name}.tsv"
    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in header})
    return len(rows)


# Tables that cannot be filled by an HTTP lookup, and precisely why.
BLOCKED = {
    "tissue_expression":
        "GTEx / HPA bulk expression matrix — a multi-GB download, and "
        "redistribution of a derived table is SO-DATA-04.",
    "repeat_expansion_loci":
        "STRipy / ExpansionHunter catalogue — needs the release pinned and "
        "its licence checked before the rows are trusted.",
    "polyadenylation_sites":
        "PolyASite / PolyA_DB bulk download; which resource is current is "
        "MUST VERIFY.",
    "alt_promoters":
        "CAGE / FANTOM promoter atlas bulk download.",
    "rbp_repressor_sites":
        "Literature curation with a PMID per row — recommended over an "
        "atlas, and not a lookup.",
    "apa_therapeutic_benefit":
        "Per-gene therapeutic judgement, not a lookup.",
    "alt_promoter_benefit":
        "Per-gene therapeutic judgement, not a lookup.",
    "intron_retention_potential":
        "Splice-site conservation plus NMD context; needs a pinned "
        "annotation release.",
    "intron_retention_benefit":
        "Per-gene therapeutic judgement, not a lookup.",
    "dominant_negative_genes":
        "Manual curation with a PMID per row. ClinGen supplies the "
        "permissive half of P6; this is the suppressing half.",
    "curated_transcript_parts":
        "Seeded automatically by real_data_cache from successful live "
        "fetches; no bulk source needed.",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="*", default=None,
                    help="Gene symbols to fetch. Defaults to the platform set.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Fetch and report without writing.")
    args = ap.parse_args()
    genes = [g.upper() for g in (args.genes or DEFAULT_GENES)]

    print("=" * 72)
    print("POPULATING REFERENCE TABLES FROM LIVE SOURCES")
    print(f"{len(genes)} genes | {'DRY RUN' if args.dry_run else 'WRITING'}")
    print("=" * 72)

    print("\nprotein_localisation  <- UniProt REST")
    loc_rows = fetch_uniprot_localisation(genes)

    print("\nclingen_dosage  <- ClinGen Dosage Sensitivity")
    clingen_rows = fetch_clingen_dosage(genes)

    if not args.dry_run:
        n1 = write_table("protein_localisation", loc_rows)
        n2 = write_table("clingen_dosage", clingen_rows)
        print(f"\nwrote protein_localisation.tsv  {n1} rows")
        print(f"wrote clingen_dosage.tsv        {n2} rows")
    else:
        print(f"\n[dry run] would write {len(loc_rows)} + {len(clingen_rows)} rows")

    print("\n" + "=" * 72)
    print("STILL EMPTY, AND WHY")
    print("=" * 72)
    for name, why in BLOCKED.items():
        print(f"  {name:28} {why}")
    print("\nEvery row written above carries its source and retrieval date. "
          "Nothing here was typed from memory.")


if __name__ == "__main__":
    sys.exit(main())
