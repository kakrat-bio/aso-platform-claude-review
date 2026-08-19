"""A32 / A33 — designs that act on pre-mRNA, not on the mature transcript.

Both mechanisms were marked design-unavailable because the isoform designer
tiles the SPLICED mRNA, which contains neither a promoter (A32) nor an
intron (A33). Neither is undesignable in principle — they just need genomic
sequence, which the mature transcript cannot supply. This module fetches it.

* **A33, targeted intron retention.** An intron is retained by blocking its
  splice sites. The intron's coordinates are the gap between consecutive
  exons of the canonical transcript, and Ensembl serves that genomic region
  directly. Oligos are tiled across the 5' splice site (donor), the branch
  point / 3' splice site (acceptor), or both, because those are the elements
  a steric blocker has to occupy to force retention.

* **A32, alternative promoter switching.** Alternative promoters are the
  distinct transcription start sites of a gene's transcripts. Ensembl gives
  a start coordinate per transcript, so the distinct TSSs are recoverable by
  clustering them. Oligos are tiled across the TSS-proximal region of the
  promoter to be suppressed, which is where a promoter-directed
  oligonucleotide acts on the nascent transcript.

Strand is handled explicitly: on the minus strand a transcript's TSS is its
genomic END and intron donor/acceptor sides swap, which is the easiest thing
to get silently wrong here.
"""

from __future__ import annotations

import logging
from typing import Any

import primer3
import RNA

from services.gene_silencing_service import ENSEMBL_REST, _ensembl_get, get_target_analysis

logger = logging.getLogger(__name__)

# How far either side of a splice junction or TSS to tile.
JUNCTION_FLANK_NT = 40
TSS_UPSTREAM_NT = 100
TSS_DOWNSTREAM_NT = 50
# Distinct TSSs closer together than this are treated as one promoter.
TSS_CLUSTER_NT = 200

_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}


def _revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def _fetch_region(chromosome: str, start: int, end: int, strand: int,
                  organism: str = "homo_sapiens") -> str:
    """Genomic sequence for a region, already oriented to the given strand."""
    if end < start:
        start, end = end, start
    url = (f"{ENSEMBL_REST}/sequence/region/{organism}/"
           f"{chromosome}:{start}..{end}:{1 if strand >= 0 else -1}"
           f"?content-type=application/json")
    resp = _ensembl_get(url, timeout=20)
    if resp is None or not getattr(resp, "ok", False):
        return ""
    payload = resp.json()
    seq = payload.get("seq", "") if isinstance(payload, dict) else ""
    return (seq or "").upper()


def _metrics(oligo: str, target: str) -> dict[str, Any]:
    dna = oligo.upper().replace("U", "T")
    tm = None
    if dna and not (set(dna) - set("ACGT")):
        tm = round(primer3.calc_tm(dna), 1)
    try:
        dg = round(RNA.duplexfold(oligo, target).energy, 2)
    except Exception:
        dg = None
    gc = round((oligo.count("G") + oligo.count("C")) / max(len(oligo), 1) * 100, 1)
    return {"meltingTempC": tm, "targetDuplexDg": dg, "gcContent": gc}


def _tile(region_seq: str, region_start: int, oligo_len: int, step: int,
          label: str, mechanism_id: str, gene: str) -> list[dict[str, Any]]:
    out = []
    for i in range(0, max(0, len(region_seq) - oligo_len + 1), step):
        window = region_seq[i:i + oligo_len].replace("T", "U")
        if set(window) - set("ACGU"):
            continue
        oligo = _revcomp(window)
        out.append({
            "constructId": f"{mechanism_id}-{gene}-{region_start + i}",
            "mechanismId": mechanism_id,
            "sequence": oligo,
            "targetSequence": window,
            "targetElement": label,
            "genomicStart": region_start + i,
            "genomicEnd": region_start + i + oligo_len,
            "length": oligo_len,
            **_metrics(oligo, window),
        })
    return out


def design_intron_retention(
    ensembl_gene_id: str,
    intron_number: int,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    oligo_length: int = 20,
    splice_element: str = "both",
    max_candidates: int = 12,
) -> dict[str, Any]:
    """A33 — block an intron's splice sites so the intron is retained."""
    target = get_target_analysis(ensembl_gene_id, gene_symbol=gene_symbol,
                                 organism=organism)
    exons = target.get("exons") or []
    canonical = target.get("canonicalTranscript") or {}
    chromosome = canonical.get("chromosome")
    strand = int(canonical.get("strand") or 1)
    if len(exons) < 2 or not chromosome:
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": (f"{gene_symbol or ensembl_gene_id} has no usable "
                            f"multi-exon canonical transcript, so it has no "
                            f"intron to retain.")}

    ordered = sorted(exons, key=lambda e: e.get("start", 0))
    n_introns = len(ordered) - 1
    if not (1 <= intron_number <= n_introns):
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": (f"Intron {intron_number} does not exist; the "
                            f"canonical transcript has {n_introns}.")}

    left, right = ordered[intron_number - 1], ordered[intron_number]
    intron_start = int(left["end"]) + 1
    intron_end = int(right["start"]) - 1
    if intron_end <= intron_start:
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": "Adjacent exons leave no intronic gap."}

    # On the minus strand the transcript reads right-to-left, so the intron's
    # donor side is the genomically HIGHER coordinate.
    donor_coord = intron_start if strand >= 0 else intron_end
    acceptor_coord = intron_end if strand >= 0 else intron_start

    wanted = []
    if splice_element in ("both", "splice_donor"):
        wanted.append(("5' splice site (donor)", donor_coord))
    if splice_element in ("both", "splice_acceptor"):
        wanted.append(("3' splice site (acceptor)", acceptor_coord))
    if not wanted:
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": f"Unknown splice_element {splice_element!r}."}

    candidates: list[dict[str, Any]] = []
    for label, coord in wanted:
        lo, hi = coord - JUNCTION_FLANK_NT, coord + JUNCTION_FLANK_NT
        seq = _fetch_region(chromosome, lo, hi, strand, organism)
        if not seq:
            continue
        candidates.extend(_tile(seq, lo, oligo_length, 4, label, "A33",
                                gene_symbol or ensembl_gene_id))

    if not candidates:
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": ("Ensembl returned no genomic sequence for the "
                            "intron's splice sites.")}

    candidates.sort(key=lambda c: (c["targetDuplexDg"] is None,
                                   c["targetDuplexDg"] or 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1
    return {
        "status": "OK",
        "mechanismId": "A33",
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "intronNumber": intron_number,
        "intronLength": intron_end - intron_start + 1,
        "strand": strand,
        "chromosome": chromosome,
        "architecture": (f"{oligo_length} nt steric blocker tiled across the "
                         f"intron's splice site(s)"),
        "ranking": {"orderedBy": "targetDuplexDg",
                    "caveat": ("Thermodynamic ordering, not a validated "
                               "retention model.")},
        "dataProvenance": {"sequence": "Ensembl genomic region (introns are "
                                       "absent from the spliced transcript)"},
        "candidates": candidates[:max_candidates],
    }


def design_alternative_promoter(
    ensembl_gene_id: str,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    oligo_length: int = 20,
    promoter_index: int = 1,
    max_candidates: int = 12,
) -> dict[str, Any]:
    """A32 — tile the TSS-proximal region of one alternative promoter."""
    resp = _ensembl_get(f"{ENSEMBL_REST}/lookup/id/{ensembl_gene_id}?expand=1",
                        timeout=20)
    if resp is None or not getattr(resp, "ok", False):
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": f"Ensembl lookup failed for {ensembl_gene_id}."}
    data = resp.json()
    chromosome = str(data.get("seq_region_name") or "")
    strand = int(data.get("strand") or 1)
    transcripts = data.get("Transcript") or []
    if not transcripts or not chromosome:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": "No transcripts annotated for this gene."}

    # A transcript's TSS is its genomic start on the plus strand and its
    # genomic end on the minus strand.
    tss_list = sorted({
        int(t["end"]) if strand < 0 else int(t["start"])
        for t in transcripts if t.get("start") and t.get("end")
    }, reverse=strand < 0)

    clusters: list[list[int]] = []
    for tss in tss_list:
        if clusters and abs(tss - clusters[-1][-1]) <= TSS_CLUSTER_NT:
            clusters[-1].append(tss)
        else:
            clusters.append([tss])
    promoters = [int(sum(c) / len(c)) for c in clusters]

    if not (1 <= promoter_index <= len(promoters)):
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": (f"Promoter {promoter_index} does not exist; "
                            f"{len(promoters)} distinct TSS cluster(s) were "
                            f"found ({promoters}).")}
    if len(promoters) < 2:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "promotersFound": promoters,
                "message": ("Only one transcription start site cluster is "
                            "annotated, so there is no alternative promoter "
                            "to switch away from.")}

    tss = promoters[promoter_index - 1]
    if strand >= 0:
        lo, hi = tss - TSS_UPSTREAM_NT, tss + TSS_DOWNSTREAM_NT
    else:
        lo, hi = tss - TSS_DOWNSTREAM_NT, tss + TSS_UPSTREAM_NT
    seq = _fetch_region(chromosome, lo, hi, strand, organism)
    if not seq:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": "Ensembl returned no genomic sequence for the promoter."}

    candidates = _tile(seq, lo, oligo_length, 4,
                       f"TSS-proximal region of promoter {promoter_index}",
                       "A32", gene_symbol or ensembl_gene_id)
    if not candidates:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": "No unambiguous window in the promoter region."}

    candidates.sort(key=lambda c: (c["targetDuplexDg"] is None,
                                   c["targetDuplexDg"] or 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1
    return {
        "status": "OK",
        "mechanismId": "A32",
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "chromosome": chromosome,
        "strand": strand,
        "promotersFound": promoters,
        "promoterIndex": promoter_index,
        "transcriptionStartSite": tss,
        "architecture": (f"{oligo_length} nt oligo tiled from "
                         f"-{TSS_UPSTREAM_NT} to +{TSS_DOWNSTREAM_NT} around "
                         f"the transcription start site"),
        "ranking": {"orderedBy": "targetDuplexDg",
                    "caveat": ("Thermodynamic ordering. Promoter-directed "
                               "oligonucleotides act on chromatin and nascent "
                               "transcript, and no model here predicts that "
                               "activity.")},
        "dataProvenance": {"sequence": "Ensembl genomic region around the "
                                       "annotated TSS"},
        "candidates": candidates[:max_candidates],
    }
