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
# Branch points cluster 18-40 nt upstream of the 3' splice site; 28 centres
# the tiling window on that range.
BRANCH_POINT_OFFSET_NT = 28
TSS_UPSTREAM_NT = 100
TSS_DOWNSTREAM_NT = 50
# Distinct TSSs closer together than this are treated as one promoter.
TSS_CLUSTER_NT = 200

# Which transcript biotypes evidence an independent transcription initiation
# event, and therefore a promoter.
#
# `retained_intron` and `protein_coding_CDS_not_defined` do NOT. They annotate
# incompletely processed or fragmentary RNA — splicing intermediates whose
# annotated 5' end is a consequence of where sequencing coverage began, not of
# where polymerase loaded. Counting them produced 13 "alternative promoters"
# for HTT, of which 10 came from retained-intron entries; the gene does not
# have 13 promoters. `nonsense_mediated_decay` transcripts ARE kept: they are
# transcribed from a real promoter and then degraded, so their 5' end is a
# genuine initiation site.
PROMOTER_EVIDENCE_BIOTYPES = {
    "protein_coding",
    "nonsense_mediated_decay",
    "lncRNA",
    "processed_transcript",
}

_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}

# Spliceosomal intron termini. The overwhelming majority of human introns are
# U2-type and begin GU / end AG; a small minority are U12-type (AU..AC). An
# intron that matches neither is a signal that the coordinates or the strand
# were resolved wrongly, which is the single easiest thing to get wrong when
# moving from transcript space to genomic space.
INTRON_CLASSES = {
    ("GU", "AG"): "U2-type (canonical GU-AG)",
    ("AU", "AC"): "U12-type (minor spliceosome, AU-AC)",
    ("GC", "AG"): "U2-type variant (GC-AG)",
}


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

    # Read the intron in transcript orientation and check its termini. A
    # genuine spliceosomal intron starts GU and ends AG; if this one does not,
    # the coordinates or the strand are wrong and every oligo tiled from them
    # would be aimed at the wrong place. Report it rather than designing on.
    intron_seq = _fetch_region(chromosome, intron_start, intron_end, strand,
                               organism)
    intron_class, splice_sites = None, {}
    if len(intron_seq) >= 4:
        rna = intron_seq.replace("T", "U")
        termini = (rna[:2], rna[-2:])
        intron_class = INTRON_CLASSES.get(termini)
        splice_sites = {
            "donorDinucleotide": termini[0],
            "acceptorDinucleotide": termini[1],
            "class": intron_class,
        }
        if intron_class is None:
            return {
                "status": "UNAVAILABLE", "mechanismId": "A33",
                "candidates": [], "spliceSites": splice_sites,
                "message": (
                    f"Intron {intron_number} reads {termini[0]}...{termini[1]}"
                    f" in transcript orientation. A spliceosomal intron begins"
                    f" GU and ends AG (or AU-AC for the minor spliceosome), so"
                    f" these coordinates do not describe an intron on this"
                    f" strand. Nothing is designed against them."
                ),
            }

    # The branch point sits ~18-40 nt upstream of the 3' splice site, with the
    # polypyrimidine tract between it and the acceptor. Blocking that region
    # is the third way to force retention, and on a long intron it is often
    # more accessible than the acceptor itself.
    bp_offset = -BRANCH_POINT_OFFSET_NT if strand >= 0 else BRANCH_POINT_OFFSET_NT
    wanted = []
    if splice_element in ("both", "splice_donor"):
        wanted.append(("5' splice site (donor)", donor_coord))
    if splice_element in ("both", "splice_acceptor"):
        wanted.append(("3' splice site (acceptor)", acceptor_coord))
    if splice_element in ("both", "branch_point"):
        wanted.append(("branch point / polypyrimidine tract",
                       acceptor_coord + bp_offset))
    if not wanted:
        return {"status": "UNAVAILABLE", "mechanismId": "A33", "candidates": [],
                "message": f"Unknown splice_element {splice_element!r}."}

    # The acceptor and branch-point windows overlap on a short intron, so the
    # same oligo can be produced twice under two labels. Keep the first and
    # record the other element it also covers.
    candidates: list[dict[str, Any]] = []
    by_sequence: dict[str, dict[str, Any]] = {}
    for label, coord in wanted:
        lo, hi = coord - JUNCTION_FLANK_NT, coord + JUNCTION_FLANK_NT
        seq = _fetch_region(chromosome, lo, hi, strand, organism)
        if not seq:
            continue
        for cand in _tile(seq, lo, oligo_length, 4, label, "A33",
                          gene_symbol or ensembl_gene_id):
            existing = by_sequence.get(cand["sequence"])
            if existing is not None:
                also = existing.setdefault("alsoCovers", [])
                if label not in also and label != existing["targetElement"]:
                    also.append(label)
                continue
            by_sequence[cand["sequence"]] = cand
            candidates.append(cand)

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
        "spliceSites": splice_sites,
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
    usable, excluded = [], {}
    for t in transcripts:
        if not (t.get("start") and t.get("end")):
            continue
        biotype = t.get("biotype") or "unknown"
        if biotype not in PROMOTER_EVIDENCE_BIOTYPES:
            excluded[biotype] = excluded.get(biotype, 0) + 1
            continue
        usable.append(int(t["end"]) if strand < 0 else int(t["start"]))
    if not usable:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "excludedBiotypes": excluded,
                "message": ("No transcript of a biotype that evidences an "
                            "independent transcription start was annotated. "
                            f"Excluded: {excluded}.")}

    tss_list = sorted(usable, reverse=strand < 0)
    clusters: list[list[int]] = []
    for tss in tss_list:
        if clusters and abs(tss - clusters[-1][-1]) <= TSS_CLUSTER_NT:
            clusters[-1].append(tss)
        else:
            clusters.append([tss])
    promoters = [int(sum(c) / len(c)) for c in clusters]
    # How many transcripts back each promoter — a promoter resting on one
    # transcript is weaker evidence than one resting on five.
    promoter_support = [len(c) for c in clusters]

    if not (1 <= promoter_index <= len(promoters)):
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "message": (f"Promoter {promoter_index} does not exist; "
                            f"{len(promoters)} distinct TSS cluster(s) were "
                            f"found ({promoters}).")}
    if len(promoters) < 2:
        return {"status": "UNAVAILABLE", "mechanismId": "A32", "candidates": [],
                "promotersFound": promoters,
                "promoterTranscriptSupport": promoter_support,
                "excludedBiotypes": excluded,
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
        "promoterTranscriptSupport": promoter_support,
        "excludedBiotypes": excluded,
        "promoterEvidenceNote": (
            "Promoters are TSS clusters from transcripts whose biotype "
            "evidences independent initiation. retained_intron and "
            "CDS-not-defined entries are excluded: they annotate incompletely "
            "processed RNA, and their 5' ends are not initiation events."
        ),
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
