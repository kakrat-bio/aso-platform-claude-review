"""A3 / A4 — the two TG02 mechanisms whose target is not the canonical mRNA.

Both were tiling the canonical transcript and labelling the result, which
cannot work for either of them:

* **A3 (TANGO, poison-exon blocking).** A poison exon is one that is SKIPPED
  in the productive transcript and INCLUDED in an NMD-destined one. It is
  therefore absent from the canonical mRNA by definition, and no amount of
  tiling that mRNA will reach it. The old code set the label "Exon junctions"
  and scanned the whole transcript, returning 891 candidates spread over all
  29 SCN1A exon junctions — none of them aimed at the poison exon the
  mechanism is named for. What is needed is the poison exon's own genomic
  sequence and the splice sites that recruit it.

* **A4 (AntagoNAT).** The target is the natural antisense transcript, which
  is a DIFFERENT GENE on the opposite strand. The old code tiled the whole
  sense transcript under the label "Full transcript (NAT complement)". The
  strand was right — an oligo antisense to the NAT reads as sense to the
  gene — but SCN1A-AS1 overlaps only part of SCN1A, so every candidate
  outside that overlap was complementary to nothing.

Both targets are recoverable from Ensembl, which is what this module does.
"""

from __future__ import annotations

import logging
from typing import Any

from services.gene_silencing_service import ENSEMBL_REST, _ensembl_get
from services.transcript_architecture_service import _fetch_region, _tile

logger = logging.getLogger(__name__)

# A poison exon is a short cassette exon. Longer "extra" exons in an NMD
# transcript are usually extended termini or retained introns, not the
# PTC-introducing cassette this mechanism targets.
MAX_POISON_EXON_NT = 300
# How far either side of the poison exon's splice sites to tile.
SPLICE_FLANK_NT = 40
MAX_CANDIDATES = 12


def find_poison_exons(ensembl_gene_id: str,
                      organism: str = "homo_sapiens") -> dict[str, Any]:
    """Exons present in an NMD transcript and absent from the canonical one.

    Ensembl annotates NMD-destined transcripts with the
    `nonsense_mediated_decay` biotype. Comparing their exon coordinates
    against the canonical transcript's yields the cassette exons whose
    inclusion sends the message to decay — the poison exons.
    """
    resp = _ensembl_get(f"{ENSEMBL_REST}/lookup/id/{ensembl_gene_id}?expand=1",
                        timeout=25)
    if resp is None or not getattr(resp, "ok", False):
        return {"status": "UNAVAILABLE", "poisonExons": [],
                "message": f"Ensembl lookup failed for {ensembl_gene_id}."}
    data = resp.json()
    transcripts = data.get("Transcript") or []
    canonical = next((t for t in transcripts if t.get("is_canonical")), None)
    if canonical is None:
        canonical = next((t for t in transcripts
                          if t.get("biotype") == "protein_coding"), None)
    if canonical is None:
        return {"status": "UNAVAILABLE", "poisonExons": [],
                "message": "No canonical protein-coding transcript annotated."}

    canonical_exons = {(int(e["start"]), int(e["end"]))
                       for e in canonical.get("Exon") or []}
    nmd = [t for t in transcripts
           if t.get("biotype") == "nonsense_mediated_decay"]
    if not nmd:
        return {
            "status": "UNAVAILABLE", "poisonExons": [],
            "message": ("No nonsense_mediated_decay transcript is annotated "
                        "for this gene, so no poison exon can be located. "
                        "The mechanism may still apply — this is the limit of "
                        "the annotation, not a negative finding."),
        }

    found: dict[tuple[int, int], dict[str, Any]] = {}
    for t in nmd:
        exons = sorted(t.get("Exon") or [], key=lambda e: int(e["start"]))
        for i, e in enumerate(exons):
            coords = (int(e["start"]), int(e["end"]))
            if coords in canonical_exons:
                continue
            # Terminal exons of the NMD transcript are extensions, not
            # cassettes; a poison exon is internal.
            if i == 0 or i == len(exons) - 1:
                continue
            length = coords[1] - coords[0] + 1
            if length > MAX_POISON_EXON_NT:
                continue
            entry = found.setdefault(coords, {
                "start": coords[0], "end": coords[1], "length": length,
                "supportingTranscripts": [],
            })
            entry["supportingTranscripts"].append(t.get("id"))

    # Most-supported first, then shortest. Length is the tie-break because a
    # poison exon is a short cassette: SCN1A's is 64 nt, and the clinically
    # targeted one (STK-001, exon 20N) is the shortest of the six this finds.
    poison = sorted(found.values(),
                    key=lambda p: (-len(p["supportingTranscripts"]), p["length"]))
    if not poison:
        return {
            "status": "UNAVAILABLE", "poisonExons": [],
            "message": (f"{len(nmd)} NMD transcript(s) are annotated, but none "
                        f"contributes a short internal exon absent from the "
                        f"canonical transcript. Their differences are terminal "
                        f"extensions or exons longer than "
                        f"{MAX_POISON_EXON_NT} nt."),
        }
    return {
        "status": "OK",
        "chromosome": str(data.get("seq_region_name") or ""),
        "strand": int(data.get("strand") or 1),
        "canonicalTranscript": canonical.get("id"),
        "nmdTranscriptCount": len(nmd),
        "poisonExons": poison,
    }


def design_poison_exon_block(
    ensembl_gene_id: str,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    poison_exon_index: int = 1,
    oligo_length: int = 20,
    splice_element: str = "both",
    max_candidates: int = MAX_CANDIDATES,
) -> dict[str, Any]:
    """A3 — block the splice sites that include a poison exon."""
    located = find_poison_exons(ensembl_gene_id, organism)
    if located["status"] != "OK":
        return {"status": "UNAVAILABLE", "mechanismId": "A3",
                "candidates": [], "message": located["message"]}

    exons = located["poisonExons"]
    if not (1 <= poison_exon_index <= len(exons)):
        return {"status": "UNAVAILABLE", "mechanismId": "A3", "candidates": [],
                "poisonExons": exons,
                "message": (f"Poison exon {poison_exon_index} does not exist; "
                            f"{len(exons)} were located.")}
    exon = exons[poison_exon_index - 1]
    chrom, strand = located["chromosome"], located["strand"]

    # On the minus strand the transcript reads right-to-left, so the exon's
    # acceptor side is the genomically higher coordinate.
    acceptor = exon["start"] if strand >= 0 else exon["end"]
    donor = exon["end"] if strand >= 0 else exon["start"]
    wanted = []
    if splice_element in ("both", "splice_acceptor"):
        wanted.append(("poison exon 3' splice site (acceptor)", acceptor))
    if splice_element in ("both", "splice_donor"):
        wanted.append(("poison exon 5' splice site (donor)", donor))

    candidates, seen = [], set()
    for label, coord in wanted:
        lo, hi = coord - SPLICE_FLANK_NT, coord + SPLICE_FLANK_NT
        seq = _fetch_region(chrom, lo, hi, strand, organism)
        if not seq:
            continue
        for cand in _tile(seq, lo, oligo_length, 4, label, "A3",
                          gene_symbol or ensembl_gene_id):
            if cand["sequence"] in seen:
                continue
            seen.add(cand["sequence"])
            candidates.append(cand)

    if not candidates:
        return {"status": "UNAVAILABLE", "mechanismId": "A3", "candidates": [],
                "message": ("Ensembl returned no genomic sequence around the "
                            "poison exon's splice sites.")}

    candidates.sort(key=lambda c: (c["targetDuplexDg"] is None,
                                   c["targetDuplexDg"] or 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1
    return {
        "status": "OK",
        "mechanismId": "A3",
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "poisonExon": exon,
        "poisonExonsFound": exons,
        "poisonExonIndex": poison_exon_index,
        "nmdTranscriptCount": located["nmdTranscriptCount"],
        "strand": strand,
        "architecture": (f"{oligo_length} nt steric blocker across the poison "
                         f"exon's splice site(s)"),
        "ranking": {"orderedBy": "targetDuplexDg",
                    "caveat": ("Thermodynamic ordering, not a validated "
                               "skipping model.")},
        "dataProvenance": {
            "poisonExon": ("Ensembl nonsense_mediated_decay transcripts, "
                           "exons absent from the canonical transcript"),
            "sequence": ("Ensembl genomic region — a poison exon is absent "
                         "from the canonical mRNA by definition"),
        },
        "candidates": candidates[:max_candidates],
    }


def find_nat(ensembl_gene_id: str,
             organism: str = "homo_sapiens") -> dict[str, Any]:
    """Genes on the opposite strand that overlap this one."""
    resp = _ensembl_get(f"{ENSEMBL_REST}/lookup/id/{ensembl_gene_id}", timeout=25)
    if resp is None or not getattr(resp, "ok", False):
        return {"status": "UNAVAILABLE", "nats": [],
                "message": f"Ensembl lookup failed for {ensembl_gene_id}."}
    gene = resp.json()
    chrom = str(gene.get("seq_region_name") or "")
    strand = int(gene.get("strand") or 1)
    start, end = int(gene["start"]), int(gene["end"])

    ov = _ensembl_get(
        f"{ENSEMBL_REST}/overlap/region/{organism}/{chrom}:{start}-{end}"
        f"?feature=gene;content-type=application/json", timeout=30)
    if ov is None or not getattr(ov, "ok", False):
        return {"status": "UNAVAILABLE", "nats": [],
                "message": "Ensembl overlap query failed."}
    nats = []
    for g in ov.json():
        if g.get("id") == gene.get("id") or int(g.get("strand", 0)) == strand:
            continue
        ov_start, ov_end = max(start, int(g["start"])), min(end, int(g["end"]))
        if ov_end <= ov_start:
            continue
        nats.append({
            "geneId": g.get("id"),
            "symbol": g.get("external_name"),
            "biotype": g.get("biotype"),
            "start": int(g["start"]), "end": int(g["end"]),
            "overlapStart": ov_start, "overlapEnd": ov_end,
            "overlapNt": ov_end - ov_start + 1,
        })
    if not nats:
        return {"status": "UNAVAILABLE", "nats": [],
                "message": ("No gene on the opposite strand overlaps this "
                            "one, so there is no natural antisense transcript "
                            "to knock down.")}
    nats.sort(key=lambda n: -n["overlapNt"])
    return {"status": "OK", "chromosome": chrom, "senseStrand": strand,
            "nats": nats}


def design_nat_knockdown(
    ensembl_gene_id: str,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    nat_index: int = 1,
    oligo_length: int = 20,
    max_candidates: int = MAX_CANDIDATES,
) -> dict[str, Any]:
    """A4 — gapmers against the NAT, tiled over the real overlap only."""
    located = find_nat(ensembl_gene_id, organism)
    if located["status"] != "OK":
        return {"status": "UNAVAILABLE", "mechanismId": "A4",
                "candidates": [], "message": located["message"]}
    nats = located["nats"]
    if not (1 <= nat_index <= len(nats)):
        return {"status": "UNAVAILABLE", "mechanismId": "A4", "candidates": [],
                "natsFound": nats,
                "message": (f"NAT {nat_index} does not exist; "
                            f"{len(nats)} overlapping antisense gene(s) found.")}
    nat = nats[nat_index - 1]

    # Read the overlap in the NAT's orientation. The oligo is the reverse
    # complement of the NAT, which — because the NAT is antisense to the gene
    # — reads as the SENSE strand of the gene. That is correct and is worth
    # stating, because it looks wrong at a glance.
    nat_strand = -located["senseStrand"]
    seq = _fetch_region(located["chromosome"], nat["overlapStart"],
                        nat["overlapEnd"], nat_strand, organism)
    if not seq:
        return {"status": "UNAVAILABLE", "mechanismId": "A4", "candidates": [],
                "message": "Ensembl returned no sequence for the NAT overlap."}

    step = max(1, (len(seq) - oligo_length) // max(max_candidates * 3, 1))
    candidates, seen = [], set()
    for cand in _tile(seq, nat["overlapStart"], oligo_length, step,
                      f"NAT {nat.get('symbol') or nat['geneId']} overlap",
                      "A4", gene_symbol or ensembl_gene_id):
        if cand["sequence"] in seen:
            continue
        seen.add(cand["sequence"])
        candidates.append(cand)
    if not candidates:
        return {"status": "UNAVAILABLE", "mechanismId": "A4", "candidates": [],
                "message": "No unambiguous window in the NAT overlap region."}

    candidates.sort(key=lambda c: (c["targetDuplexDg"] is None,
                                   c["targetDuplexDg"] or 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1
    return {
        "status": "OK",
        "mechanismId": "A4",
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "nat": nat,
        "natsFound": nats,
        "natIndex": nat_index,
        "architecture": (f"{oligo_length} nt RNase-H gapmer tiled across the "
                         f"{nat['overlapNt']} nt overlap with the NAT"),
        "strandNote": (
            "The oligo is complementary to the NAT. Because the NAT is "
            "antisense to the gene, that makes it read as the gene's sense "
            "strand — correct, though it looks inverted at a glance."
        ),
        "ranking": {"orderedBy": "targetDuplexDg",
                    "caveat": ("Thermodynamic ordering, not a validated "
                               "knockdown model.")},
        "dataProvenance": {
            "nat": "Ensembl opposite-strand gene overlap",
            "sequence": "Ensembl genomic region over the real overlap only",
        },
        "candidates": candidates[:max_candidates],
    }
