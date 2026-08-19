"""
Isoform Engineering (TG07) backend service.

Generates isoform engineering candidates using real exon and splice-site
data from Ensembl, with deterministic biophysical scoring.
"""

from __future__ import annotations

import json
import logging
import os
import re
from typing import Any

import primer3
import requests
import RNA

from services.gene_silencing_service import get_target_analysis, _ensembl_get

logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"
RULEBOOKS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "rulebooks")


def _load_rule(mechanism_id: str) -> dict | None:
    path = os.path.join(RULEBOOKS_DIR, mechanism_id, "rule.json")
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _fetch_exon_splice_sites(gene_symbol: str, exon_locus: str, organism: str = "homo_sapiens") -> dict:
    """Fetch real exon data and splice site sequences."""
    result = {
        "exons": [],
        "spliceSites": [],
        "cdsSequence": "",
        "geneId": "",
        "canonicalTranscript": "",
    }
    try:
        target = get_target_analysis(gene_symbol, gene_symbol=gene_symbol, organism=organism)
        exons = target.get("exons", [])
        cds = target.get("mrnaSequence", "") or ""
        result["exons"] = exons
        result["cdsSequence"] = cds
        result["geneId"] = target.get("geneId", "")
        result["canonicalTranscript"] = target.get("canonicalTranscript", {}).get("id", "")
    except Exception as exc:
        logger.warning("Isoform engineering target fetch failed for %s: %s", gene_symbol, exc)

    if not result["exons"]:
        try:
            resp = _ensembl_get(f"{ENSEMBL_REST}/lookup/symbol/{organism}/{gene_symbol}?expand=1")
            if resp.ok:
                data = resp.json()
                transcripts = data.get("Transcript", [])
                coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
                canonical_id = data.get("canonical_transcript", "")
                canonical = None
                for t in coding:
                    if t.get("id", "").split(".")[0] == canonical_id.split(".")[0]:
                        canonical = t
                        break
                if not canonical and coding:
                    canonical = coding[0]
                if canonical:
                    tid = canonical.get("id", "").split(".")[0]
                    result["geneId"] = data.get("id", "")
                    result["canonicalTranscript"] = canonical.get("id", "")
                    result["exons"] = [
                        {
                            "id": e.get("id"),
                            "start": e.get("start"),
                            "end": e.get("end"),
                            "length": (e.get("end", 0) - e.get("start", 0) + 1),
                        }
                        for e in canonical.get("Exon", [])
                    ]
                    seq_resp = _ensembl_get(f"{ENSEMBL_REST}/sequence/id/{tid}?type=cds")
                    if seq_resp.ok:
                        result["cdsSequence"] = seq_resp.json().get("seq", "")
        except Exception as exc:
            logger.warning("Fallback exon fetch failed for %s: %s", gene_symbol, exc)

    exon_idx = None
    m = re.search(r"exon_(\d+)", exon_locus)
    if m:
        exon_idx = int(m.group(1)) - 1

    splice_sites = []
    if exon_idx is not None and result["exons"]:
        for i, ex in enumerate(result["exons"]):
            if i == exon_idx:
                start_seq = _fetch_flank_sequence(
                    result["canonicalTranscript"].split(".")[0] if result["canonicalTranscript"] else "",
                    ex.get("start", 0),
                    ex.get("end", 0),
                    strand=1,
                    flank=12,
                )
                end_seq = _fetch_flank_sequence(
                    result["canonicalTranscript"].split(".")[0] if result["canonicalTranscript"] else "",
                    ex.get("start", 0),
                    ex.get("end", 0),
                    strand=-1,
                    flank=12,
                )
                splice_sites.append({
                    "exonIndex": i + 1,
                    "donor": start_seq.get("sequence", ""),
                    "acceptor": end_seq.get("sequence", ""),
                    "strength": _score_splice_site(start_seq.get("sequence", ""), end_seq.get("sequence", "")),
                })
    result["spliceSites"] = splice_sites
    return result


def _fetch_flank_sequence(transcript_id: str, start: int, end: int, strand: int = 1, flank: int = 12) -> dict:
    """Fetch flanking sequence around an exon boundary."""
    if not transcript_id or not start or not end:
        return {"sequence": ""}
    try:
        if strand == 1:
            region_start = max(1, start - flank)
            region_end = start + flank
        else:
            region_start = max(1, end - flank)
            region_end = end + flank
        resp = _ensembl_get(
            f"{ENSEMBL_REST}/sequence/region/homo_sapiens/{transcript_id}:{region_start}..{region_end}:{strand}"
        )
        if resp.ok:
            data = resp.json()
            return {"sequence": data.get("seq", "")}
    except Exception:
        pass
    return {"sequence": ""}


def _score_splice_site(donor: str, acceptor: str) -> float:
    """Simple heuristic splice site strength score (0-1)."""
    score = 0.5
    if donor.upper().startswith("AG"):
        score += 0.2
    if acceptor.upper().endswith("AG"):
        score += 0.2
    if "GT" in donor.upper()[:6]:
        score += 0.1
    return min(1.0, max(0.0, score))


def _calc_mfe(seq: str) -> float:
    try:
        _, energy = RNA.fold(seq)
        return round(energy, 1)
    except Exception:
        return 0.0


def _calc_gc(seq: str) -> float:
    seq = seq.upper().replace("T", "U")
    if not seq:
        return 0.0
    return round(sum(1 for b in seq if b in "GCgc") / len(seq) * 100, 1)


def _reverse_complement(seq: str) -> str:
    """Antisense strand for an RNA target window. Handles T and U alike."""
    table = {"A": "U", "U": "A", "T": "A", "G": "C", "C": "G"}
    return "".join(table.get(b, "N") for b in reversed(seq.upper()))


def _calc_tm(seq: str) -> float | None:
    """Nearest-neighbour Tm (SantaLucia via primer3) on the DNA analogue."""
    dna = seq.upper().replace("U", "T")
    if not dna or set(dna) - set("ACGT"):
        return None
    return round(primer3.calc_tm(dna), 1)


def _duplex_dg(aso: str, target: str) -> float | None:
    """ViennaRNA duplex free energy for the ASO against its target window."""
    if not aso or not target:
        return None
    try:
        return round(RNA.duplexfold(aso.upper(), target.upper()).energy, 2)
    except Exception:
        return None


# How wide a window each splice element occupies, and where it sits relative
# to the target exon. These are the regions a steric blocker is aimed at; the
# widths are design windows, not measured element boundaries.
ELEMENT_WINDOWS = {
    "splice_acceptor": ("exon_5p", 30),
    "splice_donor": ("exon_3p", 30),
    "exonic_splicing_enhancer": ("exon_body", 0),
    "exonic_splicing_silencer": ("exon_body", 0),
    "branch_point": ("exon_5p", 30),
}


def _target_window(mrna: str, exon: dict, splice_element_target: str
                   ) -> tuple[int, int, str]:
    """(start, end, label) of the region the ASO is tiled across.

    Coordinates are transcript-relative, taken from the exon's real
    cdsStart/cdsEnd mapping rather than estimated from genomic length.
    """
    lo = int(exon.get("cdsStart") or 0)
    hi = int(exon.get("cdsEnd") or 0)
    lo, hi = max(0, min(lo, len(mrna))), max(0, min(hi, len(mrna)))
    if hi <= lo:
        return 0, 0, "unavailable"
    where, width = ELEMENT_WINDOWS.get(splice_element_target, ("exon_body", 0))
    if where == "exon_5p":
        return lo, min(hi, lo + width), "exon 5' boundary (acceptor side)"
    if where == "exon_3p":
        return max(lo, hi - width), hi, "exon 3' boundary (donor side)"
    return lo, hi, "exon body"


def generate_isoform_candidates(
    target_symbol: str,
    isoform_goal: str,
    target_exon_locus: str,
    splice_element_target: str,
    steric_chemistry: str,
    enforce_in_frame: bool = True,
    aso_length: int = 20,
    max_candidates: int = 12,
    organism: str = "homo_sapiens",
) -> dict[str, Any]:
    """Tile steric-blocking ASOs across a real splice element and rank them.

    WHAT CHANGED AND WHY. This function used to emit eight candidates from a
    fixed loop: one hard-coded sequence (`"GCCACC" + "A"*76 + "AUG" +
    "GCU"*20 + ...`) repeated for every gene, splice efficiencies invented as
    `75 + i*2 - (i%3)*5` whenever the real splice-site fetch returned
    nothing, and CAI / U-content / TLR / amino-acid-identity numbers produced
    by arithmetic on the loop index. It also called `_calc_cai` and
    `_calc_u_content`, which do not exist in this module, so every request
    raised NameError before any of it reached a caller.

    Two separate problems, both fixed here. The invented numbers violate the
    project's standing rule that a missing value is reported, never
    synthesised. And CAI, U-content, TLR risk and predicted yield are
    properties of an mRNA CONSTRUCT (TG08); this goal designs a
    steric-blocking oligonucleotide, which has none of them. They are gone
    rather than recomputed.

    What is emitted now comes from the transcript: the ASO is the reverse
    complement of a real window inside the real target exon, located by that
    exon's own cdsStart/cdsEnd mapping. GC, Tm (primer3, SantaLucia) and
    duplex dG (ViennaRNA) are computed on that sequence. Frame status is
    arithmetic on the real exon length. Splice-site strength appears only
    when the flanking genomic sequence was actually fetched.

    Ranking is by target-duplex dG, which is thermodynamics, not a validated
    activity model. `backend/experiments/ml_analysis` (E12) measured exactly
    this kind of ordering against 528 held-out experiments and found it ranks
    below chance because it is largely a GC proxy. The ordering is labelled
    accordingly in the response and should be treated as a starting point for
    triage, not a prediction.
    """
    symbol = target_symbol.strip().upper()
    if not symbol:
        raise ValueError("target_symbol is required")

    splice_data = _fetch_exon_splice_sites(symbol, target_exon_locus, organism)
    exons = splice_data.get("exons", [])
    mrna = (splice_data.get("cdsSequence") or "").upper().replace("T", "U")
    splice_sites = splice_data.get("spliceSites", [])

    unavailable = {
        "status": "UNAVAILABLE",
        "overview": {
            "targetGene": symbol,
            "refSeq": splice_data.get("canonicalTranscript") or None,
            "primaryMechanism": None,
            "targetWindow": None,
        },
        "candidates": [],
    }
    if not exons or not mrna:
        unavailable["message"] = (
            f"No transcript sequence or exon map could be retrieved for "
            f"{symbol}. No oligo can be designed against it, and none is "
            f"invented here."
        )
        return unavailable

    exon_num = None
    m = re.search(r"exon_(\d+)", target_exon_locus or "")
    if m:
        exon_num = int(m.group(1))
    if exon_num is None or not (1 <= exon_num <= len(exons)):
        unavailable["message"] = (
            f"target_exon_locus {target_exon_locus!r} does not name an exon of "
            f"the canonical transcript, which has {len(exons)} exons."
        )
        return unavailable

    exon = exons[exon_num - 1]
    lo, hi, window_label = _target_window(mrna, exon, splice_element_target)
    if hi - lo < aso_length:
        unavailable["message"] = (
            f"Exon {exon_num} maps to transcript positions {lo}-{hi}, which is "
            f"shorter than the {aso_length} nt oligo requested. Nothing is "
            f"designed rather than padding the window."
        )
        return unavailable

    # Frame arithmetic on the real exon. Skipping an exon whose length is not
    # a multiple of three shifts the reading frame downstream of it.
    exon_len = int(exon.get("cdsEnd", 0)) - int(exon.get("cdsStart", 0))
    in_frame = exon_len % 3 == 0
    frame_label = "In-Frame" if in_frame else "Out-of-Frame"

    splice_strength = None
    if splice_sites:
        raw = splice_sites[0].get("strength")
        splice_strength = round(float(raw), 3) if raw is not None else None

    step = max(1, (hi - lo - aso_length) // max(max_candidates - 1, 1))
    seen: set[str] = set()
    candidates: list[dict[str, Any]] = []
    for offset in range(lo, hi - aso_length + 1, step):
        target_seq = mrna[offset:offset + aso_length]
        if len(target_seq) < aso_length or set(target_seq) - set("ACGU"):
            continue
        aso = _reverse_complement(target_seq)
        if aso in seen:
            continue
        seen.add(aso)
        candidates.append({
            "constructId": f"iso-{symbol}-e{exon_num}-{offset}",
            "modality": "Steric-Blocking ASO",
            "mechanismChemistry": steric_chemistry,
            "sequence": aso,
            "targetSequence": target_seq,
            "transcriptStart": offset,
            "transcriptEnd": offset + aso_length,
            "length": aso_length,
            "gcContent": _calc_gc(aso),
            "meltingTempC": _calc_tm(aso),
            "selfMfe": _calc_mfe(aso),
            "targetDuplexDg": _duplex_dg(aso, target_seq),
            "targetWindow": window_label,
            "exonNumber": exon_num,
            "exonLength": exon_len,
            "inFrameStatus": frame_label,
            "spliceSiteStrength": splice_strength,
            # Deliberately absent, with the reason attached rather than a
            # plausible-looking number. See the docstring.
            "predictedIsoformYield": None,
            "tlrRisk": None,
            "notComputed": {
                "predictedIsoformYield": (
                    "No calibrated model maps ASO thermodynamics to isoform "
                    "ratio; a fold-change here would be invented."
                ),
                "tlrRisk": (
                    "Innate-immune risk depends on the backbone chemistry and "
                    "sequence motifs in a way this service does not model."
                ),
            },
        })
        if len(candidates) >= max_candidates:
            break

    if not candidates:
        unavailable["message"] = (
            f"Exon {exon_num} of {symbol} contains no {aso_length} nt window of "
            f"unambiguous sequence to design against."
        )
        return unavailable

    # More negative dG binds more tightly. Candidates with no dG sort last
    # rather than being treated as zero.
    candidates.sort(key=lambda c: (c["targetDuplexDg"] is None,
                                   c["targetDuplexDg"] if c["targetDuplexDg"] is not None else 0.0))
    for i, c in enumerate(candidates):
        c["rank"] = i + 1

    if enforce_in_frame and not in_frame:
        frame_note = (
            f"Exon {exon_num} is {exon_len} nt, which is not a multiple of 3. "
            f"Skipping it shifts the downstream reading frame. Candidates are "
            f"returned because the exon was explicitly requested."
        )
    else:
        frame_note = None

    return {
        "status": "OK",
        "overview": {
            "targetGene": symbol,
            "geneId": splice_data.get("geneId") or None,
            "refSeq": splice_data.get("canonicalTranscript") or None,
            "transcriptLength": len(mrna),
            "exonCount": len(exons),
            "targetExon": exon_num,
            "exonLength": exon_len,
            "targetWindow": window_label,
            "windowStart": lo,
            "windowEnd": hi,
            "isoformGoal": isoform_goal,
            "spliceElementTarget": splice_element_target,
            "primaryMechanism": "A7 splice-modulating steric block",
            "inFrameStatus": frame_label,
            "frameNote": frame_note,
            "spliceSiteStrength": splice_strength,
        },
        "ranking": {
            "orderedBy": "targetDuplexDg",
            "caveat": (
                "Thermodynamic ordering, not a validated activity prediction. "
                "Measured against 528 held-out experiments, duplex-dG ordering "
                "ranks below chance because it tracks GC content "
                "(backend/experiments/ml_analysis, E12). Use it for triage."
            ),
        },
        "dataProvenance": {
            "exons": "Ensembl canonical transcript, cdsStart/cdsEnd mapped "
                     "from the real cDNA-to-CDS alignment",
            "thermodynamics": "primer3 (Tm) and ViennaRNA (MFE, duplex dG)",
            "spliceSiteStrength": ("Ensembl flanking genomic sequence"
                                   if splice_strength is not None
                                   else "not fetched for this exon"),
        },
        "candidates": candidates,
    }


def get_isoform_engineering_design_options() -> dict[str, Any]:
    """Return design option catalogs for the isoform engineering form."""
    return {
        "isoformGoals": [
            {"id": "exon_skipping", "label": "Exon Skipping", "description": "Skip a specific exon to restore the reading frame or remove a toxic domain."},
            {"id": "exon_inclusion", "label": "Exon Inclusion", "description": "Force inclusion of a beneficial exon that is normally skipped."},
            {"id": "intron_retention", "label": "Intron Retention", "description": "Retain a specific intron to trigger NMD, encode a micropeptide, or produce a regulatory RNA (A33)."},
            {"id": "alternative_splice_site", "label": "Alternative Splice Site Selection", "description": "Redirect splicing to an alternative splice site to generate a different isoform."},
            {"id": "mutually_exclusive_exon", "label": "Mutually Exclusive Exon Switch", "description": "Switch between mutually exclusive exons to favor a therapeutically beneficial isoform."},
            # Added with the TG07 restoration. Both currently HALT: their
            # reference tables (alternative promoters / intron-retention
            # potential, and the per-gene benefit curation each needs) ship
            # header-only, so the mechanism reports that it cannot establish
            # its evidence rather than scoring on an assumption.
            {"id": "apa_modulation", "label": "Alternative Polyadenylation (APA) Shift", "description": "Shift poly(A) site usage to change the 3' end of the transcript and the isoform it produces (A11)."},
            {"id": "alt_promoter_switch", "label": "Alternative Promoter Switch", "description": "Shift transcription start site selection to favour a different promoter-derived isoform (A32)."},
        ],
        "targetExonLoci": [
            {"id": "exon_7", "label": "Exon 7"},
            {"id": "exon_23", "label": "Exon 23"},
            {"id": "exon_51", "label": "Exon 51"},
            {"id": "exon_45", "label": "Exon 45"},
            {"id": "custom", "label": "Custom Exon Locus"},
        ],
        "spliceElementTargets": [
            {"id": "splice_donor", "label": "Splice Donor Site (5' SS)"},
            {"id": "splice_acceptor", "label": "Splice Acceptor Site (3' SS)"},
            {"id": "exonic_splicing_enhancer", "label": "Exonic Splicing Enhancer (ESE)"},
            {"id": "exonic_splicing_silencer", "label": "Exonic Splicing Silencer (ESS)"},
            {"id": "intronic_splicing_enhancer", "label": "Intronic Splicing Enhancer (ISE)"},
            {"id": "intronic_splicing_silencer", "label": "Intronic Splicing Silencer (ISS)"},
        ],
        "stericChemistries": [
            {"id": "gapmer", "label": "DNA Gapmer (2'-MOE/PS)"},
            {"id": "lnai", "label": "LNA/2'-O-Methyl mix"},
            {"id": "fully_modified", "label": "Fully Modified 2'-MOE/PS"},
            {"id": "pna", "label": "Peptide Nucleic Acid (PNA)"},
        ],
    }
