"""
Gene Silencing page backend service.

Fetches transcript/exon data from Ensembl, generates ASO candidates
for a target exon, and computes biophysical metrics (GC%, Tm, self-dimer risk).
"""

from __future__ import annotations

import json
import logging
import os
import re
import time
import requests

import RNA
import primer3

ENSEMBL_REST = "https://rest.ensembl.org"
MAX_RETRIES = 3
RETRY_BACKOFF = 1.5  # seconds

# Rulebooks live as plain rule.json files under backend/rulebooks/<MECH_ID>/.
# Same location as mechanism_service.RULEBOOKS_DIR.
RULEBOOKS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "rulebooks")

logger = logging.getLogger(__name__)


def _load_silencing_rule(mechanism_id: str) -> dict | None:
    """Load a mechanism's rule.json (None if absent). Mirror of mechanism_service."""
    path = os.path.join(RULEBOOKS_DIR, mechanism_id, "rule.json")
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _ensembl_get(url: str, timeout: int = 20, retries: int = MAX_RETRIES) -> requests.Response:
    """GET with retries and exponential backoff for transient Ensembl failures."""
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
            if resp.ok or resp.status_code == 404:
                return resp
            if resp.status_code == 429:
                # Rate-limited — honour Retry-After or back off
                retry_after = resp.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else RETRY_BACKOFF * (2 ** (attempt - 1))
                logger.warning("Ensembl rate-limited (attempt %d/%d), waiting %.1fs", attempt, retries, wait)
                time.sleep(wait)
                continue
            if resp.status_code >= 500:
                last_exc = RuntimeError(f"HTTP {resp.status_code}")
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                logger.warning("Ensembl returned %d (attempt %d/%d), retrying in %.1fs",
                               resp.status_code, attempt, retries, wait)
                time.sleep(wait)
                continue
            # Client errors other than 404/429 — return as-is
            return resp
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < retries:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                logger.warning("Ensembl request failed (attempt %d/%d), retrying in %.1fs: %s",
                               attempt, retries, wait, exc)
                time.sleep(wait)
    raise RuntimeError(f"Ensembl request failed after {retries} attempts: {last_exc}")


def _parse_exons_from_transcript(transcript: dict) -> list[dict]:
    """Extract and sort exon list from an Ensembl transcript object."""
    exons = transcript.get("Exon", [])
    strand = transcript.get("strand", 1)
    sorted_exons = sorted(exons, key=lambda e: e.get("start", 0), reverse=(strand == -1))
    return [
        {
            "id": e.get("id"),
            "index": idx + 1,
            "start": e.get("start"),
            "end": e.get("end"),
            "length": (e.get("end", 0) - e.get("start", 0) + 1),
        }
        for idx, e in enumerate(sorted_exons)
        if e.get("start") and e.get("end")
    ]


def _fetch_cdna_with_cds_offset(tid: str, cds_seq: str) -> tuple[str | None, int | None]:
    """Fetch the full spliced transcript (cdna) and locate the CDS within it.

    Ensembl's REST API exposes ``type=cdna`` (5'UTR + CDS + 3'UTR) and
    ``type=cds`` but no direct UTR type or per-exon CDS offsets, so the CDS
    is located by substring search within the cdna. Returns (cdna, cds_at)
    — cds_at is the 0-based offset of the CDS start within the cdna — or
    (None, None) when the fetch fails or the CDS cannot be located within it.
    """
    try:
        resp = _ensembl_get(f"{ENSEMBL_REST}/sequence/id/{tid}?type=cdna")
        if not resp.ok:
            return None, None
        cdna = resp.json().get("seq", "")
        if not cdna or not cds_seq:
            return None, None
        cds_at = cdna.find(cds_seq)
        if cds_at < 0:
            return None, None
        # The CDS must occur exactly once in the cDNA for its offset to be
        # unambiguous. If it occurs more than once, find() silently returns
        # the first match, which may not be the real CDS — and every exon
        # offset derived from it would then be wrong while looking fine.
        # Treat that as a failed fetch so the caller falls back to the
        # labelled proportional estimate rather than trusting a guess.
        if cdna.count(cds_seq) > 1:
            logger.warning(
                "CDS sequence occurs %d times in cDNA for %s; offset is "
                "ambiguous, falling back to estimated exon mapping",
                cdna.count(cds_seq), tid,
            )
            return None, None
        return cdna, cds_at
    except Exception as exc:
        logger.warning("cDNA sequence fetch failed for %s: %s", tid, exc)
        return None, None


def _map_exons_to_cds(exons: list[dict], cdna: str, cds_seq: str, cds_at: int) -> None:
    """Attach exact CDS-relative (cdsStart, cdsEnd) offsets to each exon, in place.

    An exon's genomic length (end - start + 1) is exactly the number of
    nucleotides it contributes to the spliced transcript — introns are
    already removed. So walking the exon list (already ordered 5'->3' by
    ``_parse_exons_from_transcript``) while accumulating cDNA offsets gives
    each exon's exact position in the cDNA, with no proportional estimate.
    Intersecting that with the CDS's [cds_at, cds_at + len(cds_seq)) window
    — in the same cDNA coordinate system — gives an exact CDS-relative
    range. Exons that fall entirely in the 5'/3' UTR get a zero-width range
    clamped to the nearest CDS boundary (0 or the CDS length), so candidate
    generation naturally skips them without needing a special case.
    """
    cds_end_offset = cds_at + len(cds_seq)
    cursor = 0
    for exon in exons:
        exon_len = exon.get("length", 0)
        exon_cdna_start = cursor
        exon_cdna_end = cursor + exon_len
        cursor = exon_cdna_end

        start = min(max(exon_cdna_start, cds_at), cds_end_offset) - cds_at
        end = min(max(exon_cdna_end, cds_at), cds_end_offset) - cds_at
        exon["cdsStart"] = start
        exon["cdsEnd"] = end


# ---------------------------------------------------------------------------
# 1. Target Analysis — transcript / exon structure for the confirmed gene
# ---------------------------------------------------------------------------

def get_target_analysis(ensembl_gene_id: str, gene_symbol: str = "", organism: str = "") -> dict:
    """Return transcript and exon structure for the gene.

    Uses the Ensembl /lookup/id endpoint with ``expand=1`` to retrieve exon
    coordinates for the canonical (or first coding) transcript.  Falls back to
    a symbol-based lookup when the ID lookup returns no exon data.
    """
    result = {
        "geneId": ensembl_gene_id,
        "canonicalTranscript": None,
        "totalCodingTranscripts": 0,
        "exons": [],
        "cdsLength": None,
        "mrnaSequence": None,
        "utr5Sequence": None,
        "utr3Sequence": None,
    }

    # --- Primary: lookup by Ensembl gene ID ---
    try:
        resp = _ensembl_get(f"{ENSEMBL_REST}/lookup/id/{ensembl_gene_id}?expand=1")
        if resp.ok:
            data = resp.json()
            transcripts = data.get("Transcript", [])
            coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
            result["totalCodingTranscripts"] = len(coding)

            canonical_id = data.get("canonical_transcript", "")
            canonical = None
            for t in coding:
                if t.get("id", "").split(".")[0] == canonical_id.split(".")[0]:
                    canonical = t
                    break
            if not canonical and coding:
                canonical = coding[0]

            if canonical:
                result["canonicalTranscript"] = {
                    "id": canonical.get("id"),
                    "biotype": canonical.get("biotype"),
                    "chromosome": canonical.get("seq_region_name"),
                    "start": canonical.get("start"),
                    "end": canonical.get("end"),
                    "strand": canonical.get("strand", 1),
                }
                result["exons"] = _parse_exons_from_transcript(canonical)

                # Fetch CDS length from the sequence endpoint
                tid = canonical.get("id", "").split(".")[0]
                try:
                    seq_resp = _ensembl_get(f"{ENSEMBL_REST}/sequence/id/{tid}?type=cds")
                    if seq_resp.ok:
                        seq = seq_resp.json().get("seq", "")
                        result["cdsLength"] = len(seq)
                        result["mrnaSequence"] = seq
                except Exception as exc:
                    logger.warning("CDS sequence fetch failed for %s: %s", tid, exc)

                # Full spliced transcript (5'UTR + CDS + 3'UTR) so
                # upregulation mechanisms can target real UTR context
                # instead of CDS approximations.
                if result.get("mrnaSequence"):
                    cdna, cds_at = _fetch_cdna_with_cds_offset(tid, result["mrnaSequence"])
                    if cdna is not None:
                        result["utr5Sequence"] = cdna[:cds_at]
                        result["utr3Sequence"] = cdna[cds_at + len(result["mrnaSequence"]):]
                        _map_exons_to_cds(result["exons"], cdna, result["mrnaSequence"], cds_at)
        else:
            logger.warning("Ensembl ID lookup returned HTTP %d for %s", resp.status_code, ensembl_gene_id)
    except Exception as exc:
        logger.warning("Ensembl ID lookup failed for %s: %s", ensembl_gene_id, exc)

    # --- Fallback: if no exons retrieved, try symbol-based lookup ---
    if not result["exons"] and gene_symbol:
        logger.info("No exons from ID lookup, trying symbol fallback for %s", gene_symbol)
        species = (organism or "homo_sapiens").lower().replace(" ", "_")
        try:
            sym_resp = _ensembl_get(
                f"{ENSEMBL_REST}/lookup/symbol/{species}/{gene_symbol}?expand=1"
            )
            if sym_resp.ok:
                sym_data = sym_resp.json()
                transcripts = sym_data.get("Transcript", [])
                coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
                result["totalCodingTranscripts"] = len(coding)

                canonical = next((t for t in coding if t.get("is_canonical")), coding[0] if coding else None)
                if canonical:
                    result["canonicalTranscript"] = {
                        "id": canonical.get("id"),
                        "biotype": canonical.get("biotype"),
                        "chromosome": canonical.get("seq_region_name"),
                        "start": canonical.get("start"),
                        "end": canonical.get("end"),
                        "strand": canonical.get("strand", 1),
                    }
                    result["exons"] = _parse_exons_from_transcript(canonical)

                    if not result["cdsLength"]:
                        tid = canonical.get("id", "").split(".")[0]
                        try:
                            seq_resp = _ensembl_get(f"{ENSEMBL_REST}/sequence/id/{tid}?type=cds")
                            if seq_resp.ok:
                                seq = seq_resp.json().get("seq", "")
                                result["cdsLength"] = len(seq)
                                result["mrnaSequence"] = seq
                        except Exception as exc:
                            logger.warning("CDS sequence fetch failed (symbol fallback) for %s: %s", tid, exc)

                    if result.get("mrnaSequence"):
                        tid = canonical.get("id", "").split(".")[0]
                        cdna, cds_at = _fetch_cdna_with_cds_offset(tid, result["mrnaSequence"])
                        if cdna is not None:
                            result["utr5Sequence"] = cdna[:cds_at]
                            result["utr3Sequence"] = cdna[cds_at + len(result["mrnaSequence"]):]
                            _map_exons_to_cds(result["exons"], cdna, result["mrnaSequence"], cds_at)
            else:
                logger.warning("Ensembl symbol lookup returned HTTP %d for %s", sym_resp.status_code, gene_symbol)
        except Exception as exc:
            logger.warning("Ensembl symbol lookup failed for %s: %s", gene_symbol, exc)

    if not result["exons"]:
        logger.error("No exon data retrieved for gene %s (ID: %s)", gene_symbol or "?", ensembl_gene_id)

    return result


# ---------------------------------------------------------------------------
# 2. ASO Candidate Generation
# ---------------------------------------------------------------------------

CHEMISTRY_OPTIONS = [
    {"id": "gapmer", "label": "DNA Gapmer (2-10-2)", "description": "Standard RNase H1-recruiting backbone; most validated.",
     "detail": "A chimeric oligonucleotide with a central DNA 'gap' of ~10 nucleotides flanked by 2\u2032-modified wings (typically 2\u2032-O-Me or LNA). The DNA gap recruits RNase H1 to cleave the target RNA, while the wings confer nuclease resistance and binding affinity. Gapmers are the most clinically validated ASO chemistry (e.g., nusinersen/Spinraza, eteplirsen)."},
    {"id": "pmo", "label": "PMO (Phosphorodiamidate Morpholino)", "description": "Steric blocker; splice-switching, no RNase H.",
     "detail": "A non-ionic backbone oligomer where each nucleoside is linked via phosphorodiamidate bonds to morpholine rings. PMOs do not recruit RNase H; instead they sterically block RNA interactions (splice junctions, ribosome binding, miRNA binding). Used in exon-skipping (e.g., eteplirsen) and translational arrest. Requires cell-penetrating peptide (CPP) conjugation for efficient uptake."},
    {"id": "lna_gapmer", "label": "LNA-enhanced Gapmer", "description": "Locked nucleic acid wings boost binding affinity and nuclease resistance.",
     "detail": "A DNA gapmer where the flanking wings contain Locked Nucleic Acids (LNA) \u2014 bicyclic RNA analogues with a methylene bridge locking the ribose in a C3\u2032-endo conformation. Each LNA substitution raises Tm by ~2\u2032-8\u2032C, dramatically increasing target affinity. LNA gapmers also have enhanced nuclease resistance. Used in miravirsen (anti-miR-122) and bepetamers. Higher off-target risk due to increased potency."},
    {"id": "2ome", "label": "2\u2032-O-Methoxyethyl (2\u2032-OMe)", "description": "Steric blocker; splicing modulation and miRNA inhibition.",
     "detail": "A ribose-modified oligonucleotide where the 2\u2032-OH is replaced with a methoxyethyl group. 2\u2032-OMe ASOs sterically block RNA interactions and are commonly used for splice-switching and miRNA inhibition. They have good nuclease resistance, low toxicity, and are often used in combination (e.g., morpholino-2\u2032-OMe mixmers). Lower binding affinity than LNA but fewer off-target effects."},
]

MODIFICATION_OPTIONS = [
    {"id": "phosphorothioate", "label": "Phosphorothioate (PS) backbone", "description": "Increases nuclease resistance and protein binding.",
     "detail": "Replaces one non-bridging oxygen in the phosphodiester backbone with sulfur. PS linkages dramatically increase nuclease resistance and promote protein binding (e.g., to albumin), extending plasma half-life from minutes to hours. However, PS backbone can increase off-target binding to unintended RNA sequences and may activate complement pathways at high doses. Most ASO drugs incorporate full or partial PS backbones."},
    {"id": "lna_wings", "label": "LNA wings (5\u2032 + 3\u2032)", "description": "Locked nucleic acids at terminal positions for higher Tm.",
     "detail": "Locked Nucleic Acid (LNA) modifications placed at the 5\u2032 and 3\u2032 ends of the oligonucleotide. Each LNA substitution raises the melting temperature (Tm) by 2\u2032-8\u2032C, increasing binding affinity to the target RNA. LNA wings also provide strong nuclease resistance at terminal positions. Typically used in gapmer wings (2-3 LNA at each end). Excessive LNA use can increase off-target activity due to hyper-stabilized binding."},
    {"id": "2omemod", "label": "2\u2032-OMe wing modifications", "description": "Ribose modification at terminal positions.",
     "detail": "2\u2032-O-Methoxyethyl or 2\u2032-O-Methyl modifications at the 5\u2032 and 3\u2032 wings. These RNA-like modifications increase nuclease resistance, improve binding affinity (Tm increase ~1\u2032-2\u2032C per substitution), and reduce immune stimulation compared to unmodified DNA. Commonly used in gapmer wings as a lower-cost alternative to LNA. Well-tolerated clinically with a favorable safety profile."},
    {"id": "pmo_core", "label": "PMO core", "description": "Non-ionic backbone; splice-switching.",
     "detail": "Phosphorodiamidate Morpholino (PMO) modifications at the central region. PMOs are non-ionic, avoiding non-specific protein interactions and reducing toxicity. They work by steric blocking rather than RNase H recruitment. The PMO core is used in splice-switching ASOs (e.g., exon-skipping for DMD) where you want to block splice junctions without degrading the RNA. Requires CPP conjugation (e.g., Vivolen/PP-PMO) for cellular uptake."},
    {"id": "pna_clamp", "label": "PNA clamp (flanking)", "description": "Peptide nucleic acid clamps to block nuclease access.",
     "detail": "Peptide Nucleic Acid (PNA) modifications flanking the ASO. PNAs have a peptide-like backbone instead of sugar-phosphate, making them extremely resistant to nucleases and proteases. PNA clamps create a 'steric shield' around the ASO core, protecting it from exonuclease degradation. They also increase binding affinity (Tm ~1\u2032C per base). PNA synthesis is expensive and delivery is challenging \u2014 typically used in research or specialized applications."},
]

LENGTH_RANGE = {"min": 12, "max": 30, "default": 18, "step": 1}

# Minimum GC fraction for a valid candidate
MIN_GC = 0.30
MAX_GC = 0.70

def _calc_gc(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def _calc_tm(seq: str) -> float:
    """Real Tm (°C) using nearest-neighbor thermodynamics (SantaLucia 1998).

    Uses primer3's calc_tm which implements the unified nearest-neighbor model
    with salt corrections — the same method used by IDT OligoAnalyzer, Primer3,
    and other professional oligo design tools.

    U is folded to T first. primer3 carries DNA nearest-neighbour parameters
    and raises ValueError on any non-ACGT base, so an RNA-alphabet oligo
    reaching it is a hard failure, not a degraded number. That is not
    hypothetical: `rna_processing_service` renders its sequences in the RNA
    alphabet (`replace("T", "U")`) before calling this, which made every
    TG04 design request 422 with "Sequence contains non-ACGT base 'U'".

    Folding U to T is the right normalisation rather than a workaround: the
    value being computed IS the DNA-analogue Tm. It does not model 2'-MOE,
    cEt, LNA or a morpholino backbone — `_effective_tm_boost` applies the
    chemistry adjustment separately, on top of this.
    """
    if not seq:
        return 0.0
    dna = seq.upper().replace("U", "T")
    if set(dna) - set("ACGT"):
        # Ambiguity codes and gaps have no nearest-neighbour parameters.
        # Returning 0.0 would read as "melts at zero"; refuse instead.
        raise ValueError(
            f"Tm is undefined for a sequence with non-ACGTU bases: "
            f"{sorted(set(dna) - set('ACGT'))}"
        )
    return round(primer3.calc_tm(dna), 1)


def _self_complement_mfe(seq: str) -> float:
    """Real MFE (minimum free energy) for self-structure prediction.

    Uses ViennaRNA's RNA.fold() which implements the Zuker algorithm for
    thermodynamic secondary structure prediction. More negative MFE indicates
    stronger tendency to form hairpins or other intramolecular structures.
    """
    if not seq or len(seq) < 4:
        return 0.0
    _, mfe = RNA.fold(seq.upper())
    return round(mfe, 2)


def _polyg_score(seq: str) -> int:
    """Number of G-tracts ≥3 in the sequence."""
    return len(re.findall(r"G{3,}", seq.upper()))


def _cpg_count(seq: str) -> int:
    """Count CpG dinucleotides (immune stimulation risk)."""
    return len(re.findall(r"CG", seq.upper()))


def _longest_homopolymer(seq: str) -> int:
    """Length of the longest homopolymer run."""
    seq = seq.upper()
    if not seq:
        return 0
    max_run = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 1
    return max_run


def _purine_content(seq: str) -> float:
    """Fraction of purines (A+G) in the sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    purines = sum(1 for b in seq if b in "AG")
    return round(purines / len(seq), 3)


def _sequence_complexity(seq: str) -> float:
    """Shannon entropy-based complexity score (0-1). Higher = more unique."""
    seq = seq.upper()
    if not seq:
        return 0.0
    from collections import Counter
    import math
    counts = Counter(seq)
    n = len(seq)
    entropy = -sum((c / n) * math.log2(c / n) for c in counts.values())
    max_entropy = math.log2(4)  # 4 nucleotides
    return round(entropy / max_entropy, 3)


def _gc_skew(seq: str) -> float:
    """GC skew: (G-C)/(G+C). Measures strand bias."""
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    if g + c == 0:
        return 0.0
    return round((g - c) / (g + c), 3)


def _molecular_weight(seq: str) -> float:
    """Approximate molecular weight (g/mol) for single-stranded DNA."""
    # Average MW of DNA nucleotides: A=331.2, T=322.2, G=347.2, C=307.2
    weights = {"A": 331.2, "T": 322.2, "G": 347.2, "C": 307.2}
    seq = seq.upper()
    mw = sum(weights.get(b, 320) for b in seq)
    # Subtract water for phosphodiester bonds
    mw -= (len(seq) - 1) * 18.0
    return round(mw, 1)


def _extinction_coefficient(seq: str) -> float:
    """UV extinction coefficient at 260nm (L/mol·cm) using nearest-neighbor method."""
    seq = seq.upper()
    # Simplified: sum of nearest-neighbor coefficients
    nn = {
        "AA": 27400, "AT": 22300, "AG": 25000, "AC": 21200,
        "TA": 23500, "TT": 16700, "TG": 20900, "TC": 16200,
        "GA": 25500, "GT": 20800, "GG": 22600, "GC": 17600,
        "CA": 21900, "CT": 16900, "CG": 18900, "CC": 15200,
    }
    total = 8800  # initiation factor
    for i in range(len(seq) - 1):
        total += nn.get(seq[i:i+2], 19000)
    return round(total, 0)


def _nuclease_resistance_score(chemistry: str, modifications: list[str]) -> float:
    """Estimated nuclease resistance (0-100). Higher = more resistant."""
    base = 20  # Unmodified DNA
    chem_scores = {"gapmer": 55, "lna_gapmer": 70, "pmo": 85, "2ome": 60}
    mod_scores = {"phosphorothioate": 25, "lna_wings": 20, "2omemod": 15, "pmo_core": 30, "pna_clamp": 35}
    score = base + chem_scores.get(chemistry, 0)
    for m in modifications:
        score += mod_scores.get(m, 0)
    return min(100, round(score, 1))


def _cellular_uptake_score(chemistry: str, length: int) -> float:
    """Estimated cellular uptake efficiency (0-100)."""
    # Shorter oligos cross membranes better
    length_factor = max(0, 100 - (length - 15) * 5)
    chem_factors = {"gapmer": 60, "lna_gapmer": 65, "pmo": 40, "2ome": 55}
    base = chem_factors.get(chemistry, 50)
    return round((base + length_factor) / 2, 1)


def _bbb_crossing_score(chemistry: str, length: int, modifications: list[str]) -> float:
    """Estimated blood-brain barrier crossing potential (0-100)."""
    base = 10  # Most ASOs don't cross BBB well
    if chemistry == "pmo":
        base += 20  # PMOs with CPP can cross
    if chemistry == "lna_gapmer":
        base += 5
    if "pna_clamp" in modifications:
        base += 15
    if length <= 18:
        base += 10  # Shorter = better BBB crossing
    return min(100, round(base, 1))


def _synthesis_difficulty(seq: str, chemistry: str, modifications: list[str]) -> float:
    """Estimated synthesis difficulty (0-100). Higher = harder to synthesize."""
    seq = seq.upper()
    difficulty = 10  # Base difficulty
    # Homopolymer runs make synthesis harder
    longest_hp = _longest_homopolymer(seq)
    difficulty += longest_hp * 5
    # GC-rich sequences are harder
    gc = _calc_gc(seq)
    if gc > 0.6:
        difficulty += 15
    elif gc > 0.55:
        difficulty += 8
    # Longer = harder
    difficulty += max(0, (len(seq) - 20) * 2)
    # Chemistry complexity
    chem_factor = {"gapmer": 5, "lna_gapmer": 15, "pmo": 10, "2ome": 8}
    difficulty += chem_factor.get(chemistry, 5)
    # Modification complexity
    mod_factor = {"phosphorothioate": 5, "lna_wings": 10, "2omemod": 5, "pmo_core": 8, "pna_clamp": 15}
    for m in modifications:
        difficulty += mod_factor.get(m, 0)
    return min(100, round(difficulty, 1))


def _off_target_risk(seq: str, complexity: float) -> float:
    """Estimated off-target binding risk (0-100). Higher = more risk."""
    risk = 20  # Base risk
    # Low complexity = higher off-target risk
    if complexity < 0.7:
        risk += 30
    elif complexity < 0.8:
        risk += 15
    # Poly-G tracts increase non-specific binding
    pg = _polyg_score(seq)
    risk += pg * 10
    # Very short sequences have more off-targets
    if len(seq) < 18:
        risk += 15
    return min(100, round(risk, 1))


def _immune_stimulation_risk(seq: str, chemistry: str) -> float:
    """Estimated immune stimulation risk (0-100)."""
    risk = 5
    cpg = _cpg_count(seq)
    risk += cpg * 15
    # Unmodified DNA is more immunostimulatory
    if chemistry == "gapmer":
        risk += 10
    elif chemistry == "pmo":
        risk -= 3  # PMOs are less immunostimulatory
    return min(100, max(0, round(risk, 1)))


def _duplex_stability(gc: float, tm: float, length: int) -> str:
    """Estimated duplex stability category."""
    # Free energy approximation
    dG = -0.36 * gc - 0.0048 * (tm + 273.15)
    dG_total = dG * length
    if dG_total < -30:
        return "Very Stable"
    elif dG_total < -20:
        return "Stable"
    elif dG_total < -10:
        return "Moderate"
    else:
        return "Weak"


def _defect_scores(defect_type: str | None, silencing_scope: str | None) -> dict:
    """Informational text for the defect type and silencing scope.

    Notes-only: no numeric bonuses vote into the ranking score. Numeric
    defect/chemistry preferences were a re-implementation of mechanism
    eligibility (already enforced by DEFECT_COMPATIBILITY upstream) with
    invented magnitudes, so they were removed rather than tuned.
    """
    defect = (defect_type or "").lower().strip().replace("_", "-").replace(" ", "-")
    scope = (silencing_scope or "").lower().strip()

    notes = "No defect-specific note available."

    if "loss-of-function" in defect or "lof" in defect:
        notes = "Loss-of-function: knockdown of the toxic/dysfunctional transcript may be appropriate; confirm eligibility in the mechanism rulebook."
    elif "haploinsufficiency" in defect:
        notes = "Haploinsufficiency: gene silencing can worsen the phenotype — consider upregulation mechanisms instead."
    elif "dominant" in defect:
        notes = "Dominant-negative: allele-specific silencing of the mutant allele is preferred; position the ASO over the variant."
    elif "gain-of-function" in defect or "gof" in defect:
        notes = "Gain-of-function: substantial reduction of the toxic protein is needed; the mechanism rulebook requires RNase H-recruiting chemistry for degradation."
    elif "toxic" in defect or "viral" in defect:
        notes = "Toxic RNA: transcript degradation is preferred; RNase H-recruiting chemistry applies where the rulebook requires it."
    elif "mirna-dysregulation" in defect or "mirna-dysreg" in defect:
        notes = "Pathogenic microRNA dysregulation: anti-miR silencing of the toxic miRNA; the A12 rulebook requires RNase H-independent steric chemistry."
    elif "overexpression" in defect or "oncogene" in defect:
        notes = "Gene overexpression: deep knockdown is needed; the A15 rulebook describes RNase H-compatible gapmer designs."
    elif "splice" in defect or "exon" in defect or "pseudoexon" in defect or "apa" in defect:
        notes = "Splice defect: splicing correction uses steric-blocking (RNase H-inactive) chemistry — a different modality than gene silencing."
    elif "nonsense" in defect or "premature" in defect:
        notes = "Nonsense mutation: exon-skipping ASOs may bypass the premature stop codon (RNA-processing mechanism)."
    elif "frameshift" in defect:
        notes = "Frameshift mutation: exon-skipping or transcript-degradation strategies may apply; eligibility is set by the mechanism rulebook."

    if scope == "total_knockdown":
        notes += " Total-transcript knockdown selected — broad targeting across the CDS."
    elif scope == "allele_specific":
        notes += " Allele-specific silencing selected — candidates are ranked to prioritise those spanning the variant coordinate."

    return {"defect_notes": notes}


def _estimated_binding_energy(gc_content: float, tm: float) -> float:
    """Estimated binding free energy (kcal/mol) from GC% and Tm."""
    # Simplified nearest-neighbor approximation
    # ΔG ≈ -RT ln(K) where K relates to Tm
    # Using empirical: ΔG ≈ -0.01 * Tm * seq_length (rough kcal/mol per bp)
    import math
    R = 1.987e-3  # kcal/(mol·K)
    Tm_K = tm + 273.15
    # Rough estimate: more negative = more stable
    dG = -0.36 * gc_content - 0.0048 * Tm_K
    return round(dG * 21, 1)  # scale to ~21-mer length


def _target_duplex_energy(candidate_seq: str, target_seq: str) -> float:
    """Real predicted binding free energy between ASO candidate and target.

    Uses ViennaRNA's duplexfold() which computes the minimum free energy
    of the duplex formed between the ASO candidate and its target mRNA region.
    More negative ΔG indicates stronger predicted target engagement.
    """
    if not candidate_seq or not target_seq:
        return 0.0
    duplex = RNA.duplexfold(candidate_seq.upper(), target_seq.upper())
    return round(duplex.energy, 2)


def _reverse_complement(seq: str) -> str:
    """Antisense oligonucleotide for a target site, in the input's alphabet.

    U USED TO PASS THROUGH UNCOMPLEMENTED. The translation table covered
    "ATGC" only, so a uracil in an RNA-alphabet target was copied into the
    oligo verbatim instead of pairing to adenine. `rna_processing_service`
    renders its splice-junction and polyadenylation windows in the RNA
    alphabet, so every TG04 candidate was non-complementary at each U of its
    target — a designed oligo that does not bind what it is aimed at, emitted
    with a composite score as if it did.

    The output alphabet follows the input: an RNA target yields an RNA-alphabet
    oligo, a DNA target a DNA-alphabet one. Mixing the two in one string, which
    is what the old behaviour produced, is not a representation of anything.
    """
    up = seq.upper()
    is_rna = "U" in up and "T" not in up
    pairs = {"A": "U" if is_rna else "T", "T": "A", "U": "A",
             "G": "C", "C": "G"}
    out = []
    for base in reversed(up):
        try:
            out.append(pairs[base])
        except KeyError:
            # An ambiguity code has no single complement. Emitting "N" would
            # look like a designed base; refuse the whole oligo instead.
            raise ValueError(
                f"Cannot reverse-complement {base!r} in {seq!r}: only "
                f"A, C, G, T and U have a defined complement."
            ) from None
    return "".join(out)


# Estimated Tm boost (°C) contributed by each chemistry's affinity-modifying
# backbone/wing design, applied to every candidate of that chemistry.
#
# The values reflect documented, cited per-base effects (not invented):
#   - LNA: ~2-8°C per substituted base (documented LNA literature; +12 for a
#     short LNA-winged gapmer is ~3°C/base — inside that cited range).
#   - 2'-OMe: ~0.5-1°C per modified base; the flat +6 for a full-length
#     2'-OMe ASO is deliberately conservative.
#   - PMO: neutral binding profile.
#   - Unmodified-DNA gapmer: primer3's Tm is already the plain DNA value, so 0.
CHEM_TM_BOOST = {
    "gapmer": 0,       # unmodified DNA wings
    "lna_gapmer": 12,  # ~2-8°C per LNA base, cited
    "pmo": 0,          # morpholino — neutral binding profile
    "2ome": 6,         # ~0.5-1°C per 2'-OMe base, conservative flat value
}

MOD_TM_BOOST = {
    "lna_wings": 8,   # additional LNA at 5'/3'
    "2omemod": 5,     # 2'-OMe wing modifications
    "pna_clamp": 10,  # PNA clamps raise affinity strongly
}

# Per-mechanism effective-Tm windows are INTERNAL DESIGN TARGETS derived from
# each rulebook's affinity description (A2 = "maintain high binding affinity"
# steric block; A12 = "optimize affinity" anti-miR; A1/A15 = moderate-affinity
# RNase H gapmer). They are mechanism-specific guidance for the designer, NOT
# cited thresholds — the _tm_fit_score term built on them is a heuristic.
OPTIMAL_TM_RANGES = {
    "A1": (50, 65),   # RNase H1 cleavage — moderate affinity
    "A2": (62, 75),   # translation block — tight binding at AUG
    "A12": (55, 68),  # anti-miR — strong complementarity
    "A15": (55, 68),  # promoter ASO — moderate
}


def _effective_tm_boost(chemistry: str, modifications: list[str]) -> float:
    """Estimated effective Tm increase from the chosen chemistry + modifications."""
    return CHEM_TM_BOOST.get(chemistry, 0) + sum(
        MOD_TM_BOOST.get(m, 0) for m in modifications
    )


def _tm_fit_score(
    tm: float, chemistry: str, modifications: list[str], mechanism_id: str
) -> float:
    """How well a candidate's chemistry-adjusted Tm fits the mechanism's
    optimal affinity window (0-100). Because the boost is chemistry-dependent
    and the base Tm is window-dependent, this term re-ranks candidates when the
    chemistry or modifications change."""
    adjusted_tm = tm + _effective_tm_boost(chemistry, modifications)
    lo, hi = OPTIMAL_TM_RANGES.get(mechanism_id, (50, 70))
    if lo <= adjusted_tm <= hi:
        dist = 0
    else:
        dist = min(abs(adjusted_tm - lo), abs(adjusted_tm - hi))
    return round(max(0.0, 100.0 - dist * 6.0), 1)


def _composite_score(dg: float, tm_fit: float) -> float:
    """Ranking score built exclusively from real, physics-based metrics:
    the ViennaRNA target duplex ΔG and the chemistry-adjusted Tm fit.

    Heuristic drug-like estimates (nuclease, uptake, BBB, off-target, immune,
    synthesis) and mechanism/defect/tissue/allele point bonuses are deliberately
    excluded from the sort order. They are surfaced to the user as labeled
    estimates in ``heuristicEstimates`` rather than silently voted into the
    ranking.
    """
    # Normalize target duplex ΔG (more negative = stronger binding).
    # Typical perfect-match ΔG for 12-30 nt oligos is roughly -8 to -40 kcal/mol.
    duplex_score = min(100.0, max(0.0, (-dg - 8.0) * 3.5))
    # Primary signal is duplex binding; Tm fit is the secondary tiebreaker.
    return round(0.65 * duplex_score + 0.35 * tm_fit, 1)


def _mechanism_design_constraints(mechanism_id: str, aso_length: int, chemistry: str) -> tuple[int, str, str]:
    """Validate mechanism-specific inputs and return effective design settings.

    Some mechanisms cannot be designed from a coding mRNA sequence alone.
    Refusing those requests is intentional: returning a CDS-derived sequence for
    a microRNA or promoter-associated-RNA mechanism would be biologically wrong.
    """
    if mechanism_id == "A1":
        return aso_length, chemistry, "mrna"
    if mechanism_id == "A2":
        return aso_length, chemistry, "translation_start"
    if mechanism_id == "A21":
        # A21 (RNAi / siRNA) is a double-stranded duplex modality, not a
        # single-stranded ASO. The ASO designer only produces single-stranded
        # candidates, so A21 requests are refused here rather than silently
        # emitting a chemically-incorrect design.
        raise ValueError(
            "A21 (RNA interference) requires a double-stranded siRNA duplex, which "
            "the single-stranded ASO designer does not support. Choose a "
            "single-stranded mechanism: A1, A2, A12, or A15."
        )
    if mechanism_id == "A12":
        # Anti-miR: cannot design from CDS — return gracefully
        return aso_length, chemistry, "mrna"
    if mechanism_id == "A15":
        # Promoter-targeting: cannot design from CDS — return gracefully
        return aso_length, chemistry, "mrna"
    raise ValueError(f"Unsupported gene-silencing mechanism: {mechanism_id}")


# Chemistry → RNase H capability, as stated in this file's own
# CHEMISTRY_OPTIONS descriptions (gapmer "recruits RNase H1"; PMO "do not
# recruit RNase H"; 2'-OMe steric blocker; LNA gapmer is a gapmer with LNA
# wings, so still RNase H-recruiting). These groupings are not invented — they
# are the mechanism of action the chemistry options themselves declare.
RNASE_H_ACTIVE_CHEMISTRIES = {"gapmer", "lna_gapmer"}
RNASE_H_INDEPENDENT_CHEMISTRIES = {"pmo", "2ome"}


def _mechanism_rnase_h_requirement(mechanism_id: str) -> str | None:
    """Classify a mechanism's chemistry requirement from its own rule.json.

    Reads the ``asoChemistry``/``limitations`` text (not an invented preference
    table). Returns ``"rnase_h"`` for RNase H-recruiting mechanisms,
    ``"steric"`` for RNase H-inactive/independent mechanisms, or None when the
    rulebook gives no signal.
    """
    rule = _load_silencing_rule(mechanism_id)
    if not rule:
        return None
    haystack = (
        str(rule.get("asoChemistry") or "")
        + " "
        + str(rule.get("limitations") or "")
    ).lower()
    if any(k in haystack for k in ("rnase h-inactive", "rnase h-independent")):
        return "steric"
    if any(k in haystack for k in ("rnase h-compatible", "rnase h-active", "requires rnase h1")):
        return "rnase_h"
    return None


def _mechanism_chemistry_compatibility(mechanism_id: str, chemistry: str) -> str | None:
    """Hard gate: is the chosen chemistry mechanically compatible with the
    mechanism, per the mechanism's own rulebook?

    A steric-blocking chemistry (PMO, 2'-OMe) cannot perform RNase
    H-mediated degradation and an RNase H-recruiting chemistry cannot act as a
    pure steric block, so a mismatch is rejected (ValueError) rather than
    scored as merely "suboptimal" — emitting a mis-designed ASO is worse than
    refusing. Returns the requirement string when a rulebook classification
    exists and the chemistry is unclassified (no gate possible).
    """
    requirement = _mechanism_rnase_h_requirement(mechanism_id)
    if requirement is None:
        return None
    if chemistry in RNASE_H_ACTIVE_CHEMISTRIES:
        capability = "rnase_h"
    elif chemistry in RNASE_H_INDEPENDENT_CHEMISTRIES:
        capability = "steric"
    else:
        return requirement
    if requirement != capability:
        rule = _load_silencing_rule(mechanism_id)
        wanted = "RNase H-recruiting (e.g. gapmer, LNA gapmer)" if requirement == "rnase_h" else "steric-blocking, RNase H-independent (e.g. PMO, 2'-OMe)"
        raise ValueError(
            f"Chemistry '{chemistry}' cannot serve mechanism {mechanism_id}: "
            f"the {mechanism_id} design rulebook specifies {wanted} chemistry "
            f"({rule.get('asoChemistry') if rule else 'no rulebook'}). Choose a compatible chemistry."
        )
    return requirement


def _mechanism_note(mechanism_id: str, chemistry: str) -> str:
    """Short per-candidate mechanism note derived from the rulebook's own text."""
    rule = _load_silencing_rule(mechanism_id)
    name = rule.get("name") if rule else mechanism_id
    requirement = _mechanism_rnase_h_requirement(mechanism_id)
    if requirement == "rnase_h":
        return f"{name} ({mechanism_id}): RNase H-recruiting chemistry required — {chemistry} is compatible."
    if requirement == "steric":
        return f"{name} ({mechanism_id}): RNase H-independent steric chemistry required — {chemistry} is compatible."
    return f"{name} ({mechanism_id})"


# ---------------------------------------------------------------------------
# HGVS coding-DNA (c.) parsing — used for allele-specific candidate design
# ---------------------------------------------------------------------------

# Regex-based subset of HGVS c. notation. This is a *simplified* parser that
# extracts a CDS position and variant type for the design pipeline; it is not a
# full HGVS-grammar implementation (no transcript/reference validation).
HGVS_C_PATTERN = re.compile(
    r"^c\.(?P<start>\d+)"
    r"(?:_(?P<end>\d+))?"
    r"(?:"
    r"(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)"                              # substitution
    r"|del(?P<del_seq>[ACGT]*)"                                       # deletion
    r"|dup(?P<dup_seq>[ACGT]*)"                                       # duplication
    r"|ins(?P<ins_seq>[ACGT]+)"                                       # insertion
    r"|del(?P<delins_del>[ACGT]*)ins(?P<delins_ins>[ACGT]+)"          # delins
    r")$",
    re.IGNORECASE,
)

HGVS_REJECT_PATTERNS = {
    "intronic_offset": re.compile(r"^c\.\d+[+-]\d+", re.IGNORECASE),
    "five_prime_utr": re.compile(r"^c\.-\d+", re.IGNORECASE),
    "three_prime_utr": re.compile(r"^c\.\*\d+", re.IGNORECASE),
}

HGVS_REJECT_REASONS = {
    "intronic_offset": "Deep intronic position — not expressible against the fetched CDS-only sequence.",
    "five_prime_utr": "5' UTR position — outside the fetched CDS.",
    "three_prime_utr": "3' UTR position — outside the fetched CDS.",
}


def parse_hgvs_c(raw: str) -> dict:
    """Parse a simplified HGVS coding-DNA (c.) variant into a CDS coordinate.

    HGVS defines c.1 as the A of the ATG start codon, which maps directly to
    index 0 of the CDS sequence fetched by ``get_target_analysis()`` — no
    separate UTR-length lookup is needed.

    Returns ``{"parsed": True, "type", "cdsStart", "cdsEnd", "length"}`` for
    recognized c. notation, or ``{"parsed": False, "reason"}`` otherwise.
    Intronic offsets, 5' UTR (c.-N), and 3' UTR (c.*N) positions cannot be
    expressed against a CDS-only sequence and are rejected with a specific
    reason. Protein (p.) notation is out of scope: it does not map to a CDS
    index without a codon-table lookup.
    """
    variant = raw.strip()
    if not variant:
        return {"parsed": False, "reason": "No variant provided."}

    lowered = variant.lower()
    for key, pattern in HGVS_REJECT_PATTERNS.items():
        if pattern.match(lowered):
            return {"parsed": False, "reason": HGVS_REJECT_REASONS[key]}

    if lowered.startswith("p."):
        return {
            "parsed": False,
            "reason": (
                "Protein-level (p.) notation is out of scope — provide the "
                "coding-DNA c. equivalent (e.g. c.1521_1523del) for CDS-based "
                "allele-specific design."
            ),
        }

    m = HGVS_C_PATTERN.match(variant)
    if not m:
        return {
            "parsed": False,
            "reason": "Not a recognized c. notation pattern. Expected e.g. c.1521C>T or c.1521_1523del.",
        }

    start = int(m.group("start")) - 1  # HGVS is 1-based; convert to 0-based index
    end = int(m.group("end")) - 1 if m.group("end") else start

    if m.group("ref"):
        return {"parsed": True, "type": "substitution", "cdsStart": start, "cdsEnd": start, "length": 1}
    if m.group("del_seq") is not None and not m.group("delins_ins"):
        return {"parsed": True, "type": "deletion", "cdsStart": start, "cdsEnd": end, "length": end - start + 1}
    if m.group("dup_seq") is not None:
        return {"parsed": True, "type": "duplication", "cdsStart": start, "cdsEnd": end, "length": end - start + 1}
    if m.group("ins_seq"):
        return {"parsed": True, "type": "insertion", "cdsStart": start, "cdsEnd": start, "length": 0}
    if m.group("delins_ins"):
        return {"parsed": True, "type": "delins", "cdsStart": start, "cdsEnd": end, "length": end - start + 1}

    return {"parsed": False, "reason": "Unhandled variant type."}


# Wing length for gap-mer chemistries. "gapmer" = 2 from this file's own
# "DNA Gapmer (2-10-2)" label. For lna_gapmer the full wing convention is not
# separately stated; MODIFICATION_OPTIONS "lna_wings" says LNA wings are
# "typically used in gapmer wings (2-3 LNA at each end)", so the lower bound of
# that stated range is used. The resulting gap boundary for lna_gapmer is still
# an estimate and is surfaced as such in the candidate note.
GAP_WING_LENGTH = {
    "gapmer": 2,
    "lna_gapmer": 2,
}


def allele_discrimination_score(variant_relative_pos: int, aso_length: int, chemistry: str) -> dict:
    """Score allele discrimination from *where the mismatch falls relative to the
    RNase H-competent gap*, not a flat per-chemistry bonus.

    For gapmer-class chemistries only the central DNA gap recruits RNase H — the
    chemically-modified wings do not — so a mismatch in the wing contributes no
    discrimination. Within the gap, discrimination is strongest when the mismatch
    sits at the gap center (where RNase H cleaves) and drops off toward the gap
    edges.

    Steric-block chemistries (PMO, 2'-OMe) have no enzymatic cleavage, so a
    single mismatch does not block binding sharply; their discrimination signal
    is capped low rather than merely slightly penalized.
    """
    wing = GAP_WING_LENGTH.get(chemistry)

    if wing is not None:  # gapmer-class chemistry — RNase H cleavage-competent
        gap_start, gap_end = wing, aso_length - wing
        gap_len = gap_end - gap_start

        if not (gap_start <= variant_relative_pos < gap_end):
            return {
                "eligibleForAlleleSpecificity": False,
                "discriminationScore": 0.0,
                "note": (
                    "Variant falls in the chemically-modified wing, not the RNase "
                    "H-competent DNA gap — a mismatch here won't meaningfully "
                    "discriminate mutant from wild-type."
                ),
            }

        gap_center = gap_start + gap_len / 2
        dist = abs(variant_relative_pos - gap_center)
        proximity = 1 - (dist / (gap_len / 2))  # 1.0 at center, 0.0 at gap edge

        return {
            "eligibleForAlleleSpecificity": True,
            "discriminationScore": round(proximity, 3),
            "note": (
                f"Mismatch is {dist:.1f} nt from gap center (RNase H discrimination "
                "is strongest at the center, weakest at the edges)."
            ),
        }

    # Steric-block chemistry — no enzymatic cleavage, discrimination is inherently
    # weaker regardless of position. Capped low, not just "penalized a bit".
    centered = 1 - abs(variant_relative_pos - aso_length / 2) / (aso_length / 2)
    proximity = centered * 0.3
    return {
        "eligibleForAlleleSpecificity": True,
        "discriminationScore": round(proximity, 3),
        "note": (
            "Steric-block chemistries discriminate single mismatches far less "
            "reliably than gapmer/RNase H cleavage — treat as a weak signal, "
            "not a strong one."
        ),
    }


def _allele_specific_scoring(
    known_variant: str | None,
    spans_variant: bool,
    variant_relative_pos: int | None,
    aso_length: int,
    chemistry: str,
) -> dict:
    """Flag allele-specific candidates based on the variant's CDS coordinate.

    An ASO can only discriminate the mutant from the wild-type allele when its
    binding window actually spans the variant's CDS coordinate (computed by the
    caller). Within a spanning window, the discriminating variable is where the
    mismatch falls relative to the RNase H gap center
    (``allele_discrimination_score``). No flat per-chemistry bonuses are
    invented. The caller applies the ``spans_variant`` and wing-eligibility hard
    gates when the user has selected allele-specific scope.
    """
    result = {
        "alleleSpecific": False,
        "alleleNotes": "",
        "alleleDiscriminationScore": None,
        "alleleDiscriminationNote": None,
    }

    if not known_variant:
        return result

    if not spans_variant or variant_relative_pos is None:
        result["alleleNotes"] = (
            f"Does not span {known_variant} at the candidate binding window; "
            "reposition the ASO over the variant for allele discrimination."
        )
        return result

    discrimination = allele_discrimination_score(variant_relative_pos, aso_length, chemistry)

    if not discrimination["eligibleForAlleleSpecificity"]:
        result["alleleNotes"] = (
            f"Binding window spans {known_variant}, but {discrimination['note']}"
        )
        return result

    note = discrimination["note"]
    if chemistry == "lna_gapmer":
        note += (
            " (LNA wing convention estimated from the 2-3 LNA/wing guidance; "
            "not separately verified.)"
        )

    result["alleleSpecific"] = True
    result["alleleDiscriminationScore"] = discrimination["discriminationScore"]
    result["alleleDiscriminationNote"] = note
    result["alleleNotes"] = f"Binding window spans {known_variant} at its CDS coordinate."
    return result


def generate_candidates(
    target_exon_indices: list[int] | None,
    aso_length: int,
    chemistry: str,
    modifications: list[str],
    mrna_sequence: str | None,
    exons: list[dict],
    mechanism_id: str,
    delivery_context: str | None = None,
    defect_type: str | None = None,
    silencing_scope: str | None = None,
    known_variant: str | None = None,
) -> list[dict]:
    """Generate candidate ASOs for the selected mechanism and target region.

    Uses each exon's exact CDS-relative offset (``cdsStart``/``cdsEnd``, as
    computed by ``_map_exons_to_cds`` from the real cDNA-to-CDS alignment) so
    a candidate labeled "Exon 5" targets that exon. If that mapping is
    unavailable — the cDNA fetch failed upstream — falls back to a
    proportional genomic-length estimate and marks every candidate with
    ``exonMappingSource: "estimated_proportional"`` so callers can tell the
    two apart; this fallback should not be trusted for position-critical
    mechanisms (e.g. splice-junction-relative targeting). When no exon list
    is supplied, candidates are generated across all exons for
    total-transcript knockdown.
    """
    candidates = []
    if not mrna_sequence or len(exons) < 2:
        return candidates

    aso_length, chemistry, targeting_mode = _mechanism_design_constraints(
        mechanism_id, aso_length, chemistry
    )

    is_allele_specific = (silencing_scope or "").lower().strip() == "allele_specific"

    # Hard gate: reject chemistry/mechanism combinations that contradict the
    # mechanism's own rulebook (e.g. a steric-blocking PMO for an RNase H
    # mechanism). Runs once, before any candidate work.
    _mechanism_chemistry_compatibility(mechanism_id, chemistry)
    mechanism_note = _mechanism_note(mechanism_id, chemistry)

    seq = mrna_sequence.upper()
    seq_len = len(seq)

    # Prefer each exon's exact CDS-relative offset, computed upstream by
    # _map_exons_to_cds from the real cDNA-to-CDS alignment. Only fall back
    # to a proportional genomic-length estimate if that mapping is missing
    # (e.g. the cDNA fetch failed) — and mark candidates accordingly, since
    # a proportional estimate wrongly assumes UTR is spread evenly across
    # all exons rather than concentrated at the 5' and 3' ends.
    has_real_cds_map = all(
        e.get("cdsStart") is not None and e.get("cdsEnd") is not None for e in exons
    )

    exon_cds_map: list[tuple[int, int]] = []
    if has_real_cds_map:
        exon_mapping_source = "ensembl_cdna"
        exon_cds_map = [(e["cdsStart"], e["cdsEnd"]) for e in exons]
    else:
        exon_mapping_source = "estimated_proportional"
        total_genomic = sum(e.get("length", 0) for e in exons)
        if total_genomic == 0:
            return candidates

        cursor = 0
        remaining_seq = seq_len
        for i, exon in enumerate(exons):
            exon_genomic_len = exon.get("length", 0)
            remaining_exons = len(exons) - i
            # Proportional CDS contribution for this exon
            if i == len(exons) - 1:
                # Last exon gets all remaining sequence
                cds_contribution = remaining_seq
            else:
                cds_contribution = round(seq_len * exon_genomic_len / total_genomic)
                # Ensure each exon gets enough room for at least one ASO binding site
                min_contribution = min(aso_length * 2, remaining_seq // remaining_exons)
                cds_contribution = max(min_contribution, cds_contribution)
                cds_contribution = min(cds_contribution, remaining_seq - (remaining_exons - 1))
            cds_start = cursor
            cds_end = cursor + cds_contribution
            exon_cds_map.append((cds_start, cds_end))
            cursor = cds_end
            remaining_seq -= cds_contribution

    exon_count = len(exons)

    # A2 sterically blocks translation, so it must cover the 5′ initiation
    # region rather than scanning arbitrary exons. The other supported
    # mechanisms use selected exon(s), or the whole transcript when none is
    # specified.
    if targeting_mode == "translation_start":
        search_ranges = [(0, min(seq_len - aso_length, 90), "5′ translation-initiation region")]
    else:
        search_ranges = []

    # A missing list means total-transcript knockdown. Invalid exon numbers
    # are ignored rather than silently using an unrelated default region.
    # A selected exon must always remain a hard targeting constraint, including
    # for allele-specific designs.
    is_total_knockdown = target_exon_indices is None
    requested_exons = (
        list(range(1, exon_count + 1))
        if is_total_knockdown
        else (target_exon_indices or [])
    )
    target_indices = sorted({index for index in requested_exons if 0 < index <= exon_count})
    if not target_indices and targeting_mode != "translation_start":
        return candidates

    # Parse the known variant's CDS coordinate once so allele-specific
    # candidates are always generated around it, even when the chosen exons
    # do not contain the variant. When the variant cannot be parsed to a CDS
    # coordinate, no positional discrimination scoring is possible; the parse
    # failure is carried into the response so the UI can surface it instead of
    # silently producing non-discriminating results.
    variant_parse = parse_hgvs_c(known_variant) if known_variant else None
    variant_pos = None
    variant_coordinate = ""
    if (
        known_variant
        and targeting_mode != "translation_start"
        and variant_parse
        and variant_parse["parsed"]
        and variant_parse["cdsStart"] < seq_len
    ):
        variant_pos = variant_parse["cdsStart"]
        variant_coordinate = f"c.{variant_parse['cdsStart'] + 1}"

    seen = set()
    # Search only windows fully contained within each selected exon. Do not
    # borrow flanking sequence from a neighbouring exon: that can make two
    # different exon choices yield overlapping candidate sets.
    step = max(1, aso_length // 3)
    if targeting_mode != "translation_start":
        for target_exon_index in target_indices:
            exon_start, exon_end = exon_cds_map[target_exon_index - 1]
            search_ranges.append((
                exon_start,
                min(seq_len - aso_length, exon_end - aso_length),
                f"Exon {target_exon_index}",
            ))

    # Keep a variant-centered window inside the requested exon. Previously it
    # was appended for every exon selection, which caused every choice to show
    # the same top-ranked, variant-centered candidates.
    if variant_pos is not None:
        variant_exon_index = next(
            (index for index, (start, end) in enumerate(exon_cds_map, 1) if start <= variant_pos < end),
            None,
        )
        if is_total_knockdown or variant_exon_index in target_indices:
            variant_exon_start, variant_exon_end = exon_cds_map[variant_exon_index - 1] if variant_exon_index else (0, seq_len)
            search_ranges.append((
                max(variant_exon_start, variant_pos - aso_length + 1),
                min(variant_pos, variant_exon_end - aso_length),
                f"Variant {variant_coordinate}",
            ))

    for search_start, search_end, target_label in search_ranges:
        if search_end < search_start:
            continue

        exon_start = exon_end = None
        if targeting_mode != "translation_start" and target_label.startswith("Exon "):
            ei = int(target_label.removeprefix("Exon "))
            exon_start, exon_end = exon_cds_map[ei - 1]

        # In allele-specific mode the variant-centered window is scanned at
        # single-base resolution so a candidate can land the mismatch exactly
        # at the RNase H gap center. The window is at most aso_length bases
        # wide, so the fine step is cheap.
        range_step = 1 if (is_allele_specific and target_label.startswith("Variant ")) else step

        for offset in range(search_start, search_end + 1, range_step):
            candidate_seq = seq[offset : offset + aso_length]
            if len(candidate_seq) < aso_length or candidate_seq in seen:
                continue
            seen.add(candidate_seq)

            spans_variant = (
                variant_pos is not None
                and offset <= variant_pos < offset + aso_length
            )

            # Hard gate for allele-specific scope: a candidate whose binding
            # window does not cover the variant's CDS coordinate cannot
            # discriminate mutant from wild-type at all — exclude it rather
            # than merely ranking it lower.
            if is_allele_specific and variant_pos is not None and not spans_variant:
                continue

            gc = _calc_gc(candidate_seq)
            if (gc < MIN_GC or gc > MAX_GC) and not spans_variant:
                continue

            tm = _calc_tm(candidate_seq)
            self_mfe = _self_complement_mfe(candidate_seq)
            pg = _polyg_score(candidate_seq)
            cpg = _cpg_count(candidate_seq)

            # Target duplex energy — the antisense ASO bound to its target
            # window. The stored ASO is the reverse complement of the mRNA
            # window, so the duplex is computed between that ASO and the target.
            aso_seq = _reverse_complement(candidate_seq)
            duplex_energy = _target_duplex_energy(aso_seq, candidate_seq)

            # Mechanism note — informational text derived from the mechanism's
            # rulebook; no numeric bonus is voted into the ranking score.
            mech_notes = mechanism_note

            # Defect-type-specific notes — same policy: informative text only,
            # excluded from the composite ranking score.
            defect_notes = _defect_scores(defect_type, silencing_scope)["defect_notes"]

            if is_total_knockdown:
                region_label = f"Full Transcript offset +{offset}"
                exon_num = None
                exon_len = None
            elif targeting_mode == "translation_start":
                region_label = f"{target_label} offset +{offset}"
                exon_num = None
                exon_len = None
            else:
                # Determine exon number for this candidate
                if target_label.startswith("Exon "):
                    exon_num = int(target_label.removeprefix("Exon "))
                else:
                    exon_num = next(
                        (ei for ei, (es, ee) in enumerate(exon_cds_map, 1) if es <= offset < ee),
                        None,
                    )
                # If the candidate falls within an exon but exon_start is not set
                # (e.g. variant-based ranges), look up the exon coordinates.
                if exon_num is not None and exon_start is None:
                    exon_start, exon_end = exon_cds_map[exon_num - 1]
                # Compute region label
                if exon_start is None:
                    region_label = f"{target_label} offset +{offset}"
                else:
                    relative_pos = offset - exon_start
                    if relative_pos < 0:
                        region_label = f"{target_label} 5' flank {relative_pos}"
                    elif offset + aso_length > exon_end:
                        region_label = f"{target_label} 3' flank +{relative_pos}"
                    else:
                        region_label = f"{target_label} offset +{relative_pos}"
                exon_len = (exon_end - exon_start) if (exon_start is not None and exon_end is not None) else None

            # Compute additional drug-like properties
            nuclease_score = _nuclease_resistance_score(chemistry, modifications)
            uptake_score = _cellular_uptake_score(chemistry, aso_length)
            bbb_score = _bbb_crossing_score(chemistry, aso_length, modifications)
            synthesis_score = _synthesis_difficulty(candidate_seq, chemistry, modifications)
            complexity_val = _sequence_complexity(candidate_seq)
            off_target = _off_target_risk(candidate_seq, complexity_val)
            immune_risk = _immune_stimulation_risk(candidate_seq, chemistry)
            stability = _duplex_stability(gc, tm, aso_length)
            mw = _molecular_weight(candidate_seq)
            ext_coeff = _extinction_coefficient(candidate_seq)

            # Allele-specific scoring — only candidates whose binding window
            # actually spans the variant's CDS coordinate are scored, and the
            # score reflects *where* the mismatch falls relative to the RNase H
            # gap center (not a flat chemistry bonus).
            allele = _allele_specific_scoring(
                known_variant=known_variant,
                spans_variant=spans_variant,
                variant_relative_pos=(variant_pos - offset) if variant_pos is not None else None,
                aso_length=aso_length,
                chemistry=chemistry,
            )

            # Second hard gate for allele-specific scope: even a spanning window
            # contributes no discrimination when the mismatch falls in the
            # chemically-modified wing (which does not recruit RNase H). Such
            # candidates are excluded from allele-specific results rather than
            # softly penalized.
            if is_allele_specific and variant_pos is not None and not allele["alleleSpecific"]:
                continue

            # Composite ranking score — real metrics only: the ViennaRNA
            # target-duplex ΔG plus the chemistry-adjusted Tm fit.
            tm_fit = _tm_fit_score(tm, chemistry, modifications, mechanism_id)
            composite_score = _composite_score(duplex_energy, tm_fit)

            candidates.append({
                "sequence": aso_seq,
                "length": aso_length,
                "compositeScore": composite_score,  # 0-100 ranking score
                "learnedEfficacy": {
                    "available": False,
                    "value": None,
                    "modelInfo": "Not yet trained",
                    "scopeCaveat": None,
                },
                # Measured / computed properties — exact physics and sequence
                # computations, the tier that drives ranking.
                "realMetrics": {
                    "targetDuplexEnergy": duplex_energy,  # Real ViennaRNA duplex ΔG (kcal/mol)
                    "meltingTempC": tm,  # Real nearest-neighbor Tm (°C) via primer3
                    "selfStructureMfe": self_mfe,  # Real ViennaRNA MFE (kcal/mol)
                    "gcContent": round(gc * 100, 1),
                    "cpgCount": cpg,
                    "longestHomopolymer": _longest_homopolymer(candidate_seq),
                    "purineContent": _purine_content(candidate_seq),
                    "gcSkew": _gc_skew(candidate_seq),
                    "sequenceComplexity": complexity_val,
                    "polyGPass": pg == 0,
                    "molecularWeight": mw,
                    "extinctionCoefficient": ext_coeff,
                    "duplexStability": stability,
                },
                # Rule-of-thumb drug-like estimates — deliberately excluded
                # from ranking, shown to the user as labeled estimates.
                "heuristicEstimates": {
                    "nucleaseResistance": {
                        "value": nuclease_score,
                        "note": "Chemistry-class rule of thumb, not measured.",
                    },
                    "cellularUptake": {
                        "value": uptake_score,
                        "note": "Length/chemistry rule of thumb, not measured.",
                    },
                    "bbbCrossing": {
                        "value": bbb_score,
                        "note": "Length/chemistry rule of thumb, not measured.",
                    },
                    "synthesisDifficulty": {
                        "value": synthesis_score,
                        "note": "Sequence/chemistry rule of thumb, not measured.",
                    },
                    "offTargetRisk": {
                        "value": off_target,
                        "note": "Length/repetitiveness heuristic — not a genome alignment check.",
                    },
                    "immuneStimulation": {
                        "value": immune_risk,
                        "note": "CpG-count heuristic, not an immunogenicity assay.",
                    },
                },
                "targetRegion": region_label,
                "mechanismId": mechanism_id,
                "chemistry": chemistry,
                "modifications": modifications,
                "exonNumber": exon_num,
                "exonLength": exon_len,
                "exonMappingSource": exon_mapping_source,
                "deliveryContext": delivery_context or "",
                "defectType": defect_type or "",
                "silencingScope": silencing_scope or "",
                "defectNotes": defect_notes,
                "mechanismNotes": mech_notes,
                "knownVariant": known_variant or "",
                "alleleSpecific": allele["alleleSpecific"],
                "alleleNotes": allele["alleleNotes"],
                "alleleDiscriminationScore": allele["alleleDiscriminationScore"],
                "alleleDiscriminationNote": allele["alleleDiscriminationNote"],
            })

    # Rank by the composite design score (higher = better). In allele-specific
    # scope with a parsed variant, all surviving candidates span the variant, so
    # the primary axis is allele discrimination (mismatch proximity to the RNase
    # H gap center), with the composite score as the tiebreaker. Otherwise a
    # supplied variant ranks spanning candidates first.
    if is_allele_specific and variant_pos is not None:
        candidates.sort(
            key=lambda c: (
                -(c["alleleDiscriminationScore"] or 0),
                -c["compositeScore"],
                c["realMetrics"]["targetDuplexEnergy"],
            )
        )
    elif variant_pos is not None:
        candidates.sort(
            key=lambda c: (
                not c["alleleSpecific"],
                -c["compositeScore"],
                c["realMetrics"]["targetDuplexEnergy"],
            )
        )
    else:
        candidates.sort(
            key=lambda c: (-c["compositeScore"], c["realMetrics"]["targetDuplexEnergy"])
        )
    return candidates[:10]
