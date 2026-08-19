"""
Gene structural feature analysis for TG02 mechanism filtering.

Queries Ensembl to determine whether a gene has the structural prerequisites
for each upregulation mechanism:
- A3 (TANGO): needs poison exons / non-productive splice variants
- A4 (NAT): needs overlapping natural antisense transcripts
- A5 (uORF): needs upstream open reading frames in 5' UTR — detected by
  fetching spliced cDNA sequences and scanning each 5' UTR for ATG start
  codons that close an in-frame stop codon before the main CDS
- A6 (miRNA site block): needs 3' UTR with miRNA binding sites (always available)
- A23 (promoter activation): always available (all genes have promoters)
"""

from __future__ import annotations

import json
import logging
import re
import threading
import time
from typing import Optional

import requests

from database.db import SessionLocal
from database.models import GeneFeatureBackup

ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_TIMEOUT = 20
# Ensembl answers most lookups in ~1 s; transient timeouts and 429s are
# common enough that a single attempt loses real data.
ENSEMBL_MAX_RETRIES = 3

# When the Ensembl REST site fails, remember when. During the cooldown window
# we skip live queries entirely and serve the stored backup / fallback, so a
# downed site never blocks Gene Function analysis for any gene.
_ENSEMBL_DOWN_SINCE: float | None = None
_ENSEMBL_DOWN_LOCK = threading.Lock()
ENSEMBL_DOWN_COOLDOWN_SECONDS = 120

logger = logging.getLogger(__name__)

# Species name mapping for Ensembl
_SPECIES_MAP = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "rattus_norvegicus": "rat",
    "danio_rerio": "zebrafish",
    "drosophila_melanogaster": "fruitfly",
    "caenorhabditis_elegans": "celegans",
    "saccharomyces_cerevisiae": "yeast",
    "bos_taurus": "cow",
    "sus_scrofa": "pig",
    "gallus_gallus": "chicken",
    "canis_lupus_familiaris": "dog",
    "felis_catus": "cat",
}


def _ensembl_request(
    method: str,
    path: str,
    params: dict | None = None,
    payload: dict | None = None,
) -> dict | list | None:
    """Call the Ensembl REST API with timeout and outage tracking.

    Records failures so a downed site is not hammered on every gene; see
    ``_ensembl_available``.
    """
    global _ENSEMBL_DOWN_SINCE

    # THIS CLIENT HAD NO RETRY, AND THAT SILENTLY BROKE TWO MECHANISMS.
    #
    # A single read timeout returned None, which `analyze_gene_features` turns
    # into `verified: False`, which `feature_service._from_annotation` declines
    # to read as evidence — so F4 (poison exon) and F6 (natural antisense
    # transcript) stayed UNRESOLVED. Both are REQUIRED features, so A3 and A4
    # halted on every real target and only scored when the user hand-asserted
    # the matching molecular defect. Ensembl answers this lookup in about 0.9 s
    # when it is not rate-limiting, so the failures were transient and a single
    # retry recovers almost all of them.
    #
    # `gene_silencing_service._ensembl_get` already retries with backoff; this
    # client simply never got the same treatment.
    last_error: Exception | str | None = None
    for attempt in range(1, ENSEMBL_MAX_RETRIES + 1):
        try:
            kwargs = {
                "headers": {"Content-Type": "application/json"},
                "timeout": ENSEMBL_TIMEOUT,
            }
            if params:
                kwargs["params"] = params
            if payload:
                kwargs["json"] = payload
            resp = requests.request(method, f"{ENSEMBL_REST}{path}", **kwargs)
            if resp.status_code == 200:
                with _ENSEMBL_DOWN_LOCK:
                    _ENSEMBL_DOWN_SINCE = None
                return resp.json()
            # 429 is rate limiting, not an outage: back off and try again.
            if resp.status_code == 429:
                last_error = "HTTP 429 (rate limited)"
                retry_after = resp.headers.get("Retry-After")
                delay = float(retry_after) if retry_after and retry_after.isdigit() else attempt * 1.5
                if attempt < ENSEMBL_MAX_RETRIES:
                    time.sleep(min(delay, 5.0))
                    continue
            elif resp.status_code >= 500:
                last_error = f"HTTP {resp.status_code}"
                if attempt < ENSEMBL_MAX_RETRIES:
                    time.sleep(attempt * 1.5)
                    continue
                # Only a persistent 5xx marks the site down.
                with _ENSEMBL_DOWN_LOCK:
                    if _ENSEMBL_DOWN_SINCE is None:
                        _ENSEMBL_DOWN_SINCE = time.time()
            else:
                last_error = f"HTTP {resp.status_code}"
            break
        except (requests.RequestException, ValueError) as e:
            last_error = e
            if attempt < ENSEMBL_MAX_RETRIES:
                logger.info("Ensembl %s %s attempt %d/%d failed (%s), retrying",
                            method, path, attempt, ENSEMBL_MAX_RETRIES, e)
                time.sleep(attempt * 1.5)
                continue
            # A timeout on one request is not evidence the site is down, and
            # marking it down degrades every subsequent gene to the permissive
            # fallback. Only connection-level failures do that.
            if isinstance(e, requests.ConnectionError):
                with _ENSEMBL_DOWN_LOCK:
                    if _ENSEMBL_DOWN_SINCE is None:
                        _ENSEMBL_DOWN_SINCE = time.time()
    logger.warning("Ensembl %s %s failed after %d attempts: %s",
                   method, path, ENSEMBL_MAX_RETRIES, last_error)
    return None


def _ensembl_get(path: str, params: dict | None = None) -> dict | list | None:
    """GET from Ensembl REST API."""
    return _ensembl_request("GET", path, params=params)


def _ensembl_post(path: str, payload: dict | None = None) -> dict | list | None:
    """POST to Ensembl REST API (used for batch sequence lookup)."""
    return _ensembl_request("POST", path, payload=payload)


def _ensembl_available() -> bool:
    """True if Ensembl is reachable (or its outage cooldown has expired)."""
    with _ENSEMBL_DOWN_LOCK:
        down_since = _ENSEMBL_DOWN_SINCE
    if down_since is None:
        return True
    return (time.time() - down_since) >= ENSEMBL_DOWN_COOLDOWN_SECONDS


def _get_transcripts(gene_id: str) -> list[dict]:
    """Fetch all transcripts for a gene from Ensembl.

    Uses the expand=1 gene lookup, which returns the transcript list under
    the (singular) "Transcript" key — each transcript carries its own
    "Exon" array. The /overlap/translation endpoint only accepts
    translation IDs (ENSP...), not gene IDs, so it cannot be used here.
    """
    data = _ensembl_get(f"/lookup/id/{gene_id}", {"expand": "1"})
    if isinstance(data, dict):
        return data.get("Transcript", []) or []
    return []


def _get_regulatory_features(region: str) -> list[dict]:
    """Fetch regulatory features (promoters, enhancers) for a genomic region."""
    data = _ensembl_get(f"/overlap/region/human/{region}", {"feature": "regulatory"})
    if isinstance(data, list):
        return data
    return []


def _check_overlapping_nats(
    gene_id: str, species: str, gene_data: dict | None = None
) -> dict:
    """
    Check for overlapping natural antisense transcripts.
    Uses Ensembl overlap API to find antisense lncRNAs.

    hasOverlappingNat is None when the gene cannot be resolved (unknown
    symbol, missing coordinates, or Ensembl unavailable) — callers treat
    that as "unverified" rather than a definitive negative, so NAT
    silencing stays available for genes that can't be checked.
    """
    result = {
        "hasOverlappingNat": False,
        "natCount": 0,
        "natGenes": [],
    }

    if not gene_id:
        result["hasOverlappingNat"] = None
        return result

    if gene_data is None:
        gene_data = _ensembl_get(f"/lookup/id/{gene_id}")
    if not gene_data or not isinstance(gene_data, dict):
        result["hasOverlappingNat"] = None
        return result

    chrom = gene_data.get("seq_region_name")
    start = gene_data.get("start")
    end = gene_data.get("end")
    strand = gene_data.get("strand")

    if not all([chrom, start, end, strand]):
        result["hasOverlappingNat"] = None
        return result

    # Query overlapping features in the region. Ensembl accepts the raw
    # species slug (e.g. "homo_sapiens") in overlap/region URLs, so genes in
    # species outside _SPECIES_MAP still resolve instead of silently being
    # treated as human.
    region = f"{chrom}:{start - 50000}-{end + 50000}"
    species_name = _SPECIES_MAP.get(species, species)

    overlap_data = _ensembl_get(
        f"/overlap/region/{species_name}/{region}",
        {"feature": "gene", "content_type": "application/json"},
    )

    if not isinstance(overlap_data, list):
        return result

    for feature in overlap_data:
        feature_type = feature.get("feature_type", "")
        feat_strand = feature.get("strand")
        feat_id = feature.get("id", "")
        feat_desc = feature.get("description", "")

        # Detect antisense: same region, opposite strand, lncRNA biotype
        if feat_strand and feat_strand != strand:
            biotype = (feature.get("biotype") or "").lower()
            if "antisense" in biotype or "lncrna" in biotype or "ncrna" in biotype:
                result["hasOverlappingNat"] = True
                result["natCount"] += 1
                result["natGenes"].append({
                    "id": feat_id,
                    "description": feat_desc or feat_id,
                })

    return result


def _scan_uorfs(utr5: str) -> list[dict]:
    """Scan a 5' UTR sequence (5'→3') for canonical upstream open reading frames.

    A uORF is an ATG start codon followed by an in-frame stop codon
    (TAA/TAG/TGA) that terminates before the main CDS — i.e. entirely inside
    the 5' UTR. Coordinates are 1-based on the transcript sequence.
    """
    uorfs = []
    pos = utr5.find("ATG")
    while pos != -1:
        for stop in range(pos + 3, len(utr5) - 2, 3):
            if utr5[stop : stop + 3] in ("TAA", "TAG", "TGA"):
                uorfs.append({
                    "start": pos + 1,
                    "end": stop + 3,
                    "length": stop + 3 - pos,
                })
                break
        pos = utr5.find("ATG", pos + 1)
    return uorfs


def _fetch_transcript_sequences(transcript_ids: list[str]) -> dict[str, str]:
    """Batch-fetch spliced cDNA sequences for transcripts (id -> sequence).

    Uses Ensembl's POST /sequence/id (up to 50 ids per call). A missing or
    failed batch just yields fewer sequences; callers decide how to treat
    partial data.
    """
    if not transcript_ids:
        return {}
    sequences: dict[str, str] = {}
    for i in range(0, len(transcript_ids), 50):
        chunk = transcript_ids[i : i + 50]
        data = _ensembl_post("/sequence/id", {"ids": chunk, "type": "cdna"})
        if isinstance(data, list):
            for item in data:
                if isinstance(item, dict) and item.get("id") and item.get("seq"):
                    sequences[item["id"]] = item["seq"]
    return sequences


def _detect_uorfs(transcripts: list[dict]) -> dict:
    """
    Detect real upstream ORFs by scanning 5' UTR sequences from Ensembl.

    For every transcript with a CDS annotation (Translation.start gives the
    1-based position of the main start codon on the spliced cDNA), the 5' UTR
    is the cDNA prefix before that codon. Each 5' UTR is scanned for ATG
    start codons that close an in-frame stop codon before the CDS begins.

    Returns:
      hasUorfPotential: True when any protein-coding transcript has a uORF,
        False when transcripts were examined and none had one, and None when
        the 5' UTR sequences could not be fetched (unverifiable).
    """
    coding = [
        t
        for t in transcripts
        if isinstance(t, dict)
        and isinstance(t.get("Translation"), dict)
        and t["Translation"].get("start")
        and t.get("id")
    ]

    if not coding:
        return {
            "hasUorfPotential": None if not transcripts else False,
            "uorfCount": 0,
            "longestUtr5": 0,
            "uorfs": [],
            "transcriptCount": len(transcripts),
        }

    sequences = _fetch_transcript_sequences([t["id"] for t in coding])

    # No 5' UTR sequence fetched at all — cannot verify.
    if not sequences:
        return {
            "hasUorfPotential": None,
            "uorfCount": 0,
            "longestUtr5": 0,
            "uorfs": [],
            "transcriptCount": len(transcripts),
        }

    uorfs = []
    longest_utr5 = 0
    for t in coding:
        seq = sequences.get(t["id"])
        if not seq:
            continue
        cds_start = int(t["Translation"]["start"])
        utr5 = seq[: cds_start - 1]
        longest_utr5 = max(longest_utr5, len(utr5))
        for u in _scan_uorfs(utr5):
            u["transcript"] = t["id"]
            uorfs.append(u)

    return {
        "hasUorfPotential": bool(uorfs),
        "uorfCount": len(uorfs),
        "longestUtr5": longest_utr5,
        "uorfs": uorfs[:5],
        "transcriptCount": len(transcripts),
    }


def _max_exon_count(transcripts: list[dict]) -> int | None:
    """Largest exon count across a gene's transcripts (from expand=1 lookup)."""
    counts = [
        len(t.get("Exon") or [])
        for t in transcripts
        if isinstance(t, dict) and t.get("Exon")
    ]
    return max(counts) if counts else None


def _check_splicing_complexity(
    transcripts: list[dict],
    exon_count: int | None = None,
    total_transcripts: int | None = None,
) -> dict:
    """
    Determine whether a gene has sufficient splicing complexity for TANGO.
    - Needs multiple transcripts (evidence of alternative splicing)
    - Needs >1 exon (can't have poison exons in single-exon genes)

    exon_count / total_transcripts fall back to pipeline-computed values
    when the live transcript fetch here is empty.
    """
    if total_transcripts is None:
        total_transcripts = len(transcripts)
    if exon_count is None:
        exon_count = _max_exon_count(transcripts)

    has_alt_splicing = total_transcripts > 1
    has_introns = (exon_count or 0) > 1

    # Look for transcripts with "retained_intron" or "nonsense_mediated_decay" biotype
    has_nmd_variants = False
    for t in transcripts:
        biotype = (t.get("biotype") or "").lower()
        if "nonsense" in biotype or "nmd" in biotype or "retained_intron" in biotype:
            has_nmd_variants = True
            break

    return {
        "hasPoisonExonPotential": has_alt_splicing and has_introns,
        "transcriptCount": total_transcripts,
        "hasNmdTranscripts": has_nmd_variants,
        "hasIntrons": has_introns,
        "exonCount": exon_count,
    }


def _save_backup(organism: str, gene_symbol: str, ensembl_id: str | None, result: dict) -> None:
    """Persist a last-known-good analysis so it survives an Ensembl outage."""
    try:
        db = SessionLocal()
        try:
            row = (
                db.query(GeneFeatureBackup)
                .filter(
                    GeneFeatureBackup.organism == organism,
                    GeneFeatureBackup.gene_symbol == gene_symbol,
                )
                .one_or_none()
            )
            now = time.time()
            if row is None:
                row = GeneFeatureBackup(
                    organism=organism,
                    gene_symbol=gene_symbol,
                    ensembl_id=ensembl_id or "",
                    result=json.dumps(result),
                    created_at=now,
                    updated_at=now,
                )
                db.add(row)
            else:
                row.ensembl_id = ensembl_id or row.ensembl_id
                row.result = json.dumps(result)
                row.updated_at = now
            db.commit()
        finally:
            db.close()
    except Exception as e:  # never let a backup write break the analysis
        logger.warning("Failed to save gene feature backup for %s: %s", gene_symbol, e)


def _load_backup(organism: str, gene_symbol: str) -> dict | None:
    """Return the last-known-good analysis for a gene, or None."""
    try:
        db = SessionLocal()
        try:
            row = (
                db.query(GeneFeatureBackup)
                .filter(
                    GeneFeatureBackup.organism == organism,
                    GeneFeatureBackup.gene_symbol == gene_symbol,
                )
                .one_or_none()
            )
            if row is None:
                return None
            result = json.loads(row.result or "{}")
            if not isinstance(result, dict):
                return None
            result["source"] = "backup"
            result["backupTimestamp"] = row.updated_at
            return result
        finally:
            db.close()
    except Exception as e:
        logger.warning("Failed to load gene feature backup for %s: %s", gene_symbol, e)
        return None


def analyze_gene_features(
    gene_symbol: str,
    organism: str = "homo_sapiens",
    ensembl_id: str | None = None,
    tissue_tpm: float | None = None,
    exon_count: int | None = None,
    total_transcripts: int | None = None,
    gene_type: str | None = None,
) -> dict:
    """
    Analyze a gene's structural features to determine TG02 mechanism availability.

    Designed to work for every gene:
      - Structural hints already computed by the main gene pipeline
        (exon_count / total_transcripts) are used first; Ensembl is only
        queried to fill in the gaps.
      - If the gene cannot be resolved or verified (unknown symbol,
        non-Ensembl species, or Ensembl unavailable), structure-dependent
        mechanisms are NOT hard-excluded. They are reported as available
        with an honest "could not verify" note, so the TG02 ranking still
        returns candidates for all genes instead of silently dropping them.
      - Resilience backup: every live Ensembl analysis is persisted
        (database.gene_feature_backups). When the Ensembl site is down or a
        gene cannot be resolved, the last-known-good analysis is replayed
        from backup so Gene Function keeps working for every gene.

    Returns a dict with:
    - features: per-mechanism availability flags
    - warnings: tissue expression / toxicity warnings
    - geneInfo: basic gene metadata used for the analysis
    - source: "live" (Ensembl), "backup" (replayed from storage), or
      "fallback" (permissive heuristics) — plus backupTimestamp for "backup".
    """
    # Resolve Ensembl ID if not provided — skipped entirely while the Ensembl
    # site is known to be down, so every gene still gets an answer quickly.
    if ensembl_id or _ensembl_available():
        if not ensembl_id:
            lookup = _ensembl_get(
                f"/lookup/symbol/{organism}/{gene_symbol}",
                {"expand": "0"},
            )
            if isinstance(lookup, dict) and lookup.get("id"):
                ensembl_id = lookup["id"]

        transcripts = _get_transcripts(ensembl_id) if ensembl_id else []
        gene_data = _ensembl_get(f"/lookup/id/{ensembl_id}") if ensembl_id else None
    else:
        transcripts = []
        gene_data = None

    # Prefer pipeline-computed structural hints; fill gaps from Ensembl.
    if exon_count is None:
        exon_count = _max_exon_count(transcripts)
    if total_transcripts is None:
        total_transcripts = len(transcripts)

    # We can make a real structural determination when we have transcript
    # evidence or the pipeline's structural counts.
    can_verify_structure = bool(transcripts) or exon_count is not None or (total_transcripts or 0) > 0
    gene_verified = bool(gene_data) or bool(transcripts)

    if can_verify_structure:
        splicing = _check_splicing_complexity(transcripts, exon_count, total_transcripts)
        uorf = _detect_uorfs(transcripts)
    else:
        splicing = {
            "hasPoisonExonPotential": None,
            "transcriptCount": 0,
            "hasNmdTranscripts": False,
            "hasIntrons": None,
            "exonCount": None,
        }
        uorf = {
            "hasUorfPotential": None,
            "uorfCount": 0,
            "longestUtr5": 0,
            "uorfs": [],
            "transcriptCount": 0,
        }

    nats = _check_overlapping_nats(ensembl_id, organism, gene_data)

    unverified_reason = (
        "Could not verify gene structure from Ensembl — treated as potentially "
        "applicable for this gene; requires experimental validation."
    )

    # Build feature availability map. A tri-state value (True/False/None) keeps
    # "verified absent" distinct from "could not verify".
    features = {
        # sourceWired=False marks the three entries below as placeholders that
        # return True for every gene. "All protein-coding genes have promoters"
        # is true and useless: it is not evidence that THIS promoter is
        # silenced, that THIS 3' UTR carries a repressive miRNA site, or that
        # THIS transcript carries a repressive RBP site. The real sources
        # (methylation atlas / TargetScan context++ / a CLIP-derived binding
        # atlas) are F8, F7 and F11 in the feature plan and none is wired, so
        # feature_service leaves them unresolved rather than reading these as
        # positive findings.
        "saRNA": {
            "available": True,
            "verified": False,
            "sourceWired": False,
            "reason": "All protein-coding genes have promoter regions that can be targeted by saRNA",
        },
        "uORF": _feature_entry(
            uorf["hasUorfPotential"],
            (
                f"Detected {uorf['uorfCount']} upstream open reading frame(s) in the 5' UTR"
                + (
                    f" (e.g. {uorf['uorfs'][0]['transcript']} at nt {uorf['uorfs'][0]['start']})"
                    if uorf["uorfs"]
                    else ""
                )
            ),
            "No upstream open reading frames found in the 5' UTR of protein-coding transcripts",
            unverified_reason,
        ),
        "TANGO": _feature_entry(
            splicing["hasPoisonExonPotential"],
            (
                f"Gene has {splicing['transcriptCount']} transcripts with introns — "
                + (
                    "including NMD-associated variants"
                    if splicing["hasNmdTranscripts"]
                    else "alternative splicing may produce poison exons"
                )
            ),
            "Single-exon gene or insufficient splicing complexity for poison exon targeting",
            unverified_reason,
        ),
        "NAT": _feature_entry(
            nats["hasOverlappingNat"],
            (
                f"Found {nats['natCount']} overlapping antisense transcript(s)"
                + (
                    f": {', '.join(g['description'][:60] for g in nats['natGenes'][:3])}"
                    if nats["natGenes"]
                    else ""
                )
            ),
            "No overlapping natural antisense transcripts detected in genomic databases",
            "Could not verify overlapping antisense transcripts from Ensembl — NAT silencing treated as potentially applicable; requires experimental validation",
        ),
        "miRNA_block": {
            "available": True,
            "verified": False,
            "sourceWired": False,
            "reason": "Most protein-coding mRNAs contain miRNA binding sites in their 3' UTR",
        },
        "RBP_block": {
            "available": True,
            "verified": False,
            "sourceWired": False,
            "reason": "Protein-coding transcripts contain RNA-binding protein (RBP) binding sites that can be masked to relieve translational repression",
        },
    }

    # Tissue expression warnings
    warnings = []
    if tissue_tpm is not None:
        if tissue_tpm > 500:
            warnings.append({
                "type": "overexpression_risk",
                "severity": "high",
                "message": (
                    f"High endogenous expression ({tissue_tpm:.0f} TPM) in target tissue — "
                    "exercise caution against overexpression toxicity. "
                    "Consider whether upregulation is appropriate for this tissue."
                ),
            })
        elif tissue_tpm > 200:
            warnings.append({
                "type": "overexpression_caution",
                "severity": "medium",
                "message": (
                    f"Moderate-high endogenous expression ({tissue_tpm:.0f} TPM) in target tissue — "
                    "monitor for potential overexpression effects."
                ),
            })

    result = {
        "features": features,
        "warnings": warnings,
        "geneInfo": {
            "ensemblId": ensembl_id,
            "transcriptCount": total_transcripts or 0,
            "exonCount": splicing["exonCount"],
            "hasIntrons": splicing["hasIntrons"],
            "hasNmdTranscripts": splicing["hasNmdTranscripts"],
            "overlappingNats": nats["natCount"] or 0,
            "uorfCount": uorf["uorfCount"],
            "verified": gene_verified,
            "geneType": gene_type,
        },
    }

    # Resilience / backup: when live Ensembl data was obtained, persist it as
    # the last-known-good analysis. When the site is down (or the gene cannot
    # be resolved), replay the stored backup so Gene Function still works for
    # every gene instead of silently dropping it. If no backup exists yet,
    # return the permissive fallback with an honest note.
    resolved_live = bool(gene_data) or bool(transcripts)
    if resolved_live:
        result["source"] = "live"
        _save_backup(organism, gene_symbol, ensembl_id, result)
    else:
        backup = _load_backup(organism, gene_symbol)
        if backup is not None:
            # Keep the analysis but refresh warnings for this request's tissue.
            backup["warnings"] = warnings
            result = backup
        else:
            result["source"] = "fallback"

    return result


def _feature_entry(
    available: bool | None,
    reason_yes: str,
    reason_no: str,
    reason_unknown: str,
) -> dict:
    """Map a tri-state availability value to a feature dict.

    None means the structural check could not be run — treat the mechanism
    as available (so it is not silently dropped from the ranking) with an
    honest note rather than a definitive negative.

    ``verified`` preserves the tri-state that ``available`` flattens. The
    permissive True-on-unknown behaviour is deliberate for the ranking UI,
    but it means ``available`` alone cannot distinguish "we checked and it is
    there" from "we could not check". feature_service reads ``verified`` to
    decide whether a feature resolves at the annotation provenance tier or
    stays unresolved, so an unverifiable gene halts a mechanism instead of
    scoring it on an assumption.
    """
    if available is None:
        return {"available": True, "verified": False, "reason": reason_unknown}
    if available:
        return {"available": True, "verified": True, "reason": reason_yes}
    return {"available": False, "verified": True, "reason": reason_no}
