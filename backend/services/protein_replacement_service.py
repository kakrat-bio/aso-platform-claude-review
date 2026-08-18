"""
Protein Replacement (TG08) backend service.

Generates ranked RNA replacement construct candidates based on real
sequence data from Ensembl, with deterministic biophysical scoring.
"""

from __future__ import annotations

import logging
import os
from typing import Any

import requests
import RNA

from services.gene_silencing_service import get_target_analysis, _ensembl_get
from services import real_data_cache as RDC

logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"

# Human relative adaptiveness for codon optimization (simplified)
_HUMAN_CODON_ADAPT: dict[str, float] = {
    "UUU": 0.52, "UUC": 0.48, "UUA": 0.07, "UUG": 0.13,
    "CUU": 0.13, "CUC": 0.20, "CUA": 0.07, "CUG": 0.40,
    "AUU": 0.36, "AUC": 0.47, "AUA": 0.18, "AUG": 1.00,
    "GUU": 0.18, "GUC": 0.24, "GUA": 0.12, "GUG": 0.46,
    "UCU": 0.19, "UCC": 0.22, "UCA": 0.15, "UCG": 0.06,
    "CCU": 0.19, "CCC": 0.20, "CCA": 0.20, "CCG": 0.06,
    "ACU": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11,
    "GCU": 0.21, "GCC": 0.27, "GCA": 0.23, "GCG": 0.09,
    "UAU": 0.44, "UAC": 0.56, "UAA": 0.30, "UAG": 0.24,
    "CAU": 0.42, "CAC": 0.58, "CAA": 0.27, "CAG": 0.73,
    "AAU": 0.47, "AAC": 0.53, "AAA": 0.43, "AAG": 0.57,
    "GAU": 0.46, "GAC": 0.54, "GAA": 0.42, "GAG": 0.58,
    "UGU": 0.45, "UGC": 0.55, "UGA": 0.26, "UGG": 1.00,
    "CGU": 0.08, "CGC": 0.19, "CGA": 0.06, "CGG": 0.21,
    "AGU": 0.15, "AGC": 0.22, "AGA": 0.21, "AGG": 0.20,
    "GGU": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}


def _codon_to_aa(codon: str) -> str:
    """Return single-letter amino acid code for a codon."""
    table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    return table.get(codon.upper(), "X")


_AMINO_ACID_TO_CODONS: dict[str, list[str]] = {}
for _codon, _adapt in sorted(_HUMAN_CODON_ADAPT.items(), key=lambda x: -x[1]):
    if _codon in ("UAA", "UAG", "UGA"):
        continue
    _aa = _codon_to_aa(_codon)
    _AMINO_ACID_TO_CODONS.setdefault(_aa, []).append(_codon)


def _optimize_codon(seq: str, strategy: str) -> str:
    """Apply codon optimization strategy to a coding sequence."""
    seq = seq.upper().replace("T", "U")
    if len(seq) % 3 != 0:
        seq = seq[: len(seq) // 3 * 3]
    if not seq.startswith("AUG"):
        idx = seq.find("AUG")
        if idx != -1:
            seq = seq[idx:]
        else:
            return seq

    optimized = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in ("UAA", "UAG", "UGA"):
            optimized.append(codon)
            continue
        aa = _codon_to_aa(codon)
        if aa == "X":
            optimized.append(codon)
            continue
        candidates = _AMINO_ACID_TO_CODONS.get(aa, [codon])
        if not candidates:
            optimized.append(codon)
            continue

        if strategy == "cai":
            best = candidates[0]
        elif strategy == "mfe":
            best = min(candidates, key=lambda c: abs(_HUMAN_CODON_ADAPT.get(c, 0.5) - 0.5))
        elif strategy == "uridine":
            best = min(candidates, key=lambda c: c.count("U"))
        else:
            best = codon
        optimized.append(best)
    return "".join(optimized)


def _calc_cai(seq: str) -> float:
    """Calculate Codon Adaptation Index for a coding sequence."""
    seq = seq.upper().replace("T", "U")
    if len(seq) < 3:
        return 0.0
    adapt_sum = 0.0
    adapt_count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in ("UAA", "UAG", "UGA"):
            continue
        adapt_sum += _HUMAN_CODON_ADAPT.get(codon, 0.5)
        adapt_count += 1
    return round(adapt_sum / adapt_count, 3) if adapt_count > 0 else 0.0


def _calc_u_content(seq: str) -> float:
    """Return uridine percentage (0-100)."""
    seq = seq.upper().replace("T", "U")
    if not seq:
        return 0.0
    return round(seq.count("U") / len(seq) * 100, 1)


def _calc_mfe(seq: str) -> float:
    """Calculate minimum free energy using ViennaRNA."""
    try:
        _, energy = RNA.fold(seq)
        return round(energy, 1)
    except Exception:
        return 0.0


def _calc_gc_content(seq: str) -> float:
    seq = seq.upper().replace("T", "U")
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return round(gc / len(seq) * 100, 1)


def _estimate_tlr_risk(u_content: float, gc_content: float) -> str:
    """Estimate TLR innate immune risk from sequence composition."""
    if u_content > 25 or gc_content < 35:
        return "High"
    if u_content > 20 or gc_content < 40:
        return "Moderate"
    if u_content > 18:
        return "Low"
    return "Very Low"


def _estimate_initiation_efficiency(cai: float, u_content: float) -> int:
    """Estimate translation initiation efficiency from CAI and U%."""
    score = 0
    score += max(0, min(50, (cai - 0.7) / 0.3 * 50))
    score += max(0, min(50, (22 - u_content) / 10 * 50))
    return int(max(55, min(98, score)))


def _estimate_protein_yield(cai: float, initiation: int, modality: str) -> str:
    """Estimate predicted protein yield category."""
    base = (cai * 0.6 + initiation / 100 * 0.4)
    if modality == "circrna":
        base += 0.15
    elif modality == "sarna":
        base += 0.25
    mult = base * 4.0
    if mult >= 3.5:
        return f"{mult:.1f}x Ultra-High"
    if mult >= 2.5:
        return f"{mult:.1f}x High"
    return f"{mult:.1f}x Medium"


def _estimate_half_life(modality: str) -> str:
    if modality == "circrna":
        return ">120 hrs"
    if modality == "sarna":
        return "72–96 hrs"
    return "24–48 hrs"


def _build_features(modality: str, utr_pair: str, ires: str | None) -> list[dict]:
    """Build transcript feature map."""
    features = []
    if modality == "circrna":
        features.append({"name": "5' Backsplicing Junction", "start": 1, "end": 1, "type": "scarsplice"})
    else:
        features.append({"name": "5' Cap-1", "start": 1, "end": 1, "type": "cap"})
    if ires:
        features.append({"name": f"{ires.upper()} IRES", "start": 2, "end": 80, "type": "ires"})
    else:
        utr_label = {
            "globin": "5' UTR (β-Globin)",
            "c3": "5' UTR (C3/CYP2E1)",
            "synthetic": "5' UTR (Synthetic ML)",
        }.get(utr_pair, "5' UTR")
        features.append({"name": utr_label, "start": 2, "end": 82, "type": "utr"})
    features.append({"name": "Kozak Consensus", "start": 83, "end": 89, "type": "kozak"})
    features.append({"name": "Codon-Optimized ORF", "start": 90, "end": 4449, "type": "orf"})
    utr3_label = {
        "globin": "3' UTR (α-Globin)",
        "c3": "3' UTR (C3/CYP2E1)",
        "synthetic": "3' UTR (Synthetic ML)",
    }.get(utr_pair, "3' UTR")
    features.append({"name": utr3_label, "start": 4450, "end": 4531, "type": "utr3"})
    if modality == "circrna":
        features.append({"name": "Splicing Scar", "start": 4532, "end": 4532, "type": "scarsplice"})
    else:
        features.append({"name": "Poly(A) Tail (120 nt)", "start": 4532, "end": 4651, "type": "polyA"})
    return features


def _get_ensembl_cds(gene_symbol: str, organism: str = "homo_sapiens") -> tuple[str, str, str]:
    """Fetch CDS, protein length, and RefSeq from Ensembl."""
    try:
        meta = get_target_analysis(gene_symbol, gene_symbol=gene_symbol, organism=organism)
        cds = meta.get("mrnaSequence") or ""
        if cds:
            return cds, f"{len(cds) // 3} aa", meta.get("canonicalTranscript", {}).get("id", "")
    except Exception as exc:
        logger.warning("CDS fetch failed for %s: %s", gene_symbol, exc)

    try:
        resp = requests.get(
            f"{ENSEMBL_REST}/lookup/symbol/{organism}/{gene_symbol}?expand=1",
            headers={"Accept": "application/json"},
            timeout=20,
        )
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
                seq_resp = requests.get(
                    f"{ENSEMBL_REST}/sequence/id/{tid}?type=cds",
                    headers={"Accept": "application/json"},
                    timeout=20,
                )
                if seq_resp.ok:
                    seq = seq_resp.json().get("seq", "")
                    if seq:
                        return seq, f"{len(seq) // 3} aa", canonical.get("id", "")
    except Exception as exc:
        logger.warning("Fallback CDS fetch failed for %s: %s", gene_symbol, exc)
    return "", "N/A", "N/A"


# Published TLR-activation reduction factors for modified nucleosides.
# m1-pseudouridine and 5-methylcytidine+pseudouridine both suppress innate
# immune sensing of synthetic mRNA; these are the reduction factors the TG08
# spec supplies, applied multiplicatively to the base TLR scores.
#
# The base scores themselves are U-content heuristics, not assays, so the
# product is a heuristic too — it is reported under `diagnostics` alongside
# the modification that produced it rather than as a measurement.
MODIFICATION_EFFECTS: dict[str, dict[str, float]] = {
    "unmodified": {"tlr_reduction": 1.0, "stability_multiplier": 1.0,
                   "translation_boost": 1.0},
    "m1psi": {"tlr_reduction": 0.1, "stability_multiplier": 1.3,
              "translation_boost": 1.2},
    "5mc_psi": {"tlr_reduction": 0.15, "stability_multiplier": 1.25,
                "translation_boost": 1.15},
}


def _apply_modification_effects(tlr_scores: dict, modification: str) -> dict:
    """Scale TLR scores by the chosen nucleotide modification.

    m1-pseudouridine reduces TLR3/7/8 activation by roughly 90%;
    5-methylcytidine with pseudouridine by roughly 85%. An unrecognised
    modification is treated as unmodified rather than assumed beneficial.
    """
    effects = MODIFICATION_EFFECTS.get(modification, MODIFICATION_EFFECTS["unmodified"])
    return {k: round(v * effects["tlr_reduction"], 3) for k, v in tlr_scores.items()}


def _fold_sequence(seq: str) -> tuple[str, float]:
    """Real ViennaRNA fold. Returns (dot_bracket, mfe).

    Returns ("", 0.0) when ViennaRNA is unavailable or the sequence is empty,
    so callers can report the structure as not computed rather than draw a
    procedural pattern that looks like one.
    """
    if not seq:
        return "", 0.0
    try:
        import RNA
    except ImportError:
        return "", 0.0
    structure, mfe = RNA.fold(seq.upper().replace("T", "U"))
    return structure, float(mfe)


def _calc_amino_acid_identity(optimized_cds: str, native_cds: str) -> float | None:
    """Protein-level identity between the optimised and native CDS.

    Synonymous codon optimisation should give 100.0; anything less means the
    optimiser changed an amino acid, which is exactly what this field exists
    to catch. Returns None when there is no native CDS to compare against —
    the previous hardcoded 100.0 asserted a match that had not been checked.
    """
    if not optimized_cds or not native_cds:
        return None
    opt_aa = "".join(_codon_to_aa(optimized_cds[i:i + 3])
                     for i in range(0, len(optimized_cds) - 2, 3))
    nat_aa = "".join(_codon_to_aa(native_cds[i:i + 3])
                     for i in range(0, len(native_cds) - 2, 3))
    if not nat_aa:
        return None
    matches = sum(1 for a, b in zip(opt_aa, nat_aa) if a == b)
    return round(100.0 * matches / max(len(opt_aa), len(nat_aa)), 1)


def _evaluate_utr_structure(utr5_sequence: str) -> dict:
    """Does 5' UTR structure obstruct ribosome scanning?

    A stable hairpin in the 5' UTR impedes the scanning 43S complex. The
    thresholds below are the ones the TG08 spec specifies; they are de-novo
    cut points on a real computed MFE, not a calibrated model, and the
    returned `mfe` is the quantity to trust.

    Returns flag "NOT_COMPUTED" when there is no UTR or no ViennaRNA, rather
    than the previous unconditional "PASSED" — which asserted a clean scan
    for every construct including ones whose UTR was never examined.
    """
    if not utr5_sequence:
        return {"flag": "NOT_COMPUTED", "mfe": None, "hairpins": None,
                "structure": None,
                "note": "No 5' UTR sequence was available to fold."}
    structure, mfe = _fold_sequence(utr5_sequence)
    if not structure:
        return {"flag": "NOT_COMPUTED", "mfe": None, "hairpins": None,
                "structure": None,
                "note": "ViennaRNA unavailable; 5' UTR structure not computed."}
    if mfe > -15:
        flag = "PASSED"
    elif mfe > -25:
        flag = "CAUTION"
    else:
        flag = "BLOCKED"
    return {
        "flag": flag,
        "mfe": round(mfe, 2),
        "hairpins": structure.count("("),
        "structure": structure,
        "note": ("Thresholds -15 / -25 kcal/mol are de-novo cut points on a "
                 "real ViennaRNA MFE, not a calibrated model."),
    }


def _fetch_real_utrs(gene_symbol: str, organism: str = "homo_sapiens") -> dict | None:
    """Real 5' and 3' UTR sequences for the canonical transcript.

    Ensembl exposes them directly via /sequence/id/{transcript}?type=utr5 and
    type=utr3, so there is no reason to construct them. Returns None when the
    lookup cannot answer, which the caller turns into an explicit
    "unavailable" rather than a substitute.
    """
    lookup = _ensembl_get(
        f"{ENSEMBL_REST}/lookup/symbol/{organism}/{gene_symbol}?expand=1",
        timeout=10,
    )
    if not lookup or not getattr(lookup, "ok", False):
        return None
    data = lookup.json()
    transcript_id = data.get("canonical_transcript") or ""
    transcript_id = transcript_id.split(".")[0]
    if not transcript_id:
        for t in data.get("Transcript", []):
            if t.get("is_canonical") or t.get("biotype") == "protein_coding":
                transcript_id = t.get("id", "")
                break
    if not transcript_id:
        return None

    out: dict = {"transcript_id": transcript_id, "assembly": data.get("assembly_name")}
    for kind, key in (("utr5", "utr5"), ("utr3", "utr3")):
        resp = _ensembl_get(
            f"{ENSEMBL_REST}/sequence/id/{transcript_id}?type={kind}",
            timeout=10,
        )
        if resp is not None and getattr(resp, "ok", False):
            payload = resp.json()
            seq = ""
            if isinstance(payload, list) and payload:
                seq = payload[0].get("seq", "")
            elif isinstance(payload, dict):
                seq = payload.get("seq", "")
            if seq:
                out[key] = seq.upper()
    # A transcript with neither UTR annotated is not a usable answer.
    if not out.get("utr5") and not out.get("utr3"):
        return None
    return out


def _resolve_transcript_parts(symbol: str, organism: str) -> dict:
    """Real CDS + UTRs: live, else a real earlier fetch, else unavailable.

    This is the one place the TG08 pipeline is allowed to obtain sequence.
    There is deliberately no synthesis branch — the previous implementation
    fell back to `"AUG" + "GCU" * 300` for the CDS and
    `"GCCACC" + "A" * n` / `"G" * n` for the UTRs, then computed CAI, GC
    content, folding MFE, U-content and a protein-yield estimate from that
    padding. Those numbers described the padding, not the gene.
    """
    def fetch() -> dict | None:
        cds, native_length, ref_seq = _get_ensembl_cds(symbol, organism)
        if not cds:
            return None
        parts: dict = {
            "cds": cds.upper(),
            "nativeLength": native_length,
            "refSeq": ref_seq,
        }
        utrs = _fetch_real_utrs(symbol, organism)
        if utrs:
            parts.update(utrs)
        return parts

    return RDC.resolve(
        "transcript_parts", f"{organism}:{symbol}", fetch,
        source="Ensembl REST", source_version="",
    )


def generate_protein_replacement_candidates(
    target_symbol: str,
    rna_modality: str,
    codon_strategy: str,
    utr_pair: str,
    ires_selection: str | None,
    nucleotide_modification: str,
    organism: str = "homo_sapiens",
) -> dict[str, Any]:
    """Generate ranked protein replacement construct candidates."""
    symbol = target_symbol.strip().upper()
    if not symbol:
        raise ValueError("target_symbol is required")

    resolved = _resolve_transcript_parts(symbol, organism)
    if resolved["status"] == RDC.UNAVAILABLE:
        # Explicitly no candidates. Every downstream metric (CAI, GC, MFE,
        # U-content, yield) is a property of the sequence, so without a real
        # sequence there is nothing honest to report.
        return {
            "targetSymbol": symbol,
            "status": RDC.UNAVAILABLE,
            "dataProvenance": resolved,
            "candidates": [],
            "message": (
                f"No real coding sequence is available for {symbol}: the "
                f"Ensembl lookup did not answer and nothing has been cached "
                f"or curated for this target. Construct metrics are computed "
                f"from the sequence, so none can be reported. Retry when the "
                f"source is reachable, or add a verified entry to "
                f"data/reference/curated_transcript_parts.tsv."
            ),
        }

    parts = resolved["data"]
    cds = parts["cds"]
    native_length = parts.get("nativeLength", "N/A")
    ref_seq = parts.get("refSeq", "N/A")
    real_utr5 = (parts.get("utr5") or "").upper().replace("T", "U")
    real_utr3 = (parts.get("utr3") or "").upper().replace("T", "U")

    cds_clean = cds.upper().replace("T", "U")
    if cds_clean.startswith("AUG"):
        cds_clean = cds_clean[cds_clean.index("AUG"):]
    if len(cds_clean) % 3 != 0:
        cds_clean = cds_clean[: len(cds_clean) // 3 * 3]

    optimized_cds = _optimize_codon(cds_clean, codon_strategy)
    cai = _calc_cai(optimized_cds)

    modalities = []
    if rna_modality in ("circrna", "any"):
        modalities.append("circrna")
    if rna_modality in ("linear", "any"):
        modalities.append("linear")
    if rna_modality in ("sarna", "any"):
        modalities.append("sarna")
    if not modalities:
        modalities = ["linear"]

    candidates = []
    rank = 1
    for modality in modalities:
        variants = []
        if modality == "linear":
            variants = [
                {
                    "label": "Cap-1/Poly(A)",
                    "five_utr": 82,
                    "three_utr": 82,
                    "poly_a": 120,
                    "orf_offset": 90,
                    "topology": "Linear Cap-1/Poly(A)",
                }
            ]
        elif modality == "circrna":
            variants = [
                {
                    "label": "Covalently Closed Loop",
                    "five_utr": 80,
                    "three_utr": 1,
                    "poly_a": 0,
                    "orf_offset": 90,
                    "topology": "Covalently Closed Loop",
                }
            ]
        else:
            variants = [
                {
                    "label": "saRNA Replicon",
                    "five_utr": 612,
                    "three_utr": 82,
                    "poly_a": 0,
                    "orf_offset": 613,
                    "topology": "saRNA Replicon",
                }
            ]

        for var in variants:
            orf_length = len(optimized_cds)
            # Real UTRs when the source gave them. The `five_utr`/`three_utr`
            # numbers on each variant describe the LENGTH a given construct
            # architecture targets; they are not sequence, and padding them
            # out to that length would make every sequence-derived metric a
            # property of the padding.
            utr5 = real_utr5
            utr3 = real_utr3
            utr_source = "ensembl_real" if (real_utr5 or real_utr3) else "absent"
            # A poly(A) tail is added post-transcriptionally and is genuinely
            # a run of A of the stated length — that is the actual molecule,
            # not a stand-in for unknown sequence.
            poly_a = "A" * var["poly_a"]

            full_seq = utr5 + optimized_cds + utr3 + poly_a
            u_content = _calc_u_content(full_seq)
            gc_content = _calc_gc_content(full_seq)
            mfe = _calc_mfe(full_seq[: min(len(full_seq), 200)])
            tlr = _estimate_tlr_risk(u_content, gc_content)
            utr_structure = _evaluate_utr_structure(utr5)
            tlr_scores = _apply_modification_effects(
                {
                    "tlr3Score": round(max(1, min(30, u_content * 0.8)), 1),
                    "tlr7Score": round(max(1, min(25, u_content * 0.6)), 1),
                    "tlr8Score": round(max(1, min(25, u_content * 0.7)), 1),
                },
                nucleotide_modification,
            )
            initiation = _estimate_initiation_efficiency(cai, u_content)
            yield_label = _estimate_protein_yield(cai, initiation, modality)

            construct_id = (
                f"cRNA-{symbol}-v{rank}"
                if modality == "circrna"
                else f"saRNA-{symbol}-v{rank}"
                if modality == "sarna"
                else f"mRNA-{symbol}-v{rank}"
            )

            candidates.append({
                "rank": rank,
                "constructId": construct_id,
                "utrSource": utr_source,
                "utr5Length": len(utr5),
                "utr3Length": len(utr3),
                "modality": {
                    "circrna": "Circular RNA (circRNA)",
                    "linear": "Linear IVT mRNA",
                    "sarna": "Self-Amplifying RNA (saRNA)",
                }[modality],
                "vectorTopology": var["topology"],
                "cai": cai,
                "uContent": u_content,
                "mfe": mfe,
                "initiationEfficiency": initiation,
                "predictedProteinYield": yield_label,
                "tlrRisk": tlr,
                "signalPeptideStatus": "N-Terminal IgK Added" if modality == "circrna" else "Native Signal Peptide",
                "secondaryStructureFlag": utr_structure["flag"],
                "secondaryStructure": utr_structure,
                "fivePrimeUtrLength": var["five_utr"],
                "orfLength": orf_length,
                "threePrimeUtrLength": var["three_utr"],
                "polyATailLength": var["poly_a"],
                "gcContent": gc_content,
                "proteinLength": native_length,
                "molecularWeight": f"~{round(len(optimized_cds) * 0.11, 1)} kDa",
                "sequence": full_seq,
                "features": _build_features(modality, utr_pair, ires_selection),
                "diagnostics": {
                    # None when there is no native CDS to compare against.
                    # The previous hardcoded 100.0 asserted a match nobody
                    # had checked.
                    "aminoAcidIdentity": _calc_amino_acid_identity(
                        optimized_cds, cds_clean),
                    **tlr_scores,
                    "tlrModification": nucleotide_modification,
                    "tlrProvenance": (
                        "U-content heuristic scaled by published modification "
                        "reduction factors. Not an assay."
                    ),
                    # Real dot-bracket from ViennaRNA over the 5' UTR, or
                    # None. The previous value was a procedural pattern
                    # (i % 7) that rendered as a structure.
                    "mfePlot": utr_structure.get("structure"),
                    "mfe": utr_structure.get("mfe"),
                    "fiveUtrHairpin": (
                        None if utr_structure.get("hairpins") is None
                        else utr_structure["hairpins"] > 0
                    ),
                },
            })
            rank += 1

    candidates.sort(key=lambda c: c["initiationEfficiency"], reverse=True)
    for i, c in enumerate(candidates):
        c["rank"] = i + 1

    primary_mech = {
        "circrna": "A20 circRNA",
        "linear": "A19 Linear mRNA",
        "sarna": "A21 saRNA Replicon",
    }.get(modalities[0], "A19 Linear mRNA")

    feasibility = int(
        max(60, min(95, cai * 70 + (candidates[0]["initiationEfficiency"] if candidates else 70) * 0.2 + (100 - u_content) * 0.1))
    )

    overview = {
        "targetGene": symbol,
        "refSeq": ref_seq or "NM_000000.1",
        "nativeLength": native_length,
        "vectorTopology": candidates[0]["vectorTopology"] if candidates else "Linear Cap-1/Poly(A)",
        "cai": cai,
        "uContent": u_content,
        "primaryMechanism": primary_mech,
        "feasibilityScore": feasibility,
        "predictedHalfLife": _estimate_half_life(modalities[0]),
    }

    return {"overview": overview, "candidates": candidates}


def get_protein_replacement_options() -> dict[str, Any]:
    """Return design option catalogs for the protein replacement form."""
    return {
        "rnaModalities": [
            {"id": "linear", "label": "Linear IVT mRNA", "description": "Gold Standard / Rapid Kinetics — Cap-1/Poly(A) architecture for transient expression."},
            {"id": "circrna", "label": "Circular RNA (circRNA)", "description": "Extended Intracellular Half-Life — Covalently closed loop eliminates exonuclease degradation."},
            {"id": "sarna", "label": "Self-Amplifying RNA (saRNA)", "description": "Dose-Sparing High Yield — Alphavirus replicon enables high protein yields at microgram doses."},
            {"id": "any", "label": "All Modalities (Comparative)", "description": "Generate candidates across all RNA architectures for comparison."},
        ],
        "codonStrategies": [
            {"id": "cai", "label": "Human Codon Adaptation Index (CAI) Maxima"},
            {"id": "mfe", "label": "Minimum Free Energy (MFE) Secondary Structure Optimization"},
            {"id": "uridine", "label": "Uridine Depletion / Rare Codon Removal"},
        ],
        "utrPairs": [
            {"id": "globin", "label": "Human β-Globin / α-Globin UTRs"},
            {"id": "c3", "label": "Complement Factor 3 (C3) / CYP2E1 UTR"},
            {"id": "synthetic", "label": "Synthetic High-Yield UTR Pair (Machine-Learning Engineered)"},
        ],
        "iresSelections": [
            {"id": "cvb3", "label": "CVB3 IRES"},
            {"id": "ev71", "label": "EV71 IRES"},
            {"id": "m6a", "label": "m6A-Mediated Translation Ring"},
        ],
        "nucleotideModifications": [
            {"id": "m1psi", "label": "100% N1-Methylpseudouridine (m1Ψ)"},
            {"id": "5mc_psi", "label": "5-Methylcytidine (5mC) / Pseudouridine (Ψ)"},
            {"id": "unmodified", "label": "Unmodified (for circRNA / low immunogenicity backbones)"},
        ],
    }
