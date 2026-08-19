"""
Gene Upregulation design pipeline backend service.

Generates ASO candidates for TG02 mechanisms (Gene Activation /
Upregulation). Mechanism IDs follow the rulebooks (backend/rulebooks/A*/rule.json):

- A3   TANGO (poison exon skipping) — splice-junction ASOs for NMD suppression
- A4   NAT silencing               — RNase H1 gapmers targeting antisense transcripts
- A5   uORF blocking               — steric-blocking ASOs at 5' UTR / start codon
- A6   miRNA site blocking         — steric-blocking ASOs masking the seed site
- A23  saRNA                       — promoter-targeted 21-mer dsRNA duplexes
- A28  RBP masking                 — steric-blocking ASOs at RBP recognition sites

Reuses biophysical scoring helpers from ``gene_silencing_service``.
"""

from __future__ import annotations

import logging
import math
from collections import Counter

from services.gene_silencing_service import (
    get_target_analysis,
    _ensembl_get,
    ENSEMBL_REST,
    _calc_gc,
    _calc_tm,
    _self_complement_mfe,
    _polyg_score,
    _cpg_count,
    _longest_homopolymer,
    _purine_content,
    _sequence_complexity,
    _gc_skew,
    _molecular_weight,
    _extinction_coefficient,
    _nuclease_resistance_score,
    _cellular_uptake_score,
    _bbb_crossing_score,
    _synthesis_difficulty,
    _off_target_risk,
    _immune_stimulation_risk,
    _duplex_stability,
    _reverse_complement,
    _target_duplex_energy,
    _tm_fit_score,
    _composite_score,
    _accessibility_profile_cached,
    _percentile_rank,
    CHEMISTRY_OPTIONS,
    MODIFICATION_OPTIONS,
    LENGTH_RANGE,
    MIN_GC,
    MAX_GC,
)

logger = logging.getLogger(__name__)

UPREGULATION_CHEMISTRY_OPTIONS = [
    {"id": "gapmer", "label": "DNA Gapmer (2-10-2)", "description": "RNase H1-recruiting; suitable for NAT silencing (A4).",
     "detail": "Central DNA gap recruits RNase H1. For NAT silencing, the ASO targets the antisense transcript. Validated in AON-based upregulation (e.g., Ataluren-class precedents)."},
    {"id": "lna_gapmer", "label": "LNA-enhanced Gapmer", "description": "High-affinity RNase H1 recruitment; NAT silencing (A4).",
     "detail": "LNA wings boost binding affinity (~2-8 C per substitution). Best for high-specificity NAT targeting where allele discrimination matters."},
    {"id": "pmo", "label": "PMO (Phosphorodiamidate Morpholino)", "description": "Steric blocker; ideal for uORF blocking (A5), splice modulation (A3), and miRNA site masking (A6).",
     "detail": "Non-ionic backbone blocks RNA interactions without degradation. Gold standard for uORF steric blocking and exon skipping."},
    {"id": "2ome", "label": "2'-O-Methoxyethyl (2'-OMe)", "description": "Steric blocker; uORF (A5), splice (A3), and miRNA site masking (A6).",
     "detail": "2'-O-Me modifications increase nuclease resistance and reduce immunostimulation. Compatible with steric-blocking and splice-switching."},
    {"id": "sirna", "label": "siRNA duplex (21-mer)", "description": "For saRNA (A23) activation — double-stranded 21-mer duplex.",
     "detail": "Small activating RNAs are 21-mer dsRNA duplexes that target promoter-associated RNAs. Delivered as a duplex guide/passenger pair."},
]

UPREGULATION_LENGTH_RANGE = {"min": 18, "max": 25, "default": 21, "step": 1}

# A design request is a shortlist, not a tiling dump. Uncapped, this tiler
# returned 69 candidates for A5 and 635 for A6 on SCN1A at 20 nt — the whole
# 3' UTR at step size 6, with nothing to distinguish the 125 that tied at
# exactly 100.0. A1 caps at 10 and the A3/A4 designer at 12; TG02 caps here.
UPREGULATION_MAX_CANDIDATES = 20

# Mechanisms whose candidates are steric blockers competing for a site on the
# MATURE mRNA. For these, the oligo has to invade local secondary structure to
# occupy its site, so ViennaRNA unpaired probability is a mechanistic
# requirement and is used as the primary ranking axis (CLAUDE.md 6).
#
# A23 is deliberately NOT in this set: it targets promoter DNA / a
# promoter-associated RNA, and folding genomic sequence as if it were the
# mature transcript would be an accessibility number with no claim behind it.
ACCESSIBILITY_RANKED_MECHANISMS = ("A5", "A6", "A28")

UPREGULATION_MECHANISM_DESIGN = {
    "A3": {
        "label": "TANGO: Poison Exon Skipping / NMD Suppression",
        "target_region": "Exon-exon junctions (poison exon splice sites)",
        "preferred_chemistry": ["pmo", "2ome", "lna_gapmer"],
        "notes": "Masks poison exon splice sites to prevent inclusion of PTC-containing exons, reducing nonsense-mediated decay. Splice-modulating ASOs at exon-exon junctions are most effective. Limited to steric-blocking chemistries only.",
        "tango_fields": [
            {"id": "target_poison_exon", "label": "Target Poison Exon", "type": "dropdown", "description": "Select the poison exon to skip"},
            {"id": "splice_element", "label": "Splice Element", "type": "dropdown", "description": "Select the splice element to mask (5'SS, 3'SS, BPS, or ISS-ISE)"},
        ],
    },
    "A4": {
        "label": "NAT Silencing (RNase H1 Gapmer)",
        "target_region": "Natural antisense transcript (overlapping lncRNA)",
        "preferred_chemistry": ["gapmer", "lna_gapmer"],
        "notes": "Degrades inhibitory antisense lncRNAs that repress the sense gene. Gapmer/LNA chemistry recruits RNase H1 to cleave the NAT transcript. CDS-derived candidates approximate the target complement.",
    },
    "A5": {
        "label": "uORF-blocking ASO",
        "target_region": "5' UTR (uAUG / uORF start site)",
        "preferred_chemistry": ["pmo", "2ome", "lna_gapmer"],
        "notes": "Targets uORF start sites in the 5' UTR to relieve translational repression. Steric-blocking chemistries (PMO/2'-OMe) are preferred — they block ribosome stalling at uORFs without cleaving the transcript. Candidates are windows across the 5' UTR and must be verified against a functionally validated inhibitory uORF.",
    },
    "A6": {
        "label": "miRNA Binding Site Blocking (Target Protector / BlockmiR)",
        "target_region": "3' UTR (miRNA seed binding site)",
        "preferred_chemistry": ["2ome", "pmo", "lna_gapmer"],
        "notes": "Masks the miRNA seed binding site on the target mRNA so the repressive miRNA cannot dock (target protection / blockmiR). miRNA sites are predominantly 3' UTR. Candidates are windows across the real 3' UTR and must be verified against a validated miRNA binding site — no miRNA-target database is integrated in this build.",
    },
    "A23": {
        "label": "saRNA (Small Activating RNA)",
        "target_region": "Promoter (-100 to -1000 bp from TSS)",
        "preferred_chemistry": ["sirna", "lna_gapmer", "gapmer"],
        "forced_length": 21,
        "notes": "saRNA targets promoter-associated RNAs to activate transcription. Requires a functional endogenous promoter. Candidates are windows in the real promoter upstream of the TSS and must be verified against validated promoter elements.",
    },
    "A28": {
        "label": "RBP Binding Site Blocking (Target Protector / RBP Masking)",
        "target_region": "RBP binding sites in 5' UTR, CDS, or 3' UTR",
        "preferred_chemistry": ["pmo", "2ome", "lna_gapmer"],
        "notes": "Masks RBP recognition elements (e.g., PTB, IRP, hnRNP, TTP sites) to relieve translational repression or prevent mRNA decay. Candidates are windows across the transcript and must be verified against validated RBP binding sites.",
    },
}


# A6 (target protection / blockmiR) fundamentally requires a validated
# miRNA binding site. As of this build no reliable machine API exists for
# the standard miRNA-target databases — tested at implementation time:
# ENCORI returns 404 on its documented endpoint, and miRDB's search CGI
# returns an empty JS-rendered template. Rather than scrape fragile HTML,
# candidates are labeled "unverified" and the requirement is stated plainly.
A6_SEED_SITE_NOTE = (
    "No validated miRNA-target database is integrated in this build "
    "(TargetScan/miRDB/ENCORI expose no reliable machine API). Each window "
    "is a putative 3' UTR seed region and must be checked against a "
    "validated miRNA binding site before use."
)


def _fetch_promoter_sequence(target: dict, organism: str) -> str | None:
    """Fetch ~1 kb of genomic sequence upstream of the transcript TSS.

    Uses the canonical transcript coordinates from get_target_analysis. For
    strand +1 the promoter is [start-1000, start-1]; for strand -1 it is
    [end+1, end+1000] requested with strand=-1 so the returned sequence is
    in 5'->3' orientation. Returns None when coordinates are unavailable.
    """
    ct = target.get("canonicalTranscript") or {}
    chromosome = ct.get("chromosome")
    start = ct.get("start")
    end = ct.get("end")
    strand = ct.get("strand", 1)
    if not chromosome or not start or not end:
        return None
    species = (organism or "homo_sapiens").lower().replace(" ", "_")
    if strand == 1:
        a, b = max(1, start - 1000), start - 1
        region_strand = 1
    else:
        a, b = end + 1, end + 1000
        region_strand = -1
    if b < a:
        return None
    try:
        resp = _ensembl_get(
            f"{ENSEMBL_REST}/sequence/region/{species}/{chromosome}:{a}-{b}:{region_strand}"
        )
        if resp.ok:
            return resp.json().get("seq", "").upper()
    except Exception as exc:
        logger.warning("Promoter sequence fetch failed for %s: %s", chromosome, exc)
    return None


def _mechanism_scoring_adjustments(
    mechanism_id: str,
    chemistry: str,
    modifications: list[str],
    gc: float,
    tm: float,
    seq: str,
) -> dict:
    """Compute upregulation-specific design notes."""
    mech_notes = ""

    if mechanism_id == "A3":
        if chemistry in ("pmo", "2ome"):
            mech_notes = "A3 (TANGO): PMO/2'-OMe sterically block splice sites at exon junctions without transcript cleavage — optimal for poison exon skipping."
        elif chemistry == "lna_gapmer":
            mech_notes = "A3 (TANGO): LNA gapmer offers high affinity for precise splice-junction targeting."
        elif chemistry == "gapmer":
            mech_notes = "A3 (TANGO): Gapmers cleave mRNA — not recommended for TANGO which requires steric blocking for precise splice control."

    elif mechanism_id == "A4":
        if chemistry in ("gapmer", "lna_gapmer"):
            mech_notes = "A4 (NAT silencing): Gapmer/LNA recruits RNase H1 to degrade the antisense lncRNA transcript."
        elif chemistry in ("pmo", "2ome"):
            mech_notes = "A4 (NAT silencing): Steric-blocking chemistries don't recruit RNase H — not optimal for NAT transcript degradation."

    elif mechanism_id == "A5":
        if chemistry in ("pmo", "2ome"):
            mech_notes = "A5 (uORF block): PMO/2'-OMe sterically block ribosome stalling at uORFs — ideal for translational upregulation."
        elif chemistry == "lna_gapmer":
            mech_notes = "A5 (uORF block): LNA gapmer can block with high affinity, though RNase H activity is secondary."
        elif chemistry == "gapmer":
            mech_notes = "A5 (uORF block): Gapmers cleave mRNA — less ideal for steric uORF blocking."

    elif mechanism_id == "A6":
        if chemistry in ("2ome", "pmo", "lna_gapmer"):
            mech_notes = "A6 (miRNA site block): Steric-blocking chemistry masks the miRNA seed site on the target mRNA without cleaving it — correct modality for target protection / blockmiR. " + A6_SEED_SITE_NOTE
        elif chemistry == "gapmer":
            mech_notes = "A6 (miRNA site block): Gapmers cleave the transcript — not appropriate for miRNA target protection."

    elif mechanism_id == "A23":
        if chemistry == "sirna":
            mech_notes = "A23 (saRNA): siRNA duplex chemistry is the native modality for transcriptional activation via promoter targeting."
        elif chemistry in ("gapmer", "lna_gapmer"):
            mech_notes = "A23 (saRNA): Gapmer/LNA can activate via RNA-mediated transcriptional activation (RNAa) with promoter proximity."
        elif chemistry in ("pmo", "2ome"):
            mech_notes = "A23 (saRNA): Steric-blocking chemistries are less suitable for promoter activation which requires RNA duplex formation."

    elif mechanism_id == "A28":
        if chemistry in ("pmo", "2ome"):
            mech_notes = "A28 (RBP masking): Steric-blocking chemistry blocks RBP recognition elements without cleaving the transcript — correct modality for translational upregulation via RBP displacement. Candidates are windows across the transcript and must be verified against validated RBP binding sites."
        elif chemistry == "lna_gapmer":
            mech_notes = "A28 (RBP masking): LNA gapmer can block with high affinity, though RNase H activity is secondary."
        elif chemistry == "gapmer":
            mech_notes = "A28 (RBP masking): Gapmers cleave mRNA — less ideal for steric RBP site blocking."

    else:
        mech_notes = "No mechanism-specific adjustments."

    return {"mechNotes": mech_notes}


def generate_upregulation_candidates(
    ensembl_gene_id: str,
    mechanism_id: str,
    aso_length: int,
    chemistry: str,
    modifications: list[str],
    defect_type: str | None = None,
    known_regulatory_element: str | None = None,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    target_poison_exon: str | None = None,
    splice_element: str | None = None,
) -> list[dict]:
    """Generate ASO candidates for upregulation mechanisms.

    Uses the CDS sequence from Ensembl and applies mechanism-specific
    targeting logic. For mechanisms that require sequence context beyond
    the CDS (promoter, 5' UTR, NAT), appropriate design notes are included.
    """
    if mechanism_id not in UPREGULATION_MECHANISM_DESIGN:
        raise ValueError(f"Unsupported upregulation mechanism: {mechanism_id}")

    target = get_target_analysis(ensembl_gene_id, gene_symbol, organism)

    candidates = []
    if not target.get("mrnaSequence") or len(target.get("exons", [])) < 1:
        return candidates

    seq = target["mrnaSequence"].upper()
    seq_len = len(seq)
    exons = target["exons"]
    # A3 and A4 do not target the canonical mRNA and cannot be served by the
    # tiler below.
    #
    # A poison exon is SKIPPED in the productive transcript, so it is absent
    # from the canonical mRNA by definition — tiling that mRNA can never reach
    # it. The old path set the label "Exon junctions" and scanned the whole
    # transcript, returning 891 candidates across all 29 SCN1A exon junctions,
    # none aimed at the poison exon.
    #
    # A NAT is a different gene on the opposite strand. The old path tiled the
    # whole sense transcript as "Full transcript (NAT complement)"; the strand
    # was right, but SCN1A-AS1 overlaps only part of SCN1A, so candidates
    # outside the overlap were complementary to nothing.
    if mechanism_id in ("A3", "A4"):
        from services.upregulation_targets_service import (
            design_nat_knockdown, design_poison_exon_block,
        )
        gene_id = target.get("geneId") or ensembl_gene_id
        if mechanism_id == "A3":
            payload = design_poison_exon_block(
                gene_id, gene_symbol=gene_symbol, organism=organism,
                oligo_length=aso_length)
        else:
            payload = design_nat_knockdown(
                gene_id, gene_symbol=gene_symbol, organism=organism,
                oligo_length=aso_length)
        # This function's contract is a list of candidates. The located target
        # (which poison exon, which NAT, and the alternatives) is attached to
        # each candidate so the caller keeps it without changing the shape.
        located = {k: v for k, v in payload.items() if k != "candidates"}
        for cand in payload["candidates"]:
            cand["mechanismNotes"] = payload.get("architecture", "")
            cand["targetLocated"] = located
        if not payload["candidates"]:
            raise ValueError(payload.get("message") or
                             f"{mechanism_id}: no target could be located.")
        return payload["candidates"]

    design = UPREGULATION_MECHANISM_DESIGN[mechanism_id]

    # Override length for saRNA
    effective_length = design.get("forced_length", aso_length)

    flank = min(10, effective_length // 2)
    step = max(1, effective_length // 3)
    # Mechanism-specific search region. A5/A6/A23 scan real sequence context
    # (5' UTR / 3' UTR / promoter upstream of TSS) when Ensembl provides it;
    # otherwise they fall back to a bounded CDS approximation with an honest
    # label so the positions are not mistaken for true UTR/promoter sites.
    scan_seq = seq
    is_utr_scan = False
    is_promoter_scan = False
    target_label = "Full transcript"
    if mechanism_id == "A3":
        # TANGO: focus on exon junctions
        target_label = "Exon junctions"
    elif mechanism_id == "A4":
        # NAT silencing: scan full transcript
        target_label = "Full transcript (NAT complement)"
    elif mechanism_id == "A5":
        utr5 = target.get("utr5Sequence")
        if utr5:
            scan_seq = utr5.upper()
            is_utr_scan = True
            target_label = "5' UTR (putative uORF windows)"
        else:
            target_label = "5' CDS region (5' UTR unavailable)"
    elif mechanism_id == "A6":
        utr3 = target.get("utr3Sequence")
        if utr3:
            scan_seq = utr3.upper()
            is_utr_scan = True
            target_label = "3' UTR (putative miRNA seed windows)"
        else:
            target_label = "3' CDS region (3' UTR unavailable)"
    elif mechanism_id == "A23":
        promoter_seq = _fetch_promoter_sequence(target, organism)
        if promoter_seq:
            scan_seq = promoter_seq
            is_promoter_scan = True
            target_label = "Promoter (upstream of TSS)"
        else:
            target_label = "5' CDS region (promoter unavailable)"
    elif mechanism_id == "A28":
        target_label = "Full transcript (putative RBP binding site windows)"

    search_start = 0
    search_end = max(0, len(scan_seq) - effective_length)
    # The 90/400 base bounds below are fallback-only heuristics applied when
    # Ensembl UTR/promoter data is unavailable. They are arbitrary and are
    # labelled as approximations; they are never claimed as sourced values.
    if mechanism_id in ("A5", "A23") and not is_utr_scan and not is_promoter_scan:
        search_end = min(search_end, 90)
    elif mechanism_id == "A6" and not is_utr_scan:
        search_start = max(0, len(scan_seq) - effective_length - 400)

    seen = set()

    # Unpaired-probability profile over the region actually scanned, in ONE
    # fold, keyed by the same offset the tiler uses. RNAplfold here is local
    # (max_bp_span 150, window 200), so folding the isolated UTR matches its
    # profile inside the full transcript except within ~200 nt of the splice
    # to the CDS — the boundary caveat is stated on every candidate.
    accessibility_profile = {}
    if mechanism_id in ACCESSIBILITY_RANKED_MECHANISMS:
        accessibility_profile = dict(
            _accessibility_profile_cached(scan_seq, effective_length)
        )

    # Build exon CDS mapping (same as gene_silencing_service)
    total_genomic = sum(e.get("length", 0) for e in exons)
    if total_genomic == 0:
        return candidates

    exon_cds_map = []
    cursor = 0
    for exon in exons:
        cds_contribution = round(seq_len * exon.get("length", 0) / total_genomic)
        exon_cds_map.append((cursor, cursor + cds_contribution))
        cursor += cds_contribution
    if exon_cds_map:
        last_start, _ = exon_cds_map[-1]
        exon_cds_map[-1] = (last_start, seq_len)

    is_poison_exon = mechanism_id == "A3"

    for offset in range(search_start, search_end + 1, step):
        candidate_seq = scan_seq[offset : offset + effective_length]
        if len(candidate_seq) < effective_length or candidate_seq in seen:
            continue
        seen.add(candidate_seq)

        gc = _calc_gc(candidate_seq)
        if gc < MIN_GC or gc > MAX_GC:
            continue

        tm = _calc_tm(candidate_seq)
        self_mfe = _self_complement_mfe(candidate_seq)
        pg = _polyg_score(candidate_seq)
        cpg = _cpg_count(candidate_seq)

        mech_adj = _mechanism_scoring_adjustments(
            mechanism_id, chemistry, modifications, gc, tm, candidate_seq
        )

        # Target duplex energy — fold the antisense ASO against its RNA target
        # window, mirroring the TG01 pipeline. The stored ASO is the reverse
        # complement of the scanned window; for A4 the RNA target is the
        # antisense transcript, so the sense window itself is the ASO.
        aso_seq = candidate_seq if mechanism_id == "A4" else _reverse_complement(candidate_seq)
        duplex_energy = _target_duplex_energy(aso_seq, _reverse_complement(aso_seq))

        # Upregulation-specific defect notes
        defect_notes = "No defect-specific adjustments applied."
        upreg_defect = (defect_type or "").lower().strip()
        if "haploinsufficiency" in upreg_defect or "loss-of-function" in upreg_defect or "lof" in upreg_defect:
            defect_notes = "Underexpression / haploinsufficiency: upregulation is the appropriate therapeutic strategy."
        elif "dominant" in upreg_defect:
            defect_notes = "Dominant-negative: consider allele-specific upregulation of the wild-type copy."

        # Determine exon number (CDS scans only; UTR/promoter scans have none)
        exon_number = None
        exon_length = None
        if not is_utr_scan and not is_promoter_scan:
            for ei, (es, ee) in enumerate(exon_cds_map):
                if es <= offset < ee:
                    exon_number = ei + 1
                    exon_length = exons[ei].get("length") if ei < len(exons) else None
                    break

        region_label = f"{target_label} offset +{offset}"
        if is_poison_exon and exon_number:
            region_label = f"Exon {exon_number} junction offset +{offset}"

        nuc_res = _nuclease_resistance_score(chemistry, modifications)
        uptake = _cellular_uptake_score(chemistry, effective_length)
        bbb = _bbb_crossing_score(chemistry, effective_length, modifications)
        synth = _synthesis_difficulty(candidate_seq, chemistry, modifications)
        off_target = _off_target_risk(candidate_seq, _sequence_complexity(candidate_seq))
        immune = _immune_stimulation_risk(candidate_seq, chemistry)

        complexity = _sequence_complexity(candidate_seq)
        skew = _gc_skew(candidate_seq)
        mw = _molecular_weight(candidate_seq)
        ec = _extinction_coefficient(candidate_seq)
        ds = _duplex_stability(gc, tm, effective_length)

        # TANGO-specific: echo the user's design inputs only. Earlier versions
        # fabricated spliceMaskingScore / predictedNmdSuppression /
        # estimatedFoldRestoration / canonicalOffSpliceHits from invented
        # constants (0.85/0.80/0.75/0.70 element bases, x0.9 NMD, x0.8 fold)
        # with no measured basis — those are removed.
        tango_fields = {}
        if is_poison_exon:
            tango_fields = {
                "targetPoisonExon": target_poison_exon or "",
                "spliceElement": splice_element or "",
            }

        # Composite ranking score — real metrics only, mirroring TG01.
        tm_fit = _tm_fit_score(tm, chemistry, modifications, mechanism_id)
        # Normalize target duplex ΔG (more negative = stronger binding), matching
        # the normalization used inside _composite_score so the two components
        # are exposed for the score-decomposition visualization.
        duplex_score = min(100.0, max(0.0, (-duplex_energy - 8.0) * 3.5))
        composite_score = _composite_score(duplex_energy, tm_fit)

        candidates.append({
            "sequence": aso_seq,
            "length": effective_length,
            "compositeScore": composite_score,  # 0-100 ranking score
            # Exposed weighted components of the composite score so the UI can
            # show exactly how each candidate was ranked.
            "scoreBreakdown": {
                "duplexScore": round(duplex_score, 1),  # normalized 0-100 binding
                "tmFitScore": tm_fit,                    # normalized 0-100 Tm fit
                "duplexRaw": duplex_energy,              # kcal/mol
                "tmRaw": tm,                             # °C
            },
            "learnedEfficacy": {
                "available": False,
                "value": None,
                "modelInfo": "Not yet trained",
                "scopeCaveat": None,
            },
            # Measured / computed properties — exact physics and sequence
            # computations, the tier that drives ranking.
            "realMetrics": {
                "siteAccessibility": accessibility_profile.get(offset),
                "targetDuplexEnergy": duplex_energy,
                "meltingTempC": tm,
                "selfStructureMfe": self_mfe,
                "gcContent": round(gc * 100, 1),
                "cpgCount": cpg,
                "longestHomopolymer": _longest_homopolymer(candidate_seq),
                "purineContent": _purine_content(candidate_seq),
                "gcSkew": skew,
                "sequenceComplexity": complexity,
                "polyGPass": pg == 0,
                "molecularWeight": mw,
                "extinctionCoefficient": ec,
                "duplexStability": ds,
            },
            # Rule-of-thumb drug-like estimates — deliberately excluded from
            # ranking, shown to the user as labeled estimates.
            "heuristicEstimates": {
                "nucleaseResistance": {
                    "value": nuc_res,
                    "note": "Chemistry-class rule of thumb, not measured.",
                },
                "cellularUptake": {
                    "value": uptake,
                    "note": "Length/chemistry rule of thumb, not measured.",
                },
                "bbbCrossing": {
                    "value": bbb,
                    "note": "Length/chemistry rule of thumb, not measured.",
                },
                "synthesisDifficulty": {
                    "value": synth,
                    "note": "Sequence/chemistry rule of thumb, not measured.",
                },
                "offTargetRisk": {
                    "value": off_target,
                    "note": "Length/repetitiveness heuristic — not a genome alignment check.",
                },
                "immuneStimulation": {
                    "value": immune,
                    "note": "CpG-count heuristic, not an immunogenicity assay.",
                },
            },
            "targetRegion": region_label,
            "mechanismId": mechanism_id,
            "chemistry": chemistry,
            "modifications": modifications,
            "exonNumber": exon_number,
            "exonLength": exon_length,
            "deliveryContext": "",
            "defectType": defect_type or "",
            "defectNotes": defect_notes,
            "mechanismNotes": mech_adj["mechNotes"],
            "knownRegulatoryElement": known_regulatory_element or "",
            # A6 candidates cannot claim to mask validated miRNA binding sites
            # (no reliable target database is integrated) — state that plainly.
            "seedSiteStatus": "unverified" if mechanism_id == "A6" else None,
            "seedSiteNote": A6_SEED_SITE_NOTE if mechanism_id == "A6" else None,
            **tango_fields,
        })

    # THE COMPOSITE SCORE SATURATES HERE FOR THE SAME REASON IT DID IN TG01.
    #
    # `_composite_score` maps duplex ΔG through `(-dg - 8) * 3.5` clipped at
    # 100. Those constants suit oligos from 12 to 30 nt, but every candidate
    # in one run is the SAME length, so the within-run spread is a few
    # kcal/mol and most of it clips. Measured on SCN1A at 20 nt before this
    # change: A5 returned 69 candidates, 19 of them tied at exactly 100.0;
    # A6 returned 635, with 125 tied at exactly 100.0.
    #
    # Percentiles are taken within THIS pool, so they discriminate by
    # construction rather than against an absolute constant that does not fit
    # the run (CLAUDE.md 5, "rank within the pool").
    acc_pool = [c["realMetrics"].get("siteAccessibility") for c in candidates]
    dg_pool = [c["realMetrics"].get("targetDuplexEnergy") for c in candidates]
    has_accessibility = any(v is not None for v in acc_pool)
    for c in candidates:
        rm = c["realMetrics"]
        c["accessibilityPercentile"] = _percentile_rank(
            rm.get("siteAccessibility"), acc_pool, higher_is_better=True)
        c["duplexEnergyPercentile"] = _percentile_rank(
            rm.get("targetDuplexEnergy"), dg_pool, higher_is_better=False)
        if has_accessibility:
            c["rankingBasis"] = {
                "primary": "siteAccessibility (ViennaRNA unpaired probability)",
                "tieBreak": "compositeScore, then duplex ΔG",
                "caveat": (
                    "A steric blocker has to occupy its site, so local "
                    "accessibility is a mechanistic requirement — it is not a "
                    "validated activity model, and no uORF or miRNA site here "
                    "is itself validated. Duplex-ΔG-first ordering was "
                    "measured below chance on 528 held-out experiments "
                    "(backend/experiments/ml_analysis, E12), so it is a "
                    "tie-break rather than the primary axis. The fold is "
                    "local to the scanned region; sites within ~200 nt of its "
                    "boundary may pair with sequence outside it."
                ),
            }
        else:
            c["rankingBasis"] = {
                "primary": "compositeScore (duplex ΔG + Tm fit)",
                "tieBreak": "duplex ΔG",
                "caveat": (
                    "No accessibility profile applies to this mechanism, so "
                    "ranking falls back to the composite score. That score "
                    "saturates when candidates share a length — read the "
                    "percentiles, not the raw value."
                ),
            }

    # Accessibility varies over orders of magnitude between sites on one
    # transcript (measured on the SCN1A 3' UTR at 20 nt: 2,823 distinct values
    # across 5,413 windows), which is why it leads and the saturating
    # composite only breaks its ties.
    if has_accessibility:
        candidates.sort(
            key=lambda c: (
                -(c["realMetrics"].get("siteAccessibility") or 0.0),
                -c["compositeScore"],
                c["realMetrics"]["targetDuplexEnergy"],
            )
        )
    else:
        candidates.sort(
            key=lambda c: (-c["compositeScore"],
                           c["realMetrics"]["targetDuplexEnergy"])
        )

    # Report how much of the tiling the shortlist came from, so a capped list
    # is never mistaken for the whole search.
    total_before_cap = len(candidates)
    candidates = candidates[:UPREGULATION_MAX_CANDIDATES]
    for c in candidates:
        c["poolSize"] = total_before_cap
        c["shortlistedFrom"] = (
            f"Top {len(candidates)} of {total_before_cap} windows tiled across "
            f"{target_label}."
        )

    return candidates


def get_upregulation_design_options() -> dict:
    """Return available chemistry, modification, and length options
    for upregulation mechanisms.
    """
    return {
        "chemistryOptions": UPREGULATION_CHEMISTRY_OPTIONS,
        "modificationOptions": MODIFICATION_OPTIONS,
        "lengthRange": UPREGULATION_LENGTH_RANGE,
        "mechanisms": UPREGULATION_MECHANISM_DESIGN,
    }
