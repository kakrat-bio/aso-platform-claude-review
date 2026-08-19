"""Sequence upload service — parsing, validation, and analysis.

Supports FASTA files, raw DNA/RNA sequences, and auto-detection of sequence type.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from services.sequence_metrics import gc_content as _gc_content

# ---------------------------------------------------------------------------
# IUPAC sets
# ---------------------------------------------------------------------------

_IUPAC_DNA = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
_IUPAC_RNA = set("ACGURYSWKMBDHVNacguryswkmbdhvn")
_IUPAC_PROTEIN = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwyBZXJbxzj")

# Common ASO modifications (simplified detection)
_ASO_PATTERNS = [
    re.compile(r"[acgtuACGTU]{15,30}"),  # typical ASO length
]


def _clean_sequence(raw: str) -> str:
    """Remove FASTA headers, whitespace, and numbering."""
    lines = raw.strip().splitlines()
    seq_lines = [l for l in lines if not l.startswith(">") and not l.startswith(";")]
    return "".join(seq_lines).replace(" ", "").replace("\t", "").replace("\n", "").upper()


def _detect_type(seq: str) -> str:
    """Auto-detect sequence type: dna, rna, protein, or unknown."""
    has_t = "T" in seq
    has_u = "U" in seq
    if has_t and not has_u:
        return "dna"
    if has_u and not has_t:
        return "rna"
    if not has_t and not has_u:
        # Could be protein or short DNA — check for protein-specific letters
        protein_letters = set("DEFHIKLMNPQRSTVWYdefhiklmnopqrstvwy")
        if any(c in protein_letters for c in seq):
            return "protein"
        return "dna"
    return "unknown"


def _find_invalid_chars(seq: str, seq_type: str) -> List[str]:
    if seq_type == "dna":
        valid = _IUPAC_DNA
    elif seq_type == "rna":
        valid = _IUPAC_RNA
    elif seq_type == "protein":
        valid = _IUPAC_PROTEIN
    else:
        valid = _IUPAC_DNA | _IUPAC_RNA | _IUPAC_PROTEIN
    return sorted(set(b for b in seq if b not in valid))


def _count_motif(seq: str, motif: str) -> int:
    return seq.count(motif.upper())


def _has_poly_g(seq: str, min_run: int = 4) -> bool:
    return bool(re.search(rf"G{{{min_run},}}", seq.upper()))


def _has_poly_a_tail(seq: str, min_run: int = 6) -> bool:
    return bool(re.search(rf"A{{{min_run},}}$", seq.upper()))


def _reverse_complement(seq: str, is_rna: bool) -> str:
    comp = {"A": "U" if is_rna else "T", "T": "A", "U": "A", "G": "C", "C": "G"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def _find_orfs(seq: str) -> List[Dict[str, Any]]:
    """
    Find open reading frames in a DNA/RNA sequence, scanning both strands
    (6 frames total). Forward-only ORF finding would miss real coding
    potential encoded on the reverse strand.
    """
    is_rna = "U" in seq and "T" not in seq
    start_codons = {"AUG"} if is_rna else {"ATG"}
    stop_codons = {"UAA", "UAG", "UGA"} if is_rna else {"TAA", "TAG", "TGA"}

    def scan(strand_seq: str, strand_label: str) -> List[Dict[str, Any]]:
        found = []
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                codon = strand_seq[i:i+3]
                if codon in start_codons:
                    start = i
                    j = i + 3
                    while j < len(strand_seq) - 2:
                        c = strand_seq[j:j+3]
                        if c in stop_codons:
                            found.append({
                                "strand": strand_label,
                                "frame": frame + 1,
                                "start": start + 1,
                                "end": j + 3,
                                "length": j + 3 - start,
                                "proteinLength": (j - start) // 3,
                            })
                            break
                        j += 3
                    i = j + 3 if j < len(strand_seq) else len(strand_seq)
                else:
                    i += 3
        return found

    orfs = scan(seq, "+")
    orfs += scan(_reverse_complement(seq, is_rna), "-")
    return orfs


def _immunostimulatory_motifs(seq: str, seq_type: str) -> List[Dict[str, Any]]:
    """
    Flags sequence patterns loosely associated with innate-immune sensing
    in the literature. This is pattern-matching against a short heuristic
    list, NOT a validated immunogenicity assay — real TLR7/8 recognition
    depends on broader GU/U-rich context than any single hexamer, and TLR9
    specifically senses unmethylated CpG in a DNA context with defined
    flanking bases, not any CG dinucleotide appearing in an RNA oligo.

    Returns every match (capped) with its position, so the frontend can
    plot hits along the sequence rather than just report a count.
    """
    motifs = []
    checks = [
        (r"[GU]{2,}U[GU]{2,}", "GU-rich stretch (literature-associated with TLR7/8 sensing; not a confirmed motif)"),
        (r"(.)\1{3,}", "Homopolymer run (4+ repeats; general repetitive-element flag)"),
    ]
    if seq_type == "dna":
        checks.append((r"[AG][AG]CG[CT][CT]", "Unmethylated CpG in a purine-purine-CG-pyrimidine-pyrimidine context (literature TLR9 motif pattern; not a confirmed assay)"))

    for pattern, label in checks:
        for m in re.finditer(pattern, seq.upper()):
            motifs.append({
                "motif": m.group(0),
                "label": label,
                "start": m.start() + 1,
                "end": m.end(),
            })
            if len(motifs) >= 40:  # cap payload size for highly repetitive sequences
                return motifs
    return motifs


def _secondary_structure_score(seq: str) -> Dict[str, Any]:
    """Simplified secondary structure prediction based on GC content and self-complementarity."""
    gc = _gc_content(seq)
    length = len(seq)

    # Estimate hairpin propensity — record actual positions, not just a count
    palindrome_positions = []
    for i in range(length - 5):
        chunk = seq[i:i+6]
        if chunk == chunk[::-1]:
            palindrome_positions.append(i + 1)

    # MFE estimate (very simplified — composition-based, not a real fold)
    gc_stability = gc / 100 * -1.5  # kcal/mol per GC pair estimate
    au_stability = (100 - gc) / 100 * -0.9  # kcal/mol per AU pair
    estimated_mfe = round((gc_stability + au_stability) * length / 2, 1)

    return {
        "estimatedMfe": estimated_mfe,
        "palindromicRegions": len(palindrome_positions),
        "palindromePositions": palindrome_positions[:50],  # cap payload
        "gcContent": gc,
        "hairpinRisk": "High" if len(palindrome_positions) > 3 else "Medium" if len(palindrome_positions) > 1 else "Low",
    }


def _gc_sliding_window(seq: str, window: int = 10, step: int = 2) -> List[Dict[str, Any]]:
    """Real per-window GC% across the sequence, for plotting GC distribution
    rather than just a single average."""
    length = len(seq)
    if length < window:
        return [{"position": 1, "gc": _gc_content(seq)}]
    points = []
    for i in range(0, length - window + 1, step):
        chunk = seq[i:i + window]
        points.append({"position": i + 1, "gc": _gc_content(chunk)})
    return points


def _nucleotide_composition(seq: str) -> Dict[str, int]:
    """Real base counts from the actual sequence."""
    seq = seq.upper()
    return {
        "A": seq.count("A"),
        "C": seq.count("C"),
        "G": seq.count("G"),
        "T": seq.count("T"),
        "U": seq.count("U"),
    }


def _amino_acid_composition(seq: str) -> Dict[str, int]:
    """Real amino acid counts from the actual sequence."""
    seq = seq.upper()
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return {a: seq.count(a) for a in aa}


def _protein_molecular_weight(seq: str) -> float:
    """Approximate molecular weight in Da for a protein sequence."""
    weights = {
        "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
        "E": 147.13, "Q": 146.15, "G": 75.07, "H": 155.16, "I": 131.17,
        "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
        "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
    }
    total = sum(weights.get(b, 110.0) for b in seq.upper())
    return round(total - (len(seq) - 1) * 18.015, 1)


def _protein_analysis(seq: str) -> Dict[str, Any]:
    """Protein-specific analysis."""
    aa_counts = _amino_acid_composition(seq)
    length = len(seq)
    mw = _protein_molecular_weight(seq)

    hydrophobic = sum(aa_counts.get(a, 0) for a in "AILMFWYV")
    hydrophilic = sum(aa_counts.get(a, 0) for a in "RNDQEKHP")
    charged = sum(aa_counts.get(a, 0) for a in "RDEKHP")
    aromatic = sum(aa_counts.get(a, 0) for a in "FWY")

    return {
        "aminoAcidComposition": aa_counts,
        "molecularWeight": mw,
        "length": length,
        "hydrophobicFraction": round(hydrophobic / length, 3) if length else 0,
        "hydrophilicFraction": round(hydrophilic / length, 3) if length else 0,
        "chargedFraction": round(charged / length, 3) if length else 0,
        "aromaticFraction": round(aromatic / length, 3) if length else 0,
    }


def validate_sequence(raw_input: str, filename: Optional[str] = None) -> Dict[str, Any]:
    """Parse and validate an uploaded sequence.

    Returns a validation report with sequence stats, detected type, and any issues.
    """
    seq = _clean_sequence(raw_input)
    if not seq:
        return {"valid": False, "error": "No valid sequence found in the input."}

    seq_type = _detect_type(seq)
    invalid = _find_invalid_chars(seq, seq_type)

    length = len(seq)

    features = []
    if seq_type == "protein":
        features.append("Protein sequence detected")
        aa_counts = _amino_acid_composition(seq)
        most_common = max(aa_counts, key=aa_counts.get)
        features.append(f"Most common residue: {most_common} ({aa_counts[most_common]}x)")
    else:
        gc = _gc_content(seq)
        features.append(f"GC content: {gc}%")
        if _has_poly_a_tail(seq):
            features.append("Poly-A tail detected")
        if _has_poly_g(seq):
            features.append("Poly-G tract detected")

        orfs = _find_orfs(seq)
        if orfs:
            best = max(orfs, key=lambda o: o["proteinLength"])
            features.append(f"Longest ORF: {best['proteinLength']} aa (frame {best['frame']})")

    return {
        "valid": len(invalid) == 0,
        "sequence": seq,
        "sequenceType": seq_type,
        "length": length,
        "gcContent": _gc_content(seq) if seq_type != "protein" else None,
        "invalidChars": invalid,
        "features": features,
        "orfs": _find_orfs(seq)[:5] if seq_type != "protein" else [],
        "filename": filename,
        "hasPolyA": _has_poly_a_tail(seq) if seq_type != "protein" else False,
        "hasPolyG": _has_poly_g(seq) if seq_type != "protein" else False,
    }


def analyze_sequence(seq: str, modality: str) -> Dict[str, Any]:
    """Run full analysis on a validated sequence for a given therapeutic modality."""
    seq_type = _detect_type(seq)
    gc = _gc_content(seq)
    length = len(seq)
    is_protein = seq_type == "protein"

    result: Dict[str, Any] = {
        "sequence": seq,
        "sequenceType": seq_type,
        "length": length,
        "gcContent": gc if not is_protein else None,
    }

    if is_protein:
        result.update({
            "proteinAnalysis": _protein_analysis(seq),
            "offTarget": {"lengthBasedRiskEstimate": "N/A", "note": "Protein sequence — off-target screening requires nucleotide-level alignment", "internalRepetitiveness": 0, "recommendedMinLength": 0, "disclaimer": "Not applicable for protein sequences."},
            "secondaryStructure": {"estimatedMfe": None, "palindromicRegions": 0, "palindromePositions": [], "gcContent": None, "hairpinRisk": "N/A"},
            "immuneScreen": [],
            "modality": _modality_analysis(seq, seq_type, modality),
            "gcCurve": [],
            "composition": {},
            "orfs": [],
            "meltingTemp": None,
            "complexity": None,
            "codonUsage": None,
            "modificationScores": None,
            "energyProfile": [],
            "grnaCandidates": [],
        })
        return result

    # Nucleic acid analysis
    off_target = _estimate_off_targets(seq)
    structure = _secondary_structure_score(seq)
    immune = _immunostimulatory_motifs(seq, seq_type)
    modality_results = _modality_analysis(seq, seq_type, modality)
    gc_curve = _gc_sliding_window(seq)
    composition = _nucleotide_composition(seq)
    orfs = _find_orfs(seq)
    tm_data = _melting_temperature(seq)
    complexity = _sequence_complexity(seq)
    codon_data = _codon_usage(seq)
    mod_scores = _modification_scorecard(seq, modality)
    energy_profile = _stacking_energy_profile(seq)
    grna_candidates = _generate_grna_candidates(seq) if modality == "sgrna" else []

    result.update({
        "offTarget": off_target,
        "secondaryStructure": structure,
        "immuneScreen": immune,
        "modality": modality_results,
        "gcCurve": gc_curve,
        "composition": composition,
        "orfs": orfs[:20],
        "meltingTemp": tm_data,
        "complexity": complexity,
        "codonUsage": codon_data,
        "modificationScores": mod_scores,
        "energyProfile": energy_profile,
        "grnaCandidates": grna_candidates,
    })
    return result


def _estimate_off_targets(seq: str) -> Dict[str, Any]:
    """
    Sequence-uniqueness heuristic — NOT a real off-target screen. This does
    not align the sequence against any genome or transcriptome, so it can't
    actually detect off-target binding sites. It only reflects length and
    local k-mer repetitiveness, which are weak, indirect proxies. Labeled
    "specificityHeuristic" rather than "offTarget" so it isn't mistaken for
    a real BLAST/alignment-based check (out of scope here — see Page 3
    candidate design notes for the same boundary).
    """
    length = len(seq)
    if length < 18:
        risk = "High"
        note = "Short sequence — generally correlates with higher off-target probability, not verified against any genome"
    elif length < 20:
        risk = "Medium"
        note = "Moderate length — not verified against any genome"
    else:
        risk = "Low"
        note = "Adequate length for specificity in general, not verified against any genome"

    k = 6
    if length >= k:
        kmers = [seq[i:i+k] for i in range(length - k + 1)]
        unique_kmers = len(set(kmers))
        repetitiveness = round(1 - unique_kmers / len(kmers), 3) if kmers else 0
    else:
        repetitiveness = 0

    return {
        "lengthBasedRiskEstimate": risk,
        "note": note,
        "internalRepetitiveness": repetitiveness,
        "recommendedMinLength": 18,
        "disclaimer": "This is a length/repetitiveness heuristic only — it does not check the sequence against any real genome or transcriptome. Use a real alignment tool (e.g. BLAST) for actual off-target screening.",
    }


def _modality_analysis(seq: str, seq_type: str, modality: str) -> Dict[str, Any]:
    """Modality-specific analysis."""
    gc = _gc_content(seq)
    length = len(seq)

    if modality == "aso":
        return _aso_analysis(seq, seq_type, gc, length)
    elif modality == "sirna":
        return _sirna_analysis(seq, gc, length)
    elif modality == "mrna":
        return _mrna_analysis(seq, seq_type, gc, length)
    elif modality == "sgrna":
        return _sgrna_analysis(seq, seq_type, gc, length)
    else:
        return {"recommendation": "Select a modality for detailed analysis"}


def _aso_analysis(seq: str, seq_type: str, gc: float, length: int) -> Dict[str, Any]:
    """ASO-specific analysis."""
    recommendations = []
    if gc < 30:
        recommendations.append("Low GC% — consider LNA or 2'-OMe modifications to boost Tm")
    elif gc > 70:
        recommendations.append("High GC% — risk of G-quadruplexes; consider shorter ASO")
    else:
        recommendations.append("GC content in optimal range for RNase H recruitment")

    if length < 15:
        recommendations.append("Very short — high off-target risk; minimum 18 nt recommended")
    elif length > 25:
        recommendations.append("Long ASO — may have reduced cellular uptake; consider gapmer design")

    chemistry = "gapmer" if gc >= 35 else "pmo"
    return {
        "recommendedChemistry": chemistry,
        "recommendations": recommendations,
        "optimalLength": "18-22 nt",
        "targetRegion": "Exon junction or mutated region recommended",
    }


def _sirna_analysis(seq: str, gc: float, length: int) -> Dict[str, Any]:
    """siRNA-specific analysis."""
    recommendations = []
    if length < 19 or length > 25:
        recommendations.append("Optimal siRNA length is 19-25 nt")
    if gc < 30 or gc > 52:
        recommendations.append("Optimal GC content for siRNA is 30-52%")

    # Check for 3' UU overhang potential
    if seq.endswith("UU") or seq.endswith("TT"):
        recommendations.append("3' UU/TT dinucleotide overhang detected — good for RISC loading")

    return {
        "strand": "Guide strand (antisense) + Passenger strand",
        "recommendations": recommendations,
        "optimalLength": "21 nt with 2-nt 3' overhangs",
        "thermodynamicBias": "Asymmetry: 5' thermodynamic instability of guide strand preferred",
    }


def _mrna_analysis(seq: str, seq_type: str, gc: float, length: int) -> Dict[str, Any]:
    """mRNA-specific analysis."""
    recommendations = []
    if seq_type != "rna":
        recommendations.append("Convert to RNA (T→U) for mRNA therapeutic design")

    orfs = _find_orfs(seq)  # handles DNA and RNA internally
    if orfs:
        best = max(orfs, key=lambda o: o["proteinLength"])
        recommendations.append(f"Longest coding ORF: {best['proteinLength']} amino acids")

    if not _has_poly_a_tail(seq):
        recommendations.append("No poly-A tail — add 100-150 nt poly(A) for stability")

    recommendations.append("Consider 5' Cap analog (Anti-Reverse Cap ARCA)")
    recommendations.append("Evaluate codon optimization for human expression")

    return {
        "recommendations": recommendations,
        "needsCodonOptimization": True,
        "needsPolyA": not _has_poly_a_tail(seq),
        "needsUTR": True,
        "nucleosideModifications": ["N1-methylpseudouridine (m1Ψ)", "5-methylcytidine (m5C)"],
    }


def _sgrna_analysis(seq: str, seq_type: str, gc: float, length: int) -> Dict[str, Any]:
    """sgRNA/CRISPR-specific analysis."""
    recommendations = []
    if length < 17 or length > 21:
        recommendations.append("Optimal sgRNA spacer length is 17-21 nt (20 nt standard)")

    # SpCas9 requires an NGG PAM adjacent to the target site
    pam = "NGG"
    recommendations.append(f"Requires {pam} PAM adjacent to target site (SpCas9)")

    if gc < 40 or gc > 80:
        recommendations.append("Optimal GC content for sgRNA is 40-80%")

    # Check for poly-T (Pol III terminator)
    if "TTTT" in seq.upper():
        recommendations.append("Poly-T tract detected — may cause premature transcription termination")

    return {
        "casProtein": "SpCas9 (NGG PAM)",
        "recommendations": recommendations,
        "optimalLength": "20 nt spacer + PAM",
        "offTargetMitigation": "Consider truncated sgRNAs (17-18 nt) for improved specificity",
    }


def _generate_grna_candidates(seq: str) -> List[Dict[str, Any]]:
    """Scan sequence for NGG PAMs and score 20 nt spacers."""
    upper = seq.upper()
    length = len(upper)
    candidates: List[Dict[str, Any]] = []

    for i in range(length - 22):
        pam = upper[i + 20 : i + 23]
        # `N` is a LITERAL N in a regex, not "any base" — `re.fullmatch(r"NGG",
        # pam)` only ever matched the three characters "NGG", which does not
        # occur in a real sequence. This scanner therefore returned zero
        # candidates for every input, on a test sequence carrying nine valid
        # NGG PAMs. Spell the degeneracy out as a character class.
        if not re.fullmatch(r"[ACGT]GG", pam):
            continue

        spacer = upper[i : i + 20]
        gc_count = sum(1 for b in spacer if b in "GC")
        gc_pct = gc_count / 20.0

        score = 0.0
        if 0.4 <= gc_pct <= 0.8:
            score += 40.0
        elif 0.3 <= gc_pct <= 0.85:
            score += 25.0
        else:
            score += 10.0

        if spacer[0] in "GC":
            score += 10.0
        if "TTTT" not in spacer:
            score += 10.0

        self_comp = _self_complementarity_score(spacer)
        if self_comp < 0.2:
            score += 15.0
        elif self_comp < 0.4:
            score += 8.0

        # Internal 6-mer repetitiveness, in [0, 1]. This used to be multiplied
        # by 12 and reported as `offTargets` — an integer COUNT of genomic
        # off-target sites, from a statistic that never touches a genome. A
        # count invites the reader to believe a search happened. Only the
        # repetitiveness it was derived from is reported now.
        repetitiveness = _kmer_repetitiveness(spacer)
        if repetitiveness <= 0.1:
            score += 15.0
        elif repetitiveness <= 0.25:
            score += 8.0
        else:
            score += 3.0

        score = max(0.0, min(100.0, score))

        color = "emerald" if score >= 70 else "amber" if score >= 40 else "rose"

        candidates.append({
            "id": f"gRNA-{len(candidates) + 1}",
            "position": i + 1,
            "sequence": spacer,
            "pam": pam,
            "strand": "+",
            "score": round(score, 1),
            "gc": round(gc_pct * 100, 1),
            "selfComplementarity": round(self_comp, 3),
            "internalRepetitiveness": round(repetitiveness, 3),
            "offTargetsNote": (
                "Not screened. No genome or transcriptome alignment is "
                "performed anywhere in this service; use BLAST or Cas-OFFinder "
                "for a real off-target search."
            ),
            "polyT": "TTTT" in spacer,
            "color": color,
        })

    return candidates


def _self_complementarity_score(seq: str) -> float:
    """Simplified self-complementarity: fraction of matching bases between sequence and its reverse complement."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N"}
    rc = "".join(comp.get(b, "N") for b in reversed(seq.upper()))
    matches = sum(1 for a, b in zip(seq.upper(), rc) if a == b)
    return matches / max(len(seq), 1)


def _kmer_repetitiveness(seq: str, k: int = 6) -> float:
    """Fraction of k-mers in the sequence that are repeats of another k-mer.

    Was `_estimate_off_target_count`, which scaled this by 12 and returned it
    as a genomic off-target count. It is not one: nothing here aligns against
    a genome. Repetitiveness is a weak proxy for a sequence being hard to
    place uniquely, and it is reported as exactly that.
    """
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    if not kmers:
        return 0.0
    return 1 - len(set(kmers)) / len(kmers)


# ---------------------------------------------------------------------------
# Melting temperature (nearest-neighbor simplified)
# ---------------------------------------------------------------------------

# DNA nearest-neighbor parameters (SantaLucia 1998) kcal/mol
_DNA_NN: Dict[str, Tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}


def _melting_temperature(seq: str) -> Dict[str, Any]:
    """
    Nearest-neighbor Tm estimate for short oligos (< 50 nt) at 50 mM Na+.
    For longer sequences, a basic GC% formula is also provided as a
    secondary reference. Neither replaces experimental Tm measurement.
    """
    seq = seq.upper().replace("U", "T")
    length = len(seq)
    if length < 6:
        return {
            "tmNearestNeighbor": 0,
            "tmBasicGC": 0,
            "length": length,
            "method": "Sequence too short for reliable Tm estimation",
            "note": "Tm < 6 nt is not meaningful for most oligonucleotide applications",
        }

    # Nearest-neighbor (simplified, assumes 50 mM Na+)
    dH = 0.0
    dS = 0.0
    n = 0
    for i in range(length - 1):
        dinuc = seq[i:i+2]
        if dinuc in _DNA_NN:
            h, s = _DNA_NN[dinuc]
            dH += h
            dS += s
            n += 1

    if n > 0:
        # Tm = dH / (dS + R * ln(C/4)) - 273.15  (simplified)
        import math
        R = 1.987  # cal/(mol·K)
        C = 250e-6  # 250 µM typical oligo concentration
        tm_nn = round(dH * 1000 / (dS + R * math.log(C / 4)) - 273.15, 1) if dS != 0 else 0
    else:
        tm_nn = 0

    # Basic GC% formula (Bolton & McCarthy, modified)
    gc = _gc_content(seq)
    tm_basic = round(64.9 + 41 * (gc - 16.4) / length, 1) if length > 0 else 0

    return {
        "tmNearestNeighbor": tm_nn,
        "tmBasicGC": tm_basic,
        "length": length,
        "gcContent": gc,
        "method": "Nearest-neighbor (SantaLucia 1998) at 50 mM Na+, 250 µM oligo",
        "note": "Estimates only — actual Tm depends on salt, DMSO, and oligo concentration. Validate experimentally.",
    }


# ---------------------------------------------------------------------------
# Sequence complexity / repeats
# ---------------------------------------------------------------------------

def _sequence_complexity(seq: str) -> Dict[str, Any]:
    """
    Detect repetitive elements, low-complexity regions, and polymeric runs.
    Used to flag problematic regions for oligonucleotide design.
    """
    seq_upper = seq.upper()
    length = len(seq)

    # Dinucleotide repeats
    dinuc_repeats = []
    for i in range(length - 3):
        dinuc = seq_upper[i:i+2]
        run = 1
        j = i + 2
        while j <= length - 2 and seq_upper[j:j+2] == dinuc:
            run += 1
            j += 2
        if run >= 3:
            dinuc_repeats.append({
                "pattern": dinuc,
                "start": i + 1,
                "end": i + run * 2,
                "repeats": run,
            })
            if len(dinuc_repeats) >= 20:
                break

    # Trinucleotide repeats
    trinuc_repeats = []
    seen = set()
    for i in range(length - 5):
        trinuc = seq_upper[i:i+3]
        if trinuc in seen:
            continue
        seen.add(trinuc)
        run = seq_upper.count(trinuc)
        if run >= 4:
            positions = [m.start() for m in re.finditer(re.escape(trinuc), seq_upper)]
            trinuc_repeats.append({
                "pattern": trinuc,
                "count": run,
                "positions": positions[:10],
            })
            if len(trinuc_repeats) >= 10:
                break

    # GC-rich / AT-rich stretches (runs of 5+)
    gc_rich = [{"start": m.start() + 1, "end": m.end(), "length": m.end() - m.start()}
               for m in re.finditer(r"[GC]{5,}", seq_upper)][:10]
    at_rich = [{"start": m.start() + 1, "end": m.end(), "length": m.end() - m.start()}
               for m in re.finditer(r"[AT]{5,}", seq_upper)][:10]

    # Self-complementarity regions (simple: reversed complement matches)
    self_comp = []
    is_rna = "U" in seq_upper and "T" not in seq_upper
    rc = _reverse_complement(seq_upper, is_rna)
    for size in [6, 8, 10]:
        for i in range(length - size + 1):
            fragment = seq_upper[i:i+size]
            if fragment in rc and i > size:
                self_comp.append({
                    "sequence": fragment,
                    "position": i + 1,
                    "size": size,
                })
                if len(self_comp) >= 10:
                    break
        if len(self_comp) >= 10:
            break

    return {
        "dinucRepeats": dinuc_repeats,
        "trinucRepeats": trinuc_repeats,
        "gcRichRegions": gc_rich,
        "atRichRegions": at_rich,
        "selfComplementarity": self_comp,
        "complexityScore": round(1 - len(dinuc_repeats) / max(length, 1), 3),
    }


# ---------------------------------------------------------------------------
# Codon usage analysis (mRNA)
# ---------------------------------------------------------------------------

# Human codon usage (relative adaptiveness, simplified from codon usage tables)
_HUMAN_CODON_ADAPT: Dict[str, float] = {
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


def _codon_usage(seq: str) -> Dict[str, Any]:
    """
    Codon usage analysis for mRNA design. Reports the codon adaptation
    index (CAI) and flags rare codons that may cause translational issues.
    """
    seq = seq.upper().replace("T", "U")
    length = len(seq)
    if length < 3:
        return {"codons": [], "cai": 0, "rareCodons": [], "totalCodons": 0}

    codons = []
    rare = []
    adapt_sum = 0
    adapt_count = 0

    # Find longest ORF first
    is_rna = "U" in seq and "T" not in seq
    start_codons = {"AUG"}
    stop_codons = {"UAA", "UAG", "UGA"}

    # Scan forward for coding region
    coding_start = None
    for i in range(length - 2):
        if seq[i:i+3] in start_codons:
            coding_start = i
            break

    if coding_start is not None:
        i = coding_start
        while i < length - 2:
            codon = seq[i:i+3]
            if codon in stop_codons:
                break
            adapt = _HUMAN_CODON_ADAPT.get(codon, 0.5)
            count = seq.count(codon)
            is_rare = adapt < 0.2
            codons.append({
                "codon": codon,
                "position": i + 1,
                "adaptiveness": round(adapt, 2),
                "isRare": is_rare,
            })
            if is_rare:
                rare.append({"codon": codon, "position": i + 1, "adaptiveness": round(adapt, 2)})
            adapt_sum += adapt
            adapt_count += 1
            i += 3

    cai = round(adapt_sum / adapt_count, 3) if adapt_count > 0 else 0

    return {
        "codons": codons[:50],  # cap payload
        "cai": cai,
        "rareCodons": rare[:20],
        "totalCodons": adapt_count,
        "note": "CAI ranges 0-1; higher = more optimized for human expression. Rare codons (<0.2) may cause ribosomal stalling.",
    }


# ---------------------------------------------------------------------------
# Modification scoring
# ---------------------------------------------------------------------------

def _modification_scorecard(seq: str, modality: str) -> Dict[str, Any]:
    """
    Chemistry-specific scoring for the given modality. Not a substitute
    for medicinal chemistry expertise — simplified heuristic flags only.
    """
    gc = _gc_content(seq)
    length = len(seq)

    scores: Dict[str, Any] = {}
    advisories: List[Dict[str, str]] = []

    if modality == "aso":
        # LNA boosting potential
        scores["lnaBoosting"] = {
            "score": max(0, min(100, round((70 - gc) * 1.5))),
            "rationale": "Low GC benefits most from LNA-mediated Tm boost" if gc < 40 else "GC already adequate — fewer LNA substitutions needed",
        }
        # Gapmer design suitability
        scores["gapmerSuitability"] = {
            "score": max(0, min(100, round(80 + (gc - 40) * 0.5 if 30 < gc < 60 else 30))),
            "rationale": "Central gap of DNA flanked by modified wings" if length >= 18 else "Too short for typical gapmer design",
        }
        # Phosphorothioate backbone
        scores["psBackbone"] = {
            "score": max(0, min(100, round(50 + length * 1.5 if length < 25 else 85))),
            "rationale": "PS bonds increase nuclease resistance and protein binding" if length >= 12 else "Very short — PS backbone may not suffice for stability",
        }
        # Uptake
        scores["cellUptake"] = {
            "score": max(0, min(100, round(90 - abs(length - 20) * 3))),
            "rationale": "18-22 nt optimal for cellular uptake of ASOs",
        }

    elif modality == "sirna":
        # Guide strand thermodynamic bias
        first_half_gc = _gc_content(seq[:len(seq)//2]) if len(seq) >= 4 else gc
        second_half_gc = _gc_content(seq[len(seq)//2:]) if len(seq) >= 4 else gc
        bias = abs(first_half_gc - second_half_gc)
        scores["thermodynamicBias"] = {
            "score": max(0, min(100, round(50 + bias * 2))),
            "rationale": f"5'-end less stable (GC {first_half_gc}%) → 3' more stable (GC {second_half_gc}%) favors guide strand loading" if first_half_gc < second_half_gc else "Consider strand polarity — guide strand should have lower 5' stability",
        }
        # RISC loading
        scores["riscLoading"] = {
            "score": max(0, min(100, round(70 + (30 if 30 <= gc <= 52 else -20)))),
            "rationale": "GC 30-52% optimal for RISC loading efficiency" if 30 <= gc <= 52 else "GC outside optimal range for RISC loading",
        }
        # Specificity
        scores["specificity"] = {
            "score": max(0, min(100, round(length * 4 if length <= 25 else 75))),
            "rationale": "19-25 nt length balances potency and specificity" if 19 <= length <= 25 else "Length outside typical siRNA range",
        }

    elif modality == "mrna":
        # Poly-A tail: an observation about the sequence, so it keeps a score.
        has_tail = seq.endswith("A" * 10)
        scores["polyAStability"] = {
            "score": 85 if has_tail else 30,
            "rationale": "Poly-A tail present — contributes to mRNA stability" if has_tail else "No poly-A tail detected — essential for mRNA stability",
        }
        # GC balance for mRNA stability
        scores["mrnaStability"] = {
            "score": max(0, min(100, round(60 + (10 if 40 <= gc <= 60 else -20)))),
            "rationale": "GC 40-60% optimal for mRNA half-life" if 40 <= gc <= 60 else "GC outside optimal range for mRNA stability",
        }
        # `capEfficiency: 70` and `nucleosideMod: 90` used to sit here as
        # scores. Neither looked at the sequence — they returned the same
        # constant for every input, and both describe how the mRNA is
        # MANUFACTURED (cap analogue, m1-pseudouridine substitution), which an
        # uploaded sequence cannot show. Averaged into overallScore they moved
        # every mRNA upload by a fixed amount and made the total look
        # sequence-derived when 40% of it was not. They are advisory notes now,
        # carry no score, and are excluded from the average.
        advisories = [
            {
                "id": "capEfficiency",
                "note": "5' cap required for ribosome recruitment; ARCA or "
                        "CleanCap are the usual choices. Not determinable from "
                        "sequence — it depends on the IVT protocol.",
            },
            {
                "id": "nucleosideMod",
                "note": "m1\u03a8 and m5C substitution reduce innate immune "
                        "activation. Not determinable from sequence — "
                        "modified bases are not represented in an A/C/G/U "
                        "upload.",
            },
        ]

    elif modality == "sgrna":
        # GC content for sgRNA
        scores["gcOptimal"] = {
            "score": max(0, min(100, round(70 + (20 if 40 <= gc <= 80 else -30)))),
            "rationale": "GC 40-80% optimal for sgRNA on-target activity" if 40 <= gc <= 80 else "GC outside optimal range for CRISPR activity",
        }
        # PAM proximity
        scores["pamProximity"] = {
            "score": 80 if "GG" in seq[-5:] else 50,
            "rationale": "NGG PAM motif detected near 3' end" if "GG" in seq[-5:] else "No obvious NGG PAM — verify PAM is adjacent to target",
        }
        # Off-target
        scores["offTargetScore"] = {
            "score": max(0, min(100, round(length * 5 if length <= 20 else 70))),
            "rationale": "20 nt spacer optimal for specificity" if 17 <= length <= 21 else "Length outside standard sgRNA range",
        }

    return {
        "modality": modality,
        "scores": scores,
        # Averaged over sequence-derived scores only. Advisories are listed
        # separately precisely so they cannot drift back into the total.
        "overallScore": round(sum(s["score"] for s in scores.values())
                              / max(len(scores), 1)),
        "advisories": advisories,
        "scoreBasis": (
            "Heuristic flags computed from length, GC and motif content of the "
            "uploaded sequence. Coefficients are rules of thumb, not fitted to "
            "any activity dataset, and overallScore is their unweighted mean — "
            "compare candidates with it, do not read it as a probability of "
            "success."
        ),
    }


# ---------------------------------------------------------------------------
# Base-stacking energy profile
# ---------------------------------------------------------------------------

# Simplified nearest-neighbor free energy (kcal/mol) for DNA at 37°C
_STACK_ENERGY: Dict[str, float] = {
    "AA": -1.0, "TT": -1.0, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}


def _stacking_energy_profile(seq: str, window: int = 10, step: int = 2) -> List[Dict[str, Any]]:
    """
    Per-window average stacking energy across the sequence. Negative
    values indicate more stable (favorable) stacking.
    """
    seq = seq.upper().replace("U", "T")
    length = len(seq)
    if length < window:
        avg = sum(_STACK_ENERGY.get(seq[i:i+2], -1.0) for i in range(length - 1)) / max(length - 1, 1)
        return [{"position": 1, "energy": round(avg, 3)}]

    points = []
    for i in range(0, length - window + 1, step):
        chunk = seq[i:i+window]
        energies = [_STACK_ENERGY.get(chunk[j:j+2], -1.0) for j in range(len(chunk) - 1)]
        avg = sum(energies) / len(energies) if energies else 0
        points.append({"position": i + 1, "energy": round(avg, 3)})
    return points


# ---------------------------------------------------------------------------
# Tm sliding window curve
# ---------------------------------------------------------------------------

def _tm_sliding_window(seq: str, window: int = 10, step: int = 2) -> List[Dict[str, Any]]:
    """Nearest-neighbor Tm at each window position along the sequence."""
    dna = seq.upper().replace("U", "T")
    length = len(dna)
    if length < window:
        tm = _melting_temperature(seq)
        return [{"position": 1, "tm": tm["tmNearestNeighbor"]}]
    points = []
    for i in range(0, length - window + 1, step):
        chunk = dna[i:i + window]
        dH = dS = 0.0
        for j in range(len(chunk) - 1):
            dinuc = chunk[j:j + 2]
            if dinuc in _DNA_NN:
                h, s = _DNA_NN[dinuc]
                dH += h; dS += s
        if dS != 0:
            import math
            R = 1.987; C = 250e-6
            tm = round(dH * 1000 / (dS + R * math.log(C / 4)) - 273.15, 1)
        else:
            tm = 0
        points.append({"position": i + 1, "tm": tm})
    return points


# ---------------------------------------------------------------------------
# Molecular weight & A260 extinction coefficient
# ---------------------------------------------------------------------------

_DNA_MW = {"A": 331.22, "T": 322.21, "G": 347.22, "C": 307.18}
_RNA_MW = {"A": 347.22, "U": 324.18, "G": 363.22, "C": 323.18}
_DNA_EC = {"A": 15400, "T": 8700, "G": 11500, "C": 7400}
_RNA_EC = {"A": 15400, "U": 10000, "G": 11500, "C": 7400}


def _molecular_weight(seq: str) -> Dict[str, Any]:
    """MW and A260 extinction coefficient for bench work."""
    seq_upper = seq.upper()
    is_rna = "U" in seq_upper and "T" not in seq_upper
    mw_table = _RNA_MW if is_rna else _DNA_MW
    ec_table = _RNA_EC if is_rna else _DNA_EC
    base_count = {b: seq_upper.count(b) for b in mw_table}
    base_mw = sum(base_count.get(b, 0) * mw_table[b] for b in mw_table)
    backbone_mw = (len(seq_upper) - 1) * 62.97
    total_mw = round(base_mw + backbone_mw, 1)
    extinction = sum(base_count.get(b, 0) * ec_table[b] for b in ec_table)
    return {
        "molecularWeight": total_mw,
        "molecularWeightKda": round(total_mw / 1000, 2),
        "extinctionCoefficient": extinction,
        "baseCounts": base_count,
        "length": len(seq_upper),
        "type": "RNA" if is_rna else "DNA",
        "note": "MW assumes standard phosphodiester backbone. PS bonds add ~16 Da per substitution.",
    }


# ---------------------------------------------------------------------------
# GC skew
# ---------------------------------------------------------------------------

def _gc_skew(seq: str, window: int = 10, step: int = 2) -> List[Dict[str, Any]]:
    """GC skew = (G − C) / (G + C) in sliding windows."""
    seq_upper = seq.upper()
    length = len(seq_upper)
    if length < window:
        g = seq_upper.count("G"); c = seq_upper.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        return [{"position": 1, "skew": round(skew, 4), "g": g, "c": c}]
    points = []
    for i in range(0, length - window + 1, step):
        chunk = seq_upper[i:i + window]
        g = chunk.count("G"); c = chunk.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        points.append({"position": i + 1, "skew": round(skew, 4), "g": g, "c": c})
    return points


# ---------------------------------------------------------------------------
# Dinucleotide frequency
# ---------------------------------------------------------------------------

_BASES_DNA = ["A", "C", "G", "T"]
_BASES_RNA = ["A", "C", "G", "U"]


def _dinucleotide_frequency(seq: str) -> Dict[str, Any]:
    """Count all 16 dinucleotide combinations for heatmap rendering."""
    seq_upper = seq.upper()
    is_rna = "U" in seq_upper and "T" not in seq_upper
    seq_upper = seq_upper.replace("T", "U") if not is_rna else seq_upper
    bases = _BASES_RNA if is_rna else ["A", "C", "G", "U"]
    if not is_rna:
        seq_upper = seq_upper.replace("U", "T")
        bases = _BASES_DNA

    counts = {}
    for b1 in bases:
        for b2 in bases:
            counts[f"{b1}{b2}"] = 0
    for i in range(len(seq_upper) - 1):
        dinuc = seq_upper[i:i + 2]
        if dinuc in counts:
            counts[dinuc] += 1
    matrix = [[counts[f"{b1}{b2}"] for b2 in bases] for b1 in bases]
    total = sum(counts.values())
    freq = {k: round(v / total, 4) if total > 0 else 0 for k, v in counts.items()}
    return {"counts": counts, "frequency": freq, "matrix": matrix, "bases": bases, "total": total}


# ---------------------------------------------------------------------------
# Self-complementarity dot plot
# ---------------------------------------------------------------------------

def _dotplot_data(seq: str, min_stem: int = 4) -> Dict[str, Any]:
    """Dot plot data: positions where complementary stretches exist."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
    seq_upper = seq.upper()
    max_len = min(len(seq_upper), 200)
    dots = []
    for i in range(max_len):
        for j in range(i + min_stem, max_len):
            match = True
            for k in range(min_stem):
                if i + k >= max_len or j + k >= max_len:
                    match = False; break
                if comp.get(seq_upper[i + k], "?") != seq_upper[j + k]:
                    match = False; break
            if match:
                dots.append({"i": i + 1, "j": j + 1, "len": min_stem})
                if len(dots) >= 500:
                    break
        if len(dots) >= 500:
            break
    return {"dots": dots, "length": max_len, "minStem": min_stem, "totalDots": len(dots)}


# ---------------------------------------------------------------------------
# Restriction enzyme sites
# ---------------------------------------------------------------------------

_RESTRICTION_ENZYMES = [
    {"name": "EcoRI", "pattern": "GAATTC", "cutAfter": 1, "overhang": "5'"},
    {"name": "BamHI", "pattern": "GGATCC", "cutAfter": 1, "overhang": "5'"},
    {"name": "HindIII", "pattern": "AAGCTT", "cutAfter": 1, "overhang": "5'"},
    {"name": "NotI", "pattern": "GCGGCCGC", "cutAfter": 2, "overhang": "5'"},
    {"name": "XhoI", "pattern": "CTCGAG", "cutAfter": 1, "overhang": "5'"},
    {"name": "SpeI", "pattern": "ACTAGT", "cutAfter": 1, "overhang": "5'"},
    {"name": "NheI", "pattern": "GCTAGC", "cutAfter": 1, "overhang": "5'"},
    {"name": "XbaI", "pattern": "TCTAGA", "cutAfter": 1, "overhang": "5'"},
    {"name": "KpnI", "pattern": "GGTACC", "cutAfter": 5, "overhang": "3'"},
    {"name": "PstI", "pattern": "CTGCAG", "cutAfter": 5, "overhang": "3'"},
    {"name": "SmaI", "pattern": "CCCGGG", "cutAfter": 3, "overhang": "blunt"},
    {"name": "MluI", "pattern": "ACGCGT", "cutAfter": 1, "overhang": "5'"},
]


def _restriction_sites(seq: str) -> List[Dict[str, Any]]:
    """Find recognition sites for common restriction enzymes."""
    seq_upper = seq.upper().replace("U", "T")
    found = []
    for enz in _RESTRICTION_ENZYMES:
        for m in re.finditer(re.escape(enz["pattern"]), seq_upper):
            found.append({
                "enzyme": enz["name"], "pattern": enz["pattern"],
                "start": m.start() + 1, "end": m.end(),
                "cutAfter": enz["cutAfter"], "overhang": enz["overhang"],
            })
    return found


# ---------------------------------------------------------------------------
# miRNA seed-match screening
# ---------------------------------------------------------------------------

_HUMAN_MIRNA_SEEDS = [
    {"name": "hsa-let-7a", "seed": "UGAGGUAG"},
    {"name": "hsa-miR-21", "seed": "AGCUUAUC"},
    {"name": "hsa-miR-155", "seed": "UUAAUGCU"},
    {"name": "hsa-miR-122", "seed": "GGAGUGUG"},
    {"name": "hsa-miR-34a", "seed": "GGCAGUGU"},
    {"name": "hsa-miR-17", "seed": "CAAAGUGC"},
    {"name": "hsa-miR-200a", "seed": "UAACACUG"},
    {"name": "hsa-miR-15a", "seed": "AGCAGCAC"},
    {"name": "hsa-miR-16", "seed": "AGCAGCAC"},
    {"name": "hsa-miR-223", "seed": "UGUCAGUU"},
    {"name": "hsa-miR-146a", "seed": "GAGAACUG"},
    {"name": "hsa-miR-125b", "seed": "UCACAAGU"},
    {"name": "hsa-miR-let-7b", "seed": "UGAGGUAG"},
]


def _mirna_seed_matches(seq: str) -> List[Dict[str, Any]]:
    """Check for complementary matches to known miRNA seed regions."""
    seq_upper = seq.upper().replace("T", "U")
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    matches = []
    for mirna in _HUMAN_MIRNA_SEEDS:
        seed = mirna["seed"]
        for i in range(len(seq_upper) - 5):
            fragment = seq_upper[i:i + min(8, len(seq_upper) - i)]
            frag_rc = "".join(comp.get(b, "N") for b in reversed(fragment))
            if len(frag_rc) >= 6 and seed[:len(frag_rc)] == frag_rc:
                matches.append({
                    "mirna": mirna["name"], "seed": seed,
                    "matchPosition": i + 1, "matchLength": len(frag_rc),
                    "matchSequence": fragment,
                })
                break
        if len(matches) >= 15:
            break
    return matches


# ---------------------------------------------------------------------------
# Approved drug comparison
# ---------------------------------------------------------------------------

_APPROVED_ASO_DRUGS = [
    {"name": "Nusinersen (Spinraza)", "length": 18, "gc": 55.6, "modality": "aso", "target": "SMN2"},
    {"name": "Eteplirsen (Exondys 51)", "length": 30, "gc": 46.7, "modality": "aso", "target": "DMD exon 51"},
    {"name": "Inotersen (Tegsedi)", "length": 18, "gc": 50.0, "modality": "aso", "target": "TTR"},
    {"name": "Tofersen (Qalsody)", "length": 20, "gc": 55.0, "modality": "aso", "target": "SOD1"},
    {"name": "Milasen (custom)", "length": 19, "gc": 47.4, "modality": "aso", "target": "Batten CLN7"},
    {"name": "Viltolarsen (Viltepso)", "length": 23, "gc": 43.5, "modality": "aso", "target": "DMD exon 53"},
    {"name": "Fomivirsen (Vitravene)", "length": 21, "gc": 47.6, "modality": "aso", "target": "CMV retinitis"},
]
_APPROVED_SIRNA_DRUGS = [
    {"name": "Patisiran (Onpattro)", "length": 21, "gc": 47.6, "modality": "sirna", "target": "TTR"},
    {"name": "Givosiran (Givlaari)", "length": 21, "gc": 42.9, "modality": "sirna", "target": "ALAS1"},
    {"name": "Lumasiran (Oxlumo)", "length": 21, "gc": 52.4, "modality": "sirna", "target": "HAO1"},
    {"name": "Inclisiran (Leqvio)", "length": 21, "gc": 47.6, "modality": "sirna", "target": "PCSK9"},
    {"name": "Vutrisiran (Amvuttra)", "length": 21, "gc": 47.6, "modality": "sirna", "target": "TTR"},
    {"name": "Nedosiran (Rivfloza)", "length": 21, "gc": 42.9, "modality": "sirna", "target": "LDHA"},
]
_APPROVED_SGRNA_DRUGS = [
    {"name": "Casgevy (exa-cel)", "length": 20, "gc": 55.0, "modality": "sgrna", "target": "BCL11A enhancer"},
]


def _approved_drug_comparison(seq: str, modality: str) -> Dict[str, Any]:
    """Compare sequence against real FDA-approved drugs of same modality."""
    gc = _gc_content(seq)
    length = len(seq)
    drugs = {"aso": _APPROVED_ASO_DRUGS, "sirna": _APPROVED_SIRNA_DRUGS, "sgrna": _APPROVED_SGRNA_DRUGS}.get(modality, [])
    my_data = {"name": "Your sequence", "length": length, "gc": gc, "modality": modality, "target": "—", "isUser": True}
    return {
        "drugs": drugs, "yourSequence": my_data,
        "lengthRange": {"min": min(d["length"] for d in drugs) if drugs else 0, "max": max(d["length"] for d in drugs) if drugs else 0},
        "gcRange": {"min": min(d["gc"] for d in drugs) if drugs else 0, "max": max(d["gc"] for d in drugs) if drugs else 0},
    }


# ---------------------------------------------------------------------------
# Design-rule checklist
# ---------------------------------------------------------------------------

def _design_rule_checklist(seq: str, modality: str) -> Dict[str, Any]:
    """Pass/fail checklist against modality-specific numeric design rules."""
    gc = _gc_content(seq)
    length = len(seq)
    rules = []
    if modality == "aso":
        rules = [
            {"rule": "Length 18–25 nt", "pass": 18 <= length <= 25, "actual": f"{length} nt"},
            {"rule": "GC content 30–70%", "pass": 30 <= gc <= 70, "actual": f"{gc}%"},
            {"rule": "No poly-G ≥ 4", "pass": "GGGG" not in seq.upper(), "actual": "Clear" if "GGGG" not in seq.upper() else "GGGG found"},
            {"rule": "No poly-C ≥ 6", "pass": "CCCCCC" not in seq.upper(), "actual": "Clear" if "CCCCCC" not in seq.upper() else "CCCCCC found"},
            {"rule": "No palindrome ≥ 8 nt", "pass": not _has_long_palindrome(seq, 8), "actual": "Clear" if not _has_long_palindrome(seq, 8) else "Found"},
        ]
    elif modality == "sirna":
        rules = [
            {"rule": "Length 19–25 nt", "pass": 19 <= length <= 25, "actual": f"{length} nt"},
            {"rule": "GC content 30–52%", "pass": 30 <= gc <= 52, "actual": f"{gc}%"},
            {"rule": "No poly-U ≥ 4 at 3' end", "pass": not seq.upper().endswith("UUUU"), "actual": "Clear" if not seq.upper().endswith("UUUU") else "Poly-U at 3'"},
        ]
    elif modality == "mrna":
        orfs = _find_orfs(seq)
        rules = [
            {"rule": "Contains ORF", "pass": len(orfs) > 0, "actual": f"{len(orfs)} found"},
            {"rule": "Poly-A tail present", "pass": _has_poly_a_tail(seq), "actual": "Yes" if _has_poly_a_tail(seq) else "No"},
            {"rule": "GC content 40–60%", "pass": 40 <= gc <= 60, "actual": f"{gc}%"},
            {"rule": "No poly-G ≥ 4", "pass": "GGGG" not in seq.upper(), "actual": "Clear" if "GGGG" not in seq.upper() else "GGGG found"},
        ]
    elif modality == "sgrna":
        rules = [
            {"rule": "Length 17–21 nt", "pass": 17 <= length <= 21, "actual": f"{length} nt"},
            {"rule": "GC content 40–80%", "pass": 40 <= gc <= 80, "actual": f"{gc}%"},
            {"rule": "NGG PAM at 3' end", "pass": bool(re.search(r"GG", seq.upper()[-5:])), "actual": "Yes" if bool(re.search(r"GG", seq.upper()[-5:])) else "Not detected"},
            {"rule": "No poly-T ≥ 4", "pass": "TTTT" not in seq.upper(), "actual": "Clear" if "TTTT" not in seq.upper() else "TTTT found"},
        ]
    passed = sum(1 for r in rules if r["pass"])
    return {"rules": rules, "passed": passed, "total": len(rules), "score": round(passed / len(rules) * 100) if rules else 0}


def _has_long_palindrome(seq: str, min_len: int = 8) -> bool:
    seq_upper = seq.upper()
    for i in range(len(seq_upper) - min_len + 1):
        chunk = seq_upper[i:i + min_len]
        if chunk == chunk[::-1]:
            return True
    return False


# ---------------------------------------------------------------------------
# siRNA thermodynamic asymmetry
# ---------------------------------------------------------------------------

def _sirna_asymmetry(seq: str) -> Dict[str, Any]:
    """Schwarz/Khvorova rule: 5' end of guide should be less stable."""
    seq_upper = seq.upper().replace("T", "U")
    length = len(seq_upper)
    half = length // 2
    if half < 5:
        return {"asymmetric": False, "firstHalfEnergy": 0, "secondHalfEnergy": 0, "note": "Too short"}
    first_half = seq_upper[:half].replace("U", "T")
    second_half = seq_upper[half:].replace("U", "T")
    def terminal_energy(subseq: str) -> float:
        e = 0.0
        for i in range(min(5, len(subseq)) - 1):
            dinuc = subseq[i:i + 2]
            if dinuc in _DNA_NN:
                e += _DNA_NN[dinuc][0]
        return round(e, 3)
    fe = terminal_energy(first_half)
    se = terminal_energy(second_half)
    return {
        "asymmetric": fe > se,
        "firstHalfEnergy": fe, "secondHalfEnergy": se,
        "firstHalfLength": half, "secondHalfLength": length - half,
        "note": "5' terminal stability of guide should be LOWER (less negative ΔG) than 3' end — favors correct RISC loading (Schwarz/Khvorova rule).",
    }


# ---------------------------------------------------------------------------
# Reverse complement view
# ---------------------------------------------------------------------------

def _reverse_complement_data(seq: str) -> Dict[str, Any]:
    """Sequence and reverse complement for side-by-side display."""
    seq_upper = seq.upper()
    is_rna = "U" in seq_upper and "T" not in seq_upper
    rc = _reverse_complement(seq_upper, is_rna)
    matches = [seq_upper[i] == rc[i] for i in range(min(len(seq_upper), len(rc)))]
    return {
        "sequence": seq_upper, "reverseComplement": rc, "matches": matches,
        "identity": round(sum(matches) / max(len(matches), 1) * 100, 1),
        "type": "RNA" if is_rna else "DNA",
    }
