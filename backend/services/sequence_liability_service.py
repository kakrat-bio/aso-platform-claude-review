"""Sequence-determined liabilities for an oligonucleotide candidate.

WHAT THIS MODULE USED TO BE, AND WHY IT ISN'T ANY MORE.

It was an "ADMET prediction" service: absorption, distribution, metabolism,
excretion and toxicity scores, a one-compartment PBPK concentration-time
curve, a Lipinski rule-of-five violation count, a charge-versus-pH profile, a
2-D "chemical space projection", and a hemolysis risk — all computed from the
bare nucleotide sequence.

Almost none of that is a property of the sequence. For an oligonucleotide:

* **Absorption, distribution and cellular uptake** are set by the backbone
  and the conjugate, not the bases. A phosphorothioate binds plasma proteins
  and is taken up broadly; a charge-neutral morpholino of the identical
  sequence is not. A GalNAc conjugate goes to hepatocytes regardless of what
  it spells. `get_admet_prediction` was never even passed the chemistry —
  `api/gene_silencing.py` called it with the sequence alone — so these scores
  could not have reflected the molecule being designed.
* **Metabolism** for this class is nuclease degradation, which is governed by
  2'-modifications and backbone linkages. An unmodified-RNA estimate is
  meaningless for a 5-10-5 MOE gapmer.
* **Excretion / half-life.** `_predict_pbpk_curve` defaulted to a 2-hour
  half-life "for unmodified RNA". A PS/2'-MOE gapmer has a tissue half-life
  of weeks. Every curve it drew was wrong by two to three orders of
  magnitude.
* **Lipinski's rule of five** describes orally absorbed small molecules. A
  6-8 kDa polyanion violates every criterion by construction; the checks that
  ran under that name (length, GC, CpG count, poly-G) are ordinary oligo
  design flags borrowing an authority they do not have.
* **Charge versus pH.** The phosphodiester backbone is fully ionised from
  pH 4 to pH 8; nucleobase pKa values sit far outside the physiological
  range. The "profile" was a flat line by construction.
* **Hemolysis** is a property of the formulation and delivery vehicle, not of
  a naked sequence.
* **`_project_to_2d`** presented hand-picked weights (0.40, 0.30, 0.20, 0.10)
  as a PCA-like projection. There was no chemical space and no fitted model.

What a sequence genuinely does determine is innate-immune recognition and a
handful of structural liabilities, and those are what remains. Every removed
endpoint is listed in the response under `notAssessed` with its reason, so a
caller sees why a field is gone rather than finding it silently absent.

Nothing in the repository consumed the old ADMET block — no frontend file
references it — so the removal changes no rendered output.
"""

import logging
import math
import re
from typing import Dict, Any, Optional, List

from services.sequence_metrics import gc_content as _compute_gc_content

logger = logging.getLogger(__name__)


def _classify_score(value: float, low_thresh: float, high_thresh: float) -> str:
    if value >= high_thresh:
        return "High"
    elif value >= low_thresh:
        return "Moderate"
    else:
        return "Low"


def _compute_sequence_descriptors(seq: str) -> Dict[str, Any]:
    """Compute physicochemical descriptors from RNA sequence."""
    if not seq:
        return {}
    seq = seq.upper()
    length = len(seq)
    gc = _compute_gc_content(seq)
    g_count = seq.count("G")
    c_count = seq.count("C")
    a_count = seq.count("A")
    u_count = seq.count("U")
    gc_skew = (g_count - c_count) / length if length else 0
    purines = g_count + a_count
    pyrimidines = c_count + u_count
    
    # Wallace rule melting temperature approximation
    tm = 2 * (a_count + u_count) + 4 * (g_count + c_count)
    
    # Motif counts
    cpg = len(re.findall(r"CG", seq))
    polyg = len(re.findall(r"G{4,}", seq))
    polyu = len(re.findall(r"U{4,}", seq))
    
    # Shannon entropy
    from collections import Counter
    counts = Counter(seq)
    entropy = 0.0
    for base, cnt in counts.items():
        p = cnt / length
        if p > 0:
            entropy -= p * math.log2(p)
    
    # Charge density (phosphates per residue, excluding 5' terminal)
    charge_density = (length - 1) / length if length > 1 else 0.0
    
    return {
        "length": length,
        "gcContent": gc,
        "gcSkew": round(gc_skew, 3),
        "atContent": round(100 * (a_count + u_count) / length, 1) if length else 0,
        "purineFraction": round(purines / length, 3) if length else 0,
        "meltingTemp": round(tm, 1),
        "cpgCount": cpg,
        "polyGCount": polyg,
        "polyUCount": polyu,
        "sequenceEntropy": round(entropy, 3),
        "chargeDensity": round(charge_density, 3),
    }


def _immune_and_structural_liabilities(seq: str) -> Dict[str, Any]:
    """Innate-immune and structural liabilities that ARE sequence-encoded.

    Each flag below is a documented, sequence-determined property of an
    oligonucleotide, and each carries the reason it is being raised:

    * **Unmethylated CpG -> TLR9.** Endosomal TLR9 recognises unmethylated CpG
      dinucleotides; this is the entire basis of CpG-oligodeoxynucleotide
      immunostimulation and a recognised liability of phosphorothioate ASOs.
      It depends on the bases, so a sequence can be screened for it.
    * **G-quadruplex (four or more consecutive G).** Runs of guanosine fold
      into quadruplexes, which drive aggregation, non-specific protein binding
      and poor synthesis yield. Purely a consequence of base order.
    * **Uridine-rich tracts -> TLR7/8.** Single-stranded uridine-rich RNA is
      the canonical TLR7/TLR8 ligand. Note this one is chemistry-modifiable:
      2'-modification or N1-methylpseudouridine substitution largely abolishes
      the response, so the flag describes the sequence, not the finished drug.

    No composite score is returned. The individual motif calls are real; the
    weights that used to be summed into a single 0-1 "immunogenicity score"
    (0.4 for CpG, 0.3 for poly-G, 0.25 for poly-U ...) were invented, and a
    number carries more authority than a flag list deserves.
    """
    if not seq:
        return {"flags": [], "note": "No sequence supplied."}

    seq = seq.upper().replace("T", "U")
    flags: List[Dict[str, str]] = []

    cpg_count = len(re.findall(r"CG", seq))
    if cpg_count:
        flags.append({
            "id": "cpg_tlr9",
            "observation": f"{cpg_count} CpG dinucleotide(s)",
            "liability": "TLR9 recognition",
            "reasoning": (
                "Endosomal TLR9 binds unmethylated CpG motifs. A synthetic "
                "oligonucleotide is unmethylated by default, so each CpG is a "
                "potential innate-immune trigger. Reduce the count or place "
                "5-methylcytosine at these positions."
            ),
        })

    g4 = re.findall(r"G{4,}", seq)
    if g4:
        flags.append({
            "id": "g_quadruplex",
            "observation": f"{len(g4)} run(s) of 4+ consecutive G: {', '.join(g4)}",
            "liability": "G-quadruplex formation",
            "reasoning": (
                "Four or more consecutive guanosines can fold into a "
                "G-quadruplex, which causes aggregation, non-specific protein "
                "binding and unreliable solid-phase synthesis. Break the run."
            ),
        })

    u_tracts = re.findall(r"U{4,}", seq)
    if u_tracts:
        flags.append({
            "id": "uridine_tract_tlr78",
            "observation": f"{len(u_tracts)} uridine tract(s) of 4+ nt",
            "liability": "TLR7/TLR8 recognition",
            "reasoning": (
                "Uridine-rich single-stranded RNA is the canonical TLR7/8 "
                "ligand. This flag describes the base sequence; 2'-sugar "
                "modification or N1-methylpseudouridine substitution largely "
                "abolishes the response, so it may not survive into the "
                "finished chemistry."
            ),
        })

    return {
        "flags": flags,
        "note": (
            "Sequence-encoded liabilities only. No composite score is given: "
            "the individual calls are documented biology, but any weighting "
            "that combined them would be invented."
        ),
    }


def _parse_loeuf_decile(value: Any) -> Optional[int]:
    """Extract numeric LOEUF decile (1-10) from a gnomAD decile string."""
    if value is None:
        return None
    m = re.match(r"Decile\s*(\d+)", str(value).strip())
    if not m:
        return None
    decile = int(m.group(1))
    return decile if 1 <= decile <= 10 else None


def _assess_gene_context_safety(gene_context: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Derive target-gene-driven safety metrics.

    These metrics describe the consequences of silencing the target gene itself
    (on-target pharmacology), not the chemistry of the molecule:

      - onTargetToxicityRisk: risk of dose-limiting on-target toxicity, driven by
        gene essentiality, gnomAD constraint (LOEUF decile) and expression in
        vital organs. Haploinsufficient/essential targets with high vital-organ
        expression leave little room between efficacy and toxicity.
      - therapeuticWindow: relative margin between effective and toxic dose,
        inferred from the same gene-level signals.
      - distributionNotes: tissue-exposure notes for biodistribution monitoring.
    """
    empty = {
        "onTargetToxicityRisk": 0.0,
        "onTargetToxicityLevel": None,
        "therapeuticWindow": {"level": None, "notes": []},
        "distributionNotes": [],
        "onTargetWarning": None,
        "onTargetStrength": None,
    }
    if not gene_context:
        return empty

    essential_raw = gene_context.get("essentialGene")
    if essential_raw is None:
        essential = None
    else:
        essential_str = str(essential_raw).strip().lower()
        if essential_str in ("essential", "yes", "true", "1"):
            essential = True
        elif essential_str in ("non-essential", "nonessential", "no", "false", "0"):
            essential = False
        else:
            essential = None

    loeuf = _parse_loeuf_decile(gene_context.get("loeufDecile"))

    vital_tpm = gene_context.get("vitalOrganTpm")
    if vital_tpm is not None:
        try:
            vital_tpm = float(vital_tpm)
        except (TypeError, ValueError):
            vital_tpm = None
    vital_tissues = gene_context.get("vitalOrganTissues") or []

    components = []
    labels = []

    if loeuf is not None:
        # Decile 1 = most constrained (haploinsufficiency-like) -> highest risk
        constraint_risk = round(max(0.15, 1.0 - (loeuf - 1) * 0.1), 2)
        components.append(constraint_risk)
        constrained = "most" if loeuf <= 3 else "moderately" if loeuf <= 6 else "least"
        labels.append(f"LOEUF decile {loeuf} ({constrained} constrained)")

    if essential is not None:
        components.append(0.85 if essential else 0.25)
        labels.append(f"{'essential' if essential else 'non-essential'} in dependency screens")

    if vital_tpm is not None and vital_tpm >= 0:
        if vital_tpm >= 50:
            expr_risk = 0.9
        elif vital_tpm >= 20:
            expr_risk = 0.6
        elif vital_tpm >= 5:
            expr_risk = 0.35
        else:
            expr_risk = 0.15
        components.append(expr_risk)
        labels.append(f"vital-organ expression TPM {vital_tpm:g}")

    if not components:
        return empty

    on_target = round(sum(components) / len(components), 2)
    if on_target >= 0.6:
        level = "High"
    elif on_target >= 0.35:
        level = "Moderate"
    else:
        level = "Low"

    notes = [f"On-target risk inferred from {', '.join(labels)}."]
    if on_target >= 0.35:
        notes.append(
            "Elevated on-target risk implies a narrow therapeutic window — partial knockdown or "
            "cell-type-restricted delivery may be needed to separate efficacy from toxicity."
        )

    if on_target >= 0.6:
        tw = "Narrow"
        tw_notes = [
            "Target is essential/constrained or highly expressed in vital organs — expect a narrow safety margin.",
            "De-risk with dose titration, knockdown-depth studies, and tissue-specific delivery.",
        ]
    elif on_target >= 0.35:
        tw = "Moderate"
        tw_notes = [
            "Moderate therapeutic window — balance knockdown depth against on-target effects.",
        ]
    else:
        tw = "Wide"
        tw_notes = [
            "Low predicted on-target toxicity — a wide therapeutic window is expected for silencing this target.",
        ]

    dist_notes = []
    if vital_tissues and vital_tpm is not None and vital_tpm > 0:
        organs = ", ".join(str(t) for t in vital_tissues[:4])
        dist_notes.append(
            f"Predominant vital-organ expression in {organs} (max TPM {vital_tpm:g}) — "
            "monitor tissue exposure and consider targeted delivery to limit on-target exposure."
        )

    if on_target >= 0.6:
        on_target_warning = (
            "High on-target toxicity risk — essential/constrained target with high "
            "vital-organ expression leaves a narrow therapeutic window"
        )
    elif on_target >= 0.35:
        on_target_warning = (
            "Moderate on-target toxicity risk — balance knockdown depth against "
            "target essentiality and vital-organ expression"
        )
    else:
        on_target_warning = None
    on_target_strength = (
        "Low on-target toxicity risk — wide therapeutic window for silencing this target"
        if 0 < on_target < 0.35
        else None
    )

    return {
        "onTargetToxicityRisk": on_target,
        "onTargetToxicityLevel": level,
        "therapeuticWindow": {"level": tw, "notes": tw_notes},
        "distributionNotes": dist_notes,
        "onTargetWarning": on_target_warning,
        "onTargetStrength": on_target_strength,
    }


# Endpoints deliberately not reported, and why. Kept in the response so a
# caller sees the reason rather than an unexplained absence.
NOT_ASSESSED: Dict[str, str] = {
    "absorption": (
        "Determined by backbone chemistry and conjugation, not by sequence. "
        "A phosphorothioate is taken up broadly through protein binding; a "
        "charge-neutral morpholino of the same sequence is not."
    ),
    "distribution": (
        "Driven by plasma-protein binding (a phosphorothioate property) and "
        "by any targeting conjugate such as GalNAc. The base sequence "
        "contributes marginally."
    ),
    "metabolism": (
        "For this class metabolism is nuclease degradation, governed by "
        "2'-modifications and backbone linkages rather than base order."
    ),
    "excretion": (
        "Renal clearance depends on molecular size and on protein binding "
        "conferred by the backbone. Length alone is too weak a proxy to "
        "report as a prediction."
    ),
    "halfLife": (
        "Requires the chemistry and a calibrated PK model. The previous "
        "1-compartment curve assumed a 2-hour half-life for unmodified RNA; "
        "a PS/2'-MOE gapmer persists in tissue for weeks."
    ),
    "lipinskiViolations": (
        "Lipinski's rule of five describes orally absorbed small molecules. "
        "A 6-8 kDa polyanion violates it by construction, so the count "
        "carries no information about this molecule class."
    ),
    "chargePhProfile": (
        "The phosphodiester backbone is fully ionised across the whole "
        "physiological pH range, so the profile is flat by construction."
    ),
    "hemolysisRisk": (
        "A property of the formulation and delivery vehicle, not of a naked "
        "sequence."
    ),
    "chemicalSpaceProjection": (
        "The previous projection used hand-picked weights presented as a "
        "PCA-like embedding. There is no fitted model behind it."
    ),
    "offTargetHybridisation": (
        "Not screened. No genome or transcriptome alignment is performed "
        "anywhere in this service; run BLAST or a dedicated off-target tool "
        "against the relevant transcriptome."
    ),
}


def get_sequence_liabilities(
    aso_sequence: Optional[str] = None,
    chemistry: Optional[str] = None,
    gene_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Sequence-determined liabilities, plus target-gene pharmacology.

    Replaces `get_admet_prediction`. Returns only what the inputs actually
    determine:

    * `sequenceDescriptors` - direct measurements of the sequence.
    * `immuneAndStructural` - CpG/TLR9, G-quadruplex and uridine-tract flags,
      each with its biological reasoning.
    * `onTargetPharmacology` - consequences of modulating the TARGET GENE
      (gnomAD LOEUF constraint, essentiality, vital-organ expression). This is
      about the gene, not the molecule, and appears only when `gene_context`
      is supplied.
    * `notAssessed` - every endpoint deliberately not reported, with why.

    `chemistry` is accepted and echoed so a caller can see which chemistry the
    candidate carries; it is NOT turned into a pharmacokinetic number, because
    no calibrated model here maps this platform's coarse chemistry labels onto
    absorption, clearance or half-life.
    """
    if not aso_sequence:
        return {
            "available": False,
            "reason": "No sequence supplied.",
            "notAssessed": NOT_ASSESSED,
        }

    descriptors = _compute_sequence_descriptors(aso_sequence.upper().replace("T", "U"))
    out: Dict[str, Any] = {
        "available": True,
        "chemistry": chemistry,
        "sequenceDescriptors": descriptors,
        "immuneAndStructural": _immune_and_structural_liabilities(aso_sequence),
        "notAssessed": NOT_ASSESSED,
        "scope": (
            "Sequence-determined properties only. Absorption, distribution, "
            "metabolism, excretion and half-life are set by backbone "
            "chemistry and conjugation and are not predicted here."
        ),
    }

    gene_safety = _assess_gene_context_safety(gene_context)
    if gene_context:
        out["onTargetPharmacology"] = {
            **gene_safety,
            "note": (
                "Consequences of modulating the target gene, not properties "
                "of the oligonucleotide. Derived from gnomAD constraint, "
                "essentiality and vital-organ expression."
            ),
        }
    return out
