"""A21 — siRNA duplex design.

A21 was the platform's one SCORED-but-DESIGN-UNAVAILABLE mechanism: it
competes for ranking, has five approved drugs behind it, and returned no
candidates because the ASO designer emits single strands and an siRNA is a
duplex. This module supplies the missing stage.

WHAT AN siRNA DESIGN ACTUALLY REQUIRES, and why the single-strand designer
could not be reused:

* **Two strands.** The guide (antisense) is the reverse complement of the
  target; the passenger (sense) matches the target. Both are emitted.
* **3' dinucleotide overhangs.** Both strands carry a 2-nt 3' overhang —
  the Tuschl architecture that Dicer products have and that RISC expects.
* **Thermodynamic asymmetry.** This is the part with no analogue in ASO
  design and the reason a good ASO site can be a bad siRNA site. Argonaute
  loads whichever strand has the less stably paired 5' end (Khvorova 2003,
  Schwarz 2003). If the passenger's 5' end is the looser one, the wrong
  strand is loaded and the design silences the wrong thing. The duplex-end
  stability difference is computed here from nearest-neighbour parameters
  and reported per candidate, and candidates are ranked by it.
* **Seed region.** Guide positions 2-8 drive miRNA-like off-target
  silencing. The seed is reported so it can be checked; its off-target load
  is NOT scored, because scoring it needs a transcriptome-wide seed match
  count and nothing here aligns against a transcriptome.

Sequences and thermodynamics are real: the target windows come from the
Ensembl transcript, duplex free energy from ViennaRNA, and the asymmetry
term from the SantaLucia nearest-neighbour parameters via primer3.
"""

from __future__ import annotations

import logging
from typing import Any

import primer3
import RNA

from services.gene_silencing_service import get_target_analysis

logger = logging.getLogger(__name__)

# Tuschl-style architecture: a 19-nt duplex region plus 2-nt 3' overhangs on
# both strands, giving 21-mers.
DUPLEX_CORE_NT = 19
OVERHANG_NT = 2
DEFAULT_OVERHANG = "dTdT"
SEED_START, SEED_END = 2, 8          # guide positions, 1-based inclusive
# How many terminal base pairs enter the 5'-end stability comparison.
ASYMMETRY_WINDOW_NT = 4
# Ago2 cleaves the target between the bases opposite guide positions 10 and 11.
AGO2_CLEAVAGE_GUIDE_POSITION = 10
# Ui-Tei's AU-rich window: guide positions 1-7.
AU_RICH_WINDOW = 7
# A GC run this long or longer impedes duplex unwinding.
MAX_GC_RUN = 9

_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


def _revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT[b] for b in reversed(seq.upper()))


def _end_stability(seq: str) -> float | None:
    """Nearest-neighbour dG of a duplex end, in kcal/mol (more negative = tighter)."""
    dna = seq.upper().replace("U", "T")
    if len(dna) < 2 or set(dna) - set("ACGT"):
        return None
    # calc_heterodimer against the perfect complement gives the duplex dG for
    # this short terminal stretch.
    comp = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    try:
        return round(primer3.calc_heterodimer(dna, comp).dg / 1000.0, 3)
    except Exception:
        return None


def _ui_tei_rules(guide: str, passenger: str) -> dict[str, Any]:
    """The four Ui-Tei asymmetry rules, checked position by position.

    Bulk duplex-end thermodynamics is only half the asymmetry story. Ui-Tei
    (2004) and Amarzguioui (2004) found that WHICH base sits at particular
    positions predicts strand selection better than the free-energy difference
    alone, because Argonaute's MID pocket reads the guide's 5' nucleotide
    directly rather than inferring it from stability:

    i.   A or U at guide position 1 (5' end). The MID domain binds a
         5'-monophosphate with an adjacent A/U far better than G/C.
    ii.  G or C at guide position 19 (the 3' end, which pairs with the
         passenger's 5' end). This makes the passenger end the tight one.
    iii. At least 4 of guide positions 1-7 are A or U, extending the loose
         5' end beyond the terminal base pair.
    iv.  No GC run of 9 nt or longer anywhere in the duplex, which would
         resist helicase unwinding.

    Returned as individually named pass/fail calls with the reason, not as a
    single score. Each rule is a documented observation; a weighting that
    combined them into one number would not be.
    """
    checks = []
    p1 = guide[0]
    checks.append({
        "rule": "guide_position_1_is_A_or_U",
        "passes": p1 in "AU",
        "observed": p1,
        "reasoning": ("Argonaute's MID pocket accommodates an A or U at the "
                      "guide 5' end far better than a G or C, and this is the "
                      "single strongest strand-selection determinant."),
    })
    p19 = guide[-1]
    checks.append({
        "rule": "guide_position_19_is_G_or_C",
        "passes": p19 in "GC",
        "observed": p19,
        "reasoning": ("The guide's 3' end pairs with the passenger's 5' end. "
                      "A G or C here makes that end the tightly bound one, so "
                      "the passenger is not the strand handed to RISC."),
    })
    au = sum(1 for b in guide[:AU_RICH_WINDOW] if b in "AU")
    checks.append({
        "rule": "guide_positions_1_to_7_AU_rich",
        "passes": au >= 4,
        "observed": f"{au}/7 A or U",
        "reasoning": ("An AU-rich 5' third keeps the whole guide 5' region "
                      "loosely paired, not just the terminal base pair."),
    })
    longest_gc = 0
    run = 0
    for b in passenger:
        run = run + 1 if b in "GC" else 0
        longest_gc = max(longest_gc, run)
    checks.append({
        "rule": "no_GC_run_of_9_or_more",
        "passes": longest_gc < MAX_GC_RUN,
        "observed": f"longest GC run {longest_gc} nt",
        "reasoning": ("A long GC stretch resists the helicase unwinding that "
                      "separates the two strands during RISC loading."),
    })
    return {
        "checks": checks,
        "passed": sum(1 for c in checks if c["passes"]),
        "total": len(checks),
        "source": "Ui-Tei et al. 2004; Amarzguioui & Prydz 2004",
    }


def _gc(seq: str) -> float:
    s = seq.upper()
    return round((s.count("G") + s.count("C")) / max(len(s), 1) * 100, 1)


def design_sirna_duplexes(
    ensembl_gene_id: str,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    max_candidates: int = 12,
    overhang: str = DEFAULT_OVERHANG,
) -> dict[str, Any]:
    """Tile 19-mer duplex cores across the transcript and rank by asymmetry."""
    target = get_target_analysis(ensembl_gene_id, gene_symbol=gene_symbol,
                                 organism=organism)
    mrna = (target.get("mrnaSequence") or "").upper().replace("T", "U")
    if not mrna:
        return {
            "status": "UNAVAILABLE",
            "mechanismId": "A21",
            "candidates": [],
            "message": (
                f"No transcript sequence available for "
                f"{gene_symbol or ensembl_gene_id}, so no duplex can be "
                f"designed against it."
            ),
        }

    span = DUPLEX_CORE_NT + OVERHANG_NT
    step = max(1, (len(mrna) - span) // max(max_candidates * 4, 1))
    seen: set[str] = set()
    candidates: list[dict[str, Any]] = []

    for start in range(0, len(mrna) - span + 1, step):
        core = mrna[start:start + DUPLEX_CORE_NT]
        if set(core) - set("ACGU"):
            continue
        guide_core = _revcomp(core)
        if guide_core in seen:
            continue
        seen.add(guide_core)

        # Guide 5' end pairs with the target's 3' end of the window, and vice
        # versa. Asymmetry compares the two duplex termini.
        guide_5p = guide_core[:ASYMMETRY_WINDOW_NT]
        passenger_5p = core[:ASYMMETRY_WINDOW_NT]
        dg_guide_5p = _end_stability(guide_5p)
        dg_pass_5p = _end_stability(passenger_5p)
        if dg_guide_5p is None or dg_pass_5p is None:
            continue
        # Positive = the guide's 5' end is LESS stable, which is what loads
        # the guide into RISC.
        asymmetry = round(dg_guide_5p - dg_pass_5p, 3)

        try:
            duplex_dg = round(RNA.duplexfold(guide_core, core).energy, 2)
        except Exception:
            duplex_dg = None

        gc = _gc(core)
        ui_tei = _ui_tei_rules(guide_core, core)
        flags = []
        for chk in ui_tei["checks"]:
            if not chk["passes"]:
                flags.append(f"Ui-Tei: {chk['rule']} not met ({chk['observed']})")
        if not (30.0 <= gc <= 52.0):
            flags.append(
                f"GC {gc}% outside the 30-52% window associated with "
                f"efficient RISC loading")
        if asymmetry <= 0:
            flags.append(
                "Passenger 5' end is the looser one — the passenger strand is "
                "the more likely to be loaded, which silences the wrong "
                "transcript")
        # The two asymmetry measures can disagree: the positional rules read
        # single bases at positions 1 and 19, while the thermodynamic term
        # averages four terminal base pairs. A guide starting U-G-C-C passes
        # rule (i) while its terminal tetramer is the stabler of the two. Say
        # so instead of letting the ranking quietly pick a winner.
        if ui_tei["passed"] == ui_tei["total"] and asymmetry <= 0:
            flags.append(
                "Ui-Tei positional rules and the duplex-end free-energy gap "
                "disagree here: the rules favour guide loading, the 4-bp "
                "terminal comparison does not. Treat the strand-selection "
                "prediction as unresolved for this site.")
        if "GGGG" in core or "GGGG" in guide_core:
            flags.append("Run of 4+ G: G-quadruplex and synthesis liability")

        candidates.append({
            "duplexId": f"siRNA-{gene_symbol or ensembl_gene_id}-{start}",
            "mechanismId": "A21",
            "guideStrand": guide_core + overhang,
            "passengerStrand": core + overhang,
            "guideCore": guide_core,
            "passengerCore": core,
            "overhang": overhang,
            "seedRegion": guide_core[SEED_START - 1:SEED_END],
            "seedPositions": f"guide {SEED_START}-{SEED_END}",
            "transcriptStart": start,
            "transcriptEnd": start + DUPLEX_CORE_NT,
            "gcContent": gc,
            "duplexDg": duplex_dg,
            "guide5pEndDg": dg_guide_5p,
            "passenger5pEndDg": dg_pass_5p,
            "asymmetryScore": asymmetry,
            "uiTeiRules": ui_tei,
            # Ago2 cleaves the target phosphodiester bond opposite guide
            # positions 10/11, so this is where the transcript is cut.
            "predictedCleavageSite": start + DUPLEX_CORE_NT - AGO2_CLEAVAGE_GUIDE_POSITION,
            "flags": flags,
            "notComputed": {
                "seedOffTargetLoad": (
                    "Counting transcriptome-wide seed matches needs an "
                    "alignment this service does not perform. Screen the seed "
                    "against the transcriptome before synthesis."
                ),
            },
        })

    if not candidates:
        return {
            "status": "UNAVAILABLE",
            "mechanismId": "A21",
            "candidates": [],
            "message": ("No unambiguous 19-nt window was available in the "
                        "transcript."),
        }

    # Rank by asymmetry: the larger the gap in favour of a loose guide 5'
    # end, the more reliably the guide strand is the one loaded.
    # Rank on the positional rules first — they read Argonaute's actual
    # preference — and use the thermodynamic gap to break ties within a rule
    # count.
    candidates.sort(key=lambda c: (-c["uiTeiRules"]["passed"],
                                   -c["asymmetryScore"]))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1

    return {
        "status": "OK",
        "mechanismId": "A21",
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "transcript": (target.get("canonicalTranscript") or {}).get("id"),
        "transcriptLength": len(mrna),
        "architecture": (
            f"{DUPLEX_CORE_NT}-nt duplex core with {OVERHANG_NT}-nt 3' "
            f"overhangs on both strands ({overhang})"
        ),
        "ranking": {
            "orderedBy": "uiTeiRules.passed, then asymmetryScore",
            "rationale": (
                "Argonaute loads the strand whose 5' end is less stably "
                "paired (Khvorova 2003, Schwarz 2003), and the Ui-Tei "
                "positional rules predict which strand that is more reliably "
                "than the bulk free-energy gap, because the MID pocket reads "
                "the guide's 5' nucleotide directly. Rules first, "
                "thermodynamic gap as the tie-break."
            ),
        },
        "dataProvenance": {
            "sequence": "Ensembl canonical transcript",
            "duplexDg": "ViennaRNA duplexfold",
            "endStability": "primer3 nearest-neighbour (SantaLucia)",
        },
        "candidates": candidates[:max_candidates],
    }
