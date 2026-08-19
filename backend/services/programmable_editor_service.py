"""A18 / A19 — guide design for protein-dependent RNA editors.

A18 (CIRTS) and A19 (REPAIR) were marked design-unavailable because the
RNA-editing designer builds linear guides that recruit an ENDOGENOUS
deaminase, and these two do not work that way: each delivers an engineered
protein and pairs it with a guide carrying a scaffold that protein
recognises. Handing back an ADAR guide labelled "A19" was the bug fixed
earlier; this module supplies the design these mechanisms actually need.

WHAT IS DESIGNED HERE, AND WHAT IS NOT.

A guide for either platform has two parts:

* **The spacer** — the target-complementary stretch. It changes with every
  target, it is the part that has to be designed, and it is fully
  determined by the transcript. That is what this module produces: spacers
  tiled across the edit site, with the editing window positioned over the
  target base, real duplex thermodynamics, and the liability checks that
  matter for these platforms.
* **The scaffold** — the direct repeat (REPAIR/Cas13) or hairpin (CIRTS)
  that the delivered protein binds. It is a FIXED sequence, identical for
  every target using that ortholog, and it is a property of the construct
  rather than of the design.

The scaffold constant is deliberately not written into this file. Getting a
direct repeat wrong by one base yields a guide that does not load, and the
project's standing rule is that a value which cannot be sourced is reported
rather than recalled. Each candidate therefore carries `scaffoldRequired`,
naming the ortholog and where the sequence must come from, and the spacer is
emitted ready to be appended to it.
"""

from __future__ import annotations

import logging
from typing import Any

import RNA

from services.gene_silencing_service import get_target_analysis

logger = logging.getLogger(__name__)

# Published guide architectures. Lengths and editing-window offsets are
# structural facts about each platform; the scaffold sequences are not
# reproduced here (see the module docstring).
PLATFORMS: dict[str, dict[str, Any]] = {
    "A19": {
        "name": "REPAIR (dCas13b-ADAR2dd)",
        "effector": "catalytically dead PspCas13b fused to the ADAR2 deaminase domain",
        "spacerNt": 50,
        # Cox et al. place the target adenosine opposite a cytidine mismatch
        # positioned toward the middle of the spacer.
        "mismatchOffsetFromSpacer5p": 26,
        "mismatchBase": "C",
        "scaffoldRole": "direct repeat recognised by PspCas13b",
        "scaffoldSource": (
            "Take the PspCas13b direct-repeat sequence from the ortholog's "
            "primary publication or the vector map you are using. It is a "
            "fixed constant for every target and is not reproduced here "
            "rather than risk a recalled base."
        ),
        "editType": "a_to_i",
    },
    "A18": {
        "name": "CIRTS (CRISPR-Cas-Inspired RNA Targeting System)",
        "effector": "engineered ssRNA-binding protein fused to an effector domain",
        "spacerNt": 30,
        "mismatchOffsetFromSpacer5p": 16,
        "mismatchBase": "C",
        "scaffoldRole": "hairpin recognised by the engineered RNA-binding protein",
        "scaffoldSource": (
            "Take the hairpin from the CIRTS construct you are building; it "
            "is specific to the RNA-binding protein in the fusion and is a "
            "fixed constant per construct."
        ),
        "editType": "a_to_i",
    },
}

_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


def _revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT[b] for b in reversed(seq.upper()))


def design_editor_guides(
    mechanism_id: str,
    ensembl_gene_id: str,
    edit_position: int,
    gene_symbol: str = "",
    organism: str = "homo_sapiens",
    max_candidates: int = 8,
) -> dict[str, Any]:
    """Spacers for a protein-dependent RNA editor, centred on the edit site.

    `edit_position` is 0-based into the transcript and names the base to be
    edited.
    """
    platform = PLATFORMS.get(mechanism_id)
    if platform is None:
        raise ValueError(
            f"{mechanism_id} is not a protein-dependent RNA editor. This "
            f"service designs guides for: {', '.join(sorted(PLATFORMS))}."
        )

    target = get_target_analysis(ensembl_gene_id, gene_symbol=gene_symbol,
                                 organism=organism)
    mrna = (target.get("mrnaSequence") or "").upper().replace("T", "U")
    if not mrna:
        return {"status": "UNAVAILABLE", "mechanismId": mechanism_id,
                "candidates": [],
                "message": f"No transcript sequence for {gene_symbol or ensembl_gene_id}."}

    if not (0 <= edit_position < len(mrna)):
        return {"status": "UNAVAILABLE", "mechanismId": mechanism_id,
                "candidates": [],
                "message": (f"edit_position {edit_position} is outside the "
                            f"{len(mrna)} nt transcript.")}

    edited_base = mrna[edit_position]
    expected = "A" if platform["editType"] == "a_to_i" else "C"
    if edited_base != expected:
        return {"status": "UNAVAILABLE", "mechanismId": mechanism_id,
                "candidates": [],
                "message": (f"Position {edit_position} is {edited_base}, but "
                            f"{platform['name']} edits {expected}. Nothing is "
                            f"designed against the wrong base.")}

    spacer_nt = platform["spacerNt"]
    ideal_offset = platform["mismatchOffsetFromSpacer5p"]
    candidates: list[dict[str, Any]] = []

    # Slide the spacer so the edited base sits at a range of positions around
    # the platform's preferred offset.
    for shift in range(-4, 5):
        offset = ideal_offset + shift
        if not (1 <= offset <= spacer_nt):
            continue
        # Spacer is antisense, so the target window runs the other way.
        win_end = edit_position + offset
        win_start = win_end - spacer_nt
        if win_start < 0 or win_end > len(mrna):
            continue
        window = mrna[win_start:win_end]
        if set(window) - set("ACGU"):
            continue
        spacer = list(_revcomp(window))
        # The deaminase acts on the adenosine opposite a cytidine mismatch.
        mismatch_index = offset - 1
        if not (0 <= mismatch_index < len(spacer)):
            continue
        spacer[mismatch_index] = platform["mismatchBase"]
        spacer_seq = "".join(spacer)

        try:
            duplex_dg = round(RNA.duplexfold(spacer_seq, window).energy, 2)
        except Exception:
            duplex_dg = None

        gc = round((spacer_seq.count("G") + spacer_seq.count("C"))
                   / len(spacer_seq) * 100, 1)
        # Adenosines elsewhere in the duplex are candidate bystander edits.
        bystanders = [i for i, b in enumerate(window)
                      if b == expected and (win_start + i) != edit_position]

        candidates.append({
            "guideId": f"{mechanism_id}-{gene_symbol or ensembl_gene_id}-{offset}",
            "mechanismId": mechanism_id,
            "platform": platform["name"],
            "effector": platform["effector"],
            "spacer": spacer_seq,
            "spacerLength": len(spacer_seq),
            "mismatchPosition": offset,
            "mismatchBase": platform["mismatchBase"],
            "targetWindow": window,
            "transcriptStart": win_start,
            "transcriptEnd": win_end,
            "editPosition": edit_position,
            "gcContent": gc,
            "duplexDg": duplex_dg,
            "bystanderCount": len(bystanders),
            "bystanderPositions": bystanders[:20],
            "scaffoldRequired": {
                "role": platform["scaffoldRole"],
                "source": platform["scaffoldSource"],
                "note": ("Append the scaffold to this spacer to obtain the "
                         "full guide. The spacer is the designed part; the "
                         "scaffold is a per-construct constant."),
            },
        })

    if not candidates:
        return {"status": "UNAVAILABLE", "mechanismId": mechanism_id,
                "candidates": [],
                "message": ("The edit site is too close to a transcript end "
                            f"to place a {spacer_nt} nt spacer around it.")}

    # Fewer bystander adenosines first, then tighter duplex.
    candidates.sort(key=lambda c: (c["bystanderCount"],
                                   c["duplexDg"] if c["duplexDg"] is not None else 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1

    return {
        "status": "OK",
        "mechanismId": mechanism_id,
        "platform": platform["name"],
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "architecture": (
            f"{spacer_nt} nt spacer with a {platform['mismatchBase']} "
            f"mismatch opposite the edited base, plus the "
            f"{platform['scaffoldRole']}"
        ),
        "ranking": {
            "orderedBy": "bystanderCount, then duplexDg",
            "rationale": (
                "Every other adenosine inside the guide-target duplex is a "
                "candidate bystander edit, so a spacer placement that "
                "encloses fewer of them is preferred."
            ),
        },
        "candidates": candidates[:max_candidates],
    }
