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

# ADAR nearest-neighbour preference, from the deaminase's own substrate
# selectivity rather than from the guide design.
#
# ADAR2 reads the bases flanking the target adenosine. The 5' neighbour is the
# dominant term and G immediately 5' is strongly disfavoured; the 3' neighbour
# is a weaker but real preference for G. Both REPAIR and CIRTS place an ADAR
# or ADAR-like deaminase on the target, so this governs whether the intended
# adenosine is edited efficiently AND how likely each bystander adenosine is
# to be hit. Ranks, not rates: the ordering is well documented, the absolute
# efficiencies are context- and construct-specific and are not asserted here.
ADAR_5P_PREFERENCE = {"U": 4, "A": 3, "C": 2, "G": 1}
ADAR_3P_PREFERENCE = {"G": 4, "C": 3, "A": 2, "U": 1}
ADAR_PREFERENCE_NOTE = (
    "ADAR2 prefers a U or A immediately 5' of the target adenosine and "
    "disfavours G there; 3' it prefers G. Reported as an ordinal rank, not a "
    "predicted editing rate."
)


def _adar_context(seq: str, index: int) -> dict[str, Any]:
    """Nearest-neighbour context of the adenosine at `index` within `seq`."""
    five_p = seq[index - 1] if index > 0 else None
    three_p = seq[index + 1] if index + 1 < len(seq) else None
    r5 = ADAR_5P_PREFERENCE.get(five_p or "", 0)
    r3 = ADAR_3P_PREFERENCE.get(three_p or "", 0)
    if r5 >= 3 and r3 >= 3:
        favourability = "favourable"
    elif r5 <= 1:
        favourability = "poor (G immediately 5' is the strongest negative)"
    else:
        favourability = "intermediate"
    return {
        "fivePrimeNeighbour": five_p,
        "threePrimeNeighbour": three_p,
        "triplet": f"{five_p or '-'}A{three_p or '-'}",
        "fivePrimeRank": r5,
        "threePrimeRank": r3,
        "favourability": favourability,
    }


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
        target_context = _adar_context(mrna, edit_position)

        # Adenosines elsewhere in the duplex are candidate bystander edits —
        # but they are not equally likely. An adenosine with a G immediately
        # 5' of it is a poor ADAR substrate and is a much smaller risk than one
        # in a UAG context. Counting them all equally overstated the risk of a
        # window full of G-preceded adenosines and understated one containing a
        # single perfect substrate.
        bystanders = []
        for i, b in enumerate(window):
            if b != expected or (win_start + i) == edit_position:
                continue
            ctx = _adar_context(mrna, win_start + i)
            bystanders.append({
                "transcriptPosition": win_start + i,
                "spacerPosition": len(window) - i,
                **ctx,
            })
        high_risk = [b for b in bystanders if b["favourability"] == "favourable"]

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
            "targetAdarContext": target_context,
            "bystanderCount": len(bystanders),
            "highRiskBystanderCount": len(high_risk),
            "bystanders": bystanders[:20],
            "adarPreferenceNote": ADAR_PREFERENCE_NOTE,
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

    # Rank on bystanders that ADAR would actually edit, then on the total,
    # then on duplex stability. A window carrying ten G-preceded adenosines is
    # a safer guide than one carrying two in UAG context.
    candidates.sort(key=lambda c: (c["highRiskBystanderCount"],
                                   c["bystanderCount"],
                                   c["duplexDg"] if c["duplexDg"] is not None else 0.0))
    for i, c in enumerate(candidates[:max_candidates]):
        c["rank"] = i + 1

    target_context = _adar_context(mrna, edit_position)
    return {
        "status": "OK",
        "mechanismId": mechanism_id,
        "platform": platform["name"],
        "targetAdarContext": target_context,
        "targetContextWarning": (
            None if target_context["fivePrimeRank"] >= 2 else
            f"The target adenosine sits in a {target_context['triplet']} "
            f"context. A guanosine immediately 5' is ADAR's least favoured "
            f"neighbour, so editing at this site is expected to be "
            f"inefficient regardless of guide placement."
        ),
        "geneSymbol": gene_symbol or ensembl_gene_id,
        "architecture": (
            f"{spacer_nt} nt spacer with a {platform['mismatchBase']} "
            f"mismatch opposite the edited base, plus the "
            f"{platform['scaffoldRole']}"
        ),
        "ranking": {
            "orderedBy": "highRiskBystanderCount, then bystanderCount, then duplexDg",
            "rationale": (
                "Every other adenosine inside the guide-target duplex is a "
                "candidate bystander edit, but not an equal one: ADAR "
                "disfavours an adenosine with a G immediately 5' of it and "
                "prefers one in a U-A-G context. Placements enclosing fewer "
                "GOOD ADAR substrates are preferred over placements enclosing "
                "fewer adenosines overall."
            ),
        },
        "candidates": candidates[:max_candidates],
    }
