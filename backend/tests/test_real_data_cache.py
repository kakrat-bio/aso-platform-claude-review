"""The no-fabrication policy, enforced.

When an external source is unreachable this platform serves real data it
fetched earlier, or it says the data is unavailable. It never synthesises a
substitute, because a substitute is indistinguishable from a measurement once
it reaches the UI.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from database.db import init_db
from services import real_data_cache as RDC

init_db()

REAL = {"cds": "AUGGCCUAA", "utr5": "GCCACCAUG", "utr3": "UGAUUU",
        "transcript_id": "ENST_TEST"}


def test_no_source_and_no_cache_is_unavailable_not_a_substitute():
    out = RDC.resolve("transcript_parts", "test:UNSEEN_GENE", lambda: None)
    assert out["status"] == RDC.UNAVAILABLE
    assert out["data"] is None
    assert "fabricated" in out["note"]


def test_a_live_fetch_is_cached_and_replayed_verbatim_on_outage():
    key = "test:REPLAY_GENE"
    live = RDC.resolve("transcript_parts", key, lambda: dict(REAL),
                       source="Ensembl REST")
    assert live["status"] == RDC.LIVE

    outage = RDC.resolve("transcript_parts", key, lambda: None)
    assert outage["status"] == RDC.CACHED
    # The same real sequence, not a regenerated approximation of it.
    assert outage["data"]["cds"] == REAL["cds"]
    assert outage["data"]["utr5"] == REAL["utr5"]
    assert "earlier fetch" in outage["note"]


def test_an_exception_is_an_outage_not_a_crash():
    key = "test:BOOM_GENE"
    RDC.resolve("transcript_parts", key, lambda: dict(REAL))

    def boom():
        raise RuntimeError("Ensembl 503")

    out = RDC.resolve("transcript_parts", key, boom)
    assert out["status"] == RDC.CACHED
    assert out["data"]["cds"] == REAL["cds"]


def test_every_outcome_declares_which_one_it_is():
    """A caller must never be able to mistake a replay for a live read."""
    key = "test:STATUS_GENE"
    statuses = {
        RDC.resolve("transcript_parts", key, lambda: dict(REAL))["status"],
        RDC.resolve("transcript_parts", key, lambda: None)["status"],
        RDC.resolve("transcript_parts", "test:NOTHING_HERE",
                    lambda: None)["status"],
    }
    assert statuses == {RDC.LIVE, RDC.CACHED, RDC.UNAVAILABLE}


def test_protein_replacement_reports_unavailable_rather_than_padding():
    """TG08 previously fell back to "AUG" + "GCU"*300 for the CDS and
    "GCCACC" + "A"*n for the UTRs, then computed CAI, GC, MFE, U-content and
    a protein-yield estimate from that padding."""
    import services.protein_replacement_service as PRS

    source = Path(PRS.__file__).read_text()
    # Only the docstring explaining the removal may mention them.
    code_lines = [
        ln for ln in source.splitlines()
        if not ln.lstrip().startswith("#") and "fell back" not in ln
        and "for the UTRs" not in ln
    ]
    code = "\n".join(code_lines)
    assert '"AUG" + "GCU" * 300' not in code
    assert '"GCCACC" + "A" *' not in code


# ---------------------------------------------------------------------------
# TG08 de-mocking (docx Part 1b-1e)
# ---------------------------------------------------------------------------

def _prs():
    import services.protein_replacement_service as PRS
    return PRS


def test_amino_acid_identity_is_computed_not_asserted():
    """It was hardcoded to 100.0, which asserted a match nobody checked."""
    P = _prs()
    assert P._calc_amino_acid_identity("AUGGCUGCUUAA", "AUGGCCGCCUAA") == 100.0
    # A non-synonymous change must show up rather than reading 100.
    assert P._calc_amino_acid_identity("AUGUGGUAA", "AUGGCCUAA") < 100.0
    # No native CDS to compare against is None, not a claimed match.
    assert P._calc_amino_acid_identity("AUGGCUUAA", "") is None


def test_utr_structure_is_folded_not_assumed_passed():
    """It was hardcoded "PASSED" for every construct, including ones whose
    UTR had never been examined."""
    import pytest
    pytest.importorskip("RNA", reason="ViennaRNA not installed")
    P = _prs()

    assert P._evaluate_utr_structure("")["flag"] == "NOT_COMPUTED"

    hairpin = P._evaluate_utr_structure("GGGGCCCCGGGGCCCC" * 4)
    assert hairpin["flag"] == "BLOCKED"
    assert hairpin["mfe"] < -25
    assert hairpin["hairpins"] > 0

    open_utr = P._evaluate_utr_structure("A" * 60)
    assert open_utr["flag"] == "PASSED"


def test_mfe_plot_is_a_real_dot_bracket_not_a_pattern():
    """The old value was ''.join('(' if i % 7 < 3 ...) — a procedural
    pattern that rendered as a structure."""
    import pytest
    pytest.importorskip("RNA", reason="ViennaRNA not installed")
    P = _prs()
    structure, mfe = P._fold_sequence("GGGCGCGCGGGCCCGGGCCCGCGCGCCC")
    assert set(structure) <= set("().")
    assert structure.count("(") == structure.count(")")
    assert mfe < 0

    source = Path(P.__file__).read_text()
    assert "i % 7 < 3" not in source


def test_modification_effects_scale_tlr_and_default_safely():
    P = _prs()
    base = {"tlr3Score": 20.0, "tlr7Score": 10.0}
    assert P._apply_modification_effects(base, "unmodified") == base
    m1 = P._apply_modification_effects(base, "m1psi")
    assert m1["tlr3Score"] < base["tlr3Score"]
    # An unrecognised modification must not be assumed beneficial.
    assert P._apply_modification_effects(base, "bogus") == base
