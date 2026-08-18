"""TG06 candidate generation.

The load-bearing assertion here is the negative one: no candidate carries a
predicted change in protein output. There is no fitted coefficient and no
calibration set for translational effect size, so a fold-change would be an
invented number rendered beside measured ones.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from services import translational_regulation_service as T


def _synthetic_target():
    """5' UTR (strong Kozak at the AUG) | CDS | 3' UTR."""
    utr5 = "GGCACGAUGCCCUAAGGCACCGCCACC" * 4
    cds = "AUGGCUGCUGCUAAAGGGCCCUUUGGGAAACCCUAA" * 6
    utr3 = "UUUAAAGGGCCCUUUAAAGGG" * 5
    mrna = utr5 + cds + utr3
    return {
        "mrnaSequence": mrna,
        "utr5": {"sequence": utr5, "start": 0, "end": len(utr5)},
        "utr3": {"sequence": utr3, "start": len(utr5 + cds), "end": len(mrna)},
        "cdsStart": len(utr5),
        "uorfs": [],
        "kozak": T._kozak_context(mrna, len(utr5)),
        "structuredElements": [],
        "ires": [],
        "polyASite": {"start": len(mrna) - 60, "end": len(mrna),
                      "sequence": mrna[-60:]},
    }


def _candidates(element, mechanism, chemistry="pmo", length=20):
    return T.generate_translational_candidates(
        element, "suppress", mechanism, length, chemistry, [],
        _synthetic_target())


# ---------------------------------------------------------------------------
# The no-fabrication rule
# ---------------------------------------------------------------------------

def test_no_candidate_claims_a_predicted_fold_change():
    """The spec asked for a log2 fold-change from a 'mechanism-specific
    sensitivity coefficient'. None has been fitted and there is no
    calibration set, so none is emitted."""
    out = _candidates("5p_utr", "A2")
    assert out["candidates"]
    banned = ("foldChange", "log2FoldChange", "translationalChangeScore",
              "predictedEfficacy", "rbpDisplacementScore")
    for cand in out["candidates"]:
        for key in banned:
            assert key not in cand, f"{key} is not a computable quantity here"


def test_real_and_heuristic_numbers_are_kept_apart():
    cand = _candidates("5p_utr", "A2")["candidates"][0]
    assert "ViennaRNA" in cand["realMetrics"]["provenance"]
    assert "not measurements" in cand["heuristicEstimates"]["provenance"]
    # The ranking signal says what it is.
    assert "NOT a predicted change" in cand["interpretation"]


def test_naming_an_rbp_yields_a_note_not_a_displacement_number():
    out = T.generate_translational_candidates(
        "5p_utr", "suppress", "A2", 20, "pmo", [], _synthetic_target(),
        target_rbp="PTB")
    assert out["rbpNote"] and "no displacement score" in out["rbpNote"]
    assert all("rbpDisplacementScore" not in c for c in out["candidates"])


# ---------------------------------------------------------------------------
# Real computation
# ---------------------------------------------------------------------------

def test_metrics_are_computed_from_the_sequence():
    pytest.importorskip("RNA", reason="ViennaRNA not installed")
    cand = _candidates("5p_utr", "A2")["candidates"][0]
    m = cand["realMetrics"]
    assert m["targetDuplexEnergy"] < 0        # a real duplex is favourable
    assert 0 < m["meltingTempC"] < 120
    assert 0 <= m["gcContent"] <= 1
    assert m["lengthNt"] == len(cand["sequence"])


def test_the_oligo_is_the_reverse_complement_of_its_target_site():
    """The shared helper translates ATGC only, so a U in an RNA target would
    pass through and produce a mixed-alphabet oligo."""
    cand = _candidates("5p_utr", "A2")["candidates"][0]
    assert set(cand["sequence"]) <= set("ACGT"), "oligo must be DNA alphabet"
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    expected = "".join(comp[b] for b in
                       reversed(cand["targetSite"].replace("U", "T")))
    assert cand["sequence"] == expected


def test_candidates_are_ranked_by_engagement():
    out = _candidates("kozak_consensus", "A30")
    scores = [c["elementEngagement"] for c in out["candidates"]]
    assert scores == sorted(scores, reverse=True)
    assert [c["rank"] for c in out["candidates"]] == list(
        range(1, len(out["candidates"]) + 1))


def test_every_candidate_actually_overlaps_the_element():
    out = _candidates("kozak_consensus", "A30")
    for c in out["candidates"]:
        assert c["elementOverlapNt"] > 0


# ---------------------------------------------------------------------------
# Gates
# ---------------------------------------------------------------------------

def test_a_cleaving_chemistry_is_refused():
    """TG06 occupies an element; a gapmer destroys the transcript instead."""
    out = _candidates("5p_utr", "A2", chemistry="gapmer")
    assert out["status"] == "incompatible_chemistry"
    assert out["candidates"] == []


def test_a_missing_element_is_reported_not_invented():
    out = _candidates("ires_element", "A29")
    assert out["status"] == "element_not_found"
    assert out["candidates"] == []
    assert "property of the transcript" in out["message"]


def test_no_transcript_means_no_candidates():
    out = T.generate_translational_candidates(
        "5p_utr", "suppress", "A2", 20, "pmo", [], {"mrnaSequence": ""})
    assert out["candidates"] == []


def test_kozak_strength_reports_consensus_matches_not_a_rate():
    t = _synthetic_target()
    k = t["kozak"]
    assert k["strength"] in ("strong", "adequate", "weak")
    assert "Not a predicted initiation rate" in k["note"]


def test_ires_windows_do_not_claim_to_be_ires_calls():
    """No validated IRES predictor is wired; these are folding observations."""
    windows = T._ires_domains("GGGCGCGCGGGCCCGGGCCCGCGCGCCC" * 6)
    for w in windows:
        assert w["predicted"] is False
        assert "Not an IRES call" in w["note"]
