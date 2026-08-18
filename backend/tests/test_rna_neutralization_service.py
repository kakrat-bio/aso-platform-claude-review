"""TG05 candidate generation.

The load-bearing assertions are negative: no off-target count and no RBP
displacement score, because neither is computed. For repeat masking the
off-target number matters more than anywhere else — a (CAG)n oligo is
complementary to every CAG-repeat transcript by construction, so a low
count would imply a selectivity that was never checked.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from services import rna_neutralization_service as N


def _a14(gene="DMPK", unit="CTG", count=">50 copies", length=17,
         chemistry="pmo", **kw):
    return N.generate_neutralization_candidates(
        gene, "A14", "steric_repeat_masking", unit, count, length,
        chemistry, **kw)


def test_no_off_target_count_is_reported():
    out = _a14()
    assert out["status"] == "ok"
    for c in out["candidates"]:
        assert "offTargetRepeatCount" not in c
        assert "offTargetRisk" not in c


def test_the_off_target_limitation_is_stated_prominently():
    """Silence here would read as 'no off-target problem'."""
    out = _a14()
    joined = " ".join(out["notes"])
    assert "OFF-TARGET RISK IS NOT QUANTIFIED" in joined
    assert "every transcript carrying that repeat" in joined


def test_no_rbp_displacement_score():
    out = _a14(target_rbp="MBNL1")
    for c in out["candidates"]:
        assert "rbpDisplacementScore" not in c
        assert "rbpBindingDeltaG" not in c
    assert any("MBNL1" in n for n in out["notes"])


def test_a_repeat_tract_yields_one_candidate_per_phase():
    """A tract is periodic: only the phases within one unit are distinct."""
    out = _a14(unit="CTG")
    assert len(out["candidates"]) == 3
    assert {c["phase"] for c in out["candidates"]} == {0, 1, 2}

    six = _a14(gene="C9orf72", unit="GGGGCC")
    assert len(six["candidates"]) == 6


def test_tract_provenance_separates_curated_from_asserted():
    """A catalogue lookup and a unit the user typed are different evidence."""
    out = _a14()
    assert out["tractProvenance"]["provenance"] == "user_asserted"
    assert "hypothesis about the target" in out["tractProvenance"]["note"]


def test_no_unit_anywhere_means_no_candidates():
    out = _a14(gene="NOSUCHGENE", unit=None, count=None)
    assert out["status"] == "target_unavailable"
    assert out["candidates"] == []


def test_an_invalid_motif_is_not_treated_as_a_unit():
    out = _a14(gene="NOSUCHGENE", unit="not a motif!", count=None)
    assert out["status"] == "target_unavailable"


def test_metrics_are_computed_from_the_sequence():
    pytest.importorskip("RNA", reason="ViennaRNA not installed")
    cand = _a14()["candidates"][0]
    m = cand["realMetrics"]
    assert m["targetDuplexEnergy"] < 0
    assert 0 < m["meltingTempC"] < 120
    assert m["lengthNt"] == len(cand["sequence"])
    assert set(cand["sequence"]) <= set("ACGT")


def test_candidates_are_ranked_by_duplex_energy():
    out = _a14()
    energies = [c["realMetrics"]["targetDuplexEnergy"] for c in out["candidates"]]
    assert energies == sorted(energies)


def test_a_cleaving_chemistry_is_refused():
    """TG05 occupies the RNA; degrading it is a different therapy (A1)."""
    out = _a14(chemistry="gapmer")
    assert out["status"] == "incompatible_chemistry"
    assert out["candidates"] == []


def test_an_antimir_needs_the_actual_mirna_sequence():
    without = N.generate_neutralization_candidates(
        "MIR21", "A12", "microrna_antagomir", None, None, 17, "moe_full_ps")
    assert without["status"] == "target_unavailable"
    assert without["candidates"] == []

    with_seq = N.generate_neutralization_candidates(
        "MIR21", "A12", "microrna_antagomir", None, None, 17, "moe_full_ps",
        mirna_sequence="UAGCUUAUCAGACUGAUGUUGA")
    assert with_seq["status"] == "ok"
    assert with_seq["candidates"]


def test_a25_is_flagged_not_designed():
    out = N.generate_neutralization_candidates(
        "VEGFA", "A25", "aptamer_decoy", None, None, 17, "pmo")
    assert out["status"] == "flagged_not_designed"
    assert out["candidates"] == []
