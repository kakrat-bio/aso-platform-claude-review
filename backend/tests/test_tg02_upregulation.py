"""TG02 regression tests — A5 (uORF blocking) and A6 (miRNA site blocking).

Each test here corresponds to a defect that was measured on this branch, on
SCN1A at 20 nt, before it was fixed:

  * A5 and A28 were REJECTED on the defect gate yet reported
    `score: 1.0, applicability: [1.0, 1.0]` — a higher number than the two
    ELIGIBLE mechanisms in the same response (A3/A4 at 0.9). A gate rejection
    is decided before any feature is collected, so the Fréchet–Hoeffding
    interval over an empty contributing set is the vacuous [1, 1].
  * `/generate` returned 69 candidates for A5 and 635 for A6 — the entire
    tiling, uncapped, where A1 caps at 10 and the A3/A4 designer at 12.
  * 19 of those 69 and 125 of those 635 tied at exactly compositeScore 100.0,
    the same `_composite_score` saturation already fixed for A1.

Note the file name and the test names: this suite is the one a session runs
as `pytest backend/tests -q -k upregulation`, which before this file matched
nothing at all and exited 5 ("no tests collected") while reading as a pass.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from services import mechanism_arbitration as A

SCN1A = "ENSG00000144285"


def _client():
    import logging
    logging.disable(logging.WARNING)
    from fastapi.testclient import TestClient
    from main import app
    return TestClient(app)


# ---------------------------------------------------------------------------
# Scoring — a mechanism that was not scored must not carry a score
# ---------------------------------------------------------------------------

def test_upregulation_gate_rejected_mechanism_reports_no_applicability():
    """The vacuous [1.0, 1.0] must not be published as a maximal score.

    A5 does not address haploinsufficiency, so it fails the defect gate. The
    gate runs before required/forbidden features are collected, leaving the
    interval vacuous. Reporting that as 1.0 tells a reader who sorts on
    `score` the exact opposite of what the status says.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SCN1A", molecular_defect="haploinsufficiency"))
    a5 = next(r for r in out["results"] if r["id"] == "A5")
    assert a5["status"] == A.REJECTED
    assert a5["score"] is None, f"gate-rejected A5 still scores {a5['score']}"
    assert a5["applicability"] is None
    assert a5["confidence"] is None
    assert a5["rationale"], "a rejection must say why"


def test_upregulation_no_unscored_mechanism_outranks_a_scored_one():
    """Whatever the status set, a number that exists must be comparable.

    This is the property the A5 bug actually violated: two REJECTED rows
    carried 1.0 while the ELIGIBLE rows carried 0.9.
    """
    # Scoped to the TG02 ranking a user actually reads, through the endpoint
    # so gene features resolve as they do in a real request. Comparing across
    # the whole unified pass hides this: some other goal's ELIGIBLE mechanism
    # sits at 1.0 and the assertion passes while TG02 is still inverted.
    resp = _client().post("/api/mechanisms/gene-upregulation", json={
        "gene_symbol": "SCN1A", "defect_type": "haploinsufficiency"})
    assert resp.status_code == 200
    results = resp.json()["results"]
    scored = [r for r in results if r["score"] is not None]
    eligible = [r["score"] for r in scored if r["status"] == "ELIGIBLE"]
    not_eligible = [r["score"] for r in scored if r["status"] != "ELIGIBLE"]
    if not eligible:
        pytest.skip("nothing eligible in this TG02 run; nothing to invert")
    assert not not_eligible or max(not_eligible) <= max(eligible), (
        "a non-eligible TG02 mechanism carries a higher number than the best "
        f"eligible one: {not_eligible} vs {eligible}")


def test_upregulation_a5_and_a6_are_scored_or_halted_with_a_reason():
    """Neither mechanism may come back with no verdict and no explanation."""
    resp = _client().post("/api/mechanisms/gene-upregulation", json={
        "gene_symbol": "SCN1A", "defect_type": "haploinsufficiency"})
    assert resp.status_code == 200
    by_id = {m["id"]: m for m in resp.json()["results"]}
    for mech in ("A5", "A6"):
        assert mech in by_id, f"{mech} missing from the TG02 ranking"
        row = by_id[mech]
        assert row["status"] in ("ELIGIBLE", "HALTED", "REJECTED", "FLAGGED")
        if row["status"] == "ELIGIBLE":
            assert row["score"] is not None
        else:
            assert row["rationale"], f"{mech} is {row['status']} with no reason"


def test_upregulation_a6_halt_names_the_missing_feature():
    """A6 halts because F7 (repressive miRNA site) has no wired source.

    "Halted" on its own is not actionable; the response has to name F7 so a
    reader knows what would unblock it.
    """
    resp = _client().post("/api/mechanisms/gene-upregulation", json={
        "gene_symbol": "SCN1A", "defect_type": "haploinsufficiency"})
    a6 = next(m for m in resp.json()["results"] if m["id"] == "A6")
    if a6["status"] != "HALTED":
        pytest.skip("A6 did not halt on this input")
    assert any("F7" in r for r in a6["rationale"])
    assert [f["id"] for f in a6["features"]["unresolved"]] == ["F7"]


# ---------------------------------------------------------------------------
# Design — candidates must come from the region the mechanism targets
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def transcript():
    from services.gene_silencing_service import get_target_analysis
    t = get_target_analysis(SCN1A, "SCN1A", "homo_sapiens")
    if not t.get("mrnaSequence"):
        pytest.skip("Ensembl unavailable; the design path cannot be exercised")
    return t


def _generate(mechanism_id, chemistry, length=20):
    resp = _client().post("/api/gene-upregulation/generate", json={
        "ensembl_gene_id": SCN1A, "mechanism_id": mechanism_id,
        "aso_length": length, "chemistry": chemistry, "modifications": [],
        "gene_symbol": "SCN1A", "organism": "homo_sapiens",
    })
    assert resp.status_code == 200, resp.text
    cands = resp.json()["candidates"]
    if not cands:
        pytest.skip(f"{mechanism_id} produced no candidates (fetch unavailable)")
    return cands


@pytest.mark.parametrize("mechanism_id,chemistry,utr_key,label", [
    ("A5", "pmo", "utr5Sequence", "5' UTR"),
    ("A6", "2ome", "utr3Sequence", "3' UTR"),
])
def test_upregulation_candidates_derive_from_the_targeted_utr(
        transcript, mechanism_id, chemistry, utr_key, label):
    """A5 blocks a uORF (5' UTR) and A6 masks a miRNA seed site (3' UTR).

    The stored `sequence` is the ASO, so its reverse complement is the target
    window and must be findable in that UTR — and NOT in the CDS, which is
    what the pre-UTR fallback path used to hand back.
    """
    from services.gene_silencing_service import _reverse_complement
    utr = (transcript.get(utr_key) or "").upper()
    if not utr:
        pytest.skip(f"no {label} in the fetched transcript")
    cds = (transcript.get("mrnaSequence") or "").upper()

    for cand in _generate(mechanism_id, chemistry):
        window = _reverse_complement(cand["sequence"])
        assert window in utr, (
            f"{mechanism_id} candidate {cand['sequence']} does not come from "
            f"the {label} it claims to target ({cand['targetRegion']})")
        assert window not in cds, (
            f"{mechanism_id} candidate {cand['sequence']} is a CDS window, "
            f"not a {label} window")


def test_upregulation_candidate_offset_matches_the_reported_region(transcript):
    """The offset in `targetRegion` must be the real index into the UTR.

    A label carrying a position that does not locate the window is worse than
    no label: it reads as a coordinate.
    """
    from services.gene_silencing_service import _reverse_complement
    utr5 = (transcript.get("utr5Sequence") or "").upper()
    if not utr5:
        pytest.skip("no 5' UTR in the fetched transcript")
    for cand in _generate("A5", "pmo"):
        offset = int(cand["targetRegion"].rsplit("+", 1)[1])
        window = _reverse_complement(cand["sequence"])
        assert utr5[offset:offset + len(window)] == window, (
            f"reported offset +{offset} does not hold {window}")


@pytest.mark.parametrize("mechanism_id,chemistry", [("A5", "pmo"), ("A6", "2ome")])
def test_upregulation_candidate_list_is_capped(mechanism_id, chemistry):
    """Uncapped this returned 69 (A5) and 635 (A6) candidates on SCN1A."""
    from services.gene_upregulation_service import UPREGULATION_MAX_CANDIDATES
    cands = _generate(mechanism_id, chemistry)
    assert len(cands) <= UPREGULATION_MAX_CANDIDATES <= 20, (
        f"{mechanism_id} returned {len(cands)} candidates")
    # A capped list must say what it was drawn from, or it reads as the whole
    # search space.
    assert cands[0]["poolSize"] >= len(cands)
    assert cands[0]["shortlistedFrom"]


@pytest.mark.parametrize("mechanism_id,chemistry", [("A5", "pmo"), ("A6", "2ome")])
def test_upregulation_composite_score_is_not_one_saturated_value(
        mechanism_id, chemistry):
    """`_composite_score` clips at 100 when candidates share a length.

    A ranking column with one value is not a ranking. This asserts the
    surviving shortlist discriminates, and that the saturated ceiling is not
    the majority of it.
    """
    cands = _generate(mechanism_id, chemistry)
    scores = [c["compositeScore"] for c in cands]
    assert len(set(scores)) > 1, (
        f"{mechanism_id}: every candidate scored {scores[0]}")
    # Pre-fix this shortlist did not exist — the ceiling held 19 of 69 (A5)
    # and 125 of 635 (A6). Within a 20-candidate shortlist, more than a
    # handful at exactly 100.0 means the axis has stopped discriminating.
    tied_at_ceiling = sum(1 for s in scores if s == 100.0)
    assert tied_at_ceiling <= max(2, len(scores) // 5), (
        f"{mechanism_id}: {tied_at_ceiling}/{len(scores)} candidates are "
        "pinned at the 100.0 ceiling")


@pytest.mark.parametrize("mechanism_id,chemistry", [("A5", "pmo"), ("A6", "2ome")])
def test_upregulation_ranks_on_site_accessibility(mechanism_id, chemistry):
    """A steric blocker has to occupy its site, so accessibility leads.

    ΔG-first ordering measured below chance on 528 held-out experiments
    (E12), and the composite saturates; accessibility varies over orders of
    magnitude between sites on one transcript.
    """
    cands = _generate(mechanism_id, chemistry)
    acc = [c["realMetrics"].get("siteAccessibility") for c in cands]
    assert all(v is not None for v in acc), (
        f"{mechanism_id}: no accessibility computed, ranking has no primary axis")
    assert all(acc[i] >= acc[i + 1] for i in range(len(acc) - 1)), (
        f"{mechanism_id}: candidates are not ordered by accessibility")
    assert cands[0]["rankingBasis"]["primary"].startswith("siteAccessibility")
    # The caveat is not decoration: neither a uORF nor a miRNA site is
    # validated here, and the ranking must not read as an activity model.
    assert cands[0]["rankingBasis"]["caveat"]


def test_upregulation_a6_candidates_are_marked_unverified():
    """No miRNA-target database is integrated, so no A6 window may claim to
    mask a validated seed site."""
    for cand in _generate("A6", "2ome"):
        assert cand["seedSiteStatus"] == "unverified"
        assert cand["seedSiteNote"]
