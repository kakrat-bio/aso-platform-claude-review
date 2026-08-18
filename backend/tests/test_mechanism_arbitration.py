"""Tests for the unified mechanism arbitration.

These lock down the decisions from
`docs/planning/therapeutic_goal_scope_plan_v3.md` (v3) that are easy to
regress: the goal being an output rather than an input, the difference
between "we looked and it is not there" and "we have no way to look", and
the three-way split between scored, halted and flagged.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from services import feature_service as F
from services import mechanism_arbitration as A
from services import mechanism_service as M


def _ids(results):
    return [r["id"] for r in results]


def _status(out, mechanism_id):
    return next(r["status"] for r in out["results"] if r["id"] == mechanism_id)


# ---------------------------------------------------------------------------
# Routing: the goal is an output
# ---------------------------------------------------------------------------

def test_nusinersen_answer_surfaces_without_being_told_the_goal():
    """The failure the whole refactor exists to remove.

    SMA's therapeutic intent is to raise SMN protein, which a user reads as
    upregulation (TG02). The mechanism that actually works is exon inclusion,
    which lives under RNA processing (TG04). Under goal-first routing a user
    who picked TG02 never saw A8 at all.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SMN2", molecular_defect="exon_inclusion_defect",
        delivery_context="cns"))

    top = next(r for r in out["results"] if r["status"] == A.ELIGIBLE)
    assert top["id"] == "A8"
    # The goal is derived from the winner, not supplied by the caller.
    assert out["therapeuticGoal"] == "TG04"


def test_goal_filter_narrows_the_result_but_not_the_scoring():
    """A filter applied after scoring must not change what won."""
    ctx = dict(gene_symbol="TTR", molecular_defect="gain_of_function")
    unfiltered = A.arbitrate(A.ArbitrationContext(**ctx))
    filtered = A.arbitrate(A.ArbitrationContext(**ctx, goal_filter=["TG01"]))

    assert filtered["goalFilterApplied"] == ["TG01"]
    assert set(_ids(filtered["results"])) <= set(_ids(unfiltered["results"]))
    # Every surviving mechanism keeps the exact score it had in the full pass.
    scores = {r["id"]: r["applicability"] for r in unfiltered["results"]}
    for r in filtered["results"]:
        assert r["applicability"] == scores[r["id"]]


def test_every_goal_route_agrees_with_the_unified_pass():
    """Mechanisms shared by TG04 and TG07 must score identically.

    TG07 is no longer a strict subset of TG04 — it has A32 and A33 of its own,
    and A11 is dual-tagged — but A7-A11 appear under both. The property that
    matters is unchanged: a mechanism cannot score differently depending on
    which page asked, because both routes filter one shared pass rather than
    running their own scorer.
    """
    processing = M.rank_rna_processing_mechanisms(
        "exon_skipping_mutation", None, None, None)
    isoform = M.rank_isoform_engineering_mechanisms("exon_skipping")

    shared = {r["id"] for r in processing} & {r["id"] for r in isoform}
    assert shared, "TG07 should overlap TG04"
    by_id = {r["id"]: r for r in processing}
    for r in isoform:
        if r["id"] in shared:
            assert r["status"] == by_id[r["id"]]["status"]
            assert r["applicability"] == by_id[r["id"]]["applicability"]


# ---------------------------------------------------------------------------
# Gates reject; missing evidence halts. These are different answers.
# ---------------------------------------------------------------------------

def test_gate_failure_rejects_rather_than_scoring_low():
    """A15 must not compete on a therapeutic-reduction case.

    Transcriptional silencing of a normal, physiologically required gene is
    not an established strategy, and A15's own rulebook scopes it to
    pathogenic overexpression. Ranking it low would still leave it in the
    running; it is out.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="APOB", molecular_defect="therapeutic_reduction"))
    assert _status(out, "A15") == A.REJECTED


def test_a15_must_not_outrank_a1_on_therapeutic_reduction():
    """The Phase 0 mipomersen finding, locked down."""
    results = M.rank_gene_silencing_mechanisms(
        "therapeutic_reduction", "total_knockdown", "liver", None)
    ranked = _ids(results)
    assert ranked.index("A1") < ranked.index("A15")


def test_mechanisms_with_no_wired_feature_source_halt():
    """A28, A11 and A14 have no way to establish their required evidence.

    F11 (repressive RBP site), F13 (polyadenylation usage) and F12 (repeat
    expansion, with nothing supplied) have no wired source. The plan is
    explicit: halt, do not guess.
    """
    for defect, mechanism in (
        ("rbp_mediated_repression", "A28"),
        ("apa_dysregulation", "A11"),
        ("toxic_rna_gain_of_function", "A14"),
    ):  # noqa: E501
        out = A.arbitrate(A.ArbitrationContext(
            gene_symbol="X", molecular_defect=defect))
        assert _status(out, mechanism) == A.HALTED, mechanism


def test_naming_an_rbp_is_not_evidence_that_a_repressive_site_exists():
    """A28 must still halt when the user names a target RBP.

    Saying which protein you have in mind is not a finding about this
    transcript, and F11's source ladder is deliberately empty.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="X", molecular_defect="rbp_mediated_repression",
        extras={"targetRbp": "PTB"}))
    assert _status(out, "A28") == A.HALTED


def test_a14_scores_once_a_repeat_is_actually_supplied():
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="DMPK", molecular_defect="toxic_rna_gain_of_function",
        repeat_unit="CTG", repeat_count=">50 copies"))
    assert _status(out, "A14") == A.ELIGIBLE


def test_subpathogenic_repeat_count_rejects_rather_than_halts():
    """"We looked and it is not there" is not "we could not look"."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="DMPK", molecular_defect="toxic_rna_gain_of_function",
        repeat_unit="CTG", repeat_count="10 copies"))
    assert _status(out, "A14") == A.REJECTED


def test_absent_feature_never_becomes_probability_zero():
    ctx = F.FeatureContext(molecular_defect="toxic_rna_gain_of_function",
                           repeat_unit="CTG", repeat_count="10 copies")
    f12 = F.resolve_features(ctx)["F12"]
    assert f12.state == F.ABSENT
    assert f12.probability is not None and f12.probability > 0


def test_unresolved_feature_carries_no_probability_at_all():
    f11 = F.resolve_features(F.FeatureContext())["F11"]
    assert f11.state == F.UNRESOLVED
    assert f11.probability is None


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def test_dropdown_derived_evidence_is_marked_as_the_users_own_input():
    """The defect dropdown standing in for SpliceAI must be visible as such.

    Otherwise a top-1 of 1.00 on TG04 reads as arbitration when it is a
    lookup: the input already contains the answer.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="DMD", molecular_defect="exon_skipping_mutation"))
    a7 = next(r for r in out["results"] if r["id"] == "A7")
    assert a7["standInOnly"] is True
    assert out["features"]["F1"]["standIn"] is True
    assert out["features"]["F1"]["provenance"] == F.USER_ASSERTED


def test_user_asserted_evidence_caps_confidence_below_annotation():
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="DMD", molecular_defect="exon_skipping_mutation"))
    a7 = next(r for r in out["results"] if r["id"] == "A7")
    assert a7["confidence"]["upper"] <= F.PROVENANCE_CAP[F.USER_ASSERTED]


def test_unverifiable_gene_features_do_not_read_as_positive_findings():
    """The gene-feature payload reports unverifiable genes as available=True
    so the UI does not silently drop them. That must not become evidence."""
    payload = {"features": {"NAT": {"available": True, "verified": False,
                                    "reason": "Could not verify"}}}
    ctx = F.FeatureContext(molecular_defect="nat_mediated_repression",
                           gene_features=payload)
    f6 = F.resolve_features(ctx)["F6"]
    # Falls through the annotation rung to the user-asserted one rather than
    # being read as a confirmed antisense transcript.
    assert f6.provenance == F.USER_ASSERTED


# ---------------------------------------------------------------------------
# The F10 split (item 10)
# ---------------------------------------------------------------------------

def test_f10_split_separates_a1_from_a2_when_a_sequence_is_available():
    """A1 needs an accessible site anywhere; A2 needs one at the 5' UTR.

    Without a transcript sequence neither query can run and the two tie, which
    is the state that produces outright_top1 = 0.545.
    """
    import random

    import pytest
    # ViennaRNA is an optional heavy dependency. feature_service leaves F10a
    # and F10b unresolved without it, which is the correct degraded behaviour
    # but makes this test vacuous.
    pytest.importorskip("RNA", reason="ViennaRNA not installed")

    random.seed(11)
    seq = "".join(random.choice("ACGU") for _ in range(1500))

    tied = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TTR", molecular_defect="gain_of_function"))
    a1, a2 = (next(r for r in tied["results"] if r["id"] == m) for m in ("A1", "A2"))
    assert a1["applicability"] == a2["applicability"]

    split = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TTR", molecular_defect="gain_of_function",
        transcript_sequence=seq, cds_start=120))
    a1, a2 = (next(r for r in split["results"] if r["id"] == m) for m in ("A1", "A2"))
    assert a1["applicability"] != a2["applicability"]


def test_accessibility_is_a_tie_break_not_a_gate():
    """A1 and A2 stay eligible when no sequence is supplied."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TTR", molecular_defect="gain_of_function"))
    assert _status(out, "A1") == A.ELIGIBLE
    assert _status(out, "A2") == A.ELIGIBLE


# ---------------------------------------------------------------------------
# Mechanism states (v3 §4)
# ---------------------------------------------------------------------------

def test_a21_is_scored_and_competes_despite_being_undesignable_here():
    """siRNA is a real alternative to a gapmer and has five approved drugs.

    It is not designable by a single-stranded designer, but hiding it removed
    an option from a decision the scientist actually makes.
    """
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TTR", molecular_defect="gain_of_function"))
    a21 = next(r for r in out["results"] if r["id"] == "A21")
    assert a21["status"] == A.ELIGIBLE
    assert a21["designAvailable"] is False
    assert a21["designUnavailableReason"]


def test_flagged_mechanisms_never_enter_the_ranking():
    """A24, A25 and A26 have no approved therapy in their indication class,
    so any applicability figure for them would be unverifiable."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="X", molecular_defect="haploinsufficiency"))
    flagged = {"A24", "A25", "A26", "A34", "A35", "A36", "A37", "A38", "A39"}
    assert not flagged & set(_ids(out["results"]))


def test_no_mechanism_was_deleted():
    """Every rulebook still resolves, including the flagged and halted ones."""
    ids = A.all_mechanism_ids()
    # 27 originals + A29-A31 (TG06) + A32/A33 (TG07) + A34-A36 (TG08 flags)
    # + A37-A39 (TG09 flags).
    assert len(ids) == 38
    assert "A22" not in ids, "A22 has never existed; see the implementation notes"
    for added in ("A29", "A30", "A31", "A32", "A33",
                  "A34", "A35", "A36", "A37", "A38", "A39"):
        assert added in ids
    for mid in ids:
        assert A.load_rule(mid)["arbitration"], mid


# ---------------------------------------------------------------------------
# The modality flag (v3 §0) — the TG02 correctness fix
# ---------------------------------------------------------------------------

def _replacement_flag(out):
    return next(f for f in out["modalityFlags"] if f["flag"] == "protein_replacement")


def test_replacement_flag_fires_when_no_boosting_mechanism_is_viable():
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SCN1A", molecular_defect="haploinsufficiency",
        tissue_tpm=0.2))
    flag = _replacement_flag(out)
    assert flag["raised"] is True
    assert flag["scored"] is False
    assert "does not evaluate or design replacement" in flag["message"]


def test_replacement_flag_is_withheld_when_transcript_is_still_boostable():
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SCN1A", molecular_defect="haploinsufficiency",
        tissue_tpm=40.0))
    flag = _replacement_flag(out)
    assert flag["raised"] is False
    assert flag["message"] is None
    assert flag["withheldBecause"]


def test_replacement_flag_is_withheld_when_expression_is_unknown():
    """Silence is not evidence. With no TPM, P2 cannot be established."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SCN1A", molecular_defect="haploinsufficiency"))
    assert _replacement_flag(out)["raised"] is False


def test_flag_never_carries_a_probability():
    """v3's central decision: emit text, not an unverifiable number."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="SCN1A", molecular_defect="haploinsufficiency",
        tissue_tpm=0.2))
    for flag in out["modalityFlags"]:
        assert flag["scored"] is False
        assert "applicability" not in flag
        assert "confidence" not in flag


def test_aptamer_flag_needs_an_extracellular_target():
    intracellular = A.arbitrate(A.ArbitrationContext(
        gene_symbol="X", molecular_defect="mirna_dysregulation",
        transcript_class="protein_coding", protein_localisation="Nucleus"))
    flag = next(f for f in intracellular["modalityFlags"] if f["flag"] == "aptamer")
    assert flag["raised"] is False


# ---------------------------------------------------------------------------
# Data quality (v3 §6.1)
# ---------------------------------------------------------------------------

def test_a24_does_not_claim_vaccines_as_protein_replacement_precedent():
    """Comirnaty and Spikevax are vaccines, a different indication class.

    Listing them gave A24 a Very High evidence rating — the same tier as A21
    with five approved siRNA drugs — and evidence rating feeds the confidence
    cap.
    """
    rule = A.load_rule("A24")
    assert rule["evidenceLevel"]["rating"] != "Very High"
    assert "No mRNA protein replacement therapy has been approved" in \
        rule["fdaApprovedDrugs"]


def test_no_rulebook_claims_an_fda_drug_it_does_not_have():
    """Every mechanism whose evidence tier is the top one must name drugs."""
    for mid in A.all_mechanism_ids():
        rule = A.load_rule(mid)
        rating = (rule.get("evidenceLevel") or {}).get("rating", "")
        if rating.strip().lower() == "very high":
            drugs = rule.get("fdaApprovedDrugs") or ""
            assert drugs and not drugs.lower().startswith("none"), mid


# ---------------------------------------------------------------------------
# Reference tables (data_sources_halted_flagged.md)
# ---------------------------------------------------------------------------

def test_populated_reference_rows_all_name_their_source():
    """Rows are allowed; recalled rows are not.

    The rule was never "these tables stay empty" — it is that nothing may be
    written from memory. Populated tables are fetched by
    data_curation/populate_reference_tables.py and stamped, so the check is
    that every row names where it came from.
    """
    import csv

    from services import reference_tables as RT

    for name, info in RT.status().items():
        if not info["populated"]:
            continue
        path = Path(RT.REFERENCE_DIR) / f"{name}.tsv"
        with open(path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                assert (row.get("source") or "").strip(), name
                assert (row.get("source_version") or "").strip(), name


def test_a11_needs_both_halves_of_f13():
    """Site presence alone would fire A11 almost everywhere.

    Most human genes carry an alternative polyadenylation site, so F13a
    without F13b is not evidence that shifting usage helps.
    """
    arb = A.load_rule("A11")["arbitration"]
    assert set(arb["requiredFeatures"]) == {"F13a", "F13b"}


def test_clingen_haploinsufficiency_score_is_not_treated_as_ordinal():
    """Code 40 means "dosage sensitivity unlikely", not "strongest evidence".

    Comparing the column as an integer would invert the meaning of the two
    highest-numbered codes.
    """
    assert F.CLINGEN_HI_SUFFICIENT == "3"
    assert F.CLINGEN_HI_UNLIKELY == "40"
    assert F.CLINGEN_HI_AUTOSOMAL_RECESSIVE == "30"
    # The codes are strings precisely so no ordinal comparison compiles.
    assert all(isinstance(c, str) for c in (
        F.CLINGEN_HI_SUFFICIENT, F.CLINGEN_HI_UNLIKELY,
        F.CLINGEN_HI_AUTOSOMAL_RECESSIVE))


def test_curated_entries_outrank_bulk_annotation():
    """A hand-curated, citation-carrying row is trusted above a bulk lookup."""
    assert F.PROVENANCE_CAP[F.CONFIRMED] > F.PROVENANCE_CAP[F.ANNOTATION]
    assert F.PROVENANCE_CAP[F.ANNOTATION] > F.PROVENANCE_CAP[F.PREDICTED]
    assert F.PROVENANCE_CAP[F.PREDICTED] > F.PROVENANCE_CAP[F.USER_ASSERTED]


# ---------------------------------------------------------------------------
# SpliceAI-backed F1/F2/F3
# ---------------------------------------------------------------------------

def _synthetic_pre_mrna(seed=7, ex1=150, intron_body=300, ex2=150):
    """Exon | intron with canonical GT..AG | exon."""
    import random
    rng = random.Random(seed)
    e1 = "".join(rng.choice("ACGT") for _ in range(ex1))
    intron = ("GTAAGT" + "".join(rng.choice("ACGT") for _ in range(intron_body))
              + "TTTTTTTTTTGCAG")
    e2 = "".join(rng.choice("ACGT") for _ in range(ex2))
    return e1, intron, e2


def _require_spliceai():
    import pytest
    from services import spliceai_service as SAI
    ok, why = SAI.available()
    if not ok:
        pytest.skip(f"SpliceAI unavailable: {why}")
    return SAI


def test_spliceai_finds_the_acceptor_where_it_was_placed():
    """Ground truth we control: the model must locate a canonical AG.

    Without this, every downstream F1/F2/F3 assertion is testing plumbing
    rather than predictions.
    """
    import numpy as np
    SAI = _require_spliceai()
    e1, intron, e2 = _synthetic_pre_mrna()
    pre = e1 + intron + e2
    probs = SAI.predict(pre)
    assert probs is not None

    acceptor = probs[:, 1]
    true_pos = len(e1) + len(intron)   # 3' end of the intron

    # Not "is it the global maximum" — random filler contains plenty of AG
    # dinucleotides and a spurious site can outscore the planted one on some
    # seeds. The robust claim is that the planted acceptor is among the very
    # top positions and sits orders of magnitude above the background.
    rank = int((acceptor > acceptor[true_pos]).sum())
    assert rank < 5, f"planted acceptor ranked {rank} of {len(acceptor)}"
    assert acceptor[true_pos] > 100 * float(np.median(acceptor))


def test_f1_resolves_from_spliceai_not_the_user_dropdown():
    """With a pre-mRNA supplied, the answer stops resting on the user's input.

    This is the whole point of wiring SpliceAI: standInOnly must go False and
    the confidence cap must rise from user_asserted to predicted.
    """
    _require_spliceai()
    e1, intron, e2 = _synthetic_pre_mrna()
    pre = e1 + intron + e2
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TESTG", molecular_defect="exon_skipping_mutation",
        pre_mrna_sequence=pre,
        exon_start=len(e1 + intron), exon_end=len(pre) - 1))

    f1 = out["features"]["F1"]
    assert f1["provenance"] == F.PREDICTED
    assert f1["standIn"] is False
    a7 = next(r for r in out["results"] if r["id"] == "A7")
    assert a7["standInOnly"] is False
    assert a7["confidence"]["upper"] <= F.PROVENANCE_CAP[F.PREDICTED]


def test_f1_falls_back_to_the_stand_in_without_a_pre_mrna():
    """No sequence, no prediction — the documented stand-in still applies."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="DMD", molecular_defect="exon_skipping_mutation"))
    assert out["features"]["F1"]["standIn"] is True


def test_f1_does_not_score_a_boundary_that_does_not_exist():
    """A terminal exon has no downstream donor.

    Scoring the missing side reads 0.00 and would mark every first and
    terminal exon "weakly recognised" — a prediction about a splice site
    that is not there.
    """
    _require_spliceai()
    e1, intron, e2 = _synthetic_pre_mrna()
    pre = e1 + intron + e2
    first = F.resolve_features(F.FeatureContext(
        pre_mrna_sequence=pre, exon_start=0, exon_end=len(e1) - 1))["F1"]
    assert "not scored" in (first.detail or "")
    # The real donor at the end of exon 1 is strong, so it must NOT be
    # reported as weakly recognised on the strength of an absent acceptor.
    assert first.state == F.ABSENT


def test_f2_requires_a_gain_not_merely_a_change():
    """Destroying an acceptor is a loss; F2 asks whether one was CREATED."""
    _require_spliceai()
    e1, intron, e2 = _synthetic_pre_mrna()
    pre = e1 + intron + e2
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="TESTG", molecular_defect="cryptic_splice_site",
        pre_mrna_sequence=pre,
        variant_offset=len(e1) + len(intron) - 2, variant_alt="T"))
    f2 = out["features"]["F2"]
    assert f2["state"] == F.ABSENT
    assert f2["provenance"] == F.PREDICTED


def test_delta_scores_refuse_length_changing_variants():
    """An indel shifts the coordinate frame, so a positionwise diff would
    report an artefact of the shift rather than a splicing change."""
    SAI = _require_spliceai()
    e1, intron, e2 = _synthetic_pre_mrna()
    pre = e1 + intron + e2
    assert SAI.delta_scores(pre, pre + "A") is None


def test_aptamer_candidates_carry_no_fabricated_numbers():
    """TG09 returns guidance strings, never invented measurements.

    No aptamer has been selected, so there is no sequence to fold, no Tm, no
    dG and no measured Kd. The endpoint previously derived all of those from
    hash() of the form inputs and rendered them as measurements.
    """
    from api.mechanisms import _aptamer_candidate_from_rule

    for mid in ("A25", "A37", "A38", "A39"):
        cand = _aptamer_candidate_from_rule(mid, 1)
        for numeric in ("tm", "deltaGFolding", "targetSpecificityScore",
                        "sequence", "dotBracket", "foldingScore", "tHalfScore"):
            assert cand[numeric] is None, f"{mid}.{numeric} must not be invented"
        assert cand["scored"] is False
        # Guidance is a string about what SELEX would determine, not a number.
        assert isinstance(cand["kdPrediction"], str)


def test_flagged_mechanisms_are_always_listed_even_when_no_flag_fires():
    """Unscorable is not a reason to be unlistable.

    A24/A26/A34-A36 and A25/A37-A39 were reachable only through a raised
    modality flag, and a flag needs tissue expression (P2) or subcellular
    localisation (B1) — inputs the page does not collect. On a
    haploinsufficiency case, where every transcript-acting mechanism halts and
    protein replacement is the obvious move, all nine were invisible.
    """
    # A gene with no reference-table coverage, so no flag can fire.
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="NOT_A_REAL_GENE", molecular_defect="haploinsufficiency"))

    # No flag fires ...
    assert all(not f["raised"] for f in out["modalityFlags"])
    # ... but the mechanisms are listed regardless. That is the point: their
    # visibility must not depend on data the page never collects.
    listed = {m["id"] for m in out["flaggedMechanisms"]}
    assert listed == {"A24", "A25", "A26", "A34", "A35", "A36",
                      "A37", "A38", "A39"}


def test_flagged_mechanisms_carry_no_score_shaped_numbers():
    """An empty applicability interval degenerates to [1.0, 1.0], which reads
    as 'perfectly applicable' — the opposite of what FLAGGED means."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="CFTR", molecular_defect="haploinsufficiency"))
    for m in out["flaggedMechanisms"]:
        assert m["status"] == A.FLAGGED
        assert m["applicability"] is None
        assert m["confidence"] is None
        assert m["score"] is None
        assert m["rationale"], "a flagged mechanism must say why it is flagged"


def test_flagged_mechanisms_stay_out_of_the_ranking():
    """Listed alongside, never ranked against a scored mechanism."""
    out = A.arbitrate(A.ArbitrationContext(
        gene_symbol="CFTR", molecular_defect="haploinsufficiency"))
    ranked = {r["id"] for r in out["results"]}
    flagged = {m["id"] for m in out["flaggedMechanisms"]}
    assert not (ranked & flagged)
    assert all(r["score"] is not None for r in out["results"])
