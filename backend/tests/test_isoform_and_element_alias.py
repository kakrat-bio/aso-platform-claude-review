"""Regressions for two defects found by probing the live goal endpoints.

1. `/api/isoform-engineering/generate` raised NameError on every request
   (`_calc_cai` / `_calc_u_content` were never defined in that module), so
   TG07 candidate generation had never returned anything. When it did run it
   emitted a hard-coded sequence and index-derived numbers.

2. `/api/mechanisms/options` advertises the uORF target element as `uorf`
   while the TG06 design service calls it `5p_uorf`. The frontend reads the
   options list and posts the value straight back, so the default element
   400'd.
"""

import ast
import json
import pathlib

import pytest

from api.translational_regulation import (
    ELEMENT_ALIASES, VALID_ELEMENTS, normalise_element,
)
from services import isoform_engineering_service as iso

SERVICE = pathlib.Path(iso.__file__)


def test_isoform_service_calls_no_undefined_names():
    """The NameError that made TG07 unreachable must not come back."""
    tree = ast.parse(SERVICE.read_text())
    defined = {n.name for n in ast.walk(tree)
               if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            defined |= {a.asname or a.name for a in n.names}
        elif isinstance(n, ast.Import):
            defined |= {(a.asname or a.name).split(".")[0] for a in n.names}
    import builtins
    called = {n.func.id for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    missing = sorted(c for c in called
                     if c not in defined and not hasattr(builtins, c))
    assert not missing, f"undefined names called: {missing}"


def _nucleotide_literals_in_code(path: pathlib.Path) -> list[str]:
    """Nucleotide-only string constants that reach executable code.

    Checked over the AST with docstrings removed, not over raw text — the
    module's own docstring quotes the literals it used to emit, and a
    substring search would flag the explanation as the offence.
    """
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if (isinstance(first, ast.Expr)
                    and isinstance(first.value, ast.Constant)
                    and isinstance(first.value.value, str)):
                node.body = node.body[1:] or [ast.Pass()]
    def is_nt(node) -> bool:
        return (isinstance(node, ast.Constant) and isinstance(node.value, str)
                and node.value
                and set(node.value.upper()) <= set("ACGTU"))

    found = []
    for node in ast.walk(tree):
        # A long nucleotide literal is a sequence, not an alphabet. Short
        # constants like "ACGU" or "GCgc" are membership sets and are fine.
        if is_nt(node) and len(node.value) >= 8:
            found.append(node.value)
        # `"GCU" * 20` is the repetition pattern that produced the old
        # placeholder transcript, at any length.
        if (isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult)
                and (is_nt(node.left) or is_nt(node.right))):
            found.append(ast.unparse(node) if hasattr(ast, "unparse")
                         else "<nucleotide literal * n>")
    return found


def test_isoform_service_has_no_hardcoded_sequence():
    """No placeholder transcript may be substituted for a real fetch."""
    literals = _nucleotide_literals_in_code(SERVICE)
    assert not literals, (
        f"nucleotide string literals reach executable code in "
        f"{SERVICE.name}: {literals}. The standing rule is that a missing "
        f"transcript is reported, never synthesised."
    )


def test_isoform_refuses_when_transcript_unavailable(monkeypatch):
    """No exons means no candidates and a stated reason — not invented ones."""
    monkeypatch.setattr(iso, "_fetch_exon_splice_sites",
                        lambda *a, **k: {"exons": [], "cdsSequence": "",
                                         "spliceSites": [],
                                         "canonicalTranscript": ""})
    out = iso.generate_isoform_candidates(
        target_symbol="NOSUCHGENE", isoform_goal="exon_skipping",
        target_exon_locus="exon_7", splice_element_target="splice_donor",
        steric_chemistry="pmo")
    assert out["status"] == "UNAVAILABLE"
    assert out["candidates"] == []
    assert "invented" in out["message"] or "No transcript" in out["message"]


def test_isoform_candidates_are_reverse_complements_of_the_real_window(monkeypatch):
    """Every emitted ASO must derive from the supplied transcript."""
    mrna = "ACGU" * 40                      # 160 nt of unambiguous sequence
    monkeypatch.setattr(iso, "_fetch_exon_splice_sites", lambda *a, **k: {
        "exons": [{"index": 1, "cdsStart": 0, "cdsEnd": 90}],
        "cdsSequence": mrna, "spliceSites": [],
        "canonicalTranscript": "ENST00000000001", "geneId": "ENSG00000000001",
    })
    out = iso.generate_isoform_candidates(
        target_symbol="TESTGENE", isoform_goal="exon_skipping",
        target_exon_locus="exon_1",
        splice_element_target="exonic_splicing_enhancer",
        steric_chemistry="pmo", aso_length=20)
    assert out["status"] == "OK"
    assert out["candidates"]
    for c in out["candidates"]:
        window = mrna[c["transcriptStart"]:c["transcriptEnd"]]
        assert c["targetSequence"] == window
        assert c["sequence"] == iso._reverse_complement(window)
        assert len(c["sequence"]) == 20
        # Nothing is invented for the two quantities with no model behind them.
        assert c["predictedIsoformYield"] is None
        assert c["tlrRisk"] is None
        assert set(c["notComputed"]) >= {"predictedIsoformYield", "tlrRisk"}


def test_isoform_frame_status_is_real_arithmetic(monkeypatch):
    """In-frame is exon length % 3, not a loop index."""
    mrna = "ACGU" * 40
    for cds_end, expected in ((90, "In-Frame"), (91, "Out-of-Frame")):
        monkeypatch.setattr(iso, "_fetch_exon_splice_sites", lambda *a, _e=cds_end, **k: {
            "exons": [{"index": 1, "cdsStart": 0, "cdsEnd": _e}],
            "cdsSequence": mrna, "spliceSites": [],
            "canonicalTranscript": "ENST00000000001",
        })
        out = iso.generate_isoform_candidates(
            target_symbol="TESTGENE", isoform_goal="exon_skipping",
            target_exon_locus="exon_1",
            splice_element_target="exonic_splicing_enhancer",
            steric_chemistry="pmo", aso_length=20)
        assert out["overview"]["inFrameStatus"] == expected, cds_end


@pytest.mark.parametrize("advertised,expected", [
    ("uorf", "5p_uorf"),
    ("5p_uorf", "5p_uorf"),
    ("3p_utr_mirna", "3p_utr_mirna"),
    ("5p_utr", "5p_utr"),
])
def test_translational_element_aliases_resolve(advertised, expected):
    assert normalise_element(advertised) == expected
    assert expected in VALID_ELEMENTS


def test_every_advertised_target_element_is_accepted():
    """Whatever /api/mechanisms/options offers must be postable to TG06."""
    from services.mechanism_service import TRANSLATIONAL_TARGET_ELEMENTS
    for advertised in TRANSLATIONAL_TARGET_ELEMENTS:
        assert normalise_element(advertised) in VALID_ELEMENTS, (
            f"/api/mechanisms/options advertises target element "
            f"{advertised!r}, which the translational-regulation design "
            f"endpoint rejects. Add it to ELEMENT_ALIASES "
            f"(currently {sorted(ELEMENT_ALIASES)})."
        )


# ---------------------------------------------------------------------------
# TG04 designer and the upload path
# ---------------------------------------------------------------------------

def test_reverse_complement_pairs_uracil():
    """U used to pass through uncomplemented, so TG04 oligos did not bind."""
    from services.gene_silencing_service import _reverse_complement as rc
    pair = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"}
    for target in ("AAUUAAUAGGAUUAUUAG", "UUACCAACAGGACCACCAG",
                   "AACCAACAGGACCACCAG", "ATGCATGCATGC"):
        got = rc(target)
        expected = "".join(pair[b] for b in reversed(target))
        assert got.replace("T", "U") == expected.replace("T", "U"), target
        # The alphabet must not be mixed: T and U cannot both appear.
        assert not ("T" in got and "U" in got), f"mixed alphabet: {got}"


def test_reverse_complement_refuses_ambiguity_codes():
    from services.gene_silencing_service import _reverse_complement as rc
    with pytest.raises(ValueError):
        rc("ACGN")


def test_calc_tm_accepts_rna_alphabet():
    """primer3 raises on U; TG04 renders its windows in the RNA alphabet."""
    from services.gene_silencing_service import _calc_tm
    assert _calc_tm("GAAAUAUUCCUUAUAGCC") == _calc_tm("GAAATATTCCTTATAGCC")
    with pytest.raises(ValueError):
        _calc_tm("GAAANATTCC")


def test_grna_scanner_finds_real_pams():
    """/^NGG$/ matched the literal characters "NGG", so this found nothing."""
    from services.upload_service import _generate_grna_candidates
    seq = ("GACGTTGCAGGTACCATGGCTAGCTAGGTACCGGTAGCTAGCTAGCTAGGTTACGATCGATCGG"
           "ATCGATCGGCTAGCTAGCTAGCTAAGGCTTGCATGCATGCAGGTACGT")
    candidates = _generate_grna_candidates(seq)
    assert candidates, "no gRNA found in a sequence containing NGG PAMs"
    for c in candidates:
        assert len(c["sequence"]) == 20
        assert c["pam"][1:] == "GG" and c["pam"][0] in "ACGT"
        # No fabricated off-target count may reappear.
        assert "offTargets" not in c
        assert 0.0 <= c["internalRepetitiveness"] <= 1.0


def test_scorecard_excludes_sequence_independent_constants():
    """capEfficiency=70 and nucleosideMod=90 ignored the sequence entirely."""
    from services.upload_service import _modification_scorecard
    a = _modification_scorecard("AUGGCUAGCUAGCUAGC" + "A" * 12, "mrna")
    b = _modification_scorecard("GGGGGGCCCCCCGGGGGG", "mrna")
    assert "capEfficiency" not in a["scores"]
    assert "nucleosideMod" not in a["scores"]
    assert {adv["id"] for adv in a["advisories"]} == {"capEfficiency",
                                                     "nucleosideMod"}
    # overallScore must now move with the sequence.
    assert a["overallScore"] != b["overallScore"]


# ---------------------------------------------------------------------------
# TG08 protein replacement
# ---------------------------------------------------------------------------

def test_protein_replacement_rejects_unknown_choices():
    """Unknown enum values used to fall through and be designed as `linear`."""
    from services.protein_replacement_service import (
        VALID_RNA_MODALITIES, generate_protein_replacement_candidates,
    )
    base = dict(target_symbol="CFTR", rna_modality="linear",
                codon_strategy="cai", utr_pair="globin",
                ires_selection=None, nucleotide_modification="m1psi")
    for field, bad in (("rna_modality", "NONSENSE"),
                       ("codon_strategy", "BOGUS"),
                       ("utr_pair", "FAKE"),
                       ("nucleotide_modification", "INVALID")):
        with pytest.raises(ValueError) as exc:
            generate_protein_replacement_candidates(**{**base, field: bad})
        assert field in str(exc.value)
    assert "any" in VALID_RNA_MODALITIES


def test_protein_replacement_utrs_come_from_the_working_endpoint(monkeypatch):
    """`?type=utr5` is a 400 from Ensembl; the UTRs must come via cDNA/CDS.

    The old implementation asked for `/sequence/id/{tx}?type=utr5`, which
    Ensembl answers with `{"error":"The type 'utr5' is not understood by this
    service"}`. It therefore returned None for every gene and every construct
    was emitted with `utrSource: "absent"`.
    """
    from services import protein_replacement_service as prs

    calls: list[str] = []

    def fake_target(gene_id, gene_symbol="", organism="homo_sapiens"):
        calls.append(gene_symbol or gene_id)
        return {
            "utr5Sequence": "GCCACCATGG" * 7,
            "utr3Sequence": "TTTATTTAAA" * 20,
            "canonicalTranscript": {"id": "ENST00000003084.11"},
        }

    monkeypatch.setattr(prs, "get_target_analysis", fake_target)
    out = prs._fetch_real_utrs("CFTR")
    assert out is not None, "UTR lookup returned nothing"
    assert len(out["utr5"]) == 70
    assert len(out["utr3"]) == 200
    assert out["transcript_id"] == "ENST00000003084"
    assert calls, "did not go through get_target_analysis"


def _string_constants_in_code(path: pathlib.Path) -> list[str]:
    """Every string constant that reaches executable code, docstrings removed."""
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if (isinstance(first, ast.Expr)
                    and isinstance(first.value, ast.Constant)
                    and isinstance(first.value.value, str)):
                node.body = node.body[1:] or [ast.Pass()]
    return [n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)]


def test_protein_replacement_source_has_no_utr_type_request():
    """Guard against the unsupported Ensembl parameter coming back.

    Checked over string constants with docstrings stripped — the docstring
    explaining why `type=utr5` is wrong necessarily contains it.
    """
    from services import protein_replacement_service as prs
    constants = _string_constants_in_code(pathlib.Path(prs.__file__))
    for bad in ("type=utr5", "type=utr3"):
        offenders = [c for c in constants if bad in c]
        assert not offenders, (
            f"{bad} is an unsupported Ensembl sequence type and returns "
            f"HTTP 400; use the cDNA/CDS alignment instead. Found in: "
            f"{offenders}"
        )


# ---------------------------------------------------------------------------
# Mechanism-to-designer coverage, and the ADMET withdrawal
# ---------------------------------------------------------------------------

def test_editing_designer_refuses_mechanisms_it_cannot_build():
    """`mechanism_id` used to be echoed into the output and change nothing.

    Asking for A19 returned an ADAR-recruiting guide labelled "A19" — a REPAIR
    crRNA in name only, since REPAIR guides carry a Cas13b direct-repeat
    scaffold this service never emits.
    """
    from services.rna_editing_service import _validate_editing_mechanism
    _validate_editing_mechanism("A13", "a_to_i")
    _validate_editing_mechanism("A17", "a_to_i")
    _validate_editing_mechanism("A16", "c_to_u")
    _validate_editing_mechanism("A20", "trans_splicing")
    for mech in ("A18", "A19"):
        with pytest.raises(ValueError, match="not designable"):
            _validate_editing_mechanism(mech, "a_to_i")
    with pytest.raises(ValueError, match="does not perform"):
        _validate_editing_mechanism("A13", "c_to_u")


def test_every_designable_mechanism_has_a_designer():
    """A mechanism advertising designAvailable must be built by something."""
    import glob
    from services.rna_processing_service import RNA_PROCESSING_MECHANISMS
    from services.translational_regulation_service import (
        TRANSLATIONAL_MECHANISM_CHEMISTRY,
    )
    from services.gene_upregulation_service import UPREGULATION_MECHANISM_DESIGN
    from services.rna_neutralization_service import (
        NEUTRALIZATION_MECHANISM_CHEMISTRY,
    )
    from services.rna_editing_service import EDITING_MECHANISM_EDIT_TYPES
    from services.programmable_editor_service import PLATFORMS as EDITOR_PLATFORMS

    covered = (set(RNA_PROCESSING_MECHANISMS)
               | set(TRANSLATIONAL_MECHANISM_CHEMISTRY)
               | set(UPREGULATION_MECHANISM_DESIGN)
               | set(NEUTRALIZATION_MECHANISM_CHEMISTRY)
               | set(EDITING_MECHANISM_EDIT_TYPES)
               # gene_silencing_service._mechanism_design_constraints
               | {"A1", "A2", "A12", "A15"}
               # Dedicated designers added for mechanisms the goal services
               # cannot build: siRNA duplexes, protein-dependent RNA editors,
               # and the two pre-mRNA mechanisms that need genomic sequence.
               | {"A21"}
               | set(EDITOR_PLATFORMS)
               | {"A32", "A33"})

    orphans = []
    for path in glob.glob("rulebooks/A*/rule.json"):
        rule = json.loads(pathlib.Path(path).read_text())
        arb = rule.get("arbitration", {})
        mech = rule.get("mechanismId") or pathlib.Path(path).parent.name
        # A mechanism is covered either by a designer whitelist above, or by
        # an explicit designRoute naming the endpoint and the parameters that
        # build it (A24/A26 route into the protein-replacement designer's
        # linear and circRNA architectures).
        route = arb.get("designRoute")
        if route:
            assert route.get("endpoint"), f"{mech} designRoute names no endpoint"
            continue
        if arb.get("designAvailable") and mech not in covered:
            orphans.append(mech)
    assert not orphans, (
        f"these mechanisms advertise designAvailable=True but no designer "
        f"accepts them, so a request for one silently returns another "
        f"mechanism's candidates: {sorted(orphans)}"
    )


def test_admet_endpoints_are_withdrawn_with_reasons():
    """Sequence-independent ADMET endpoints must not come back silently."""
    from services.sequence_liability_service import (
        NOT_ASSESSED, get_sequence_liabilities,
    )
    out = get_sequence_liabilities("GCCGCGGGTTTTCCCGGAAA", chemistry="gapmer")
    assert out["available"] is True
    for gone in ("absorptionScore", "distributionScore", "metabolismScore",
                 "excretionScore", "pbpkTimeSeries", "lipinskiViolations",
                 "chargePhProfile", "hemolysisRisk",
                 "chemicalSpaceProjection", "cellUptake", "renalClearance"):
        assert gone not in out, f"{gone} is back in the payload"
    # Every withdrawal carries a stated biological reason.
    for field, reason in NOT_ASSESSED.items():
        assert len(reason) > 40, f"{field} has no real reason attached"
    assert {"absorption", "distribution", "metabolism", "excretion",
            "halfLife", "lipinskiViolations"} <= set(NOT_ASSESSED)


def test_sequence_liabilities_keeps_what_sequence_determines():
    from services.sequence_liability_service import get_sequence_liabilities
    # CpG-rich, G-quadruplex-forming, uridine-tract-carrying.
    out = get_sequence_liabilities("CGCGCGGGGGAAUUUUCGCG")
    ids = {f["id"] for f in out["immuneAndStructural"]["flags"]}
    assert {"cpg_tlr9", "g_quadruplex", "uridine_tract_tlr78"} <= ids
    for flag in out["immuneAndStructural"]["flags"]:
        assert flag["reasoning"], f"{flag['id']} has no biological reasoning"
    # A clean sequence raises nothing.
    quiet = get_sequence_liabilities("AACAACAACAACAACAACAA")
    assert quiet["immuneAndStructural"]["flags"] == []


def test_sirna_duplex_strands_are_complementary():
    """A21 emits two strands; the guide must be the target's reverse complement."""
    from services.sirna_duplex_service import design_sirna_duplexes, _revcomp
    out = design_sirna_duplexes("ENSG00000197386", "HTT", max_candidates=4)
    assert out["status"] == "OK"
    assert out["candidates"]
    for c in out["candidates"]:
        assert _revcomp(c["passengerCore"]) == c["guideCore"]
        assert len(c["guideCore"]) == 19
        assert c["guideStrand"].endswith(c["overhang"])
        assert c["passengerStrand"].endswith(c["overhang"])
        assert c["seedRegion"] == c["guideCore"][1:8]
    # Ranking is Ui-Tei positional rules first, thermodynamic gap as the
    # tie-break — the rules read Argonaute's MID-pocket preference directly,
    # so they outrank the bulk free-energy difference. This test used to
    # assert a pure asymmetry ordering, which described the earlier design.
    keys = [(-c["uiTeiRules"]["passed"], -c["asymmetryScore"])
            for c in out["candidates"]]
    assert keys == sorted(keys)
    # When the two measures disagree the candidate must say so rather than
    # letting the ranking silently pick a side.
    for c in out["candidates"]:
        if c["uiTeiRules"]["passed"] == c["uiTeiRules"]["total"] and c["asymmetryScore"] <= 0:
            assert any("disagree" in f for f in c["flags"])


def test_editor_guides_refuse_the_wrong_target_base():
    """An A-to-I editor must not be pointed at a non-adenosine."""
    from services.programmable_editor_service import design_editor_guides
    with pytest.raises(ValueError, match="not a protein-dependent"):
        design_editor_guides("A13", "ENSG00000197249", 1096)
    out = design_editor_guides("A19", "ENSG00000197249", edit_position=1096,
                               gene_symbol="SERPINA1", max_candidates=3)
    if out["status"] == "OK":
        for c in out["candidates"]:
            assert c["spacer"][c["mismatchPosition"] - 1] == c["mismatchBase"]
            assert len(c["spacer"]) == c["spacerLength"]
            # The scaffold is named, never invented.
            assert "source" in c["scaffoldRequired"]


def test_editor_scaffold_sequence_is_never_invented():
    """No literal direct-repeat / hairpin sequence may appear in the module."""
    from services import programmable_editor_service as pes
    literals = _nucleotide_literals_in_code(pathlib.Path(pes.__file__))
    assert not literals, (
        f"a scaffold sequence appears to be hard-coded: {literals}. Getting a "
        f"direct repeat wrong by one base yields a guide that does not load, "
        f"so it must come from the construct, not from recall."
    )


def test_no_mechanism_is_silently_unavailable():
    """Every mechanism either has a designer or states why it has none."""
    import glob
    silent = []
    for path in glob.glob("rulebooks/A*/rule.json"):
        rule = json.loads(pathlib.Path(path).read_text())
        arb = rule.get("arbitration", {})
        mech = rule.get("mechanismId") or pathlib.Path(path).parent.name
        if not arb.get("designAvailable") and not arb.get("designUnavailableReason"):
            silent.append(mech)
    assert not silent, (
        f"these mechanisms produce no candidates and give no reason, so a "
        f"user sees an empty result with no explanation: {sorted(silent)}"
    )


# ---------------------------------------------------------------------------
# A1 / A2 — the two mechanisms the platform leans on hardest
# ---------------------------------------------------------------------------

def test_a1_ranking_is_not_a_saturated_tie(monkeypatch):
    """compositeScore clipped at 100 for every candidate of one length.

    `_composite_score` maps duplex dG through `(-dg - 8) * 3.5` capped at 100.
    Those constants suit oligos from 12 to 30 nt, but every candidate in one
    run has the SAME length, so the within-run spread is a few kcal/mol and
    all of it clips. Measured on HTT exons 1-3 at 20 nt: 10 candidates, one
    distinct score, a ten-way tie with the displayed order decided by nothing.
    """
    from services.gene_silencing_service import (
        get_target_analysis, generate_candidates,
    )
    target = get_target_analysis("ENSG00000197386", gene_symbol="HTT")
    cands = generate_candidates(
        [1, 2, 3], 20, "gapmer", [], target["mrnaSequence"], target["exons"],
        "A1", defect_type="gain_of_function", silencing_scope="total_knockdown")
    assert len(cands) >= 5

    accs = [c["realMetrics"].get("siteAccessibility") for c in cands]
    assert all(a is not None for a in accs), "site accessibility not computed"
    # The primary ranking axis must actually separate the candidates.
    assert len(set(accs)) > 1, (
        "every candidate has identical site accessibility — the ranking is "
        "degenerate again")
    # And it must be the axis they are sorted on.
    assert accs == sorted(accs, reverse=True), (
        f"candidates are not ordered by accessibility: {accs}")
    for c in cands:
        assert c["rankingBasis"]["primary"].startswith("siteAccessibility")
        assert "below chance" in c["rankingBasis"]["caveat"]
        assert c["accessibilityPercentile"] is not None


def test_a1_a2_are_separated_by_f10_when_a_transcript_is_supplied():
    """The F10a/F10b split exists to break the A1 vs A2 tie.

    Without a transcript both sit at rulebook evidence. With one, A1 asks
    whether ANY accessible cleavable site exists and A2 asks specifically
    about the translation-initiation region, so a transcript that is open
    elsewhere but structured at the AUG must separate them.
    """
    from services.gene_silencing_service import get_target_analysis
    from services import mechanism_arbitration as MA
    target = get_target_analysis("ENSG00000197386", gene_symbol="HTT")
    out = MA.arbitrate(MA.ArbitrationContext(
        gene_symbol="HTT", molecular_defect="gain_of_function",
        transcript_sequence=target["mrnaSequence"],
        cds_start=len(target.get("utr5Sequence") or ""), oligo_length=20))
    by_id = {r["id"]: r for r in out["results"]}
    assert "A1" in by_id and "A2" in by_id
    assert by_id["A1"]["score"] != by_id["A2"]["score"], (
        "A1 and A2 scored identically with a real transcript supplied; the "
        "F10a/F10b split is not discriminating")


# ---------------------------------------------------------------------------
# A3 / A4 — the two TG02 mechanisms that depend on fetched gene features
# ---------------------------------------------------------------------------

def test_halted_mechanisms_carry_no_score():
    """A halt means a required feature is unresolved — so it has no score.

    HALTED was leaking `score: 1.0, applicability: [1.0, 1.0]` from the
    vacuous interval an empty feature list produces. That asserts perfect
    applicability for the one mechanism the system just said it cannot assess.
    """
    from services import mechanism_arbitration as MA
    out = MA.arbitrate(MA.ArbitrationContext(
        gene_symbol="NOSUCHGENE", molecular_defect="haploinsufficiency"))
    halted = [r for r in out["results"] if r["status"] == MA.HALTED]
    assert halted, "expected at least one halt with no gene features supplied"
    for r in halted:
        assert r["score"] is None, f"{r['id']} halted but reports a score"
        assert r["applicability"] is None
        assert r["confidence"] is None


def test_gene_feature_client_retries():
    """A single transient timeout used to make A3/A4 halt on every gene.

    `_ensembl_request` had no retry: one read timeout returned None, which
    marked TANGO/NAT unverified, which left F4/F6 unresolved, which halted
    both mechanisms. Ensembl answers this lookup in about a second, so the
    failures were transient.
    """
    import services.gene_feature_service as G
    assert G.ENSEMBL_MAX_RETRIES >= 2, "the client must retry"
    calls = {"n": 0}

    class _Resp:
        status_code = 200
        @staticmethod
        def json():
            return {"id": "ENSG00000000001"}

    def flaky(method, url, **kwargs):
        calls["n"] += 1
        if calls["n"] < 2:
            import requests
            raise requests.Timeout("read timed out")
        return _Resp()

    import requests as _requests
    original = _requests.request
    G.requests.request = flaky
    try:
        out = G._ensembl_get("/lookup/id/ENSG00000000001")
    finally:
        G.requests.request = original
    assert out == {"id": "ENSG00000000001"}, "a retryable timeout was not retried"
    assert calls["n"] == 2


def test_a3_a4_score_from_fetched_features_not_an_asserted_defect():
    """SCN1A has a poison exon and an antisense transcript; both are fetchable.

    The endpoint used to pass `payload.gene_features` straight through, so a
    caller that sent none left F4 and F6 unresolved and both mechanisms
    halted — the platform never checked whether the gene actually has these
    features, it only believed the user's defect selection.
    """
    import logging
    logging.disable(logging.WARNING)
    from fastapi.testclient import TestClient
    from main import app
    client = TestClient(app)
    resp = client.post("/api/mechanisms/gene-upregulation", json={
        "gene_symbol": "SCN1A", "defect_type": "haploinsufficiency"})
    assert resp.status_code == 200
    by_id = {m["id"]: m for m in resp.json()["results"]}
    for mech in ("A3", "A4"):
        assert mech in by_id, f"{mech} missing from the TG02 ranking"
        if by_id[mech]["status"] == "HALTED":
            pytest.skip("Ensembl unavailable; the fetched-feature path cannot "
                        "be exercised offline")
        assert by_id[mech]["score"] is not None


def test_gene_feature_payload_shape_is_normalised():
    """`gene_feature()` reads payload["features"][key]; the inner dict alone
    silently resolves nothing, which is an easy call-site mistake."""
    from api.mechanisms import _resolve_gene_features
    inner = {"NAT": {"available": True, "verified": True, "reason": "x"}}
    assert _resolve_gene_features("X", inner)["features"] == inner
    whole = {"features": inner, "source": "live"}
    assert _resolve_gene_features("X", whole) == whole


def test_a3_targets_a_poison_exon_not_every_junction():
    """A poison exon is absent from the canonical mRNA by definition.

    The old path set the label "Exon junctions" and tiled the whole
    transcript, returning 891 candidates across all 29 SCN1A exon junctions —
    none aimed at the exon the mechanism is named for. It could not work from
    the canonical mRNA at all, because the exon is skipped there.
    """
    from services.upregulation_targets_service import (
        MAX_POISON_EXON_NT, design_poison_exon_block, find_poison_exons,
    )
    located = find_poison_exons("ENSG00000144285")
    if located["status"] != "OK":
        pytest.skip("Ensembl unavailable")
    assert located["nmdTranscriptCount"] > 0
    for exon in located["poisonExons"]:
        assert exon["length"] <= MAX_POISON_EXON_NT
        assert exon["supportingTranscripts"]
    # SCN1A's clinically targeted poison exon (20N, the STK-001 target) is
    # 64 nt and must be among those located.
    lengths = [e["length"] for e in located["poisonExons"]]
    assert 64 in lengths, f"SCN1A's 64 nt poison exon not located: {lengths}"

    out = design_poison_exon_block("ENSG00000144285", gene_symbol="SCN1A",
                                   oligo_length=20)
    assert out["status"] == "OK"
    assert 0 < len(out["candidates"]) <= 12, "candidate count must be capped"
    exon = out["poisonExon"]
    for c in out["candidates"]:
        assert "poison exon" in c["targetElement"]
        assert (abs(c["genomicStart"] - exon["start"]) < 5000
                or abs(c["genomicStart"] - exon["end"]) < 5000)


def test_a4_targets_the_nat_overlap_not_the_whole_transcript():
    """The NAT is a different gene, and it overlaps only part of this one."""
    from services.upregulation_targets_service import (
        design_nat_knockdown, find_nat,
    )
    located = find_nat("ENSG00000144285")
    if located["status"] != "OK":
        pytest.skip("Ensembl unavailable")
    symbols = {n.get("symbol") for n in located["nats"]}
    assert "SCN1A-AS1" in symbols, f"SCN1A's antisense gene not found: {symbols}"
    for nat in located["nats"]:
        assert nat["overlapNt"] > 0

    out = design_nat_knockdown("ENSG00000144285", gene_symbol="SCN1A",
                               oligo_length=20)
    assert out["status"] == "OK"
    assert 0 < len(out["candidates"]) <= 12
    nat = out["nat"]
    for c in out["candidates"]:
        # Confined to the real overlap — the old code tiled the entire
        # transcript, so candidates outside it were complementary to nothing.
        assert nat["overlapStart"] <= c["genomicStart"] <= nat["overlapEnd"], (
            f"candidate at {c['genomicStart']} is outside the "
            f"{nat['overlapStart']}-{nat['overlapEnd']} NAT overlap")
    assert "antisense to the gene" in out["strandNote"]
