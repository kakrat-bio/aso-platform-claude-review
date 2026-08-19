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
