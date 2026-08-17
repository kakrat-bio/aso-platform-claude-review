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
