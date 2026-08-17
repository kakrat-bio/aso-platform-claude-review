"""Last-known-good cache for real external data. Never a synthesiser.

THE RULE THIS ENFORCES
----------------------
When an external source is unreachable, this platform serves **real data it
fetched earlier** or it says the data is unavailable. It does not fabricate a
substitute.

The distinction matters because a fabricated substitute is indistinguishable
from a measurement once it reaches the UI. `protein_replacement_service` built
5' UTRs as `"GCCACC" + "A" * n` and 3' UTRs as `"G" * n`, then computed GC
content, folding MFE, U-content and a translation-risk estimate from that
padding. Every one of those numbers described the padding rather than the
gene, and nothing in the response said so.

THREE OUTCOMES, ALWAYS DISTINGUISHABLE
--------------------------------------
    live        fetched from the source just now
    cached      a real earlier fetch, replayed, with its age
    unavailable no source and no cache entry — say so, return nothing

There is deliberately no fourth outcome. A caller that wants to render
something regardless must decide what to show for `unavailable`; it cannot be
handed a plausible-looking sequence.

SEEDING
-------
The cache fills itself from successful live fetches, so it warms up in normal
use. It can also be pre-seeded from a curated TSV of verified sequences for
common therapeutic targets — see `data/reference/curated_transcript_parts.tsv`
and the README beside it. Curated rows carry their own provenance and are
never overwritten by a live fetch.
"""

from __future__ import annotations

import json
import logging
import time

from database.db import SessionLocal
from database.models import RealDataCache

from . import reference_tables as RT

logger = logging.getLogger(__name__)

LIVE = "live"
CACHED = "cached"
CURATED = "curated"
UNAVAILABLE = "unavailable"

# How old a cached entry may be before the caller is warned. It is still
# served — real but stale beats invented and fresh — just flagged.
STALE_AFTER_SECONDS = 30 * 24 * 3600


def _now() -> float:
    return time.time()


def store(namespace: str, key: str, payload: dict,
          source: str = "", source_version: str = "") -> None:
    """Persist a successful live fetch as the new last-known-good.

    Curated rows are never overwritten: a hand-verified sequence outranks
    whatever an API returned today.
    """
    try:
        db = SessionLocal()
        try:
            row = (
                db.query(RealDataCache)
                .filter(RealDataCache.namespace == namespace,
                        RealDataCache.cache_key == key)
                .one_or_none()
            )
            if row is not None and row.origin == CURATED:
                return
            now = _now()
            if row is None:
                db.add(RealDataCache(
                    namespace=namespace, cache_key=key,
                    payload=json.dumps(payload), origin=LIVE,
                    source=source, source_version=source_version,
                    created_at=now, updated_at=now,
                ))
            else:
                row.payload = json.dumps(payload)
                row.origin = LIVE
                row.source = source
                row.source_version = source_version
                row.updated_at = now
            db.commit()
        finally:
            db.close()
    except Exception as exc:  # noqa: BLE001 — caching must never break a request
        logger.warning("real_data_cache store failed (%s/%s): %s",
                       namespace, key, exc)


def load(namespace: str, key: str) -> dict | None:
    """Replay the last-known-good entry, or None."""
    try:
        db = SessionLocal()
        try:
            row = (
                db.query(RealDataCache)
                .filter(RealDataCache.namespace == namespace,
                        RealDataCache.cache_key == key)
                .one_or_none()
            )
            if row is None:
                return None
            age = _now() - (row.updated_at or 0)
            return {
                "payload": json.loads(row.payload or "{}"),
                "origin": row.origin or CACHED,
                "source": row.source or None,
                "sourceVersion": row.source_version or None,
                "fetchedAt": row.updated_at,
                "ageSeconds": age,
                "stale": age > STALE_AFTER_SECONDS,
            }
        finally:
            db.close()
    except Exception as exc:  # noqa: BLE001
        logger.warning("real_data_cache load failed (%s/%s): %s",
                       namespace, key, exc)
        return None


def resolve(namespace: str, key: str, fetch, source: str = "",
            source_version: str = "") -> dict:
    """Live, else cached, else explicitly unavailable.

    `fetch` is a zero-argument callable returning the real payload, or None /
    raising when the source cannot answer. Its result is cached on success.

    The return is always the same shape, and `status` always says which of
    the three happened, so no caller can mistake a replay for a live read or
    an absence for a value.
    """
    seeded = _from_curated(namespace, key)
    try:
        payload = fetch()
    except Exception as exc:  # noqa: BLE001 — an outage is not a crash
        logger.info("%s/%s live fetch failed: %s", namespace, key, exc)
        payload = None

    if payload:
        store(namespace, key, payload, source, source_version)
        return {
            "status": LIVE, "data": payload, "source": source or None,
            "sourceVersion": source_version or None,
            "fetchedAt": _now(), "ageSeconds": 0.0, "stale": False,
            "note": None,
        }

    cached = load(namespace, key) or seeded
    if cached:
        origin = cached["origin"]
        age_days = (cached.get("ageSeconds") or 0) / 86400
        return {
            "status": CURATED if origin == CURATED else CACHED,
            "data": cached["payload"],
            "source": cached.get("source"),
            "sourceVersion": cached.get("sourceVersion"),
            "fetchedAt": cached.get("fetchedAt"),
            "ageSeconds": cached.get("ageSeconds"),
            "stale": bool(cached.get("stale")),
            "note": (
                "Curated verified entry — the live source was not reachable."
                if origin == CURATED else
                f"Replayed from a real earlier fetch ({age_days:.0f} days old)"
                f"{' — STALE' if cached.get('stale') else ''}. The live source "
                f"was not reachable."
            ),
        }

    return {
        "status": UNAVAILABLE, "data": None, "source": None,
        "sourceVersion": None, "fetchedAt": None, "ageSeconds": None,
        "stale": False,
        "note": (
            "The live source was not reachable and nothing real has been "
            "cached or curated for this target. No substitute is generated: a "
            "fabricated sequence would be indistinguishable from a measured "
            "one once rendered."
        ),
    }


def _from_curated(namespace: str, key: str) -> dict | None:
    """Pre-seeded verified rows, from the curated reference table."""
    if namespace != "transcript_parts":
        return None
    row = RT.row_for("curated_transcript_parts", key)
    if not row:
        return None
    payload = {
        k: (row.get(k) or "")
        for k in ("utr5", "cds", "utr3", "transcript_id", "assembly")
        if row.get(k)
    }
    if not payload:
        return None
    return {
        "payload": payload,
        "origin": CURATED,
        "source": RT.provenance_of(row) or "curated transcript parts",
        "sourceVersion": row.get("source_version"),
        "fetchedAt": None,
        "ageSeconds": None,
        "stale": False,
    }


def status() -> dict:
    """Counts per namespace, for diagnostics and the /scope endpoint."""
    out: dict = {"available": True, "namespaces": {}}
    try:
        db = SessionLocal()
        try:
            for (ns, origin, count) in (
                db.query(RealDataCache.namespace, RealDataCache.origin,
                         __import__("sqlalchemy").func.count())
                .group_by(RealDataCache.namespace, RealDataCache.origin)
                .all()
            ):
                out["namespaces"].setdefault(ns, {})[origin or CACHED] = int(count)
        finally:
            db.close()
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "reason": str(exc), "namespaces": {}}
    return out
