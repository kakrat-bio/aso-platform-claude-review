"""SpliceAI predictions for the splice-related features F1, F2 and F3.

WHAT THIS IS AND IS NOT
-----------------------
SpliceAI ships two things: a command-line variant annotator, and the five
trained networks underneath it. The annotator (`spliceai.utils.Annotator` +
`get_delta_scores`) needs a VCF record and a whole-genome reference FASTA —
several gigabytes, and genomic coordinates this platform does not carry.

The networks themselves need neither. Each takes a one-hot sequence
`(batch, length, 4)` and returns `(batch, length - 2*CONTEXT, 3)`: per
position, the probability of being [neither, acceptor, donor]. That is
exactly the quantity F1/F2/F3 are defined in terms of, so this module drives
the networks directly on a supplied pre-mRNA sequence.

The consequence, stated plainly: these features resolve **only when a
pre-mRNA sequence is supplied**. A mature mRNA is not enough — splice-site
recognition is a statement about exon/intron boundaries, and a spliced
transcript has none left. Without one, F1/F2/F3 fall through to the
documented user-asserted stand-in exactly as before.

ENSEMBLE
--------
SpliceAI averages five independently trained models. Using one is not
SpliceAI, so all five are loaded and averaged. First call pays the load
cost; the ensemble is then cached for the process.

CALIBRATION
-----------
These are raw network outputs. The field routinely thresholds them at
0.2 / 0.5 / 0.8 as though they were probabilities, and `model_training_specs.md`
M3 exists precisely because they are not calibrated. Everything emitted here
therefore carries the PREDICTED provenance tier and its confidence cap, and
M3 remains open. Do not read a delta of 0.4 as "40% chance".
"""

from __future__ import annotations

import logging
import threading

logger = logging.getLogger(__name__)

# SpliceAI's trained context: 5,000 nt either side of every scored position.
# Sequence shorter than 2*CONTEXT + 1 cannot be scored at all, and a position
# nearer than CONTEXT to either end is scored with padding rather than real
# sequence, which the network was not trained on.
CONTEXT = 5000
N_MODELS = 5

_ALPHABET = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}

_models = None
_load_lock = threading.Lock()
_load_error: str | None = None


def _probe_import() -> tuple[bool, str | None]:
    """Import the dependencies once, at module load, on the main thread.

    `spliceai/__init__.py` installs a SIGINT handler at import time, and
    `signal.signal` raises ValueError anywhere but the main thread. FastAPI
    runs sync endpoints in a worker threadpool, so importing it lazily inside
    a request handler raises and 500s the endpoint. Importing here — while
    the module graph is still being built on the main thread — installs the
    handler once and makes every later call a cached lookup.

    Catches Exception rather than ImportError: a broken native dependency
    raises all sorts of things, and none of them should take down a request.
    """
    try:
        import spliceai  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        if isinstance(exc, ImportError):
            return False, "the spliceai package is not installed"
        return False, f"spliceai failed to import: {type(exc).__name__}: {exc}"
    try:
        import keras  # noqa: F401
    except Exception:  # noqa: BLE001
        return False, "keras is not installed"
    try:
        import tensorflow  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        return False, (
            f"no usable tensorflow backend — keras cannot load the bundled "
            f".h5 models ({type(exc).__name__})"
        )
    return True, None


_AVAILABLE, _UNAVAILABLE_REASON = _probe_import()


def available() -> tuple[bool, str | None]:
    """Can SpliceAI actually run here? Returns (ok, reason_if_not)."""
    return _AVAILABLE, _UNAVAILABLE_REASON


def _load_models():
    """Load and cache the five-model ensemble.

    A failure is cached too. Retrying a broken TensorFlow install on every
    feature resolution would turn one bad environment into a slow one.
    """
    global _models, _load_error
    if _models is not None or _load_error is not None:
        return _models
    with _load_lock:
        if _models is not None or _load_error is not None:
            return _models
        ok, why = available()
        if not ok:
            _load_error = why
            return None
        try:
            from keras.models import load_model
            from pkg_resources import resource_filename

            loaded = []
            for i in range(1, N_MODELS + 1):
                path = resource_filename("spliceai", f"models/spliceai{i}.h5")
                loaded.append(load_model(path, compile=False))
            _models = loaded
            logger.info("SpliceAI ensemble loaded (%d models)", len(loaded))
        except Exception as exc:  # noqa: BLE001 — cache any load failure
            _load_error = f"{type(exc).__name__}: {exc}"
            logger.warning("SpliceAI models failed to load: %s", _load_error)
    return _models


def load_error() -> str | None:
    return _load_error


def _one_hot(seq: str):
    import numpy as np

    x = np.zeros((len(seq), 4), dtype=np.float32)
    for i, base in enumerate(seq.upper()):
        idx = _ALPHABET.get(base)
        if idx is not None:
            x[i, idx] = 1.0
    # Anything not A/C/G/T/U stays all-zero, which is how SpliceAI encodes N.
    return x


def predict(sequence: str):
    """Per-position [neither, acceptor, donor] probabilities.

    The sequence is padded with CONTEXT of N on both sides, so the returned
    array aligns 1:1 with the input. Positions within CONTEXT of an end are
    scored against padding and should be treated as less reliable.

    Returns None when SpliceAI cannot run.
    """
    import numpy as np

    models = _load_models()
    if not models or not sequence:
        return None

    padded = ("N" * CONTEXT) + sequence + ("N" * CONTEXT)
    x = _one_hot(padded)[np.newaxis, :, :]

    total = None
    for model in models:
        y = model.predict(x, verbose=0)
        total = y if total is None else total + y
    return (total / len(models))[0]


def splice_site_scores(sequence: str) -> dict | None:
    """Summarise how strongly a sequence is recognised as spliced.

    Returns max and mean acceptor/donor probability, and the positions where
    the maxima sit — enough for F1 to say whether a boundary is recognised
    weakly or strongly.
    """
    import numpy as np

    probs = predict(sequence)
    if probs is None:
        return None
    acceptor, donor = probs[:, 1], probs[:, 2]
    return {
        "maxAcceptor": float(acceptor.max()),
        "maxDonor": float(donor.max()),
        "meanAcceptor": float(acceptor.mean()),
        "meanDonor": float(donor.mean()),
        "argmaxAcceptor": int(acceptor.argmax()),
        "argmaxDonor": int(donor.argmax()),
        "length": len(sequence),
    }


def delta_scores(ref_sequence: str, alt_sequence: str,
                 window: int | None = None) -> dict | None:
    """SpliceAI delta scores: what the variant does to splicing.

    The four published quantities — acceptor gain, acceptor loss, donor gain,
    donor loss — each the largest change in the corresponding probability
    between the reference and alternate sequences.

    `window` restricts the comparison to that many bases either side of the
    point where the two sequences first differ, matching SpliceAI's own
    distance cutoff. None compares the full length.

    Both sequences must be the same length: an indel changes the coordinate
    frame, and comparing shifted positions would report a change that is an
    artefact of the shift rather than of splicing.
    """
    import numpy as np

    if not ref_sequence or not alt_sequence:
        return None
    if len(ref_sequence) != len(alt_sequence):
        return None

    ref = predict(ref_sequence)
    alt = predict(alt_sequence)
    if ref is None or alt is None:
        return None

    if window is not None:
        diff = np.nonzero(
            np.frombuffer(ref_sequence.encode(), dtype=np.uint8)
            != np.frombuffer(alt_sequence.encode(), dtype=np.uint8)
        )[0]
        if len(diff):
            centre = int(diff[0])
            lo = max(0, centre - window)
            hi = min(len(ref_sequence), centre + window + 1)
            ref, alt = ref[lo:hi], alt[lo:hi]

    d_acceptor = alt[:, 1] - ref[:, 1]
    d_donor = alt[:, 2] - ref[:, 2]
    return {
        "acceptorGain": float(max(0.0, d_acceptor.max())),
        "acceptorLoss": float(max(0.0, (-d_acceptor).max())),
        "donorGain": float(max(0.0, d_donor.max())),
        "donorLoss": float(max(0.0, (-d_donor).max())),
        "maxDelta": float(
            max(abs(d_acceptor).max(), abs(d_donor).max())
        ),
    }


def status() -> dict:
    """For the /scope endpoint: is SpliceAI usable, and if not, why."""
    ok, why = available()
    return {
        "installed": ok,
        "reason": why,
        "loaded": _models is not None,
        "loadError": _load_error,
        "contextNt": CONTEXT,
        "nModels": N_MODELS,
        "calibrated": False,
        "calibrationNote": (
            "Raw network outputs, not calibrated probabilities. Features "
            "derived from them carry the PREDICTED provenance tier. See "
            "model_training_specs.md M3."
        ),
    }
