"""Run the whole suite, in dependency order, and write a summary index.

Several experiments share trained models through `ranker.train_or_load`'s
cache, so order matters for wall time (not for correctness): E1 trains the
pooled model that E2, E4, E7, E8 and E12 reuse.

Each experiment is isolated — a failure is recorded with its traceback and
the run continues, because a suite that stops at the first missing
dependency tells you less than one that reports eleven results and one
blocker.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402

ORDER = [
    ("exp01_class_imbalance", "epochs"),
    ("exp02_per_mechanism", "epochs"),
    ("exp04_ranking_metrics", "epochs"),
    ("exp06_significance", "epochs"),
    ("exp07_error_analysis", "epochs"),
    ("exp08_attribution", "epochs"),
    ("exp12_ml_vs_heuristic", "epochs"),
    ("exp03_within_mechanism_split", "epochs"),
    ("exp05_cross_chemistry", "epochs"),
    ("exp09_motif_analysis", None),
    ("exp10_thermodynamics", None),
    ("exp11_rnafm_mlp", "epochs_conv"),
]


def main(epochs: int, only: list[str] | None = None) -> None:
    index = {}
    for name, kwarg in ORDER:
        if only and not any(o in name for o in only):
            continue
        mod = importlib.import_module(
            f"backend.experiments.ml_analysis.{name}")
        t0 = time.time()
        print(f"\n{'=' * 70}\n{name}\n{'=' * 70}", flush=True)
        try:
            mod.main(**({kwarg: epochs} if kwarg else {}))
            index[name] = {"status": "ok",
                           "seconds": round(time.time() - t0, 1)}
        except Exception:
            tb = traceback.format_exc()
            print(tb, flush=True)
            index[name] = {"status": "failed",
                           "seconds": round(time.time() - t0, 1),
                           "traceback": tb.strip().splitlines()[-3:]}
    C.write_result("_index", {"epochs": epochs, "runs": index})
    print("\n" + json.dumps(index, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--only", nargs="*", default=None)
    main(**vars(ap.parse_args()))
