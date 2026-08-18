"""Generate the design sets that E9, E10 and E12 analyse.

Loads the trained mechanism-conditioned CVAE from
`backend/results/benchmark/generative_v3/generator.pt` and samples, for each
mechanism, a design set plus two controls:

    generated   CVAE samples conditioned on the mechanism
    shuffled    each generated sequence dinucleotide-shuffled. This is the
                right control for motif analysis: it preserves length,
                composition AND dinucleotide frequency, so any motif that
                survives comparison to it is not an artefact of base
                composition.
    random      uniform i.i.d. nucleotides at matched lengths. The weaker
                control, kept because it is what the existing benchmark
                compares against.

Written once to `results/generated_sequences.csv` and reused, so E9, E10
and E12 all analyse the same designs rather than three separate draws.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from backend.experiments.ml_analysis import common as C  # noqa: E402
from backend.experiments.benchmark.generative_design import (  # noqa: E402
    CHEM_OOV, NUCLEOTIDES, gc_mean, generate, load_model,
)

CKPT = C.BACKEND / "results" / "benchmark" / "generative_v3" / "generator.pt"
OUT = C.OUT_DIR / "generated_sequences.csv"
N_PER_MECHANISM = 500


def dinucleotide_shuffle(seq: str, rng: np.random.Generator,
                         tries: int = 20) -> str:
    """Altschul-Erickson style shuffle: preserve dinucleotide counts.

    Implemented as a random Eulerian walk on the dinucleotide graph. A walk
    can strand vertices, so it retries; on repeated failure it returns a
    plain mononucleotide shuffle and the caller can see the difference in
    the recorded dinucleotide distance.
    """
    if len(seq) < 4:
        return seq
    for _ in range(tries):
        edges: dict[str, list[str]] = {}
        for a, b in zip(seq, seq[1:]):
            edges.setdefault(a, []).append(b)
        for k in edges:
            rng.shuffle(edges[k])
        out, cur = [seq[0]], seq[0]
        ok = True
        for _ in range(len(seq) - 1):
            if not edges.get(cur):
                ok = False
                break
            cur = edges[cur].pop()
            out.append(cur)
        if ok and len(out) == len(seq):
            return "".join(out)
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def main(n: int = N_PER_MECHANISM) -> Path:
    if not CKPT.exists():
        raise FileNotFoundError(
            f"generator checkpoint not found at {CKPT}. E9/E10/E12 analyse "
            f"generated designs and cannot run without it; train one with "
            f"`python -m backend.experiments.benchmark.generative_design "
            f"--mode train`."
        )
    df = C.load_benchmark()
    model, mechs, chems = load_model(CKPT)
    rng = np.random.default_rng(C.SEED)

    rows = []
    for mech in sorted(df["modality"].unique()):
        target_gc = gc_mean(df, mech)
        seqs = generate(df, model, mechs, chems, mech, CHEM_OOV, n,
                        gc_target=target_gc)
        for s in seqs:
            rows.append({"mechanism": mech, "kind": "generated", "seq": s,
                         "gc_target": target_gc})
            rows.append({"mechanism": mech, "kind": "shuffled",
                         "seq": dinucleotide_shuffle(s, rng),
                         "gc_target": target_gc})
            rows.append({"mechanism": mech, "kind": "random",
                         "seq": "".join(rng.choice(list(NUCLEOTIDES), len(s))),
                         "gc_target": target_gc})
        print(f"[generate] {mech}: {len(seqs)} designs "
              f"(gc_target={target_gc:.3f})", flush=True)

    out = pd.DataFrame(rows)
    C.OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(out)} rows)")
    return OUT


def load_generated() -> pd.DataFrame:
    if not OUT.exists():
        main()
    return pd.read_csv(OUT)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else N_PER_MECHANISM)
