# ASO therapeutic-design platform — working agreement

FastAPI + SQLAlchemy backend (`backend/`), Next.js 15 frontend (`frontend/`).
Nine therapeutic goals (TG01–TG09), 38 mechanisms (A1–A39; **A22 does not
exist** — a known gap, not an oversight to fix by renumbering).

---

## 1. Never fabricate a value

This is a therapeutic-design tool. A plausible-looking invented number is
worse than a missing one, because downstream nothing distinguishes it from a
measurement.

**Resolution order for every external fact: live → cached → curated →
`UNAVAILABLE`.** There is no synthesis branch. `services/real_data_cache.py`
implements this; use it.

Concretely, this means:

- No placeholder sequences. Not `"AUG" + "GCU" * 300`, not `"GCCACC" + "A" * 76`.
  If the transcript fetch fails, return `status: "UNAVAILABLE"` with the reason.
- No constants dressed as measurements. A score that returns the same value
  for every input is not a score.
- No derived quantity without the thing it derives from. If a claim needs a
  fitted coefficient or a calibration set that does not exist, do not ship the
  claim — emit `null` plus a `notComputed` entry saying what is missing.
- No counts from proxies. "7 off-target sites" from k-mer repetitiveness is a
  fabrication; nothing in this repo aligns against a genome.
- Numbers that cannot be computed are reported as blocked **with the specific
  blocker**, never filled with a plausible value. Same rule applies to
  experiments in `backend/experiments/`.

Three separate quantities — applicability, confidence, evidence — are never
blended into one number.

## 2. "Done" means you ran it

Editing files is not completion. Before reporting any mechanism, goal or
endpoint as working, **issue the request and paste the output.**

Nearly every serious defect found in this codebase survived because someone
(including Claude) marked work complete from the diff rather than from an
executed probe. TG07's designer raised `NameError` on every call and had never
returned a candidate; it was recorded as done. Each of these took under two
minutes to find once actually invoked.

`uvicorn` is **not installed**. Probe in-process:

```python
import sys, logging; sys.path.insert(0, 'backend'); logging.disable(logging.WARNING)
from fastapi.testclient import TestClient
from main import app
c = TestClient(app)
print(c.post('/api/mechanisms/gene-silencing', json={...}).json())
```

A mechanism is working when it (a) scores, (b) returns candidates with real
sequences, and (c) those sequences are derived from the fetched transcript —
check one by hand.

## 3. Commands

```bash
python3 -m pytest backend/tests -q          # from the REPO ROOT, not backend/
```

Several tests import `backend.*`; running from `backend/` fails collection.
Use `-k <pattern>` while working, full suite before committing (~70 s).

```bash
cd frontend && npx --no-install tsc --noEmit
```

Python 3.9. Installed: fastapi, ViennaRNA (`import RNA`), primer3, lightgbm,
torch, spliceai, tensorflow<2.17, multimolecule.

`backend/database/notifications.json` picks up rows whenever an endpoint is
exercised. **Revert it before committing** — it is runtime state, not a change.

## 4. Architecture

`POST /api/mechanisms/arbitrate` is the primary path: it scores **every**
mechanism in one pass, and the therapeutic goal is an *output label*, not an
input. This exists to fix the nusinersen case — intent is upregulation (TG02)
while the mechanism is splice modulation (TG04), so goal-first routing hides
the right answer. Per-goal endpoints are filtered views over that one ranking.

Mechanism states: SCORED+DESIGNABLE / SCORED+DESIGN-UNAVAILABLE / HALTED /
FLAGGED. A **required** feature that is UNRESOLVED halts; a **discriminating**
one is skipped. An absent feature returns ABSENT, never probability zero.

Provenance caps confidence: `measured 1.00 > confirmed 0.95 > annotation 0.90
> predicted 0.75 > user_asserted 0.60`. A computed feature (e.g. a ViennaRNA
fold) is PREDICTED — expect confidence to drop when one resolves.

Designers, by mechanism:

| mechanisms | service |
| --- | --- |
| A1, A2, A12, A15 | `gene_silencing_service` |
| A3–A6, A23, A28 | `gene_upregulation_service` |
| A7–A11 | `rna_processing_service` |
| A13, A16, A17, A20 | `rna_editing_service` |
| A14, A12 | `rna_neutralization_service` |
| A2, A5, A6, A27, A29–A31 | `translational_regulation_service` |
| A21 | `sirna_duplex_service` |
| A18, A19 | `programmable_editor_service` |
| A32, A33 | `transcript_architecture_service` |
| A24, A26 | `protein_replacement_service` via `arbitration.designRoute` |

31 of 38 are designable. The other 7 (A25, A34–A39) each carry a specific
`designUnavailableReason`. A test fails if any mechanism produces no
candidates **and** gives no reason.

## 5. Traps that have already cost time

- **`arbitration.flagReason` is the field the API serves**, not a top-level
  `flagReason`. Writing the top-level one creates a duplicate nothing reads.
- **Ensembl has no `type=utr5` / `type=utr3`.** Those return HTTP 400. Valid
  types: genomic, cds, cdna, protein. Get UTRs from
  `gene_silencing_service.get_target_analysis`, which aligns cDNA against CDS.
- **`N` in a regex is a literal N**, not a degenerate base. `/^NGG$/` matched
  nothing for years. Use `[ACGT]GG`.
- **primer3 is DNA-only** and raises on `U`. Fold U→T before calling it —
  the value being computed is the DNA-analogue Tm; chemistry adjustment is a
  separate term.
- **Reverse complement must map U→A.** A table covering only `ATGC` lets U
  pass through uncomplemented, producing an oligo that does not bind its
  target at every U position.
- **`np.argsort` is stable**, so ties break by row order — which in the
  benchmark is tiling order along the transcript. Shuffle each group before
  sorting or the ranking silently encodes position.
- **Gene symbols are ambiguous in free text.** `HTT` matches "hyalinizing
  trabecular tumours"; UniProt returns SLC6A4 for `gene_exact:HTT`. Use
  record-linked routes (NCBI Gene→PubMed elink) and verify the symbol matches.
- **`.gitignore` trailing comments are part of the pattern.** Put comments on
  their own line.
- **Score normalisation against absolute constants saturates.** Every
  candidate in one run shares a length, so a range tuned for 12–30 nt clips
  entirely. Rank within the pool.

## 6. Ranking

`gene_silencing_service._composite_score` is largely a GC proxy: it correlates
0.67 with GC and 0.033 with measured activity. Experiment E12 in
`backend/experiments/ml_analysis` measured ΔG-first ordering against 528
held-out experiments and found it selects **below chance** (top-5 lands at the
46.8th activity percentile; random gets 51.3).

So: rank on site accessibility where the mechanism justifies it (an RNase-H
heteroduplex cannot form inside a hairpin), keep ΔG as a tie-break, and attach
the caveat. Do not present a thermodynamic ordering as a validated activity
model.

`backend/experiments/ml_analysis/FINDINGS.md` holds the measured results —
consult it before adding a ranking term.

## 7. Data reality

`backend/data/benchmark/unified_benchmark.parquet`, 165,237 rows:

| modality | rows | genes | rows/gene | split key |
| --- | ---: | ---: | ---: | --- |
| `rnase_h` | 159,109 | 339 | 469 | `target_gene` |
| `splice_switching` | 2,181 | 5 | 436 | `target_gene` (held-out set = 1 gene) |
| `sirna` | 3,947 | 3,947 | 1.00 | `experiment_id` — no gene split exists |

siRNA's `target_gene` holds the target *site*, one per row. A "gene split"
there is a random row split wearing a label; `common.gene_split()` raises
rather than run it. **Effective n is genes, not rows** — bootstrap over genes.

Any siRNA or splice-switching result rests on 2–4 held-out groups. Label those
UNDERPOWERED and say a null result means "no effect large enough to see".

ML models: M2 and M6 have run. M1, M3, M4, M5, M7 are blocked on external data
or unmade research decisions — see `docs/planning/model_training_results.md`.
SpliceAI is installed and wired but **not calibrated** (`calibrated: False`).

## 8. Working style

- One mechanism (or pair) per session. State the acceptance check up front.
- Write the probe before the fix, so completion is mechanical.
- Add a regression test for every defect found — those tests are the durable
  memory across sessions, and cost nothing per session.
- Commit messages: what was broken, the evidence it was broken, what changed.
  Include the measured numbers.
- If a concern about the request is raised and the user reaffirms it, that is
  their decision — proceed with the full request.
