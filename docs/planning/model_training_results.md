# Model training specs — what was run, and what it found

Reports execution of `model_training_specs.md` M1–M7. One section per model:
what the spec asked for, what happened, and what is blocked.

Two models were runnable here (M2, M6). Two data defects the spec marks
"VERIFIED as broken" were confirmed against the committed data. The rest are
blocked on external data or model weights that are not in this repository,
and are listed with the specific thing that blocks them rather than skipped.

Environment: Python 3.11 was unavailable, so everything ran on 3.9 with
lightgbm 4.6.0, numpy 2.0.2, pandas 2.3.3, torch 2.8.0, ViennaRNA 2.7.0.
Seed 42 throughout.

---

## Summary

| Model | Status | Outcome |
|---|---|---|
| **M1** token cross-attention | **PARTLY UNBLOCKED** | `Hu.csv` supplied and loading (2,361 rows). Still blocked on RNA-FM weights (download 403s) and the gene re-annotation. |
| **M2** LightGBM lambdarank | **RUN** | Acceptance PASSED. Best model is `regress-rank`, not lambdarank. |
| **M3** SpliceAI calibration | **BLOCKED** | SpliceAI is reported as installed but is not present on this machine; reference set is MUST VERIFY. |
| **M4** NMD-exon classifier | **BLOCKED** | Needs a GENCODE release; MUST VERIFY. |
| **M5** uORF regressor | **BLOCKED** | Needs ribosome profiling data; MUST VERIFY. |
| **M6** Mondrian conformal | **RUN** | **Found and fixed an off-by-one that broke the coverage guarantee.** |
| **M7** OOD scorer | **BLOCKED** | Needs RNA-FM embeddings; weight download 403s. |

---

## Data defects — both confirmed

The spec marks two defects "VERIFIED as broken" and says to fix them before
retraining. Both were checked against `data/benchmark/unified_benchmark.parquet`
rather than taken on trust, and both are real.

### `target_gene` on the siRNA rows

Confirmed, and **more specific than the spec states**. The spec says the
column "contains the siRNA sequence". It actually holds the **mRNA target
site** — the reverse complement of the guide, in DNA letters. Verified for
all 3,947 of 3,947 rows: reverse-complement match 3,947, direct match 0.

The distinction matters to whoever fixes it: you are looking for the gene the
target site belongs to, not the gene the guide came from.

3,947 rows, 3,947 distinct values. A split on that column is a random row
split.

### 106 sequences in two modalities

Confirmed exactly. 106 sequences appear in both `rnase_h` and
`splice_switching` — leakage along the axis a modality comparison tests.

### Effective n

Confirmed exactly as stated: 159,215 rows, 1,941 experiments, 339 genes,
67 cell lines.

### A guard that was passing on an average

`unified_gbm_baseline.py` carried a Phase-0 guard refusing a gene split when
rows-per-gene drops below 2.0. It was computing that **globally**:

| modality | rows | genes | rows/gene |
|---|---|---|---|
| rnase_h | 159,215 | 339 | 469.66 |
| **sirna** | **3,947** | **3,947** | **1.00** |
| splice_switching | 2,287 | 6 | 381.17 |
| **global** | 165,449 | 4,292 | **38.58** |

The global mean is 38.58, so the guard passed — on the strength of the
healthy arm, while the broken one sat at exactly 1.00. Now computed per
modality, and it correctly refuses.

---

## M2 — LightGBM lambdarank (RNase-H)

**Ran.** `backend/experiments/benchmark/m2_rnase_h_ranker.py`, 32 s.
Results in `backend/results/benchmark/m2_rnase_h_ranker.json`.

All three specified fixes applied: 106 dual-modality sequences dropped,
restricted to the RNase-H arm (the only one with real gene symbols), split by
gene with effective n reported as 339 genes.

| model | top-10 | 95% CI (over genes) | pooled Pearson | per-experiment Pearson |
|---|---|---|---|---|
| random guessing (computed) | **0.174** | — | — | — |
| lambdarank-rank | 0.296 | [0.268, 0.320] | 0.278 | 0.285 |
| regress-raw | 0.298 | [0.268, 0.323] | 0.270 | 0.297 |
| **regress-rank** | **0.302** | [0.272, 0.325] | 0.309 | 0.312 |

**Acceptance: PASS.** The CI lower bound (0.272) clears random guessing
(0.174).

### The random baseline was computed, and it landed on 0.174

The spec quotes 0.174. Simulating a random ranker over the actual test
groups gives **0.174** (sd 0.004, 50 trials) — independent agreement to
three decimals.

One correction along the way: the first implementation permuted a constant
vector, which leaves every value tied and hands `nlargest` the first k rows
in index order. That is not random guessing. Drawing fresh uniform scores per
trial gives the number above.

### Lambdarank is not the best model

`regress-rank` beats it, 0.302 against 0.296. The spec's own reference table
shows the same ordering (0.299 against 0.292), so this reproduces rather than
contradicts it — but M2 is titled "LightGBM lambdarank" and the honest
reading of its own cross-cutting rule ("if a model doesn't beat its baseline,
ship the baseline") is that **plain regression on the rank label should be
the shipped model**. It is simpler and it wins.

### The metric warning, partly resolved

The spec flags that pooled and per-experiment Pearson are not comparable, and
that pooled is inflated by between-experiment variance. Both are computed
here. On this split the gap is small and runs in **both** directions:

- lambdarank: pooled 0.278, per-experiment 0.285 (pooled *lower*)
- regress-rank: pooled 0.309, per-experiment 0.312 (pooled *lower*)

So on this data the inflation the spec warns about does not appear. That does
not make the two interchangeable — they answer different questions, and
per-experiment is the one the product asks — but the choice is not currently
worth a large number either way. Reported so the decision can be made on
evidence.

Bootstrap CIs resample **unique genes**, not rows. With 339 genes behind
159k rows a row-level bootstrap would report an interval roughly an order of
magnitude too narrow.

---

## M6 — Mondrian conformal predictor

**Ran.** `backend/experiments/benchmark/m6_conformal_audit.py`, ~5 s.
Results in `backend/results/benchmark/m6_conformal_audit.json`.

The spec describes this implementation as "VERIFIED as already implemented,
with caveats" and reports empirical coverage 0.887 at α = 0.10. Verifying it
by simulation rather than accepting that found a bug.

### The off-by-one

`q_hat` was the **ceil**(α·(n_cal+1))-th smallest calibration tau. The
guarantee needs the **floor**.

The calibration taus and the test tau are n+1 exchangeable draws, so the test
tau's rank among all n+1 is uniform. Coverage holds exactly when
`tau_test >= q_hat`, which for the m-th smallest of the other n means rank
≥ m+1. So

```
P(cover) = 1 - m/(n_cal + 1)
```

and 1−α needs `m <= α·(n_cal+1)` — the **largest** such integer, the floor.
Taking the ceil overshoots by one whenever α·(n_cal+1) is not an integer,
raising the threshold and dropping coverage below nominal.

Simulation at n_cal = 30, 200 runs per α:

| α | nominal | ceil predicts | **measured before fix** | floor predicts | **measured after fix** |
|---|---|---|---|---|---|
| 0.05 | 0.95 | 0.935 | 0.941 | 0.968 | **0.973** |
| 0.10 | 0.90 | 0.871 | 0.867 | 0.903 | **0.898** |
| 0.20 | 0.80 | 0.774 | 0.770 | 0.806 | **0.807** |

Every pre-fix measurement matches the ceil prediction and sits below nominal.
Every post-fix measurement matches the floor prediction and clears it. The
spec's own 0.887 is this same failure at a different n_cal, not a passing
result.

The audit tests against the exact finite-sample expectation rather than a
hand-picked tolerance. That mattered: the first version used an invented
slack and the α = 0.10 verdict landed 0.001 from flipping — a coin toss
deciding whether a guarantee holds.

### The guarantee is unavailable for two of three classes

| class | n_groups | calibration groups | guarantee |
|---|---|---|---|
| sirna | 6 | 3 | **unavailable** |
| splice_switching | 12 | 6 | **unavailable** |
| rnase_h | 100 | 50 | valid |

At α = 0.10 a non-trivial threshold needs at least 1/α − 1 = 9 calibration
groups. `conformal_topk` now returns `guarantee: "unavailable"` with the
reason instead of a coverage number, because a number there looks like a
guarantee and is not one.

**Do not fix this by pooling.** Pooling into RNase-H would produce RNase-H
guarantees wearing another label, which is precisely what Mondrian
calibration exists to prevent. Either gather more experiments per class,
raise α and say so, or report no guarantee.

### Stored results are stale, and cannot be regenerated here

Confirmed as the spec says:

| file | stored coverage | n_groups |
|---|---|---|
| `final_gc_auto/rnase_h/pipeline_result.json` | 0.04 | 100 |
| `final_gc_auto/sirna/pipeline_result.json` | 0.167 | 6 |
| `final_gc_auto/splice_switching/pipeline_result.json` | 0.0 | 12 |
| `ranker_v2/pipeline_result.json` | 0.0 | 12 |

These predate the fix and must not be quoted. Regenerating them means
retraining the ranker that produced them — a full pipeline run, not something
this audit can do. **Left in place rather than deleted**, because deleting
them would hide that published numbers were wrong; the audit output records
that they are void.

### Also fixed

`_bootstrap_ci` ran a 10,000-iteration Python loop per call, twice per
`conformal_topk`. Vectorised — same estimator, ~100× faster. It was the
reason the audit could not complete inside two minutes.

---

## Blocked models, and precisely what blocks them

None of these is skipped for convenience. Each needs something that is not in
this repository and that the specs mark MUST VERIFY.

### M1 — token cross-attention
- `OligoFormer/data/Hu.csv` is not in the repo.
- RNA-FM weights (`backend/pretrained/RNA-FM_pretrained.pth`) are not in the
  repo; `pretrained/` is gitignored.
- The embedding cache `backend/data/hu_embeddings.pt` does not exist.
- The gene split it must be evaluated under depends on the re-annotation
  below.

`backend/data_curation/annotate_sirna_genes.py` already exists and documents
this exactly. It needs either NCBI BLAST+ with a downloaded human transcript
FASTA, or hours of rate-limited remote NCBI BLAST with an identifying email.
**Which backend, which database release, and whose email are your calls** —
the remote route submits ~4,000 queries under your identity.

Its acceptance criteria are already strict (≥98% identity, ≥50 nt alignment,
best hit must beat the second-best different gene by ≥10 bitscore, otherwise
left NA). Nothing there should be loosened to raise the annotation rate.

### M3 — SpliceAI calibration
No SpliceAI installed. The calibration reference set is MUST VERIFY and the
spec warns specifically about picking a set SpliceAI was trained on
(GTEx/GENCODE overlap). Choosing it is a research decision.

### M4 — NMD-exon classifier
Needs a specific GENCODE release, marked MUST VERIFY with an explicit
instruction not to accept a recalled version number.

### M5 — uORF repression regressor
Needs public ribosome profiling paired with RNA-seq. The spec also requires
deciding *before starting* whether the target is main-ORF translation
efficiency or the uORF/main-ORF ribosome density ratio — they answer
different questions. That decision is not made.

### M7 — OOD scorer
Needs RNA-FM embeddings of all training transcripts, so it is blocked behind
the same missing weights as M1.

---

## Cross-cutting compliance

| Requirement | Status |
|---|---|
| Split by gene, always | M2 yes. M1 cannot until re-annotation lands; the guard now refuses rather than mislabelling. |
| Bootstrap CI over unique genes | M2 yes, over 339 genes. |
| Calibration for anything entering §3.2 | Not applicable to M2 (ranking) or M6 (already a coverage method). Outstanding for M3–M5 when they run. |
| Beat the stated baseline | M2 yes (0.302 vs 0.174, computed). |
| Seed 42, versions recorded | Yes — every result file carries package versions, git commit and full hyperparameters. |

---

## Open items this raised

- **SO-ML-01** — M2's shipped model should be `regress-rank`, not
  lambdarank. It wins on the spec's own metric and is simpler.
- **SO-ML-02** — pooled vs per-experiment Pearson: pick one for publication.
  The gap is small here and does not favour either, so the choice should be
  made on which question you are answering, not on which number is larger.
- **SO-ML-03** — the conformal guarantee is unavailable for siRNA and
  splice-switching. Decide between gathering more experiments, raising α, or
  publishing "no guarantee" for those classes. Pooling is not an option.
- **SO-ML-04** — the stale `pipeline_result.json` files need a full pipeline
  rerun. Until then any coverage figure in the paper drawn from them is void.
- **SO-ML-05** — the siRNA re-annotation backend, database release and
  contact email need deciding before M1 can be evaluated honestly.

---

## Round 2 — unblocking pass

Following the ML section supplied with `organism tier and other qns.docx`.
Environment claims in that section were checked rather than assumed; two did
not hold here.

### 1. Hu.csv — LANDED

Placed at `OligoFormer/data/Hu.csv`, which is the path every experiment
script hardcodes. Loads through `HueskenDataset` at **2,361 rows**, matching
the M1 spec exactly. Columns `siRNA, mRNA, label, y, td` line up with what
the loader expects.

`test_dataset.py`: 2 pass, 2 skip. The two that skip need the RNA-FM
embedding cache, which needs weights (see item 3).

**Also fixed a false green.** Those two tests previously did
`print("Skipping"); return` on a missing cache, which pytest reports as
PASSED. A test that executes nothing and reports green is worse than one
that fails. Now `pytest.skip` with the reason and the command to fix it.

### 2. SpliceAI — NOT INSTALLED, F1/F2/F3 still on stand-ins

The doc states "spliceai is now installed". It is not present in any
interpreter or anywhere on disk on this machine:

```
python3           -> ModuleNotFoundError: No module named 'spliceai'
scratch venv      -> ModuleNotFoundError: No module named 'spliceai'
find / -name 'spliceai*'  -> nothing outside Downloads
```

So F1, F2 and F3 remain on the user-asserted stand-in rung and everything
said about that in the arbitration notes still holds. **Nothing was wired**,
because wiring an import that does not resolve would turn a documented
stand-in into a broken one.

Also note the interface named in the doc (`from spliceai import SpliceAI`,
`SpliceAI(spliceai_models_path)`) should be confirmed against the package
that actually gets installed — the pip distribution and the Illumina
repository do not expose the same API.

### 3. RNA-FM — package installs, weights do not download

`pip install rna-fm` succeeds and `from fm.pretrained import rna_fm_t12`
imports. The auto-download then fails:

```
Downloading: "https://proj.cse.cuhk.edu.hk/rnafm/api/download?filename=RNA-FM_pretrained.pth"
HTTP Error 403: Forbidden
```

Two corrections to the verification command in the doc:
- the symbol is `rna_fm_t12` (from `fm.pretrained`), not `rnafm_t12`
- it is not importable from `features.rnafm` as written

So M1, M7 and the two skipped dataset tests stay blocked on weights. Fetching
them needs either a working mirror or a manual download to
`backend/pretrained/RNA-FM_pretrained.pth`, which `RNAFMEmbedder` already
prefers when present.

### 4. test_main.py — one stale test fixed, two specify fabricated data

The doc describes all three as "stale refs". One is. Two are not.

**Fixed (genuinely stale).** `monkeypatch.setattr(main_module,
"get_top_dbsnp_id", ...)` in two tests. `main.py` has never defined or called
that function — no commit in the history touches it — so the patch raised
before the test reached anything real. Removed, along with the unreachable
stub. `test_initialize_target_preserves_aso_metrics_when_analysis_is_slow`
now passes for real.

**Left failing, deliberately, as `xfail(strict=True)`:**

`test_get_aso_analysis_uses_gene_metadata_without_network` asks for a
`gene_metadata` fast path that this codebase has never had. Two of its four
assertions cannot be met honestly:

- `structuralAccessibility` is computed from CDS **sequence** GC content in
  `_compute_aso_metrics_from_sequence`. `gene_metadata` carries counts and
  coordinates, no sequence. With the network forbidden there is nothing to
  compute it from.
- `spliceSwitches == 94` implies `totalTranscripts - 1`, while the network
  path defines it as `(distinct exon counts) - 1`. Satisfying the test gives
  one output field two different meanings depending on which path ran.

`test_initialize_target_uses_generic_aso_fallback_for_any_gene_when_analysis_times_out`
asks the endpoint to emit `activeIsoforms=1`, `spliceSwitches=0` and a truthy
`structuralAccessibility` **when the analysis times out**. Those are invented
values delivered in the same fields as measured ones, with nothing marking
them as a fallback — the same defect class as the `hash()`-derived binding
affinities and the always-True feature placeholders removed earlier in this
review. On timeout `_run_sync` correctly returns `{}` and the fields stay
None.

Making either green requires inventing numbers, so both are marked with the
reasoning inline. **This is a decision to confirm, not a fix to apply.**

### 5. Full suite

```
48 passed, 2 skipped, 2 xfailed
```

| file | result |
|---|---|
| test_benchmark.py | 5 passed |
| test_dataset.py | 2 passed, 2 skipped (no RNA-FM weights) |
| test_features.py | 1 passed |
| test_gene_service.py | 2 passed |
| test_gene_silencing_service.py | 2 passed |
| test_main.py | 1 passed, 2 xfailed (see above) |
| test_mechanism_arbitration.py | 30 passed |
| test_token_model.py | 5 passed |

Before this pass the suite could not be collected at all in a clean
environment: `fastapi`, `sqlalchemy`, `aiohttp`, `torch`, `rna-fm`,
`lightgbm`, `pandas` and `ViennaRNA` were all missing. Installing them is
what makes the numbers above reproducible; `backend/requirements.txt` lists
them but nothing pins versions.

### SO-ML decisions — recorded, with one discrepancy

| item | decision | status |
|---|---|---|
| SO-ML-01 | Ship regress-rank as M2 | Matches the measurement (top-10 0.302). Applied. |
| SO-ML-02 | "Use pooled Pearson (0.362) for publication" | **Number does not match any measurement here.** See below. |
| SO-ML-03 | Publish "no guarantee" for siRNA and splice-switching | Matches; `conformal_topk` already returns that, and the calibration-group counts (3 and 6) are as stated. |
| SO-ML-04 | pipeline_result.json stale, regenerate after F1-F3 | Agreed; blocked on SpliceAI. Old coverage numbers remain void. |
| SO-ML-05 | siRNA re-annotation skipped, not on critical path | Applied. The gene-split guard still refuses rather than mislabelling. |

**SO-ML-02 needs checking before it is published.** The measured pooled
Pearson values on the gene split are:

| model | pooled | per-experiment |
|---|---|---|
| lambdarank-rank | 0.278 | 0.285 |
| regress-raw | 0.270 | 0.297 |
| **regress-rank** | **0.309** | 0.312 |

There is no 0.362 in this run. It may come from a different split, a
different subset, or the pre-fix pipeline. Until its origin is identified it
should not go into a paper — every number in this project is supposed to
trace to a run. `m2_rnase_h_ranker.json` carries the full provenance for the
figures above.
