# ML analysis suite

Twelve experiments on `backend/data/benchmark/unified_benchmark.parquet`.
Results land in `results/<experiment>.json`; the narrative is in
[`FINDINGS.md`](FINDINGS.md).

```bash
python -m backend.experiments.ml_analysis.run_all --epochs 8
python -m backend.experiments.ml_analysis.run_all --only exp05 exp09
```

Run from the repository root. Trained models are cached under
`results/models/` on a hash of (config, training rows, chemistry vocabulary),
so a rerun is cheap and two experiments that say "the pooled model" mean the
same weights.

## The experiments

| file | question |
| --- | --- |
| `exp01_class_imbalance.py` | Does upweighting siRNA / splice-switching pairs change held-out ranking? |
| `exp02_per_mechanism.py` | Is the pooled ranker an RNase-H ranker? Do the scarce arms do better alone? |
| `exp03_within_mechanism_split.py` | Within one mechanism, how much of the score is memorising the target gene? |
| `exp04_ranking_metrics.py` | NDCG / MAP / MRR / top-k at k ∈ {1,3,5,10,20}, floor and ceiling included |
| `exp05_cross_chemistry.py` | Train on chemistry X, test on Y — the transfer matrix |
| `exp06_significance.py` | Is conditioned significantly better than seq-only or invariant? |
| `exp07_error_analysis.py` | Hardest genes, mechanisms, chemistries; what the misses look like |
| `exp08_attribution.py` | Which nucleotide positions drive the score (integrated gradients + occlusion) |
| `exp09_motif_analysis.py` | Gapmer architecture, Reynolds rules, recurrent k-mers |
| `exp10_thermodynamics.py` | Tm / ΔG / MFE of generated designs vs the training distribution |
| `exp11_rnafm_mlp.py` | An MLP on frozen RNA-FM embeddings |
| `exp12_ml_vs_heuristic.py` | The learned ranker against the platform's own composite heuristic |

Support modules: `common.py` (loading, splits, metrics), `ranker.py`
(shared trainer with pair weighting and a model cache), `rnafm.py` (RNA-FM
embeddings), `generate_sequences.py` (the design sets E9/E10/E12 analyse).

## Conventions that apply to every number here

**Effective n is genes, not rows.** 159,109 RNase-H rows come from 339
genes. Every confidence interval resamples genes; a row-level bootstrap on
this data would be roughly 20× too narrow.

**Metrics are computed within an experiment group, then averaged over
groups.** Ranking two oligos measured in different assays is not a question
the product asks.

**The cross-modality overlap is dropped by default.** 106 sequences appear
in both `rnase_h` and `splice_switching`; `load_benchmark()` removes them, so
any cross-modality comparison is not measuring that overlap.

**Each modality is split on the strongest key it actually supports.**

| modality | rows | genes | rows/gene | split key |
| --- | ---: | ---: | ---: | --- |
| `rnase_h` | 159,109 | 339 | 469 | `target_gene` |
| `splice_switching` | 2,181 | 5 | 436 | `target_gene` (held-out set is **one gene**) |
| `sirna` | 3,947 | 3,947 | 1.00 | `experiment_id` — see below |

The siRNA arm's `target_gene` column holds the target *site*, one value per
row, so a "gene split" there is a random row split wearing a label.
`common.gene_split()` raises rather than run it; E3 reports that refusal as
its result for that arm.

**Nothing is filled in when it cannot be measured.** An experiment that
cannot run records its blocker. Cells with too few held-out groups record
`unavailable` and the count, not a number.

## Reading the metrics

`rank_label` is a within-experiment activity percentile in [0, 100] — dense
and positive. NDCG's gain is therefore never zero and NDCG cannot fall near
zero even for a bad scorer; E4 carries a drawn-random floor and a label
oracle ceiling in the same table so every NDCG has something to be read
against. Top-k overlap and MRR are the more honest headline for "did we
surface the good ones".

The random floor is **drawn**, not shuffled. Permuting a constant score
vector leaves every tie intact and scores far too well.
