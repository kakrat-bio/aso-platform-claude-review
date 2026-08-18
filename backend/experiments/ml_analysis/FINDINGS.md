# ML analysis suite — findings

Twelve experiments on `unified_benchmark.parquet` (165,237 rows after the
cross-modality overlap is dropped). Every number below is in
`results/<experiment>.json`; `run_all.py` regenerates all twelve in about
11 minutes from cached models, or roughly 4 hours from scratch on 8 CPU
cores.

Read the [conventions section of the README](README.md#conventions-that-apply-to-every-number-here)
first — in particular, effective n is **genes**, not rows, and the split
each modality gets is the strongest one its columns can actually support.

---

## The six findings that matter

### 1. The platform's current heuristic ranks below chance

`gene_silencing_service._composite_score` sorts candidates by
`0.65 * duplex-ΔG + 0.35 * Tm-fit`. Scored against real held-out
measurements over 528 experiment groups:

| policy | NDCG@10 | top-10 | mean activity percentile of its top 5 | best oligo recovered |
| --- | ---: | ---: | ---: | ---: |
| `_composite_score` | 0.504 | 0.144 | **46.8** | 7.4% |
| learned ranker | 0.725 | 0.340 | 69.1 | 27.5% |
| drawn random | 0.541 | 0.166 | 51.3 | 10.0% |

Taking the heuristic's top 5 gets you oligos at the 46.8th activity
percentile. Picking at random gets you the 50th. The gap to random is small
but consistent and it is in the wrong direction (NDCG@10 −0.038 vs random,
p = 0.000; top-10 −0.022, p = 0.003).

**Why**, decomposed on 4,000 held-out RNase-H oligos:

* `_composite_score` correlates **0.67 with GC content** and **0.033 with
  measured activity**.
* 54% of oligos sit at exactly `tm_fit == 100` — the plateau inside the
  A1 (50–65 °C) window — so for most of the pool that term contributes
  nothing but ties. A further 15% of duplex scores saturate at 100.
* The two remaining discriminating signals, duplex ΔG and Tm, are both
  monotone in GC at fixed length. The composite is a GC ranker.

E8 shows the learned model going the other way: C carries **negative**
attribution and U the most positive. E7 finds group difficulty correlates
−0.24 with mean GC. High-GC-first is the wrong sort order for this data.

Two coverage gaps found by running it, worth fixing regardless of the
ranking question:

* `CHEM_TM_BOOST` has keys for `gapmer`, `lna_gapmer`, `pmo`, `2ome` — and
  **none for MOE or cEt**, which are essentially the entire RNase-H arm and
  all of splice-switching. Every benchmark oligo gets a boost of 0 and its
  Tm-fit is computed from the unmodified-DNA Tm.
* `OPTIMAL_TM_RANGES` has no siRNA key, so that arm falls through to the
  default (50, 70) window. Its heuristic row is the heuristic applied
  outside its intended scope and should be read that way.

> A confound that had to be removed before this result could be believed:
> the heuristic rounds to one decimal and ties **27%** of oligos within an
> experiment, with a tied top-5 cutoff in **58%** of groups. `np.argsort` is
> stable, so those ties were being resolved by row order — which in this
> benchmark is tiling order along the transcript. Every metric in the suite
> now shuffles each group before sorting.

### 2. A random row split roughly doubles the top-10 number

Same mechanism, same model, same metrics; only the split changes.

| modality | rows/gene | random top-10 | gene top-10 | random NDCG@10 | gene NDCG@10 | random Pearson | gene Pearson |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `rnase_h` | 469 | 0.667 | 0.345 | 0.794 | 0.727 | 0.394 | 0.380 |
| `splice_switching` | 436 | 0.605 | 0.100 | 0.758 | 0.454 | 0.323 | **−0.158** |
| `sirna` | 1.00 | 0.699 | *refused* | 0.699 | *refused* | 0.325 | *refused* |

The important part is which metric hides it. For RNase-H, NDCG@10 moves
0.067 and Pearson moves 0.015 — you could report either and sound
reproducible — while **top-10 overlap moves 0.32 and MRR moves 0.19**. Any
result quoted from a random split should be quoted as top-k or not at all.

For splice-switching, holding out a whole gene takes Pearson to −0.158: on
the one gene it has never seen, the model is anti-correlated with activity.
That is four experiment groups from a single held-out gene — the only test
this arm can offer, and not a basis for anything stronger than a warning.

The siRNA gene split was **refused** by `common.gene_split()`: 3,947
distinct `target_gene` values for 3,947 rows. The column holds the target
site, not a gene. That refusal is the result for that arm.

### 3. The chemistry-invariant model — "the method centerpiece" — collapses

Three model classes, identical data, split and budget, paired over the same
528 groups, Holm-corrected across all 12 comparisons:

| model | NDCG@10 | top-10 | Pearson | MRR |
| --- | ---: | ---: | ---: | ---: |
| conditioned | 0.725 | 0.340 | 0.366 | 0.188 |
| sequence-only | 0.712 | 0.329 | 0.340 | 0.191 |
| invariant (adversarial) | **0.600** | **0.231** | **0.134** | **0.123** |

* conditioned − invariant: −0.126 NDCG@10 (Holm p ≈ 1.6e−55), −0.232
  Pearson (p ≈ 9.6e−73). The gradient-reversal model is worse than feeding
  chemistry in *and* worse than ignoring chemistry entirely.
* conditioned − sequence-only: +0.013 NDCG@10 (Holm p = 0.012) and +0.026
  Pearson (p = 2.2e−5) are significant; **top-10 (+0.011, p = 0.095) and
  MRR (−0.002, p = 0.22) are not**. Conditioning improves the shape of the
  ranking without measurably improving what a user takes off the top.

`invariant_ranker.py` documents the domain-adversarial model as the method
centerpiece. Under a grouped split on this benchmark it is not defensible as
one.

### 4. The imbalance is worse than the row counts say, and correcting it costs

RNase-H is 95.5% of training rows but **98.03% of training pairs**
(88,358 of ~90,131 per epoch) — the pairwise loss concentrates it further,
because pairs are capped per experiment and RNase-H owns 1,410 of 1,438
groups. siRNA contributes 512 pairs; splice-switching 1,261.

Four weightings, one protocol:

| scheme | rnase_h NDCG@10 | Δ vs none | p | sirna NDCG@10 | splice NDCG@10 |
| --- | ---: | ---: | ---: | ---: | ---: |
| none | 0.725 | — | — | 0.816 | 0.536 |
| `inverse_rows` (siRNA ×34, splice ×55) | 0.702 | −0.024 | 0.000 | 0.776 | 0.601 |
| `inverse_exps` | 0.687 | −0.039 | 0.000 | 0.847 | 0.602 |
| `sqrt_inverse` | 0.727 | +0.001 | 0.72 | 0.878 | 0.488 |

Aggressive inverse weighting costs RNase-H a real, well-measured amount.
`sqrt_inverse` is free on RNase-H.

**But the columns this experiment is about cannot be read.** The held-out
set contains **3 siRNA groups and 2 splice-switching groups**, because the
benchmark contains only 11 siRNA experiments and 5 splice genes in total. A
paired bootstrap over 3 numbers resamples 3 numbers; its interval can come
out narrow by accident. Every p-value in the siRNA and splice columns is
indicative only, and no conclusion about whether upweighting *helps* the
scarce arms is supportable from this data. `results/exp01_class_imbalance.json`
records this in an `evaluation_power` block rather than leaving it to be
noticed.

### 5. The specialists don't rescue the scarce arms either

| modality | pair share | pooled NDCG@10 | specialist | Δ | p |
| --- | ---: | ---: | ---: | ---: | ---: |
| `rnase_h` | 98.03% | 0.725 | 0.719 | −0.007 | 0.10 |
| `sirna` | 0.57% | 0.816 | 0.812 | −0.004 | 0.78 |
| `splice_switching` | 1.40% | 0.536 | 0.628 | +0.093 | 0.51 (n = 2) |

The pooled ranker *is* an RNase-H ranker by composition, but training the
scarce arms alone does not beat it — the siRNA specialist has 8 training
groups and 512 pairs to learn from. Taken with finding 4, the bottleneck is
minority **data volume**, not the pooling and not the loss weighting. No
reweighting scheme creates siRNA experiments.

### 6. The generator has collapsed onto the wrong gapmer length

The benchmark's chemistry strings are position-level annotations
(`L20 MOE|sugar|1,2,3,4,5,16,17,18,19,20 …`), so the architecture is parsed,
not inferred. The dominant layout is **5-10-5 MOE at 20 nt (88,966 rows)**,
then 3-10-3 cEt at 16 nt (51,625 rows).

What the CVAE emits for `rnase_h`:

| length | 15 | **16** | 17 | **18** | 19 | **20** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| designs | 5 | **331** | 5 | 157 | 1 | **1** |

**0.2% of RNase-H designs are 20 nt** — the length a 5-10-5 gapmer requires.
66% are 16 nt. A 16-mer cannot be built as the architecture that accounts for
57% of the training rows, at all. (Splice-switching is fine: 75% land on the
18 nt uniform-MOE length that arm actually uses.)

E10 shows the same collapse thermodynamically: RNase-H designs bind
**3.4 kcal/mol weaker** and melt **2.9 °C lower** than the training oligos
(Cohen's d 0.70 and 0.54; KS p ≈ 2e−49 and 7e−45). siRNA designs, whose modal
length of 19 nt matches, are indistinguishable in effect size (d = 0.05).

---

## The rest

### E4 — ranking metrics at every k

| scorer | NDCG@10 | MAP@10 | top-10 | MRR |
| --- | ---: | ---: | ---: | ---: |
| model | 0.725 [0.701, 0.743] | 0.548 | 0.340 [0.314, 0.366] | 0.188 |
| drawn random | 0.541 [0.531, 0.551] | 0.293 | 0.166 [0.145, 0.192] | 0.079 |
| label oracle | 1.000 | 1.000 | 1.000 | 0.944 |

CIs resample genes. Two things to read carefully: NDCG cannot fall near zero
here because `rank_label` is a dense percentile — the random row, not 0, is
the floor. And the oracle's MRR is 0.944 rather than 1.0 because
`rank_label` has ties; that is the data, not a bug.

### E5 — cross-chemistry transfer costs 0.035 NDCG@10

Sequence-only models, so an off-diagonal cell is not measuring an untrained
chemistry embedding. Restricted to RNase-H, because across the whole
benchmark chemistry and mechanism are confounded (every splice-switching row
is MOE, every siRNA row unmodified).

top-10 overlap, rows = train, columns = test:

| | MOE | MOE+cEt | cEt |
| --- | ---: | ---: | ---: |
| **MOE** (94,818) | **0.318** | 0.311 | 0.260 |
| **MOE+cEt** (9,921) | 0.237 | 0.307 | 0.275 |
| **cEt** (52,141) | 0.249 | 0.311 | **0.321** |

Mean diagonal 0.704 NDCG@10, off-diagonal 0.670. MOE and cEt are each best
predicted by their own model — but **MOE+cEt, the smallest class, is
predicted slightly better by both large-class models (0.311) than by its own
(0.307)**. Below some data volume, borrowing beats specialising.

### E8 — attribution is method-robust and points 5'

Integrated gradients (midpoint rule, 64 steps, zero baseline) and occlusion
agree at r = 0.71–0.75 per position. IG completeness holds (median relative
error 0.005–0.014), which is checked and reported because a violated
completeness check means the attributions should not be read at all.

Peak attribution sits at position 3–5 from the 5' end in **every** modality
and every fixed length, decaying monotonically toward the 3' end. For a
20-mer that is the 5' wing of the 5-10-5. Per-nucleotide mean IG:

| | A | C | G | U |
| --- | ---: | ---: | ---: | ---: |
| rnase_h | +0.0005 | **−0.0013** | +0.0149 | **+0.0175** |
| sirna | +0.0102 | **−0.0115** | −0.0041 | **+0.0255** |
| splice_switching | +0.0094 | +0.0019 | +0.0037 | **+0.0176** |

U is the most favoured base everywhere and C is penalised in two of three
arms — the opposite of the GC-first ordering in finding 1.

### E9 — Reynolds rules don't separate active from inactive here

Scored on the 3,947 19-mers of the siRNA arm. Mean 3.79 of 8 criteria.
Comparing the top 20% of measured activity against the bottom 20%:

**3.919 of 8 versus 3.910 — a difference of 0.009, 95% CI [−0.149, +0.165],
p = 0.92** (unpaired two-sample bootstrap over 836 and 790 sequences; they
are different oligos, so pairing them would invent a correspondence).

The rules do not distinguish active from inactive siRNAs in this benchmark.
A generator that matches them is matching a convention, not predicted
activity — and indeed the designs (3.76) sit between the dinucleotide-shuffled
control (3.62) and uniform random (3.98).

Recurrent 4-mers, against a dinucleotide-shuffled control that holds length,
composition and dinucleotide frequency fixed: **no k-mer appears in the top 15
of more than one mechanism**. Pairwise Jaccard 0.034–0.071. The mechanisms'
sequence preferences are disjoint at this resolution.

### E10 — thermodynamics

Tm and hairpin ΔG from primer3 (SantaLucia), MFE and duplex ΔG from
ViennaRNA 2.7.0. Everything is computed on the **unmodified backbone** —
neither tool models MOE, cEt or PS, and those shift real duplex stability
substantially. The numbers compare sequence sets to each other; they are not
predictions of a modified oligo's Tm.

Generated vs training, mean difference (Cohen's d):

| | Tm (°C) | duplex ΔG (kcal/mol) | MFE |
| --- | ---: | ---: | ---: |
| rnase_h | **−2.92 (0.54)** | **+3.38 (0.70)** | +0.42 (0.30) |
| sirna | −0.23 (0.05) | +0.18 (0.05) | +0.32 (0.21) |
| splice_switching | −0.84 (0.15) | +0.65 (0.14) | +0.43 (0.28) |

KS rejects equality on nearly every cell, including siRNA's Tm at d = 0.05
(p ≈ 2e−24). At n ≈ 2,000 per set it detects differences far too small to
matter, which is why the effect size is the column to read.

Matching a training distribution shows the generator learned the training
composition. It is not evidence of activity, and it is not claimed as such.

### E11 — frozen RNA-FM buys nothing on 12–28 nt oligos

Subsampled to 36,882 rows (whole genes) so the encoder could be run; every
scorer sees the identical subsample, so the comparison is fair even though
the absolute numbers are not comparable to E1–E7.

| scorer | NDCG@10 | top-10 | Pearson | Δ NDCG@10 vs conv | p |
| --- | ---: | ---: | ---: | ---: | ---: |
| RNA-FM + MLP | 0.672 | 0.279 | 0.287 | +0.001 | 0.87 |
| conv ranker | 0.670 | 0.264 | 0.276 | — | — |
| RNA-FM linear probe | 0.665 | 0.276 | 0.262 | −0.006 | 0.53 |
| 4-mer ridge | 0.661 | 0.259 | 0.234 | −0.010 | 0.25 |
| drawn random | 0.540 | 0.134 | 0.019 | −0.130 | 0.000 |

A frozen RNA-FM representation is indistinguishable from a small supervised
conv net, and a *linear probe* on it is too. RNA-FM was pretrained on ncRNA
hundreds of nucleotides long; this is evidence about transfer to very short
sequences, not about RNA-FM generally. The encoder was not fine-tuned.

> **A silent failure worth knowing about.**
> `RnaFmModel.from_pretrained("multimolecule/rnafm")` returns a model with
> **all 205 tensors randomly initialised**. The published checkpoint is an
> `RnaFmForPreTraining` export storing the base model under a `model.`
> prefix, while `RnaFmModel.base_model_prefix` is `rnafm`, so nothing
> matches and transformers initialises everything from scratch — announcing
> it only in a warning inside a wall of text. `rnafm.py` strips the prefix,
> loads explicitly, and `assert_loaded()` raises unless every non-pooler
> tensor arrived. Anything embedded via the obvious route and labelled
> "RNA-FM" would be a randomly initialised transformer.

### E7 — where it fails

527 groups analysed (one excluded for having no label spread — that
hypothesis is not the story here).

Hardest genes: CTGF 0.599, STAT3 0.641, F7 0.645, SPDEF 0.658, SCN2A 0.701.
Easiest: PRNP 0.808, DNM2 0.801, TMPRSS6 0.794, ATXN2 0.783.
By chemistry class the spread is narrow: MOE 0.730, cEt 0.717,
MOE+cEt 0.715, unmodified 0.710.

Group difficulty correlates with mean GC (r = −0.24) and with nothing else
(|r| < 0.10 for group size, label spread, label range, length spread).

**The misses have no signature.** Oligos the model wrongly promoted into its
top 10 versus the ones that belonged there: GC 0.406 vs 0.413, length 18.51
vs 18.50. Whatever the model gets wrong is not visible in composition or
length, which rules out the cheap fixes.

---

## What this suite says about the platform

1. **Replace `_composite_score` as a sort order.** It is a GC ranker and it
   selects below chance. The learned ranker moves the mean activity
   percentile of a top-5 selection from 46.8 to 69.1 and more than triples
   the rate at which the single best oligo is recovered (7.4% → 27.5%).
2. **Add MOE and cEt to `CHEM_TM_BOOST`, or stop applying the Tm-fit term to
   them.** Right now they silently take a boost of 0.
3. **Quote top-k, not NDCG, and only from a gene split.** Finding 2 shows
   NDCG is the metric that best conceals a leaky split.
4. **Retire or re-justify the adversarial invariant ranker.** Finding 3.
5. **Fix the generator's length prior for RNase-H before its designs are
   used.** Finding 6 — 0.2% buildable as the dominant architecture.
6. **The scarce-modality questions need more data, not better methods.**
   Findings 4 and 5. With 11 siRNA experiments and 5 splice genes, the
   benchmark cannot currently answer them, and the honest move is to say so
   rather than quote a p-value from 3 groups.

## What is not claimed

* No causal claim about *why* a sequence works. Every number here is
  correlational, on a benchmark whose labels are within-experiment activity
  percentiles from heterogeneous assays.
* No claim that the generated designs are active. E10 tests distribution
  matching, which is a weaker thing.
* No transfer claim beyond the three chemistry classes with enough data;
  `unmodified` (1,578 rows), `PS` (628), `OMe+cEt` (22) and `F+cEt` (1) were
  excluded from E5 rather than reported thinly.
* Nothing about siRNA or splice-switching generalisation to unseen genes.
  One arm cannot be gene-split at all and the other has five genes.
