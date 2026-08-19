# Repair status

Branch `tg02-arbitration`. 125 passed, 2 skipped, 2 xfailed.
Read [CLAUDE.md](CLAUDE.md) first — it holds the working rules this file
assumes.

**How to read the verification column.** "Probed" means the endpoint was
actually called and its output inspected this session. "Inherited" means the
code was touched but the mechanism was never individually exercised — treat
those exactly as unverified, because every mechanism probed so far has turned
out to have at least one defect.

---

## Mechanisms

31 of 38 designable. A22 does not exist (known gap, not a renumbering task).

### Probed end to end — scoring, design, and a regression test

| | verdict | what was wrong |
| --- | --- | --- |
| **A1** RNase-H gapmer | fixed | Ranking was a ten-way tie at exactly 100.0 — `_composite_score` normalises against constants tuned for 12–30 nt while every candidate in a run shares one length, so all of it clipped. Now ranks on ViennaRNA site accessibility (1000-fold spread), ΔG demoted to tie-break. |
| **A2** steric translation block | no fix needed | Correctly separated from A1 by F10a/F10b; targets the 5′ initiation region; honest that no approved drug has this primary mechanism. |
| **A3** poison-exon blocking | fixed | Could never have worked: a poison exon is absent from the canonical mRNA, which is what it was tiling. Now locates it from NMD-transcript exon diffs — finds SCN1A exon 20N (64 nt), the STK-001 target. |
| **A4** AntagoNAT | fixed | Tiled the whole sense transcript; SCN1A-AS1 overlaps only part of it, so candidates outside the overlap bound nothing. Now confined to the real 198,166 nt overlap. |

### Fixed but never individually probed — assume defects remain

| mechanisms | last change | what has NOT been checked |
| --- | --- | --- |
| A7–A11 (TG04) | U-alphabet fixes in `_calc_tm` / `_reverse_complement` | whether each mechanism targets its own geometry, or all five share one tiling |
| A13, A16, A17, A20 (TG03) | mechanism/edit-type gate added | guide geometry per platform; bystander handling beyond A18/A19 |
| A5, A6, A23, A28 (TG02) | none | **still on the generic tiler, so the 100.0 saturation that hit A3/A4 almost certainly applies** (inferred from the shared code path, not measured) |
| A12, A14 (TG05) | none | repeat-tract phase logic beyond the 3-candidate check |
| A27, A29–A31 (TG06) | element alias fix | per-mechanism target regions |
| A21, A18/A19, A32/A33 | built this session | only build-time probes; no independent review |
| A24/A26 (TG08) | routed to protein-replacement | that the linear/circRNA split matches what each mechanism means |

### Not designable — each carries a stated reason

A25, A37–A39 (aptamers: sequence is a SELEX output).
A34 (needs a pre-trans-splicing molecule), A35 (a small molecule acting on the
ribosome — no complementary sequence exists to design), A36 (a tRNA body
requiring a curated per-isoacceptor table).

A test fails if any mechanism yields no candidates **and** no reason.

---

## Pending

### Known defects, not yet fixed

1. **Score saturation in `gene_upregulation_service`** for A5/A6/A23/A28. The
   same `_composite_score` clipping fixed in `gene_silencing_service`. A3/A4
   escaped it by being routed to a dedicated designer, not by the bug being
   fixed.
2. **`exon_cds_map` uses a proportional estimate** in
   `gene_upregulation_service` (`seq_len * exon.length / total_genomic`)
   rather than the real `cdsStart`/`cdsEnd`. `gene_silencing_service`
   documents why that is wrong — UTR is not spread evenly across exons — and
   prefers the real mapping. Exon numbers in TG02 output are therefore
   approximate.
3. **Four new endpoints are not wired to the frontend**: `A21/sirna-duplex`,
   `editor-guide`, `A33/intron-retention`, `A32/alternative-promoter`. They
   work over the API and are unreachable from the UI.

### Blocked on decisions or data that are yours to make

- **F7 (miRNA site)** needs TargetScan; **F8 (promoter methylation)** needs a
  methylation atlas. Until wired, A6 and A23 halt — correctly.
- **SpliceAI is installed and wired but not calibrated** (`calibrated: False`).
  M3 needs a calibration reference set, and the spec warns against picking one
  SpliceAI was trained on.
- **M1, M4, M5, M7 unrun.** M1 needs the siRNA gene re-annotation (BLAST
  backend, database release and contact email are your calls — the remote
  route submits ~4,000 queries under your identity). M4 needs a specific
  GENCODE release. M5 needs ribosome profiling plus a decision about what the
  target actually is.
- **11 reference tables empty**, blocked on the SO-DATA-04 licence review.
- The benchmark **cannot answer siRNA or splice-switching generalisation** —
  11 siRNA experiments, 5 splice genes. No method fixes that.

---

## Starting the next session

Paste this, adjusted for the mechanism pair:

> Read CLAUDE.md and STATUS.md. Working on **A5 and A6** only.
>
> Done means, and show me the output of each:
> 1. `POST /api/mechanisms/gene-upregulation` returns them scored or halted
>    with a stated reason.
> 2. `POST /api/gene-upregulation/generate` returns ≤ 20 candidates whose
>    sequences derive from the region the mechanism actually targets — verify
>    one by hand against the fetched transcript.
> 3. `compositeScore` has more than one distinct value across the candidates.
> 4. `python3 -m pytest backend/tests -q -k upregulation` passes.
>
> Write the probe before the fix. Add a regression test for anything you find.

Suggested order — earlier pairs are the ones the platform leans on hardest:
**A5/A6 → A7/A8 → A9/A10 → A11/A12 → A13/A14 → A15/A16 → A17/A18 →
A19/A20 → A21/A23 → A27/A28 → A29/A30 → A31/A32 → A33**.

---

## Parallel sessions

Useful, but not on arbitrary pairs — the mechanisms share files.

**Safe to run at the same time** (disjoint services, disjoint tests):

| session | files it owns |
| --- | --- |
| TG02 — A5, A6, A23, A28 | `gene_upregulation_service.py`, `upregulation_targets_service.py` |
| TG03 — A13, A16, A17, A20 | `rna_editing_service.py`, `programmable_editor_service.py` |
| TG05 — A12, A14 | `rna_neutralization_service.py` |
| TG06 — A27, A29–A31 | `translational_regulation_service.py` |
| frontend wiring | `frontend/**` only |

**Do not parallelise:**

- Anything touching `gene_silencing_service.py` — A1, A2, A12, A15 *and* the
  TG04 dispatch all route through it.
- Two sessions appending to `backend/tests/test_isoform_and_element_alias.py`.
  That file is append-only by habit and will conflict every time. Give each
  session its own test file, e.g. `test_tg02_mechanisms.py`.
- Anything touching `mechanism_arbitration.py` or `feature_service.py` —
  changes there affect every mechanism at once.

**Two practical limits.** Ensembl rate-limits: a single session already hit
429s and read timeouts often enough to need a retry layer, and parallel
sessions will trip each other's. Two or three at once is realistic; six is
not. And each session should work on **its own branch off `tg02-arbitration`**,
merging one at a time — same-branch parallel work will conflict on
`rulebooks/*/rule.json` and the shared test file.

If you only want one thing running, do the frontend wiring in parallel with
mechanism work: it touches no backend file.
