"""L1 evidence-feature layer for mechanism arbitration.

Implements the fixed feature vocabulary from
`docs/planning/scoring_and_ml_plan.md` §3.1 and the gap list in
`docs/planning/therapeutic_goal_scope_plan_v3.md` §6.

WHAT A FEATURE IS
-----------------
A feature is a statement about the *target transcript*, not about the user's
intent. "This transcript contains an NMD-inducing exon" is a feature. "The
user wants to upregulate this gene" is not — that is the therapeutic goal,
which under inverted routing is an output, not an input.

Every feature resolves to a `Feature` carrying:

  state       PRESENT / ABSENT / UNRESOLVED
  probability a number in [0,1], or None when UNRESOLVED
  provenance  how we know it — this is what caps confidence downstream
  source      the specific thing consulted
  stand_in    True when the "source" is really the user's own form input

The two non-negotiables carried over from the plan:

  * An absent feature returns **ABSENT**, never probability zero. ABSENT
    means "we looked and it is not there"; UNRESOLVED means "we have no way
    to look". They are different facts and they must not collapse into the
    same number.
  * A predicted feature and a literature-confirmed feature never enter the
    score identically. That is what `provenance` is for.

THE SOURCE LADDER
-----------------
Each feature declares an ordered ladder of sources. The first rung that
fires wins, and its provenance tier is recorded. Most ladders end in a
`user_asserted` rung: the user's own form input, echoed back. That rung is
marked `stand_in=True` and capped hard, because a mechanism whose evidence
is the dropdown the user just picked has not been *arbitrated* — it has been
looked up. Making that visible in the output is the point; the companion
plan's whole critique of the current TG04 numbers is that the input already
contains the answer.

Two features have a deliberately EMPTY ladder:

  F11 repressive RBP site      → blocks A28
  F13 polyadenylation usage    → blocks A11

Neither has a wired source and neither has a user input that constitutes
evidence about the transcript, so mechanisms requiring them halt rather than
score (plan §6.2, §6.4, checklist item 9). Naming a target RBP in a form
field says which protein you have in mind; it is not evidence that a
repressive site for it exists in this transcript.

F12 (repeat expansion) keeps a user rung because the plan documents the
current input as user-supplied free text that *should become* an annotation
lookup (§5) — so A14 halts only when nothing is supplied at all.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Callable

from . import reference_tables as RT
from . import spliceai_service as SAI

# ---------------------------------------------------------------------------
# States and provenance
# ---------------------------------------------------------------------------

PRESENT = "PRESENT"
ABSENT = "ABSENT"
UNRESOLVED = "UNRESOLVED"

# Provenance tiers, strongest first. The cap is the ceiling a mechanism's
# confidence can reach when this is the weakest feature supporting it
# (plan §3.3: quality and reliability stay on separate axes).
#
# The companion plan describes each feature as a triple carrying both an
# evidence tier and a provenance tier. In this implementation the two
# collapse: every rung of every ladder below is distinguished by *where the
# statement came from*, and no second axis had a distinct rule attached to
# it. Carrying an `evidence_tier` field with no rule that reads it would be
# decoration. If a genuine second axis appears later, it belongs here.
MEASURED = "measured"
CONFIRMED = "confirmed"
ANNOTATION = "annotation"
PREDICTED = "predicted"
USER_ASSERTED = "user_asserted"

PROVENANCE_CAP: dict[str, float] = {
    MEASURED: 1.00,
    CONFIRMED: 0.95,
    ANNOTATION: 0.90,
    PREDICTED: 0.75,
    USER_ASSERTED: 0.60,
}

PROVENANCE_LABEL: dict[str, str] = {
    MEASURED: "Experimentally validated for this transcript",
    CONFIRMED: "Curated catalogue entry, carrying a citation",
    ANNOTATION: "Genome annotation lookup",
    PREDICTED: "Model prediction",
    USER_ASSERTED: "Supplied by the user on the input form",
}


@dataclass(frozen=True)
class Feature:
    """One resolved feature observation about the target transcript."""

    id: str
    state: str
    probability: float | None = None
    provenance: str | None = None
    source: str | None = None
    stand_in: bool = False
    detail: str | None = None
    # Multi-way classification for features whose answer is a category rather
    # than a yes/no — P2 (ABUNDANT / LOW / ABSENT_IN_TISSUE) and B1
    # (secreted / membrane / intracellular / unknown). `state` still carries
    # the gating answer; `call` carries the distinction underneath it.
    call: str | None = None

    @property
    def resolved(self) -> bool:
        return self.state != UNRESOLVED

    @property
    def cap(self) -> float:
        """Confidence ceiling this observation imposes."""
        if self.provenance is None:
            return 0.0
        return PROVENANCE_CAP.get(self.provenance, 0.0)

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "label": FEATURE_CATALOG.get(self.id, {}).get("label", self.id),
            "state": self.state,
            "probability": self.probability,
            "provenance": self.provenance,
            "provenanceLabel": PROVENANCE_LABEL.get(self.provenance or ""),
            "source": self.source,
            "call": self.call,
            "standIn": self.stand_in,
            "detail": self.detail,
        }


def _unresolved(fid: str, why: str) -> Feature:
    return Feature(id=fid, state=UNRESOLVED, detail=why)


# ---------------------------------------------------------------------------
# The vocabulary
# ---------------------------------------------------------------------------

FEATURE_CATALOG: dict[str, dict[str, Any]] = {
    "F1": {
        "label": "Exon weakly recognised by spliceosome",
        "intendedSource": "SpliceAI",
        "wired": True,
    },
    "F2": {
        "label": "Variant creates a cryptic splice site",
        "intendedSource": "SpliceAI + MaxEntScan",
        "wired": True,
    },
    "F3": {
        "label": "Deep-intronic pseudoexon activated",
        "intendedSource": "SpliceAI",
        "wired": True,
    },
    "F4": {
        "label": "Transcript contains an NMD-inducing (poison) exon",
        "intendedSource": "GENCODE nonsense_mediated_decay biotype",
        "wired": True,
    },
    "F5": {
        "label": "Repressive uORF in the 5' UTR",
        "intendedSource": "literature-validated uORF list; 5' UTR ORF scan",
        "wired": True,
    },
    "F6": {
        "label": "Overlapping natural antisense transcript",
        "intendedSource": "annotation lookup",
        "wired": True,
    },
    "F7": {
        "label": "Repressive miRNA site in the 3' UTR",
        "intendedSource": "TargetScan context++",
        "wired": False,
    },
    "F8": {
        "label": "Promoter methylation / accessibility in the target tissue",
        "intendedSource": "methylation atlas",
        "wired": False,
        # Serves A23 (activate a silenced promoter) and A15 (silence an
        # accessible one) — the same measurement, opposite intents.
        "serves": ["A23", "A15"],
    },
    "F9": {
        "label": "Allele-distinguishing variant in the transcript",
        "intendedSource": "dbSNP / ClinVar",
        "wired": False,
    },
    "F10a": {
        "label": "Density of accessible designable sites, transcript-wide",
        "intendedSource": "ViennaRNA RNAplfold accessibility",
        "wired": True,
    },
    "F10b": {
        "label": "Density of accessible designable sites, 5' UTR / start codon",
        "intendedSource": "ViennaRNA RNAplfold accessibility",
        "wired": True,
    },
    "F11": {
        "label": "Repressive RNA-binding-protein site in the transcript",
        # Deliberately the curated route, not the atlas route. An atlas
        # answers "an RBP binds here", not "a repressor binds here and
        # masking it raises protein output" — binding is not repression, and
        # most RNA-binding proteins are context-dependent. A short list of
        # validated cases is worth more than genome-wide coverage of a
        # quantity that cannot be interpreted.
        "intendedSource": "curated list of validated repressive RBP sites",
        "wired": False,
        "blocks": ["A28"],
    },
    "F12": {
        "label": "Pathogenic repeat expansion, with unit and length",
        "intendedSource": "annotation lookup against known repeat-expansion loci",
        "wired": False,
        "blocks": ["A14"],
    },
    # F13 is split. Site LOCATION is annotatable; site usage SHIFTING being
    # therapeutic is not. Most human genes carry alternative polyadenylation
    # sites, so F13a alone would fire A11 almost everywhere — A11 requires
    # both halves.
    "F13a": {
        "label": "Alternative polyadenylation site present",
        "intendedSource": "PolyASite / PolyA_DB / APAatlas (MUST VERIFY)",
        "wired": False,
        "blocks": ["A11"],
    },
    "F13b": {
        "label": "Shifting polyadenylation usage is therapeutically beneficial",
        "intendedSource": "per-gene disease-specific curation",
        "wired": False,
        "blocks": ["A11"],
    },
    "F14a": {
        "label": "Alternative promoter exists for this gene",
        "intendedSource": "CAGE / FANTOM promoter atlas (MUST VERIFY)",
        "wired": False,
        "blocks": ["A32"],
    },
    "F14b": {
        "label": "Shifting promoter usage is therapeutically beneficial",
        "intendedSource": "per-gene disease-specific curation",
        "wired": False,
        "blocks": ["A32"],
    },
    "F15a": {
        "label": "Intron with therapeutic retention potential exists",
        "intendedSource": "splice-site conservation + NMD context (MUST VERIFY)",
        "wired": False,
        "blocks": ["A33"],
    },
    "F15b": {
        "label": "Retention of this intron is therapeutically beneficial",
        "intendedSource": "per-gene disease-specific curation",
        "wired": False,
        "blocks": ["A33"],
    },
    # --- modality-flag features (NOT part of the scored family) -------------
    # These three feed the qualitative modality flag and nothing else. They
    # never enter an applicability interval and never contribute to a score.
    # Family T features answer "how well does this mechanism fit this
    # transcript"; these answer "is a transcript-acting mechanism the right
    # kind of answer at all".
    "P2": {
        "label": "Residual endogenous transcript present and boostable",
        "intendedSource": "tissue expression atlas",
        "wired": False,
        "family": "modality_flag",
    },
    "P6": {
        "label": "Dominant-negative allele present",
        "intendedSource": "ClinVar + literature",
        "wired": False,
        "family": "modality_flag",
    },
    "B1": {
        "label": "Target protein is extracellular or cell-surface",
        "intendedSource": "UniProt subcellular localisation",
        "wired": False,
        "family": "modality_flag",
    },
}

MODALITY_FLAG_FEATURES = ("P2", "P6", "B1")

# P2 returns a three-way call, not a number: there is a real difference
# between "expressed but low" and "not expressed in this tissue at all", and
# the second is a much stronger signal that boosting has nothing to work on.
P2_ABUNDANT = "ABUNDANT"
P2_LOW = "LOW"
P2_ABSENT_IN_TISSUE = "ABSENT_IN_TISSUE"

# DE-NOVO PARAMETERS — SO-TG-09 / SO-DATA-01, on the critical path for TG02.
# Neither cut point is derived from any measurement; both need calibration
# register entries before they are relied on. They decide only whether a
# qualitative signpost is shown — no score is multiplied by either.
P2_ABUNDANT_TPM = 5.0
P2_ABSENT_TPM = 0.5

# B1 returns a localisation class. The aptamer flag fires only on the first
# two: pegaptanib is delivered intravitreally against a secreted growth
# factor, and an intracellular target is out of an aptamer's reach.
B1_SECRETED = "secreted"
B1_MEMBRANE = "membrane"
B1_INTRACELLULAR = "intracellular"
B1_UNKNOWN = "unknown"
B1_APTAMER_ACCESSIBLE = (B1_SECRETED, B1_MEMBRANE)


# ---------------------------------------------------------------------------
# Resolution context
# ---------------------------------------------------------------------------

@dataclass
class FeatureContext:
    """Everything the feature layer is allowed to look at.

    Deliberately does NOT carry a therapeutic goal. Goal is an output of
    arbitration, not an input to it.
    """

    gene_symbol: str = ""
    molecular_defect: str | None = None
    known_variant: str | None = None
    variant_hgvs: str | None = None
    transcript_class: str | None = None
    allele_selective: bool | None = None
    gene_features: dict | None = None
    transcript_sequence: str | None = None
    cds_start: int | None = None
    repeat_unit: str | None = None
    repeat_count: str | None = None
    oligo_length: int = 18
    # SpliceAI inputs. A pre-mRNA (genomic, introns retained) sequence is
    # required: splice-site recognition is a statement about exon/intron
    # boundaries, and a mature transcript has none left.
    pre_mrna_sequence: str | None = None
    variant_offset: int | None = None      # 0-based, into pre_mrna_sequence
    variant_alt: str | None = None         # single substituted base
    exon_start: int | None = None          # 0-based offset of the exon of interest
    exon_end: int | None = None            # 0-based, inclusive
    # Modality-flag inputs.
    tissue_tpm: float | None = None
    protein_localisation: str | None = None
    extras: dict = field(default_factory=dict)

    def gene_feature(self, key: str) -> dict | None:
        feats = (self.gene_features or {}).get("features")
        if isinstance(feats, dict):
            entry = feats.get(key)
            if isinstance(entry, dict):
                return entry
        return None


# ---------------------------------------------------------------------------
# Accessibility (F10a / F10b)
# ---------------------------------------------------------------------------

# RNAplfold parameters. Local folding keeps this linear in transcript length —
# a full partition function over a multi-kb mRNA is not affordable in a
# request path.
PLFOLD_WINDOW = 80
PLFOLD_MAX_BP_SPAN = 40

# A window counts as "designable" when the probability that its whole
# oligo-length stretch is unpaired clears this bar.
#
# NOT CALIBRATED. This threshold is a placeholder chosen to be permissive
# enough to separate structured from unstructured transcripts; it has not
# been fitted against measured ASO activity. It shifts how sharply F10a and
# F10b discriminate, never whether a mechanism is eligible, because both are
# discriminating features rather than gates. Calibrating it is a sign-off
# item (see docs/planning/therapeutic_goal_scope_plan_implementation.md).
ACCESSIBLE_SITE_THRESHOLD = 0.05

# How many distinct designable sites count as "enough". Above this, more
# sites do not make the mechanism more applicable — they make the downstream
# tiling step easier, which is feasibility (plan §3.4) and is reported
# separately rather than folded in here. Also uncalibrated; see the note on
# ACCESSIBLE_SITE_THRESHOLD.
DESIGNABLE_SITE_TARGET = 5

# How far past the start codon the 5'-UTR / start-codon window extends.
# Steric-block translation inhibition (A2) acts on the ribosome scanning and
# initiation region, not the whole transcript.
START_CODON_WINDOW_NT = 30

_VALID_RNA = set("ACGU")


def _clean_rna(seq: str) -> str:
    return "".join(b for b in seq.upper().replace("T", "U") if b in _VALID_RNA)


def _unpaired_profile(seq: str, oligo_length: int) -> list[float] | None:
    """Per-position probability that an oligo-length stretch is unpaired.

    Index i holds the probability for the stretch ENDING at 1-based position
    i+1, matching ViennaRNA's pfl_fold_up layout. Returns None when
    ViennaRNA is unavailable or the transcript is shorter than the oligo.
    """
    if len(seq) < oligo_length:
        return None
    try:
        import RNA  # noqa: PLC0415 — optional heavy dependency
    except ImportError:
        return None

    up = RNA.pfl_fold_up(seq, oligo_length, PLFOLD_WINDOW, PLFOLD_MAX_BP_SPAN)
    # up is 1-indexed with a dummy row 0; each row is indexed by stretch length.
    return [up[i][oligo_length] for i in range(oligo_length, len(seq) + 1)]


def _distinct_sites(profile: list[float], oligo_length: int) -> int:
    """Count NON-OVERLAPPING designable windows.

    Adjacent windows share almost all of their sequence, so counting every
    position that clears the bar counts one accessible region many times over.
    Stepping a full oligo length past each hit counts distinct candidate
    sites, which is what "how many oligos could I actually design here"
    means.
    """
    count = 0
    i = 0
    while i < len(profile):
        if profile[i] >= ACCESSIBLE_SITE_THRESHOLD:
            count += 1
            i += oligo_length
        else:
            i += 1
    return count


def _site_sufficiency(count: int) -> float:
    """Turn a site count into a probability that a design target exists.

    Deliberately NOT a density. A1 needs an accessible cleavable site
    anywhere in the transcript; A2 needs one specifically at the 5' UTR or
    start codon. Those windows differ in length by an order of magnitude, so
    comparing densities would penalise A1 for the transcript being long
    rather than for being inaccessible. A saturating count is comparable
    across window sizes: once there are enough distinct sites to design
    against, more sites do not make the mechanism more applicable.
    """
    return min(1.0, count / DESIGNABLE_SITE_TARGET)


def _resolve_accessibility(ctx: FeatureContext) -> tuple[Feature, Feature]:
    """Resolve F10a and F10b from one fold of the transcript.

    Item 8 of the plan: A1 needs an accessible cleavable site ANYWHERE in the
    transcript; A2 needs one specifically at the 5' UTR / start codon. Those
    are different queries over the same accessibility profile, and running
    both is what separates two mechanisms that otherwise tie on half of all
    inputs.
    """
    why = "No transcript sequence supplied — accessibility not computed."
    if not ctx.transcript_sequence:
        return _unresolved("F10a", why), _unresolved("F10b", why)

    seq = _clean_rna(ctx.transcript_sequence)
    profile = _unpaired_profile(seq, ctx.oligo_length)
    if profile is None:
        why = (
            "Transcript shorter than the oligo, or ViennaRNA unavailable — "
            "accessibility not computed."
        )
        return _unresolved("F10a", why), _unresolved("F10b", why)

    whole_n = _distinct_sites(profile, ctx.oligo_length)
    f10a = Feature(
        id="F10a",
        state=PRESENT if whole_n > 0 else ABSENT,
        probability=_site_sufficiency(whole_n),
        provenance=PREDICTED,
        source=f"ViennaRNA RNAplfold over {len(seq)} nt",
        detail=(
            f"{whole_n} distinct accessible {ctx.oligo_length} nt sites across "
            f"the whole transcript"
        ),
    )

    # The 5'-UTR / start-codon window. Without a CDS start we cannot say where
    # it is, so F10b stays unresolved rather than silently reusing F10a.
    if ctx.cds_start is None:
        f10b = _unresolved(
            "F10b",
            "No CDS start position supplied — the 5' UTR / start-codon window "
            "cannot be located, so the initiation-region query was not run.",
        )
        return f10a, f10b

    end = min(len(profile), max(0, ctx.cds_start + START_CODON_WINDOW_NT))
    window = profile[:end]
    if not window:
        f10b = _unresolved(
            "F10b",
            "The 5' UTR / start-codon window is shorter than one oligo length.",
        )
        return f10a, f10b

    local_n = _distinct_sites(window, ctx.oligo_length)
    f10b = Feature(
        id="F10b",
        state=PRESENT if local_n > 0 else ABSENT,
        probability=_site_sufficiency(local_n),
        provenance=PREDICTED,
        source=f"ViennaRNA RNAplfold over the first {end} nt",
        detail=(
            f"{local_n} distinct accessible {ctx.oligo_length} nt sites in the "
            f"5' UTR and first {START_CODON_WINDOW_NT} nt of CDS"
        ),
    )
    return f10a, f10b


# ---------------------------------------------------------------------------
# Ladder rungs
# ---------------------------------------------------------------------------

def _from_annotation(ctx: FeatureContext, fid: str, key: str,
                     source: str) -> Feature | None:
    """Read an Ensembl-derived structural check from the gene-feature payload.

    Only a *verified* entry counts. The payload reports unverifiable genes as
    available=True so the ranking UI does not silently drop them; treating
    that as a positive finding here would be exactly the guess the plan
    forbids, so an unverified entry falls through to the next rung.
    """
    entry = ctx.gene_feature(key)
    if not entry or not entry.get("verified"):
        return None
    available = bool(entry.get("available"))
    return Feature(
        id=fid,
        state=PRESENT if available else ABSENT,
        # Annotation is a yes/no lookup, not a calibrated probability. ABSENT
        # keeps a small non-zero value so a single annotation miss can never
        # drive a product of features to a hard zero — the state is what
        # gates, the number only ranks.
        probability=0.9 if available else 0.05,
        provenance=ANNOTATION,
        source=source,
        detail=entry.get("reason"),
    )


def _from_defect(ctx: FeatureContext, fid: str, defects: set[str],
                 label: str) -> Feature | None:
    """The user's own molecular-defect selection, echoed back as evidence.

    This is the rung that makes today's top-1 numbers look better than the
    system is: the user asserts the defect, the mechanism is gated on that
    defect, and the same assertion then satisfies the mechanism's required
    feature. It is marked `stand_in` and capped so the output says so.
    """
    if not ctx.molecular_defect:
        return None
    if ctx.molecular_defect not in defects:
        return None
    return Feature(
        id=fid,
        state=PRESENT,
        probability=1.0,
        provenance=USER_ASSERTED,
        source=f"user-selected molecular defect '{ctx.molecular_defect}'",
        stand_in=True,
        detail=(
            f"{label} is taken from your own input, not from "
            f"{FEATURE_CATALOG[fid]['intendedSource']}. Recovery numbers on "
            "this path measure lookup, not arbitration."
        ),
    )


def _from_repeat_text(ctx: FeatureContext) -> Feature | None:
    """F12 from the free-text repeat unit / count fields."""
    unit = _normalize_repeat_unit(ctx.repeat_unit)
    count = _extract_repeat_count(ctx.repeat_count)
    if not unit and count is None:
        return None
    if ctx.repeat_unit and ctx.repeat_unit.strip() and not unit:
        return Feature(
            id="F12",
            state=ABSENT,
            probability=0.02,
            provenance=USER_ASSERTED,
            source="user-supplied repeat unit",
            stand_in=True,
            detail=(
                f"'{ctx.repeat_unit}' is not a valid nucleotide repeat motif."
            ),
        )
    if count is not None and count < PATHOGENIC_REPEAT_THRESHOLD:
        return Feature(
            id="F12",
            state=ABSENT,
            probability=0.05,
            provenance=USER_ASSERTED,
            source="user-supplied repeat count",
            stand_in=True,
            detail=(
                f"~{count} copies is below the pathogenic expansion threshold "
                f"(~{PATHOGENIC_REPEAT_THRESHOLD})."
            ),
        )
    known = KNOWN_REPEAT_UNITS.get(unit or "")
    return Feature(
        id="F12",
        state=PRESENT,
        probability=0.9 if known else 0.7,
        provenance=USER_ASSERTED,
        source="user-supplied repeat unit / count",
        stand_in=True,
        detail=(
            f"Repeat unit {unit} recognised ({known})"
            if known
            else f"Repeat unit {unit} accepted as a nucleotide repeat motif "
                 "(not in the curated reference list)"
        ),
    )


def _from_variant_text(ctx: FeatureContext) -> Feature | None:
    """F9 from a user-supplied variant description.

    Only resolved when allele-selective silencing is actually requested. F9
    asks whether there is a variant the oligo can discriminate on; for a
    total-knockdown design there is nothing to discriminate, so a variant
    description in the free-text box is not evidence about anything the
    design depends on.
    """
    if not ctx.allele_selective:
        return None
    text = (ctx.variant_hgvs or ctx.known_variant or "").strip()
    if not text:
        return None
    return Feature(
        id="F9",
        state=PRESENT,
        probability=0.8,
        provenance=USER_ASSERTED,
        source="user-supplied variant description",
        stand_in=True,
        detail=(
            f"'{text}' taken as an allele-distinguishing variant. Not checked "
            "against dbSNP or ClinVar — no variant database is wired."
        ),
    )


def _repeat_from_catalogue(ctx: FeatureContext) -> Feature | None:
    """F12 from the curated repeat-expansion catalogue.

    The set of pathogenic repeat-expansion loci is small and closed — roughly
    sixty diseases — so this is a lookup table, not a predictor. A gene that
    is NOT in the catalogue leaves this rung unfired and falls through to the
    user-supplied text; it is never read as "no repeat here".
    """
    row = RT.row_for("repeat_expansion_loci", ctx.gene_symbol)
    if not row:
        return None
    unit = _normalize_repeat_unit(row.get("repeat_unit"))
    region = (row.get("transcript_region") or "").strip()
    return Feature(
        id="F12",
        state=PRESENT,
        probability=0.95,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "repeat-expansion catalogue",
        call=region or None,
        detail=(
            f"{ctx.gene_symbol} is a curated repeat-expansion locus: "
            f"{unit or row.get('repeat_unit')}"
            + (f" in the {region}" if region else "")
            + ". Repeat location is a design input, not just metadata — it "
              "decides whether and where the transcript is targetable."
        ),
    )


def _localisation_from_table(ctx: FeatureContext) -> Feature | None:
    """B1 from the protein-localisation table."""
    row = RT.row_for("protein_localisation", ctx.gene_symbol)
    if not row:
        return None
    call = (row.get("localisation_class") or "").strip().lower()
    if call not in (B1_SECRETED, B1_MEMBRANE, B1_INTRACELLULAR):
        return None
    accessible = call in B1_APTAMER_ACCESSIBLE
    return Feature(
        id="B1",
        state=PRESENT if accessible else ABSENT,
        probability=0.9 if accessible else 0.05,
        provenance=ANNOTATION,
        source=RT.provenance_of(row) or "protein localisation table",
        call=call,
        detail=(
            f"{ctx.gene_symbol} classified as {call}"
            + (" (signal peptide annotated)"
               if (row.get("has_signal_peptide") or "").strip().lower()
               in ("1", "true", "yes") else "")
        ),
    )


def _expression_from_table(ctx: FeatureContext) -> Feature | None:
    """P2 from the tissue-expression table, when a TPM was not passed in."""
    row = RT.row_for("tissue_expression", ctx.gene_symbol)
    if not row:
        return None
    try:
        tpm = float(row.get("median_tpm"))
    except (TypeError, ValueError):
        return None
    return _classify_expression(
        tpm, RT.provenance_of(row) or "tissue expression table")


def _dominant_negative_from_table(ctx: FeatureContext) -> Feature | None:
    """P6 from the curated dominant-negative gene list.

    Presence in this list SUPPRESSES the replacement flag: supplying more
    wild-type protein does not fix a disease driven by a mutant product that
    actively interferes.
    """
    row = RT.row_for("dominant_negative_genes", ctx.gene_symbol)
    if not row:
        return None
    return Feature(
        id="P6",
        state=PRESENT,
        probability=0.9,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "curated dominant-negative list",
        detail=(
            f"{ctx.gene_symbol} is curated as having a documented "
            f"dominant-negative mechanism"
            + (f" (PMID {row.get('pmid')})" if row.get("pmid") else "")
            + ". " + (row.get("evidence_summary") or "")
        ).strip(),
    )


def _classify_expression(tpm: float, source: str) -> Feature:
    """Turn a TPM into the three-way P2 call.

    Three-way rather than binary because "expressed but low" and "not
    expressed in this tissue at all" are different findings, and the second
    is a much stronger signal that a boosting mechanism has nothing to act on.
    """
    if tpm >= P2_ABUNDANT_TPM:
        call, state, prob, note = (
            P2_ABUNDANT, PRESENT, 0.9,
            "enough endogenous transcript for a boosting mechanism to act on",
        )
    elif tpm >= P2_ABSENT_TPM:
        call, state, prob, note = (
            P2_LOW, ABSENT, 0.2,
            "expressed, but too low for boosting to have much to work with",
        )
    else:
        call, state, prob, note = (
            P2_ABSENT_IN_TISSUE, ABSENT, 0.02,
            "effectively not expressed in this tissue — there is no "
            "transcript to boost",
        )
    return Feature(
        id="P2",
        state=state,
        probability=prob,
        provenance=ANNOTATION,
        source=source,
        call=call,
        detail=(
            f"{tpm:.1f} TPM in the target tissue: {call} — {note} "
            f"(cut points {P2_ABSENT_TPM} / {P2_ABUNDANT_TPM} TPM are de-novo "
            f"parameters, SO-TG-09 / SO-DATA-01)"
        ),
    )


def _rbp_site_from_curation(ctx: FeatureContext) -> Feature | None:
    """F11 from the curated validated-repressive-site list."""
    row = RT.row_for("rbp_repressor_sites", ctx.gene_symbol)
    if not row:
        return None
    return Feature(
        id="F11",
        state=PRESENT,
        probability=0.9,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "curated repressive RBP site list",
        call=(row.get("transcript_region") or "").strip() or None,
        detail=(
            f"{row.get('rbp')} has a validated repressive site in "
            f"{ctx.gene_symbol}"
            + (f" ({row.get('transcript_region')})" if row.get("transcript_region") else "")
            + (f", PMID {row.get('pmid')}" if row.get("pmid") else "")
        ),
    )


def _apa_site_from_table(ctx: FeatureContext) -> Feature | None:
    """F13a — is there an alternative polyadenylation site at all?"""
    row = RT.row_for("polyadenylation_sites", ctx.gene_symbol)
    if not row:
        return None
    n = (row.get("alternative_site_count") or "").strip()
    try:
        count = int(n)
    except (TypeError, ValueError):
        return None
    return Feature(
        id="F13a",
        state=PRESENT if count > 0 else ABSENT,
        probability=0.9 if count > 0 else 0.05,
        provenance=ANNOTATION,
        source=RT.provenance_of(row) or "polyadenylation site atlas",
        detail=(
            f"{count} alternative polyadenylation site(s) annotated for "
            f"{ctx.gene_symbol}. Presence alone is not sufficient for A11 — "
            "most human genes have them."
        ),
    )


def _apa_benefit_from_curation(ctx: FeatureContext) -> Feature | None:
    """F13b — is shifting the site actually therapeutic for this gene?"""
    row = RT.row_for("apa_therapeutic_benefit", ctx.gene_symbol)
    if not row:
        return None
    return Feature(
        id="F13b",
        state=PRESENT,
        probability=0.9,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "curated APA benefit list",
        detail=(
            (row.get("evidence_summary") or
             f"Shifting polyadenylation usage is curated as therapeutic in "
             f"{ctx.gene_symbol}")
            + (f" (PMID {row.get('pmid')})" if row.get("pmid") else "")
        ),
    )


def _alt_promoter_from_table(ctx: FeatureContext) -> Feature | None:
    """F14a — does this gene have an annotated alternative promoter?"""
    row = RT.row_for("alt_promoters", ctx.gene_symbol)
    if not row:
        return None
    try:
        count = int((row.get("alternative_promoter_count") or "").strip())
    except (TypeError, ValueError):
        return None
    return Feature(
        id="F14a",
        state=PRESENT if count > 0 else ABSENT,
        probability=0.9 if count > 0 else 0.05,
        provenance=ANNOTATION,
        source=RT.provenance_of(row) or "alternative promoter atlas",
        detail=(
            f"{count} alternative promoter(s) annotated for {ctx.gene_symbol}. "
            "Presence alone is not sufficient for A32 — many genes have them."
        ),
    )


def _alt_promoter_benefit(ctx: FeatureContext) -> Feature | None:
    """F14b — is shifting promoter usage therapeutic in this gene?"""
    row = RT.row_for("alt_promoter_benefit", ctx.gene_symbol)
    if not row:
        return None
    return Feature(
        id="F14b",
        state=PRESENT,
        probability=0.9,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "curated promoter-switch benefit list",
        detail=(
            (row.get("evidence_summary")
             or f"Promoter switching curated as therapeutic in {ctx.gene_symbol}")
            + (f" (PMID {row.get('pmid')})" if row.get("pmid") else "")
        ),
    )


def _intron_retention_from_table(ctx: FeatureContext) -> Feature | None:
    """F15a — is there an intron whose retention could be therapeutic?"""
    row = RT.row_for("intron_retention_potential", ctx.gene_symbol)
    if not row:
        return None
    target = (row.get("target_intron") or "").strip()
    return Feature(
        id="F15a",
        state=PRESENT if target else ABSENT,
        probability=0.9 if target else 0.05,
        provenance=ANNOTATION,
        source=RT.provenance_of(row) or "intron retention potential table",
        call=target or None,
        detail=(
            f"Intron {target} annotated as a retention candidate in "
            f"{ctx.gene_symbol}"
            + (f"; mechanism {row.get('retention_mechanism')}"
               if row.get("retention_mechanism") else "")
            if target else
            f"No retention-candidate intron annotated for {ctx.gene_symbol}"
        ),
    )


def _intron_retention_benefit(ctx: FeatureContext) -> Feature | None:
    """F15b — is retaining that intron therapeutic in this gene?"""
    row = RT.row_for("intron_retention_benefit", ctx.gene_symbol)
    if not row:
        return None
    return Feature(
        id="F15b",
        state=PRESENT,
        probability=0.9,
        provenance=CONFIRMED,
        source=RT.provenance_of(row) or "curated intron-retention benefit list",
        detail=(
            (row.get("evidence_summary")
             or f"Intron retention curated as therapeutic in {ctx.gene_symbol}")
            + (f" (PMID {row.get('pmid')})" if row.get("pmid") else "")
        ),
    )


# SpliceAI's own reporting convention: 0.2 high-recall, 0.5 recommended,
# 0.8 high-precision. 0.5 is used here.
#
# DE-NOVO IN THIS CONTEXT — these cut points come from SpliceAI's variant
# interpretation guidance, not from a calibration against this platform's
# outcome. The scores are uncalibrated network outputs (M3 is open), so
# everything derived from them carries the PREDICTED tier.
SPLICEAI_SITE_THRESHOLD = 0.5
SPLICEAI_DELTA_THRESHOLD = 0.5
# Beyond this distance from an exon boundary a gained site is treated as
# deep-intronic, i.e. pseudoexon territory (F3) rather than a cryptic site
# adjacent to a real junction (F2).
DEEP_INTRONIC_NT = 100


def _alt_sequence(ctx: FeatureContext) -> str | None:
    """Build the variant sequence by substituting one base.

    Substitutions only. An indel changes the coordinate frame, and comparing
    shifted positions reports a difference that is an artefact of the shift
    rather than of splicing.
    """
    seq, off, alt = ctx.pre_mrna_sequence, ctx.variant_offset, ctx.variant_alt
    if not seq or off is None or not alt or len(alt) != 1:
        return None
    if not (0 <= off < len(seq)):
        return None
    return seq[:off] + alt.upper() + seq[off + 1:]


def _spliceai_exon_recognition(ctx: FeatureContext) -> Feature | None:
    """F1 — is the exon weakly recognised by the spliceosome?

    Scores the acceptor at the exon start and the donor at the exon end.
    Weak recognition of EITHER boundary is what makes an exon skippable or
    in need of inclusion help, so the weaker of the two governs.

    Without exon coordinates the strongest site anywhere in the sequence is
    used instead, which answers a weaker question and says so.
    """
    if not ctx.pre_mrna_sequence:
        return None
    scores = SAI.splice_site_scores(ctx.pre_mrna_sequence)
    if scores is None:
        return None

    if ctx.exon_start is not None and ctx.exon_end is not None:
        probs = SAI.predict(ctx.pre_mrna_sequence)
        if probs is None:
            return None
        n = len(probs)
        # Only score boundaries that exist. A first exon has no upstream
        # acceptor and a terminal exon has no downstream donor, so scoring
        # the missing side reads 0.00 and would mark every terminal exon
        # "weakly recognised" — a prediction about a splice site that is not
        # there. A boundary sitting at the edge of the supplied sequence is
        # treated the same way, because there is no flanking context to
        # score it against.
        edge = 1
        scored: list[tuple[str, float]] = []
        if ctx.exon_start is not None and edge <= ctx.exon_start < n:
            scored.append(("acceptor", float(probs[ctx.exon_start][1])))
        if ctx.exon_end is not None and 0 <= ctx.exon_end < n - edge:
            scored.append(("donor", float(probs[ctx.exon_end][2])))

        if not scored:
            return _unresolved(
                "F1",
                "Both exon boundaries sit at the edge of the supplied "
                "sequence, so neither can be scored. Supply a pre-mRNA "
                "window with flanking intron on the side(s) that matter.",
            )

        weakest = min(v for _, v in scored)
        weak = weakest < SPLICEAI_SITE_THRESHOLD
        shown = ", ".join(f"{name} {v:.2f}" for name, v in scored)
        omitted = (
            "" if len(scored) == 2
            else " (the other boundary is at the sequence edge — a first or "
                 "terminal exon has no site there — so it was not scored)"
        )
        detail = (
            f"SpliceAI {shown}; the weaker scored boundary governs "
            f"({weakest:.2f} vs threshold {SPLICEAI_SITE_THRESHOLD})"
            f"{omitted}"
        )
    else:
        weakest = max(scores["maxAcceptor"], scores["maxDonor"])
        weak = weakest < SPLICEAI_SITE_THRESHOLD
        detail = (
            f"No exon coordinates supplied, so the strongest site anywhere in "
            f"the {scores['length']} nt window was used: acceptor "
            f"{scores['maxAcceptor']:.2f}, donor {scores['maxDonor']:.2f}. "
            "Supply exon_start/exon_end to score the boundaries that matter."
        )

    # PRESENT means "weakly recognised" — that is what the feature asserts.
    return Feature(
        id="F1",
        state=PRESENT if weak else ABSENT,
        probability=round(1.0 - weakest, 4) if weak else round(1.0 - weakest, 4),
        provenance=PREDICTED,
        source="SpliceAI ensemble (5 models)",
        detail=detail,
    )


def _spliceai_cryptic_site(ctx: FeatureContext) -> Feature | None:
    """F2 — does the variant create a cryptic splice site?"""
    alt = _alt_sequence(ctx)
    if alt is None:
        return None
    d = SAI.delta_scores(ctx.pre_mrna_sequence, alt, window=DEEP_INTRONIC_NT)
    if d is None:
        return None
    gain = max(d["acceptorGain"], d["donorGain"])
    return Feature(
        id="F2",
        state=PRESENT if gain >= SPLICEAI_DELTA_THRESHOLD else ABSENT,
        probability=round(gain, 4) if gain > 0 else 0.02,
        provenance=PREDICTED,
        source="SpliceAI ensemble delta scores",
        detail=(
            f"Within {DEEP_INTRONIC_NT} nt of the variant: acceptor gain "
            f"{d['acceptorGain']:.2f}, donor gain {d['donorGain']:.2f} "
            f"(threshold {SPLICEAI_DELTA_THRESHOLD}). Raw SpliceAI deltas, "
            "not calibrated probabilities."
        ),
    )


def _spliceai_pseudoexon(ctx: FeatureContext) -> Feature | None:
    """F3 — does the variant activate a deep-intronic pseudoexon?

    A pseudoexon needs a gained ACCEPTOR and a gained DONOR bracketing a new
    exon, far enough from any real junction to be intronic. Requiring both
    is what separates this from F2's single cryptic site.
    """
    alt = _alt_sequence(ctx)
    if alt is None:
        return None
    full = SAI.delta_scores(ctx.pre_mrna_sequence, alt, window=None)
    if full is None:
        return None

    both = min(full["acceptorGain"], full["donorGain"])
    near = SAI.delta_scores(
        ctx.pre_mrna_sequence, alt, window=DEEP_INTRONIC_NT)
    deep = both >= SPLICEAI_DELTA_THRESHOLD and (
        near is None or max(near["acceptorGain"], near["donorGain"]) < both
    )
    return Feature(
        id="F3",
        state=PRESENT if deep else ABSENT,
        probability=round(both, 4) if both > 0 else 0.02,
        provenance=PREDICTED,
        source="SpliceAI ensemble delta scores",
        detail=(
            f"Pseudoexon activation needs a gained acceptor AND a gained "
            f"donor bracketing a new exon: acceptor gain "
            f"{full['acceptorGain']:.2f}, donor gain {full['donorGain']:.2f}, "
            f"weaker of the two {both:.2f} (threshold "
            f"{SPLICEAI_DELTA_THRESHOLD}), and the gain must sit further than "
            f"{DEEP_INTRONIC_NT} nt from the variant to be deep-intronic."
        ),
    )


def _residual_transcript(ctx: FeatureContext) -> Feature | None:
    """P2 — is there endogenous transcript left for a boosting mechanism?

    This is the load-bearing feature for the gene-activation correctness fix.
    Without it, upregulation recommends transcript-boosting unconditionally,
    including for genes where no boostable transcript exists.
    """
    if ctx.tissue_tpm is None:
        return None
    return _classify_expression(
        ctx.tissue_tpm, "tissue expression (TPM) supplied by the gene pipeline")


# ClinGen haploinsufficiency codes. NOT an ordinal 0-3 scale: 30 and 40 are
# categorical codes that happen to be written as numbers, and treating the
# field as an integer to compare would read "dosage sensitivity unlikely" as
# the strongest possible evidence. SO-DATA-02.
CLINGEN_HI_SUFFICIENT = "3"          # sufficient evidence for dosage sensitivity
CLINGEN_HI_AUTOSOMAL_RECESSIVE = "30"
CLINGEN_HI_UNLIKELY = "40"


def _dominant_negative_from_clingen(ctx: FeatureContext) -> Feature | None:
    """P6 from ClinGen Dosage Sensitivity curations.

    The logic, from the data-sources document:

      HI score 3   loss of function through DOSAGE is the established
                   mechanism, so dominant-negative is unlikely -> P6 ABSENT,
                   replacement flag permitted
      HI 0/1/2     dosage insufficiency is not established, so the mechanism
                   may be something else including dominant-negative ->
                   unresolved, flag suppressed pending curation
      30           autosomal recessive phenotype -> not a dominant-negative
                   question; unresolved
      40           dosage sensitivity unlikely -> unresolved
      not curated  unresolved, flag suppressed

    Only score 3 resolves. Everything else withholds the flag, which is the
    safe direction: a missing suggestion costs a prompt, a wrong suggestion
    to replace a protein in a dominant-negative disease is harmful.
    """
    row = RT.row_for("clingen_dosage", ctx.gene_symbol)
    if not row:
        return None
    score = (row.get("haploinsufficiency_score") or "").strip()
    if score != CLINGEN_HI_SUFFICIENT:
        return None
    return Feature(
        id="P6",
        state=ABSENT,
        probability=0.1,
        provenance=ANNOTATION,
        source=RT.provenance_of(row) or "ClinGen Dosage Sensitivity",
        call=f"HI={score}",
        detail=(
            f"ClinGen records sufficient evidence that {ctx.gene_symbol} is "
            "haploinsufficient, so loss of function through dosage is the "
            "established mechanism and a dominant-negative one is unlikely."
        ),
    )


def _dominant_negative(ctx: FeatureContext) -> Feature | None:
    """P6 — is the disease driven by a dominant-negative allele?

    Hard suppressor of the replacement flag: supplying more wild-type protein
    does not fix a dominant-negative disease, because the mutant product
    actively interferes.

    No ClinVar lookup is wired. The one inference available is from the
    user's own defect selection: haploinsufficiency and dominant-negative are
    the two standard, mutually exclusive models for a dominant disorder, so
    a user who has classified the defect as haploinsufficiency has said the
    disease is a dosage problem rather than an interference one. That is a
    user assertion, recorded as such, and it is the tier a real ClinVar
    lookup should displace.

    Any other defect leaves P6 unresolved, which withholds the flag. That is
    the safe direction: the cost of not showing a signpost is much lower than
    the cost of recommending replacement for a dominant-negative disease.
    """
    if ctx.molecular_defect != "haploinsufficiency":
        return None
    return Feature(
        id="P6",
        state=ABSENT,
        probability=0.1,
        provenance=USER_ASSERTED,
        source="user-selected molecular defect 'haploinsufficiency'",
        stand_in=True,
        detail=(
            "Defect classified as haploinsufficiency, a dosage model rather "
            "than a dominant-negative one. Not checked against ClinVar — no "
            "variant database is wired."
        ),
    )


# Localisation strings that place a protein where a circulating aptamer can
# reach it. Pegaptanib is intravitreal against a secreted factor (VEGF165);
# an aptamer has no route to a cytoplasmic or nuclear target.
_ACCESSIBLE_LOCALISATIONS = (
    "secreted",
    "extracellular",
    "cell surface",
    "cell_surface",
    "membrane",
    "plasma membrane",
)


def _extracellular_target(ctx: FeatureContext) -> Feature | None:
    """B1 — is the target protein reachable by an aptamer?

    Returns a localisation class, not a yes/no. The aptamer flag fires only
    on `secreted` or `membrane`.
    """
    loc = (ctx.protein_localisation or "").strip().lower()
    if not loc:
        return None

    if any(k in loc for k in ("secreted", "extracellular", "signal peptide")):
        call = B1_SECRETED
    elif any(k in loc for k in ("membrane", "cell surface", "cell_surface")):
        call = B1_MEMBRANE
    elif any(k in loc for k in ("cytoplasm", "nucleus", "nuclear",
                                "mitochond", "intracellular", "cytosol")):
        call = B1_INTRACELLULAR
    else:
        call = B1_UNKNOWN

    if call == B1_UNKNOWN:
        # An unrecognised localisation string is not a negative finding; it
        # is an unread one. Leave it unresolved so the flag is withheld.
        return _unresolved(
            "B1",
            f"Localisation '{ctx.protein_localisation}' could not be "
            "classified, so whether an aptamer could reach the target is "
            "unknown.",
        )

    accessible = call in B1_APTAMER_ACCESSIBLE
    return Feature(
        id="B1",
        state=PRESENT if accessible else ABSENT,
        probability=0.85 if accessible else 0.05,
        provenance=ANNOTATION,
        source="UniProt subcellular localisation",
        call=call,
        detail=(
            f"Localisation '{ctx.protein_localisation}' classified as {call} — "
            + (
                "reachable by a circulating or locally delivered aptamer"
                if accessible
                else "not reachable by an aptamer, which cannot cross into "
                     "the cytoplasm or nucleus"
            )
        ),
    )


# Repeat-expansion reference data, moved here from mechanism_service because
# it is feature evidence (F12), not ranking logic.
KNOWN_REPEAT_UNITS = {
    "CUG": "DMPK (Myotonic Dystrophy Type 1)",
    "CTG": "DMPK (Myotonic Dystrophy Type 1)",
    "CAG": "HTT / ATXN1 / ATXN2 / ATN1 (polyglutamine disorders)",
    "GGGGCC": "C9orf72 (ALS / FTD)",
    "G4C2": "C9orf72 (ALS / FTD)",
    "CCUG": "CNBP (Myotonic Dystrophy Type 2)",
    "CGG": "FMR1 / FXN (FXTAS, Fragile X)",
    "GAA": "FXN (Friedreich Ataxia)",
    "TTC": "FXN (Friedreich Ataxia)",
}

PATHOGENIC_REPEAT_THRESHOLD = 30


def _extract_repeat_count(repeat_text: str | None) -> int | None:
    """Pull the largest number out of free text like '>50 copies' or '55–200'."""
    if not repeat_text or not repeat_text.strip():
        return None
    numbers = [int(n.replace(",", "")) for n in re.findall(r"\d[\d,]*", repeat_text)]
    return max(numbers) if numbers else None


def _normalize_repeat_unit(repeat_unit: str | None) -> str | None:
    """Strip punctuation / non-nucleotide characters and uppercase the unit."""
    if not repeat_unit:
        return None
    cleaned = re.sub(r"[^ACGTUacgtu]", "", repeat_unit).upper()
    return cleaned or None


# ---------------------------------------------------------------------------
# The ladders
# ---------------------------------------------------------------------------

# Each entry is an ordered list of callables. The first to return a Feature
# wins. An empty list means the feature has no source at all and always
# resolves UNRESOLVED — the mechanisms requiring it halt.
_LADDERS: dict[str, list[Callable[[FeatureContext], Feature | None]]] = {
    # SpliceAI first; the user's defect dropdown remains as the documented
    # stand-in beneath it, for callers with no pre-mRNA sequence to give.
    "F1": [_spliceai_exon_recognition, lambda c: _from_defect(
        c, "F1", {"exon_skipping_mutation", "exon_inclusion_defect"},
        "Weak exon recognition")],
    "F2": [_spliceai_cryptic_site, lambda c: _from_defect(
        c, "F2", {"cryptic_splice_site"}, "Cryptic splice-site creation")],
    "F3": [_spliceai_pseudoexon, lambda c: _from_defect(
        c, "F3", {"pseudoexon_activation"}, "Pseudoexon activation")],
    "F4": [
        lambda c: _from_annotation(
            c, "F4", "TANGO", "Ensembl transcript biotypes / splicing complexity"),
        lambda c: _from_defect(
            c, "F4", {"poison_exon_inclusion"}, "Poison-exon presence"),
    ],
    "F5": [
        lambda c: _from_annotation(
            c, "F5", "uORF", "Ensembl 5' UTR open-reading-frame scan"),
        lambda c: _from_defect(
            c, "F5", {"uorf_mediated_repression"}, "Repressive uORF presence"),
    ],
    "F6": [
        lambda c: _from_annotation(
            c, "F6", "NAT", "Ensembl overlapping-transcript lookup"),
        lambda c: _from_defect(
            c, "F6", {"nat_mediated_repression"},
            "Overlapping antisense transcript"),
    ],
    "F7": [lambda c: _from_defect(
        c, "F7", {"mirna_mediated_repression"}, "Repressive miRNA site")],
    "F8": [lambda c: _from_defect(
        c, "F8", {"epigenetic_promoter_silencing"}, "Promoter silencing")],
    "F9": [_from_variant_text],
    # F10a / F10b are resolved together from a single fold; see
    # _resolve_accessibility. They are listed here for completeness only.
    "F10a": [],
    "F10b": [],
    # Deliberately empty — plan §6.2 and §6.4.
    "F11": [_rbp_site_from_curation],
    "F13a": [_apa_site_from_table],
    "F13b": [_apa_benefit_from_curation],
    "F14a": [_alt_promoter_from_table],
    "F14b": [_alt_promoter_benefit],
    "F15a": [_intron_retention_from_table],
    "F15b": [_intron_retention_benefit],
    "F12": [_repeat_from_catalogue, _from_repeat_text],
    # Modality-flag family. Unresolved simply withholds the flag.
    "P2": [_residual_transcript, _expression_from_table],
    # A curated dominant-negative finding (which SUPPRESSES the flag) outranks
    # a ClinGen haploinsufficiency call (which PERMITS it), which outranks the
    # user's own defect classification.
    "P6": [_dominant_negative_from_table, _dominant_negative_from_clingen,
           _dominant_negative],
    "B1": [_localisation_from_table, _extracellular_target],
}

_NO_SOURCE_REASON = {
    "F11": (
        "No curated list of validated repressive RBP sites has been "
        "populated. An atlas would answer 'an RBP binds here', not 'a "
        "repressor binds here and masking it helps', so A28 halts rather "
        "than scoring on a quantity that cannot be interpreted."
    ),
    "F13a": (
        "No polyadenylation-site table has been populated, so whether this "
        "transcript has an alternative poly(A) site is not established."
    ),
    "F14a": (
        "No alternative-promoter table has been populated, so whether this "
        "gene has one is not established."
    ),
    "F14b": (
        "Whether shifting promoter usage is therapeutically useful in this "
        "gene is a disease-specific judgement and has not been curated. A32 "
        "halts rather than guessing — alternative promoters are common, so "
        "their presence alone would fire it almost everywhere."
    ),
    "F15a": (
        "No intron-retention-potential table has been populated, so whether "
        "this gene has a retainable intron is not established."
    ),
    "F15b": (
        "Whether retaining that intron is therapeutically useful in this gene "
        "is a disease-specific judgement and has not been curated. A33 halts "
        "rather than guessing."
    ),
    "F13b": (
        "Whether shifting polyadenylation usage is therapeutically useful in "
        "this gene is a disease-specific judgement and has not been curated. "
        "A11 halts rather than guessing — an alternative site exists in most "
        "human genes, so site presence alone would fire it almost everywhere."
    ),
}


def resolve_features(ctx: FeatureContext) -> dict[str, Feature]:
    """Resolve the whole vocabulary for one target transcript."""
    out: dict[str, Feature] = {}

    f10a, f10b = _resolve_accessibility(ctx)
    out["F10a"] = f10a
    out["F10b"] = f10b

    for fid, ladder in _LADDERS.items():
        if fid in out:
            continue
        resolved: Feature | None = None
        for rung in ladder:
            resolved = rung(ctx)
            if resolved is not None:
                break
        if resolved is None:
            reason = _NO_SOURCE_REASON.get(fid)
            if reason is None:
                intended = FEATURE_CATALOG[fid]["intendedSource"]
                reason = (
                    f"Not established for this transcript. Intended source "
                    f"({intended}) is not wired, and no user input stands in "
                    f"for it here."
                )
            resolved = _unresolved(fid, reason)
        out[fid] = resolved

    return out


def unwired_features() -> list[str]:
    """Feature IDs whose intended source is not installed. For reporting."""
    return sorted(f for f, spec in FEATURE_CATALOG.items() if not spec["wired"])
