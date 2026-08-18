"use client";

import { useState } from "react";
import Sidebar from "@/components/Sidebar";
import Topbar from "@/components/Topbar";
import { BookOpen, ChevronDown, ChevronRight } from "lucide-react";

interface Mechanism {
  code: string;
  name: string;
  description: string;
  targetRegion: string;
  modality: string;
  clinicalExample?: string;
  tone: string;
}

interface TherapeuticGoalSection {
  id: string;
  title: string;
  subtitle: string;
  badge: string;
  badgeTone: string;
  mechanisms: Mechanism[];
}

const DESIGN_PARAMS = [
  {
    id: "chemistry",
    title: "Chemistry",
    description: "The oligonucleotide backbone chemistry determines the mechanism of action, nuclease resistance, and binding properties.",
    fields: [
      {
        name: "DNA Gapmer (2-10-2)",
        detail: "A chimeric oligonucleotide with a central DNA 'gap' of ~10 nucleotides flanked by 2'-modified wings (typically 2'-O-Me or LNA). The DNA gap recruits RNase H1 to cleave the target RNA, while the wings confer nuclease resistance and binding affinity. Gapmers are the most clinically validated ASO chemistry (e.g., nusinersen/Spinraza, eteplirsen).",
      },
      {
        name: "PMO (Phosphorodiamidate Morpholino)",
        detail: "A non-ionic backbone oligomer where each nucleoside is linked via phosphorodiamidate bonds to morpholine rings. PMOs do not recruit RNase H; instead they sterically block RNA interactions (splice junctions, ribosome binding, miRNA binding). Used in exon-skipping (e.g., eteplirsen) and translational arrest. Requires cell-penetrating peptide (CPP) conjugation for efficient uptake.",
      },
      {
        name: "LNA-enhanced Gapmer",
        detail: "A DNA gapmer where the flanking wings contain Locked Nucleic Acids (LNA) — bicyclic RNA analogues with a methylene bridge locking the ribose in a C3'-endo conformation. Each LNA substitution raises Tm by ~2-8°C, dramatically increasing target affinity. LNA gapmers also have enhanced nuclease resistance. Used in miravirsen (anti-miR-122) and bepetamers. Higher off-target risk due to increased potency.",
      },
      {
        name: "2'-O-Methoxyethyl (2'-OMe)",
        detail: "A ribose-modified oligonucleotide where the 2'-OH is replaced with a methoxyethyl group. 2'-OMe ASOs sterically block RNA interactions and are commonly used for splice-switching and miRNA inhibition. They have good nuclease resistance, low toxicity, and are often used in combination. Lower binding affinity than LNA but fewer off-target effects.",
      },
      {
        name: "siRNA duplex",
        detail: "A 21-23 nt duplex in which the antisense guide strand is loaded into RISC to direct AGO2-mediated cleavage of the complementary mRNA target.",
      },
    ],
  },
  {
    id: "length",
    title: "Oligo Length",
    description: "The length of the oligonucleotide in nucleotides. Default is 20 nt, with a supported range of 12-30 nt.",
    fields: [
      {
        name: "Short (12-16 nt)",
        detail: "Higher cell penetration and lower synthesis cost, but reduced target specificity and binding affinity. May be suitable for highly abundant targets.",
      },
      {
        name: "Standard (18-22 nt)",
        detail: "Optimal balance of specificity, affinity, and pharmacokinetics. Most FDA-approved ASOs fall in this range.",
      },
      {
        name: "Long (23-30 nt)",
        detail: "Higher specificity and binding affinity, but reduced cell penetration and higher synthesis cost. Used when targeting low-abundance transcripts or requiring high potency.",
      },
    ],
  },
  {
    id: "modifications",
    title: "Modifications",
    description: "Chemical modifications applied to the oligonucleotide backbone or wings to enhance stability, binding, and pharmacokinetic properties.",
    fields: [
      {
        name: "Phosphorothioate (PS) backbone",
        detail: "Replaces one non-bridging oxygen in the phosphodiester backbone with sulfur. PS linkages dramatically increase nuclease resistance and promote protein binding (e.g., to albumin), extending plasma half-life from minutes to hours. However, PS backbone can increase off-target binding to unintended RNA sequences and may activate complement pathways at high doses. Most ASO drugs incorporate full or partial PS backbones.",
      },
      {
        name: "LNA wings (5' + 3')",
        detail: "Locked Nucleic Acid (LNA) modifications placed at the 5' and 3' ends of the oligonucleotide. Each LNA substitution raises the melting temperature (Tm) by 2-8°C, increasing binding affinity to the target RNA. LNA wings also provide strong nuclease resistance at terminal positions. Typically used in gapmer wings (2-3 LNA at each end). Excessive LNA use can increase off-target activity due to hyper-stabilized binding.",
      },
      {
        name: "2'-OMe wing modifications",
        detail: "2'-O-Methyl modifications at the 5' and 3' wing positions of a gapmer. Provides moderate nuclease resistance and binding affinity improvement without the potency increase of LNA. Lower off-target risk than LNA wings; often used as a safer alternative for initial screening.",
      },
      {
        name: "PMO core",
        detail: "Phosphorodiamidate Morpholino backbone in the central region. Non-ionic, does not recruit RNase H. Used for splice-switching ASOs where RNA cleavage is not desired.",
      },
      {
        name: "PNA clamp (flanking)",
        detail: "Peptide Nucleic Acid clamps at flanking positions. PNA has an uncharged pseudopeptide backbone, providing extremely strong binding affinity and nuclease resistance. PNA clamps block nuclease access to the oligo termini, extending half-life.",
      },
    ],
  },
];

const DESIGN_TONE_MAP: Record<string, { bg: string; text: string; pill: string }> = {
  chemistry: { bg: "bg-violet-50/60", text: "text-violet-500", pill: "bg-violet-100 text-violet-700" },
  length: { bg: "bg-blue-50/60", text: "text-blue-500", pill: "bg-blue-100 text-blue-700" },
  modifications: { bg: "bg-emerald-50/60", text: "text-emerald-500", pill: "bg-emerald-100 text-emerald-700" },
};

function DesignParamSection({ param }: { param: typeof DESIGN_PARAMS[number] }) {
  const [expanded, setExpanded] = useState(false);
  const t = DESIGN_TONE_MAP[param.id];

  return (
    <section className="rounded-xl border border-[#E5E7EB] bg-white">
      <button
        onClick={() => setExpanded((e) => !e)}
        className="flex w-full items-center gap-3 px-5 py-4 text-left"
      >
        <span className={`rounded-full px-2.5 py-0.5 text-[11px] font-semibold ${t.pill}`}>
          {param.id === "tissue" ? "OPTIONAL" : param.id.toUpperCase()}
        </span>
        <h2 className="text-[15px] font-bold text-slate-900">{param.title}</h2>
        {expanded ? (
          <ChevronDown className="ml-auto h-4 w-4 text-slate-400" />
        ) : (
          <ChevronRight className="ml-auto h-4 w-4 text-slate-400" />
        )}
      </button>
      {expanded && (
        <div className="border-t border-slate-100 px-5 pb-5 pt-4">
          <p className="mb-4 text-[13px] leading-relaxed text-slate-600">{param.description}</p>
          <div className="grid gap-2">
            {param.fields.map((field) => (
              <div
                key={field.name}
                className={`rounded-lg border border-[#E5E7EB] ${t.bg} p-3`}
              >
                <p className="text-[12.5px] font-semibold text-slate-800">{field.name}</p>
                <p className="mt-1 text-[11.5px] leading-relaxed text-slate-600">{field.detail}</p>
              </div>
            ))}
          </div>
        </div>
      )}
    </section>
  );
}

const GOALS: TherapeuticGoalSection[] = [
  {
    id: "tg01",
    title: "TG01 — Gene Silencing",
    subtitle: "Reduce expression of a pathogenic gene or transcript through transcriptional or post-transcriptional mechanisms.",
    badge: "TG01",
    badgeTone: "bg-rose-100 text-rose-700",
    mechanisms: [
      {
        code: "A1",
        name: "RNase H-mediated Gapmer Knockdown",
        description:
          "DNA gap flanked by chemically modified nucleotides (2'-MOE, LNA or cEt) with a phosphorothioate backbone. The ASO binds target mRNA and recruits endogenous RNase H1 to cleave the RNA strand, catalytically degrading the transcript. Effective in both nucleus and cytoplasm on coding and non-coding RNAs.",
        targetRegion: "Mature mRNA (CDS, exons, UTRs)",
        modality: "DNA gapmer ASO (13-20 nt, PS backbone)",
        clinicalExample: "Mipomersen (Kynamro); Inotersen (Tegsedi); Tofersen (Qalsody); Olezarsen (Tryngolza)",
        tone: "rose",
      },
      {
        code: "A2",
        name: "Steric-Blocking Translation Inhibition",
        description:
          "Fully modified RNase H-inactive ASOs (PMO, PNA, 2'-OMe, 2'-MOE) bind the 5' UTR or translation initiation codon (AUG), physically blocking ribosome access or progression. Reduces protein synthesis without degrading the target RNA, preserving transcript integrity.",
        targetRegion: "5' UTR / AUG start codon",
        modality: "Steric-block ASO (15-25 nt, RNase H-inactive)",
        tone: "red",
      },
      {
        code: "A12",
        name: "microRNA Inhibition (Anti-miR / AntagomiR)",
        description:
          "Fully modified LNA, 2'-O-Me, or 2'-MOE oligonucleotides bind the mature microRNA (including seed region), preventing it from loading into RISC and silencing its downstream target mRNAs. One anti-miR can simultaneously derepress multiple disease-relevant genes regulated by a single pathogenic miRNA.",
        targetRegion: "Mature miRNA sequence (seed region)",
        modality: "LNA / 2'-O-Me / 2'-MOE anti-miR (15-23 nt)",
        clinicalExample: "Miravirsen (anti-miR-122, Hep C); Cobomarsen (anti-miR-155, oncology)",
        tone: "pink",
      },
      {
        code: "A15",
        name: "Transcriptional Gene Silencing (Promoter-Targeting ASOs)",
        description:
          "ASOs bind promoter-associated RNAs (paRNAs) or promoter-overlapping transcripts, recruiting chromatin-modifying complexes that alter histone marks and DNA methylation at the promoter, reducing transcription initiation. Regulates gene expression at the DNA level for potentially durable suppression.",
        targetRegion: "Promoter-associated RNA (paRNA) / TSS",
        modality: "Phosphorothioate gapmer ASO (~18-20 nt)",
        tone: "violet",
      },
      {
        code: "A21",
        name: "RNA Interference (siRNA-mediated Gene Silencing)",
        description:
          "Double-stranded 21-23 bp siRNA duplex is loaded into the RISC complex. The guide strand directs Argonaute-2 (AGO2) to complementary mRNA, which is cleaved and degraded. Highly potent catalytic mechanism with efficient liver delivery via GalNAc conjugation.",
        targetRegion: "Mature mRNA (CDS or UTRs)",
        modality: "siRNA duplex (21-23 bp, GalNAc-conjugated)",
        clinicalExample: "Patisiran; Givosiran; Lumasiran; Inclisiran; Vutrisiran",
        tone: "fuchsia",
      },
    ],
  },
  {
    id: "tg02",
    title: "TG02 — Gene Upregulation",
    subtitle: "Increase expression of a therapeutic gene through transcriptional activation, translational enhancement, or suppression of negative regulators.",
    badge: "TG02",
    badgeTone: "bg-indigo-100 text-indigo-700",
    mechanisms: [
      {
        code: "saRNA",
        name: "Small Activating RNA (saRNA / RNAa)",
        description:
          "21-mer dsRNA duplexes target promoter regions (−100 to −1000 bp from TSS). They load into AGO2 and recruit transcriptional co-factors and RNA Polymerase II to upregulate gene transcription at the DNA level.",
        targetRegion: "Promoter (−100 to −1000 bp from TSS)",
        modality: "21-mer dsRNA Duplex",
        clinicalExample: "MiNA Therapeutics — MTL-CEBPA for liver cancer",
        tone: "indigo",
      },
      {
        code: "uORF_block",
        name: "uORF-Targeting ASOs",
        description:
          "Steric-blocking ASOs bind to uORF start sites (uAUG) in the 5' UTR, blocking translational repressors and directing ribosomes straight to the main protein-coding region (pAUG) to increase translation.",
        targetRegion: "5' UTR (uAUG locus)",
        modality: "Steric Block (2'-MOE / PMO)",
        tone: "blue",
      },
      {
        code: "poison_exon",
        name: "NMD Suppression / Poison Exon Skipping (TANGO)",
        description:
          'Steric-blocking ASOs target splice junctions of "poison exons" (containing premature termination codons) to prevent their inclusion. The cell then produces more full-length, functional mRNA instead of NMD-degraded transcripts.',
        targetRegion: "Cryptic Intronic Splice Sites",
        modality: "Splice-modulating ASO",
        clinicalExample: "Stoke Therapeutics — Zorevunersen (STK-001) for Dravet Syndrome",
        tone: "emerald",
      },
      {
        code: "NAT_silencing",
        name: "Natural Antisense Transcript (NAT) Silencing",
        description:
          "RNase H1 gapmer ASOs degrade inhibitory long antisense RNAs (lncRNAs) transcribed from the opposite strand, derepressing the sense gene and increasing its protein production.",
        targetRegion: "Overlapping Antisense Transcript",
        modality: "RNase H1 Gapmer",
        clinicalExample: "BDNF-AS / UBE3A-ATS for Angelman syndrome",
        tone: "amber",
      },
    ],
  },
  {
    id: "tg04",
    title: "TG04 — RNA Processing Modulation",
    subtitle: "Modify RNA maturation, including splicing, polyadenylation, transcript processing, and RNA stability.",
    badge: "TG04",
    badgeTone: "bg-cyan-100 text-cyan-700",
    mechanisms: [
      {
        code: "A7",
        name: "Exon Skipping",
        description:
          "Steric-blocking ASOs (PMO or 2'-OMe PS) mask splice-regulatory sequences — splice donor/acceptor sites or exonic splicing enhancers (ESEs) — preventing spliceosome recognition of a target exon and promoting its exclusion from mature mRNA. Restores reading frame when a frameshift/nonsense mutation is in the skipped exon.",
        targetRegion: "ESEs / Splice donor / Splice acceptor",
        modality: "PMO or 2'-OMe PS steric-block ASO (20-30 nt)",
        clinicalExample: "Eteplirsen (Exondys 51); Golodirsen; Casimersen; Viltolarsen — all for DMD",
        tone: "cyan",
      },
      {
        code: "A8",
        name: "Exon Inclusion (Splice Correction)",
        description:
          "Steric-blocking ASOs bind intronic or exonic splicing silencers (ISS/ESS), preventing binding of negative splicing factors and promoting exon recognition by the spliceosome. Restores inclusion of a therapeutically important exon.",
        targetRegion: "Intronic Splicing Silencers (ISS) / ESS",
        modality: "2'-MOE PS or PMO steric-block ASO (18-25 nt)",
        clinicalExample: "Spinraza (nusinersen) — SMN2 exon 7 inclusion for SMA",
        tone: "sky",
      },
      {
        code: "A9",
        name: "Pseudoexon Suppression (Pseudoexon Skipping)",
        description:
          "ASOs mask activated pseudoexons or cryptic splice sites generated by deep intronic variants, preventing their inclusion in mature mRNA. Directly corrects deep intronic splicing defects without genome editing.",
        targetRegion: "Activated Pseudoexon / Cryptic Splice Sites",
        modality: "2'-OMe PS or PMO steric-block ASO (18-25 nt)",
        clinicalExample: "QR-1011 (ABCA4, Stargardt disease); Sepofarsen (CEP290, LCA10)",
        tone: "teal",
      },
      {
        code: "A10",
        name: "Cryptic Splice Site Blocking",
        description:
          "ASOs bind and mask pathogenic cryptic splice donor or acceptor sites created by mutations, preventing aberrant splicing and restoring use of the authentic splice site. Highly mutation-specific.",
        targetRegion: "Cryptic Splice Donor / Acceptor Site",
        modality: "2'-OMe PS or PMO steric-block ASO (18-25 nt)",
        tone: "emerald",
      },
      {
        code: "A11",
        name: "Alternative Polyadenylation (APA) Modulation",
        description:
          "ASOs block aberrant polyadenylation signals (AAUAAA and variants) or regulatory elements to redirect 3' end processing toward the desired poly(A) site. Alters transcript isoform expression, stability, and protein output while preserving genomic sequence.",
        targetRegion: "Poly(A) Signal / Cleavage Site / DSE / USE",
        modality: "Steric-blocking ASO (18-25 nt, PS / 2'-MOE / PMO)",
        tone: "slate",
      },
    ],
  },
  {
    id: "tg03",
    title: "TG03 — RNA Editing / Correction",
    subtitle: "Repair pathogenic RNA sequences or correct disease-causing RNA alterations without modifying genomic DNA.",
    badge: "TG03",
    badgeTone: "bg-purple-100 text-purple-700",
    mechanisms: [
      {
        code: "A13",
        name: "ADAR-Mediated Site-Directed RNA Editing (SDRE)",
        description:
          "Correct pathogenic RNA sequences by recruiting endogenous ADAR enzymes to catalyze adenosine-to-inosine (A→I) editing. A chemically modified antisense guide RNA hybridizes with the target transcript, forming a duplex that recruits endogenous ADAR to deaminate a specific adenosine to inosine.",
        targetRegion: "Target adenosine within pre-mRNA or mature mRNA",
        modality: "ADAR-recruiting antisense oligonucleotide (20–30 nt)",
        tone: "violet",
      },
      {
        code: "A16",
        name: "C-to-U RNA Editing (APOBEC / RESCUE)",
        description:
          "Correct pathogenic RNA sequences by recruiting or engineering cytidine deaminases (APOBEC) to catalyze cytidine-to-uridine (C→U) editing. A guide RNA targets a specific cytidine and recruits a deaminase, converting C→U and restoring protein function.",
        targetRegion: "Target cytidine within pre-mRNA or mature mRNA",
        modality: "Guide RNA (30–120 nt) + deaminase effector",
        tone: "fuchsia",
      },
      {
        code: "A17",
        name: "LEAPER (Leveraging Endogenous ADAR for Programmable Editing of RNA)",
        description:
          "Correct pathogenic RNA mutations through endogenous ADAR-mediated A→I (G) editing without introducing exogenous editing proteins. An engineered ADAR-recruiting RNA hybridizes to the target transcript, forming a dsRNA structure that recruits endogenous ADAR enzymes.",
        targetRegion: "Exonic or splice-site regions containing an editable adenosine",
        modality: "Engineered ADAR-recruiting RNA (70–200 nt)",
        tone: "purple",
      },
      {
        code: "A18",
        name: "CIRTS (CRISPR–Cas-Inspired RNA Targeting System)",
        description:
          "Recruit programmable RNA effectors to specific transcripts for RNA editing, degradation, translation control, or imaging. A guide RNA directs the CIRTS effector complex to the target RNA, where an effector domain performs site-specific modifications.",
        targetRegion: "Exonic or splice-site regions complementary to the guide RNA",
        modality: "Guide RNA + programmable effector",
        tone: "indigo",
      },
      {
        code: "A19",
        name: "REPAIR (RNA Editing for Programmable A-to-I Replacement)",
        description:
          "Correct pathogenic RNA mutations by programmable adenosine-to-inosine (A→I) editing using a catalytically inactive Cas13b fused to ADAR2. A guide RNA directs the fusion protein to the target RNA for site-specific editing.",
        targetRegion: "Exonic or splice-site regions complementary to the guide RNA",
        modality: "Guide RNA (~30–80 nt) + dCas13b–ADAR2 fusion",
        tone: "blue",
      },
      {
        code: "A20",
        name: "Spliceosome-Mediated RNA Trans-splicing (SMaRT)",
        description:
          "Replace the mutant portion of a pre-mRNA with a corrected RNA sequence through spliceosome-mediated trans-splicing. An engineered RNA trans-splicing molecule binds the target pre-mRNA and recruits the endogenous spliceosome to perform trans-splicing.",
        targetRegion: "Intronic sequences adjacent to the intended trans-splicing site",
        modality: "Engineered RNA trans-splicing molecule (>1 kb)",
        tone: "cyan",
      },
    ],
  },
  {
    id: "tg05",
    title: "TG05 — RNA Neutralization",
    subtitle: "Neutralize toxic RNA molecules or block pathogenic RNA–protein interactions without necessarily degrading the RNA.",
    badge: "TG05",
    badgeTone: "bg-slate-100 text-slate-700",
    mechanisms: [
      {
        code: "A12",
        name: "microRNA Inhibition (Anti-miR / AntagomiR)",
        description:
          "Inhibit a pathogenic microRNA to restore expression of its downstream target genes. Fully modified oligonucleotides bind the mature miRNA, preventing it from loading into RISC and silencing its downstream target mRNAs.",
        targetRegion: "Mature microRNA sequence (seed region)",
        modality: "LNA / 2'-O-Me / 2'-MOE anti-miR (15–23 nt)",
        clinicalExample: "Miravirsen (anti-miR-122, Hep C); Cobomarsen (anti-miR-155, oncology)",
        tone: "pink",
      },
      {
        code: "A14",
        name: "RNA Toxicity Neutralization (RNA Foci Disruption)",
        description:
          "Restore normal RNA-binding protein function by disrupting pathogenic RNA-protein complexes. RNase H-independent steric blockers bind the repeat tract, displace sequestered proteins like MBNL, and restore splicing homeostasis without degrading the transcript.",
        targetRegion: "Expanded repeat tract within pre-mRNA/mRNA (e.g. CUG, CCUG, GGGGCC)",
        modality: "PMO / 2'-MOE / LNA-DNA mixmer steric blocker (15–30 nt)",
        tone: "slate",
      },
      {
        code: "A25",
        name: "RNA Aptamer Therapeutics",
        description:
          "Modulate disease by specifically binding and inhibiting target molecules without altering gene expression. RNA aptamers fold into defined 3D structures that bind target proteins with high affinity and specificity, preventing ligand binding or protein function.",
        targetRegion: "Extracellular or intracellular target protein",
        modality: "RNA aptamer (20–100 nt)",
        clinicalExample: "Pegaptanib (AMD)",
        tone: "amber",
      },
    ],
  },
  {
    id: "tg06",
    title: "TG06 — Translational Regulation",
    subtitle: "Increase or decrease protein synthesis by regulating mRNA translation without altering RNA abundance.",
    badge: "TG06",
    badgeTone: "bg-yellow-100 text-yellow-700",
    mechanisms: [
      {
        code: "A2",
        name: "Steric-Blocking Translation Inhibition",
        description:
          "Reduce protein synthesis by sterically blocking ribosome initiation or progression. Fully modified RNase H-inactive ASOs bind the 5' UTR or translation initiation codon, physically blocking ribosome access without degrading the target RNA.",
        targetRegion: "5' UTR / AUG start codon",
        modality: "Steric-block ASO (15–25 nt, RNase H-inactive)",
        tone: "red",
      },
      {
        code: "A5",
        name: "Upstream Open Reading Frame (uORF) Blocking",
        description:
          "Increase translation of the primary open reading frame by preventing translation initiation at inhibitory upstream open reading frames. Steric-blocking ASOs mask the uORF start site to favor translation of the primary coding region.",
        targetRegion: "5' UTR, specifically inhibitory uORFs and their AUG initiation codons",
        modality: "RNase H-inactive steric-blocking ASO (~16–20 nt)",
        tone: "amber",
      },
      {
        code: "A6",
        name: "miRNA Binding Site Blocking (Target Protector / BlockmiR)",
        description:
          "Increase expression of a specific target gene by preventing binding of a regulatory miRNA to its binding site on the target mRNA. Blocking only that interaction restores expression of the specific gene without inhibiting the miRNA globally.",
        targetRegion: "miRNA recognition element (MRE) within the 3' UTR",
        modality: "Morpholino / LNA / 2'-O-Me steric-block ASO (15–25 nt)",
        tone: "yellow",
      },
      {
        code: "A27",
        name: "Riboswitch / RNA Structure Targeting",
        description:
          "Modulate translation by stabilizing or destabilizing structured RNA elements (IRES, G-quadruplexes, riboswitches) that control ribosomal entry or scanning. ASOs bind structured elements to alter their fold, redirecting ribosome access and restoring normal translation levels.",
        targetRegion: "Structured RNA elements: IRES, G-quadruplexes, riboswitch aptamer domains",
        modality: "PMO / 2'-O-MOE / LNA steric-block ASO (18–25 nt)",
        tone: "orange",
      },
    ],
  },
  {
    id: "tg07",
    title: "TG07 — Isoform Engineering",
    subtitle: "Generate therapeutically beneficial transcript or protein isoforms by modifying RNA processing or transcript architecture.",
    badge: "TG07",
    badgeTone: "bg-green-100 text-green-700",
    mechanisms: [
      {
        code: "A7",
        name: "Exon Skipping",
        description:
          "Restore production of a functional protein by preventing inclusion of a disease-associated exon. Steric-blocking ASOs mask splice-regulatory sequences, promoting exon exclusion from mature mRNA and restoring reading frame.",
        targetRegion: "ESEs / Splice donor / Splice acceptor",
        modality: "PMO or 2'-OMe PS steric-block ASO (20–30 nt)",
        clinicalExample: "Eteplirsen (Exondys 51); Golodirsen; Casimersen; Viltolarsen — all for DMD",
        tone: "emerald",
      },
      {
        code: "A8",
        name: "Exon Inclusion (Splice Correction)",
        description:
          "Restore production of a functional protein by promoting inclusion of a therapeutically important exon. Steric-blocking ASOs bind splice-silencing elements, preventing negative splicing factors and promoting exon recognition by the spliceosome.",
        targetRegion: "Intronic Splicing Silencers (ISS) / ESS",
        modality: "2'-MOE PS or PMO steric-block ASO (18–25 nt)",
        clinicalExample: "Spinraza (nusinersen) — SMN2 exon 7 inclusion for SMA",
        tone: "green",
      },
      {
        code: "A9",
        name: "Pseudoexon Suppression (Pseudoexon Skipping)",
        description:
          "Restore normal mRNA splicing by preventing inclusion of a pathogenic pseudoexon generated by a cryptic splice site. ASOs mask activated pseudoexons or cryptic splice sites, preventing their inclusion in mature mRNA.",
        targetRegion: "Activated Pseudoexon / Cryptic Splice Sites",
        modality: "2'-OMe PS or PMO steric-block ASO (18–25 nt)",
        clinicalExample: "QR-1011 (ABCA4, Stargardt disease); Sepofarsen (CEP290, LCA10)",
        tone: "teal",
      },
      {
        code: "A10",
        name: "Cryptic Splice Site Blocking",
        description:
          "Restore normal pre-mRNA splicing by masking a pathogenic cryptic splice donor or acceptor site and promoting use of the authentic splice site. Highly mutation-specific approach for splice correction.",
        targetRegion: "Cryptic Splice Donor / Acceptor Site",
        modality: "2'-OMe PS or PMO steric-block ASO (18–25 nt)",
        tone: "cyan",
      },
      {
        code: "A11",
        name: "Alternative Polyadenylation (APA) Modulation",
        description:
          "Modify polyadenylation site selection to restore normal transcript processing or alter mRNA stability, localization, or protein expression. ASOs block aberrant polyadenylation signals to redirect processing toward the desired poly(A) site.",
        targetRegion: "Poly(A) Signal / Cleavage Site / DSE / USE",
        modality: "Steric-blocking ASO (18–25 nt, PS / 2'-MOE / PMO)",
        tone: "slate",
      },
    ],
  },
  {
    id: "tg08",
    title: "TG08 — Protein Replacement",
    subtitle: "Restore protein function by delivering therapeutic RNA molecules that encode the missing or defective protein.",
    badge: "TG08",
    badgeTone: "bg-sky-100 text-sky-700",
    mechanisms: [
      {
        code: "A24",
        name: "Messenger RNA (mRNA) Replacement Therapy",
        description:
          "Restore functional protein expression by delivering synthetic mRNA encoding the therapeutic protein. Exogenous mRNA is delivered into the cytoplasm, translated by host ribosomes into the therapeutic protein, thereby restoring or supplementing endogenous protein function.",
        targetRegion: "N/A (delivered transcript)",
        modality: "Synthetic linear mRNA (1–5 kb, modified nucleosides)",
        clinicalExample: "Comirnaty; Spikevax",
        tone: "sky",
      },
      {
        code: "A26",
        name: "Circular RNA (circRNA)-Mediated Protein Replacement",
        description:
          "Produce sustained therapeutic protein expression using circular RNA molecules with enhanced intracellular stability. Synthetic circular RNAs are translated through cap-independent translation, with the circular topology conferring resistance to exonuclease degradation.",
        targetRegion: "N/A (delivered circRNA)",
        modality: "Circular RNA (1–5 kb, IRES-mediated translation)",
        tone: "teal",
      },
    ],
  },
  {
    id: "tg09",
    title: "TG09 — Protein Function Modulation",
    subtitle: "Modulate protein activity directly using RNA molecules that bind and regulate protein function rather than gene expression.",
    badge: "TG09",
    badgeTone: "bg-orange-100 text-orange-700",
    mechanisms: [
      {
        code: "A25",
        name: "RNA Aptamer Therapeutics",
        description:
          "RNA aptamers fold into defined three-dimensional structures that bind target proteins with high affinity and specificity, preventing ligand binding, receptor activation, or protein function. High specificity; reversible activity; relatively low immunogenicity.",
        targetRegion: "Extracellular or cell-surface target protein",
        modality: "RNA aptamer (20–100 nt)",
        clinicalExample: "Pegaptanib (AMD)",
        tone: "amber",
      },
      {
        code: "A37",
        name: "Allosteric RNA Aptamer",
        description:
          "An aptamer selected against an allosteric site rather than the active site, modulating protein activity by shifting conformational equilibrium instead of blocking the catalytic pocket directly.",
        targetRegion: "Allosteric site on the target protein",
        modality: "RNA aptamer (20–100 nt)",
        clinicalExample: "None — flagged, not scored",
        tone: "indigo",
      },
      {
        code: "A38",
        name: "Aptamer–Drug Conjugate (ApDC)",
        description:
          "An aptamer against an internalising cell-surface receptor carries a cytotoxic or oligonucleotide payload, giving cell-type-selective delivery analogous to an antibody–drug conjugate.",
        targetRegion: "Internalising cell-surface receptor",
        modality: "RNA aptamer + conjugated payload",
        clinicalExample: "None — flagged, not scored",
        tone: "emerald",
      },
      {
        code: "A39",
        name: "Bifunctional RNA Aptamer",
        description:
          "A single construct engaging two targets — for example bridging an effector cell to a tumour antigen, or simultaneously blocking two receptors on one pathway.",
        targetRegion: "Two distinct target proteins",
        modality: "Chimeric RNA aptamer",
        clinicalExample: "None — flagged, not scored",
        tone: "purple",
      },
    ],
  },
];

const TONE_MAP: Record<string, { bg: string; text: string; pill: string }> = {
  rose: { bg: "bg-rose-50/60", text: "text-rose-500", pill: "bg-rose-100 text-rose-700" },
  red: { bg: "bg-red-50/60", text: "text-red-500", pill: "bg-red-100 text-red-700" },
  pink: { bg: "bg-pink-50/60", text: "text-pink-500", pill: "bg-pink-100 text-pink-700" },
  violet: { bg: "bg-violet-50/60", text: "text-violet-500", pill: "bg-violet-100 text-violet-700" },
  fuchsia: { bg: "bg-fuchsia-50/60", text: "text-fuchsia-500", pill: "bg-fuchsia-100 text-fuchsia-700" },
  indigo: { bg: "bg-indigo-50/60", text: "text-indigo-500", pill: "bg-indigo-100 text-indigo-700" },
  blue: { bg: "bg-blue-50/60", text: "text-blue-500", pill: "bg-blue-100 text-blue-700" },
  emerald: { bg: "bg-emerald-50/60", text: "text-emerald-500", pill: "bg-emerald-100 text-emerald-700" },
  amber: { bg: "bg-amber-50/60", text: "text-amber-500", pill: "bg-amber-100 text-amber-700" },
  cyan: { bg: "bg-cyan-50/60", text: "text-cyan-500", pill: "bg-cyan-100 text-cyan-700" },
  sky: { bg: "bg-sky-50/60", text: "text-sky-500", pill: "bg-sky-100 text-sky-700" },
  teal: { bg: "bg-teal-50/60", text: "text-teal-500", pill: "bg-teal-100 text-teal-700" },
  slate: { bg: "bg-slate-50/60", text: "text-slate-500", pill: "bg-slate-100 text-slate-700" },
  orange: { bg: "bg-orange-50/60", text: "text-orange-500", pill: "bg-orange-100 text-orange-700" },
  purple: { bg: "bg-purple-50/60", text: "text-purple-500", pill: "bg-purple-100 text-purple-700" },
  yellow: { bg: "bg-yellow-50/60", text: "text-yellow-500", pill: "bg-yellow-100 text-yellow-700" },
  green: { bg: "bg-green-50/60", text: "text-green-500", pill: "bg-green-100 text-green-700" },
};

function GoalSection({ goal }: { goal: TherapeuticGoalSection }) {
  const [expanded, setExpanded] = useState(false);
  const t = TONE_MAP[goal.badgeTone.match(/bg-(\w+)-\d+/)?.[1] ?? "slate"];

  return (
    <div className="rounded-xl border border-[#E5E7EB] bg-white shadow-sm hover:shadow-md transition-all">
      <button
        onClick={() => setExpanded((e) => !e)}
        className="flex w-full items-center gap-3 px-5 py-4 text-left"
      >
        <span className={`rounded-full px-2.5 py-0.5 text-[11px] font-semibold ${goal.badgeTone}`}>
          {goal.badge}
        </span>
        <div className="flex-1 min-w-0">
          <h2 className="text-[15px] font-bold text-slate-900 truncate">{goal.title}</h2>
          <p className="text-[11px] text-slate-500 truncate">{goal.subtitle}</p>
        </div>
        <span className="text-[10px] font-medium text-slate-400">
          {goal.mechanisms.length} mechanism{goal.mechanisms.length !== 1 ? "s" : ""}
        </span>
        {expanded ? (
          <ChevronDown className="h-4 w-4 text-slate-400 shrink-0" />
        ) : (
          <ChevronRight className="h-4 w-4 text-slate-400 shrink-0" />
        )}
      </button>
      {expanded && (
        <div className="border-t border-slate-100 px-5 pb-5 pt-4">
          <p className="mb-4 text-[13px] leading-relaxed text-slate-600">{goal.subtitle}</p>
          <div className="grid gap-3">
            {goal.mechanisms.map((mech) => {
              const mt = TONE_MAP[mech.tone];
              return (
                <div
                  key={mech.code}
                  className="rounded-lg border border-[#E5E7EB] bg-slate-50/50 p-4"
                >
                  <div className="flex items-center gap-2">
                    <span className={`rounded-full px-2.5 py-0.5 text-[11px] font-semibold ${mt.pill}`}>
                      {mech.code}
                    </span>
                    <h3 className="text-[13px] font-semibold text-slate-800">{mech.name}</h3>
                  </div>
                  <p className="mt-2 text-[12px] leading-relaxed text-slate-600">
                    {mech.description}
                  </p>
                  <div className="mt-3 grid gap-2 sm:grid-cols-3">
                    <div className={`rounded-md ${mt.bg} px-2.5 py-1.5`}>
                      <p className={`text-[10px] font-medium uppercase tracking-wider ${mt.text}`}>
                        Target Region
                      </p>
                      <p className="mt-0.5 text-[11px] text-slate-700">{mech.targetRegion}</p>
                    </div>
                    <div className="rounded-md bg-white px-2.5 py-1.5 border border-slate-100">
                      <p className="text-[10px] font-medium uppercase tracking-wider text-slate-400">
                        Modality
                      </p>
                      <p className="mt-0.5 text-[11px] text-slate-700">{mech.modality}</p>
                    </div>
                    {mech.clinicalExample && (
                      <div className="rounded-md bg-emerald-50/60 px-2.5 py-1.5">
                        <p className="text-[10px] font-medium uppercase tracking-wider text-emerald-500">
                          Clinical Example
                        </p>
                        <p className="mt-0.5 text-[11px] text-slate-700">{mech.clinicalExample}</p>
                      </div>
                    )}
                  </div>
                </div>
              );
            })}
          </div>

          {/* Summary table for this goal */}
          <div className="mt-4 overflow-x-auto rounded-lg border border-[#E5E7EB]">
            <table className="w-full text-left text-[11px]">
              <thead className="bg-slate-100 text-slate-600">
                <tr>
                  <th className="px-3 py-2 font-semibold">Code</th>
                  <th className="px-3 py-2 font-semibold">Name</th>
                  <th className="px-3 py-2 font-semibold">Target Region</th>
                  <th className="px-3 py-2 font-semibold">Modality</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-100 text-slate-700">
                {goal.mechanisms.map((m) => (
                  <tr key={m.code}>
                    <td className="px-3 py-2 font-medium">{m.code}</td>
                    <td className="px-3 py-2">{m.name.split("(")[0].trim()}</td>
                    <td className="px-3 py-2">{m.targetRegion}</td>
                    <td className="px-3 py-2">{m.modality}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
}

export default function DocumentationPage() {
  return (
    <div className="flex h-screen overflow-hidden bg-[#F8FAFC]">
      <Sidebar />
      <div className="flex min-w-0 flex-1 flex-col">
        <Topbar />
        <main className="flex-1 overflow-y-auto px-6 py-5">
          <div className="mx-auto max-w-4xl space-y-5">
            <div className="flex items-center gap-3">
              <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-indigo-100">
                <BookOpen className="h-5 w-5 text-indigo-600" />
              </div>
              <div>
                <h1 className="text-[18px] font-bold text-slate-900">
                  Therapeutic Mechanism Documentation
                </h1>
                <p className="text-[12px] text-slate-500">
                  Reference documentation for RNA therapeutics mechanism design — all 9 therapeutic goals
                </p>
              </div>
            </div>

            <div className="grid gap-5">
              {GOALS.map((goal) => (
                <GoalSection key={goal.id} goal={goal} />
              ))}
            </div>

            <section className="rounded-xl border border-[#E5E7EB] bg-white">
              <div className="px-5 py-4">
                <h2 className="text-[15px] font-bold text-slate-900">Understanding Mechanism Rankings</h2>
                <p className="mt-1 text-[12.5px] text-slate-500">
                  The ranking page shows two independent signals to help distinguish general scientific support from fit for the selected design context.
                </p>
              </div>
              <div className="grid gap-3 border-t border-slate-100 px-5 py-4 sm:grid-cols-2">
                <div className="rounded-lg bg-indigo-50/60 p-3">
                  <h3 className="text-[12px] font-semibold text-indigo-800">Contextual fit score</h3>
                  <p className="mt-1 text-[11.5px] leading-relaxed text-slate-600">
                    A rule-based score calculated from the inputs selected for the current ranking, including therapeutic goal, defect type, scope, delivery context, and compatible chemistry. It orders mechanisms for this use case; it is not a percentage, probability, or clinical-confidence measure.
                  </p>
                </div>
                <div className="rounded-lg bg-emerald-50/60 p-3">
                  <h3 className="text-[12px] font-semibold text-emerald-800">Evidence rating</h3>
                  <p className="mt-1 text-[11.5px] leading-relaxed text-slate-600">
                    A rulebook label describing the level of published and clinical support for the mechanism in general. It is shown separately and is not added to the contextual fit score. A well-established mechanism can therefore have a low fit for a particular target or delivery setting.
                  </p>
                </div>
              </div>
            </section>

            {/* Tissue Expression & Delivery */}
            <section className="rounded-xl border border-[#E5E7EB] bg-white">
              <div className="px-5 py-4">
                <h2 className="text-[15px] font-bold text-slate-900">Tissue Expression &amp; Delivery</h2>
                <p className="mt-1 text-[12.5px] text-slate-500">
                  How the platform evaluates tissue-specific ASO delivery, expression data sources, and scoring modifiers.
                </p>
              </div>
              <div className="border-t border-slate-100 px-5 pb-5 pt-4 space-y-5">

                {/* Data Sources */}
                <div>
                  <h3 className="mb-2 text-[13px] font-semibold text-slate-800">Expression Data Sources</h3>
                  <div className="grid gap-2 sm:grid-cols-3">
                    <div className="rounded-lg border border-[#E5E7EB] bg-indigo-50/40 p-3">
                      <p className="text-[11px] font-bold text-indigo-700">GTEx Portal v8</p>
                      <p className="mt-1 text-[11px] leading-relaxed text-slate-600">
                        Primary source. Median gene expression (TPM) across 54 human tissues. Provides tissue-level and transcript-level expression data via API v2.
                      </p>
                    </div>
                    <div className="rounded-lg border border-[#E5E7EB] bg-blue-50/40 p-3">
                      <p className="text-[11px] font-bold text-blue-700">Human Protein Atlas</p>
                      <p className="mt-1 text-[11px] leading-relaxed text-slate-600">
                        Fallback source. RNA tissue specific nTPM, tissue specificity index, and single-cell type enrichment data from HPA JSON API.
                      </p>
                    </div>
                    <div className="rounded-lg border border-[#E5E7EB] bg-emerald-50/40 p-3">
                      <p className="text-[11px] font-bold text-emerald-700">UniProt</p>
                      <p className="mt-1 text-[11px] leading-relaxed text-slate-600">
                        Last-resort source. Tissue specificity comments from protein annotations when GTEx and HPA are unavailable.
                      </p>
                    </div>
                  </div>
                </div>

                {/* Expression Classification */}
                <div>
                  <h3 className="mb-2 text-[13px] font-semibold text-slate-800">Expression Level Classification</h3>
                  <div className="grid gap-2 sm:grid-cols-3">
                    <div className="rounded-lg bg-emerald-50/60 px-3 py-2">
                      <p className="text-[11px] font-bold text-emerald-700">High</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">&gt; 25 TPM in target tissue</p>
                    </div>
                    <div className="rounded-lg bg-amber-50/60 px-3 py-2">
                      <p className="text-[11px] font-bold text-amber-700">Medium</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">5 – 25 TPM in target tissue</p>
                    </div>
                    <div className="rounded-lg bg-rose-50/60 px-3 py-2">
                      <p className="text-[11px] font-bold text-rose-700">Low</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">&lt; 5 TPM in target tissue</p>
                    </div>
                  </div>
                </div>

                {/* Delivery Contexts */}
                <div>
                  <h3 className="mb-2 text-[13px] font-semibold text-slate-800">Delivery Contexts</h3>
                  <div className="grid gap-2 sm:grid-cols-2 lg:grid-cols-3">
                    {[
                      { key: "cns", label: "CNS / Intrathecal", desc: "Direct injection into cerebrospinal fluid. Required for brain/spinal cord targets." },
                      { key: "systemic", label: "Systemic / Subcutaneous", desc: "IV or SC administration. Distributes throughout body; liver and kidney accumulate highest doses." },
                      { key: "liver", label: "Liver-Targeted", desc: "GalNAc conjugation for hepatocyte-specific uptake via ASGPR. Most validated delivery pathway." },
                      { key: "local_intramuscular", label: "Local / Intramuscular", desc: "Direct injection into target muscle. Used for DMD exon-skipping (e.g., eteplirsen)." },
                      { key: "ocular", label: "Ocular", desc: "Intravitreal injection. Immune-privileged site with minimal systemic exposure." },
                      { key: "other", label: "Other / TBD", desc: "Not yet determined or experimental delivery route." },
                    ].map((ctx) => (
                      <div key={ctx.key} className="rounded-lg border border-[#E5E7EB] bg-slate-50/60 px-3 py-2">
                        <p className="text-[11px] font-bold text-slate-700">{ctx.label}</p>
                        <p className="mt-0.5 text-[11px] text-slate-500">{ctx.desc}</p>
                      </div>
                    ))}
                  </div>
                </div>

                {/* Chemistry-Tissue Bonuses */}
                <div>
                  <h3 className="mb-2 text-[13px] font-semibold text-slate-800">Chemistry–Tissue Interaction Bonuses</h3>
                  <div className="grid gap-2 sm:grid-cols-2">
                    <div className="rounded-lg border border-[#E5E7EB] bg-indigo-50/40 px-3 py-2">
                      <p className="text-[11px] font-bold text-indigo-700">CNS / Brain / Spinal Cord</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">PMO or LNA gapmer: <span className="font-semibold text-emerald-600">+8</span> bonus. Standard gapmer: <span className="font-semibold text-rose-500">-3</span> penalty.</p>
                    </div>
                    <div className="rounded-lg border border-[#E5E7EB] bg-amber-50/40 px-3 py-2">
                      <p className="text-[11px] font-bold text-amber-700">Liver</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">Gapmer chemistry: <span className="font-semibold text-emerald-600">+5</span> bonus for hepatic uptake.</p>
                    </div>
                    <div className="rounded-lg border border-[#E5E7EB] bg-emerald-50/40 px-3 py-2">
                      <p className="text-[11px] font-bold text-emerald-700">Eye / Retina</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">PMO chemistry: <span className="font-semibold text-emerald-600">+5</span> bonus for intravitreal delivery.</p>
                    </div>
                    <div className="rounded-lg border border-[#E5E7EB] bg-rose-50/40 px-3 py-2">
                      <p className="text-[11px] font-bold text-rose-700">Length Penalty (CNS)</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">ASO &gt; 20 nt targeting CNS/Brain/Spinal Cord: <span className="font-semibold text-rose-500">-5</span> penalty for reduced BBB permeability.</p>
                    </div>
                  </div>
                </div>

                {/* Overexpression Warnings */}
                <div>
                  <h3 className="mb-2 text-[13px] font-semibold text-slate-800">Overexpression Risk Warnings</h3>
                  <div className="grid gap-2 sm:grid-cols-2">
                    <div className="rounded-lg border border-amber-200 bg-amber-50/60 px-3 py-2">
                      <p className="text-[11px] font-bold text-amber-700">High Risk (&gt; 500 TPM)</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">High endogenous expression in target tissue — exercise caution against overexpression toxicity. Consider whether upregulation is appropriate.</p>
                    </div>
                    <div className="rounded-lg border border-yellow-200 bg-yellow-50/60 px-3 py-2">
                      <p className="text-[11px] font-bold text-yellow-700">Caution (&gt; 200 TPM)</p>
                      <p className="mt-0.5 text-[11px] text-slate-600">Moderate-high endogenous expression — monitor for potential overexpression effects during upregulation therapy.</p>
                    </div>
                  </div>
                </div>

              </div>
            </section>

            {/* ASO Design Parameters */}
            <div className="mt-8">
              <div className="mb-4">
                <h2 className="text-[16px] font-bold text-slate-900">ASO Design Parameters</h2>
                <p className="mt-1 text-[12.5px] text-slate-500">
                  Reference for all configurable fields in the Gene Silencing ASO Design form. These parameters control how candidate oligonucleotides are generated, scored, and ranked.
                </p>
              </div>
              <div className="grid gap-4">
                {DESIGN_PARAMS.map((param) => (
                  <DesignParamSection key={param.id} param={param} />
                ))}
              </div>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
}
