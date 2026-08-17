export type OrganismStatus = "live" | "curated" | "comingSoon";

export interface OrganismCapabilities {
  /** Gene lookup from Ensembl / NCBI. Available for every organism. */
  geneInfo: boolean;
  /** Therapeutic-goal + mechanism arbitration flow. */
  mechanisms: boolean;
}

export interface Organism {
  id: string; // stable key used in UI state + API calls
  commonName: string;
  scientificName: string;
  tier: 1 | 2 | 3 | 4 | 5 | 6;
  status: OrganismStatus;
  primaryUse?: string;
  /** Ensembl production species name, used by /lookup/symbol/:species/:symbol. Live tiers only. */
  ensemblSpecies?: string;
  /** NCBI/UniProt taxonomy ID, used for UniProt organism_id filtering. */
  taxonId?: number;
  /**
   * Per-organism override of the tier default. Omit to inherit
   * `defaultCapabilitiesForTier`.
   */
  capabilities?: Partial<OrganismCapabilities>;
}

/**
 * Tier defaults for the mechanism flow.
 *
 * Tier 1-3 (clinical, model and veterinary species) are on by default: they
 * resolve through Ensembl with full transcript structure, which is what the
 * mechanism arbitration needs.
 *
 * Tier 4-6 (plants, viruses, bacteria) are OPT-IN rather than blocked. The
 * rulebooks are written around mammalian RNA biology — NMD, ADAR editing,
 * RNase H1, the spliceosome — so applying them to a bacterium or a plant is
 * a judgement the user makes explicitly, not a default the product makes for
 * them. Nothing is hidden; the flow is simply not pre-enabled.
 */
export function defaultCapabilitiesForTier(tier: number): OrganismCapabilities {
  return { geneInfo: true, mechanisms: tier >= 1 && tier <= 3 };
}

/** Resolved capabilities: tier default, with any per-organism override. */
export function organismCapabilities(
  organism: Organism | undefined,
): OrganismCapabilities {
  if (!organism) return { geneInfo: false, mechanisms: false };
  return { ...defaultCapabilitiesForTier(organism.tier), ...organism.capabilities };
}

/**
 * Does this organism support the mechanism flow, accounting for an explicit
 * opt-in? `optedIn` carries the user's own choice for a Tier 4-6 organism.
 */
export function supportsMechanisms(
  organismId: string,
  optedIn = false,
): boolean {
  const org = getOrganism(organismId);
  if (!org) return false;
  if (organismCapabilities(org).mechanisms) return true;
  return optedIn && canOptIntoMechanisms(org);
}

/** Tier 4-6 may be opted into; Tier 1-3 are already on. */
export function canOptIntoMechanisms(organism: Organism | undefined): boolean {
  if (!organism) return false;
  return !organismCapabilities(organism).mechanisms;
}

/**
 * Show the advisory banner for anything outside Tier 1. Tier 2-3 run the full
 * flow but are not clinical species; Tier 4-6 need an opt-in. Tier 1 gets no
 * banner because there is nothing to caveat.
 */
export function needsOrganismBanner(organismId: string): boolean {
  const org = getOrganism(organismId);
  return org ? org.tier !== 1 : false;
}

export function organismBannerMessage(organismId: string): string | null {
  const org = getOrganism(organismId);
  if (!org || org.tier === 1) return null;
  if (org.tier <= 3) {
    return (
      `${org.commonName} is a Tier ${org.tier} species. Gene data and the ` +
      `mechanism flow are both available, but the rulebooks and the ` +
      `delivery-precedent tables are built on human and clinical-species ` +
      `evidence — read mechanism rankings as indicative for this organism.`
    );
  }
  return (
    `${org.commonName} is a Tier ${org.tier} organism. Gene information is ` +
    `available. The mechanism flow is off by default here because the ` +
    `rulebooks describe mammalian RNA biology — NMD, ADAR editing, RNase H1, ` +
    `the spliceosome — which may not apply. You can enable it explicitly if ` +
    `that is appropriate for your target.`
  );
}

export const TIER_LABELS: Record<number, { title: string; subtitle: string }> = {
  1: {
    title: "Clinical Species",
    subtitle: "~90–95% of current RNA therapeutic research",
  },
  2: {
    title: "Model Organisms",
    subtitle: "Mechanistic studies and basic biological research",
  },
  3: {
    title: "Veterinary Species",
    subtitle: "Companion and agricultural animal health",
  },
  4: {
    title: "Plants",
    subtitle: "Future expansion — plant RNA therapeutics",
  },
  5: {
    title: "Viruses",
    subtitle: "Curated reference genes — viral RNA therapeutic targets",
  },
  6: {
    title: "Bacteria",
    subtitle: "Optional — future antimicrobial RNA therapeutics",
  },
};

export const ORGANISMS: Organism[] = [
  // Tier 1 — Clinical species (default, live via Ensembl)
  { id: "human", commonName: "Human", scientificName: "Homo sapiens", tier: 1, status: "live", ensemblSpecies: "homo_sapiens", taxonId: 9606 },
  { id: "mouse", commonName: "Mouse", scientificName: "Mus musculus", tier: 1, status: "live", ensemblSpecies: "mus_musculus", taxonId: 10090 },
  { id: "rat", commonName: "Rat", scientificName: "Rattus norvegicus", tier: 1, status: "live", ensemblSpecies: "rattus_norvegicus", taxonId: 10116 },
  { id: "cynomolgus", commonName: "Cynomolgus monkey", scientificName: "Macaca fascicularis", tier: 1, status: "live", ensemblSpecies: "macaca_fascicularis", taxonId: 9541 },
  { id: "rhesus", commonName: "Rhesus macaque", scientificName: "Macaca mulatta", tier: 1, status: "live", ensemblSpecies: "macaca_mulatta", taxonId: 9544 },

  // Tier 2 — Model organisms (live via Ensembl)
  { id: "zebrafish", commonName: "Zebrafish", scientificName: "Danio rerio", tier: 2, status: "live", primaryUse: "Developmental biology", ensemblSpecies: "danio_rerio", taxonId: 7955 },
  { id: "fruitfly", commonName: "Fruit fly", scientificName: "Drosophila melanogaster", tier: 2, status: "live", primaryUse: "Genetics", ensemblSpecies: "drosophila_melanogaster", taxonId: 7227 },
  { id: "celegans", commonName: "C. elegans", scientificName: "Caenorhabditis elegans", tier: 2, status: "live", primaryUse: "RNA interference, aging, neurobiology", ensemblSpecies: "caenorhabditis_elegans", taxonId: 6239 },
  { id: "yeast", commonName: "Yeast", scientificName: "Saccharomyces cerevisiae", tier: 2, status: "live", primaryUse: "Molecular biology", ensemblSpecies: "saccharomyces_cerevisiae", taxonId: 4932 },
  { id: "fissionyeast", commonName: "Fission yeast", scientificName: "Schizosaccharomyces pombe", tier: 2, status: "live", primaryUse: "Cell cycle studies", ensemblSpecies: "schizosaccharomyces_pombe", taxonId: 4896 },

  // Tier 3 — Veterinary species (live via Ensembl)
  { id: "dog", commonName: "Dog", scientificName: "Canis lupus familiaris", tier: 3, status: "live", ensemblSpecies: "canis_lupus_familiaris", taxonId: 9615 },
  { id: "cat", commonName: "Cat", scientificName: "Felis catus", tier: 3, status: "live", ensemblSpecies: "felis_catus", taxonId: 9685 },
  { id: "pig", commonName: "Pig", scientificName: "Sus scrofa", tier: 3, status: "live", ensemblSpecies: "sus_scrofa", taxonId: 9823 },
  { id: "cow", commonName: "Cow", scientificName: "Bos taurus", tier: 3, status: "live", ensemblSpecies: "bos_taurus", taxonId: 9913 },
  { id: "horse", commonName: "Horse", scientificName: "Equus caballus", tier: 3, status: "live", ensemblSpecies: "equus_caballus", taxonId: 9796 },
  { id: "sheep", commonName: "Sheep", scientificName: "Ovis aries", tier: 3, status: "live", ensemblSpecies: "ovis_aries", taxonId: 9940 },
  { id: "goat", commonName: "Goat", scientificName: "Capra hircus", tier: 3, status: "live", ensemblSpecies: "capra_hircus", taxonId: 9925 },
  { id: "chicken", commonName: "Chicken", scientificName: "Gallus gallus", tier: 3, status: "live", ensemblSpecies: "gallus_gallus", taxonId: 9031 },

  // Tier 4 — Plants (live via Ensembl Plants API)
  { id: "arabidopsis", commonName: "Arabidopsis", scientificName: "Arabidopsis thaliana", tier: 4, status: "live", ensemblSpecies: "arabidopsis_thaliana", taxonId: 3702 },
  { id: "rice", commonName: "Rice", scientificName: "Oryza sativa", tier: 4, status: "live", ensemblSpecies: "oryza_sativa", taxonId: 39947 },
  { id: "maize", commonName: "Maize", scientificName: "Zea mays", tier: 4, status: "live", ensemblSpecies: "zea_mays", taxonId: 4577 },
  { id: "wheat", commonName: "Wheat", scientificName: "Triticum aestivum", tier: 4, status: "live", ensemblSpecies: "triticum_aestivum", taxonId: 4565 },
  { id: "tomato", commonName: "Tomato", scientificName: "Solanum lycopersicum", tier: 4, status: "live", ensemblSpecies: "solanum_lycopersicum", taxonId: 4081 },

  // Tier 5 — Viruses (curated reference gene sets, not a live connector)
  { id: "sars-cov-2", commonName: "SARS-CoV-2", scientificName: "Severe acute respiratory syndrome coronavirus 2", tier: 5, status: "curated" },
  { id: "influenza-a", commonName: "Influenza A", scientificName: "Influenza A virus", tier: 5, status: "curated" },
  { id: "hiv-1", commonName: "HIV-1", scientificName: "Human immunodeficiency virus 1", tier: 5, status: "curated" },
  { id: "hbv", commonName: "HBV", scientificName: "Hepatitis B virus", tier: 5, status: "curated" },
  { id: "hcv", commonName: "HCV", scientificName: "Hepatitis C virus", tier: 5, status: "curated" },
  { id: "rsv", commonName: "RSV", scientificName: "Respiratory syncytial virus", tier: 5, status: "curated" },

  // Tier 6 — Bacteria (live via NCBI Gene API fallback)
  { id: "ecoli", commonName: "Escherichia coli", scientificName: "Escherichia coli", tier: 6, status: "live", ensemblSpecies: "escherichia_coli", taxonId: 511145 },
  { id: "saureus", commonName: "Staphylococcus aureus", scientificName: "Staphylococcus aureus", tier: 6, status: "live", ensemblSpecies: "staphylococcus_aureus", taxonId: 1280 },
  { id: "mtuberculosis", commonName: "Mycobacterium tuberculosis", scientificName: "Mycobacterium tuberculosis", tier: 6, status: "live", ensemblSpecies: "mycobacterium_tuberculosis", taxonId: 83332 },
  { id: "paeruginosa", commonName: "Pseudomonas aeruginosa", scientificName: "Pseudomonas aeruginosa", tier: 6, status: "live", ensemblSpecies: "pseudomonas_aeruginosa", taxonId: 208964 },
];

export function getOrganism(id: string): Organism | undefined {
  return ORGANISMS.find((o) => o.id === id);
}

/**
 * Resolve an organism from whatever identifier is to hand.
 *
 * The confirmed-target object carries `organism` as the Ensembl production
 * species name ("homo_sapiens"), not the UI id ("human"), so a lookup by id
 * alone silently misses and every capability check falls open.
 */
export function resolveOrganism(
  identifier: string | undefined | null,
): Organism | undefined {
  if (!identifier) return undefined;
  const needle = identifier.trim().toLowerCase();
  return (
    ORGANISMS.find((o) => o.id.toLowerCase() === needle) ??
    ORGANISMS.find((o) => o.ensemblSpecies?.toLowerCase() === needle) ??
    ORGANISMS.find((o) => o.scientificName.toLowerCase() === needle) ??
    ORGANISMS.find((o) => o.commonName.toLowerCase() === needle)
  );
}

export function organismsByTier(tier: number): Organism[] {
  return ORGANISMS.filter((o) => o.tier === tier);
}

/** Every tier that has organisms, ascending. For grouped rendering. */
export function populatedTiers(): number[] {
  return Array.from(new Set(ORGANISMS.map((o) => o.tier))).sort((a, b) => a - b);
}

/** Organisms grouped by tier, for a grouped selector. Nothing is filtered out. */
export function organismsGroupedByTier(): {
  tier: number;
  label: { title: string; subtitle: string };
  mechanismsByDefault: boolean;
  organisms: Organism[];
}[] {
  return populatedTiers().map((tier) => ({
    tier,
    label: TIER_LABELS[tier],
    mechanismsByDefault: defaultCapabilitiesForTier(tier).mechanisms,
    organisms: organismsByTier(tier),
  }));
}

/** Returns true if the organism is Tier 1-3 (Ensembl-based, full enrichment data). */
export function isEnsemblOrganism(organismId: string): boolean {
  const org = getOrganism(organismId);
  return org ? org.tier >= 1 && org.tier <= 3 : false;
}

/** Returns true if the organism is a virus (Tier 5, curated data). */
export function isViralOrganism(organismId: string): boolean {
  const org = getOrganism(organismId);
  return org ? org.tier === 5 : false;
}

/** Returns true if the organism is a plant (Tier 4) or bacteria (Tier 6). */
export function isNonModelOrganism(organismId: string): boolean {
  const org = getOrganism(organismId);
  return org ? org.tier === 4 || org.tier === 6 : false;
}
