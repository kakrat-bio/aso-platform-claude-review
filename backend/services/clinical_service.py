"""Clinical disease details from NCBI E-utilities.

This module fetches comprehensive clinical information for genes including:
- Disease mechanisms
- Diagnostic biomarkers and tests
- Clinical symptoms and manifestations
- Therapeutic options and treatments
- Carrier manifestations (for X-linked disorders)

The service searches PubMed for clinical literature and extracts structured
information using keyword-based pattern matching and sentence extraction.
"""

import logging
import re
import xml.etree.ElementTree as ET
from typing import Optional, List

import requests

logger = logging.getLogger(__name__)

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HUMAN_FILTER = "humans[MeSH Terms] OR human[All Fields] OR patient[Title/Abstract]"

NOISE_PATTERNS = [
    r"keywords?\s+included", r"boolean\s+operators", r"combined\s+with",
    r"this\s+review\s+(aims?|seeks?|attempts?|provides?)\s+to",
    r"systematic\s+(search|review)", r"meta-analysis",
    r"we\s+(sought|aimed|attempted)\s+to",
    r"this\s+study\s+(aims?|seeks?|attempts?)",
    r"the\s+purpose\s+of\s+this", r"in\s+this\s+(review|article|study)",
    r"a\s+literature\s+review", r"to\s+our\s+knowledge",
    r"there\s+is\s+still\s+a\s+lack\s+of",
    r"the\s+genetic\s+counseling\s+capacity",
    r"with\s+the\s+advancement\s+of", r"more\s+and\s+more",
    r"current\s+treatments.*offer\s+limited",
    r"recent\s+research\s+has\s+advanced",
    r"further\s+(research|studies)", r"future\s+studies",
    r"more\s+studies\s+are\s+needed",
    r"in\s+conclusion", r"to\s+conclude",
    r"these\s+findings\s+should\s+be\s+interpreted",
    r"retrospective\s+(design|cohort|study)",
    r"the\s+advent\s+of", r"we\s+conducted\s+a\s+retrospective",
    r"we\s+further\s+explore", r"candidate\s+measures\s+across",
    r"their\s+potential\s+multidimensional",
    r"the\s+identification\s+that",
    r"designed\s+with\s+strong\s+translational",
    r"functional\s+studies\s+showed\s+that",
    r"we\s+developed\s+an?\s+\w+\s+gene\s+therapy",
    r"we\s+describe\s+\d+\s+\w+\s+patients",
    r"both\s+patients\s+lack", r"three\s+cases\s+were\s+subsidiary",
    r"due\s+to\s+their\s+rarity",
    r"evidence-based\s+treatment\s+guidelines",
    r"further\s+study\s+is\s+needed",
    r"using\s+this\s+framework", r"we\s+propose\s+possible",
    r"increased\s+genetic\s+testing\s+is\s+identifying",
    r"conflicting\s+(interpretations|results)",
    r"these\s+findings\s+support",
    r"body\s+mass\s+index", r"z-scores",
    r"a\s+better\s+understanding\s+and\s+systematic",
    r"to\s+evaluate\s+changes\s+in\s+peripheral",
    r"our\s+findings\s+therefore\s+provide",
    r"the\s+extent\s+of\s+surgical\s+resection",
    r"preliminary\s+evidence\s+suggests",
    r"single-cell\s+whole\s+genome\s+sequencing",
    r"artificial\s+intelligence\s+has\s+emerged",
    r"definite\s+AIDs",
    r"homozygosity-by-descent",
    r"analysis\s+of\s+\d+\s+fistulae",
    r"the\s+proposed\s+XNNLM",
    r"the\s+viral\s+vector\s+was\s+tested",
    r"results\s+show\s+that\s+the\s+proposed",
    # Experimental or association-study prose is not a patient symptom,
    # diagnostic test, or treatment recommendation.
    r"\b(mice|mouse|murine|rat|rodent|zebrafish|xenograft|animal\s+model|in\s+vitro|in\s+vivo|cell\s+line|cell\s+culture|primary\s+cells|transgenic|knockout|irradiation|organoid|nonhuman|monkey|macaque|rabbit|dog|pig)\b",
    r"^(here|our\s+(results|findings)|to\s+address\s+this|since\s+these|infection\s+of|in\s+fresh)",
    r"\b(odds\s+ratio|p\s*[=<]|p\(|fdr|genotyp|resequencing)\b",
    # Case-report-specific detail (individual patient, not general knowledge)
    r"\d+[- ]year[- ]old",                    # "43-year-old", "11-year-old girl"
    r"we\s+report\s+(the\s+)?(case|here)",    # "we report the case of"
    r"this\s+case\s+report",                   # "this case report"
    r"the\s+patient\s+presented\s+with",       # individual patient presentation
    r"we\s+describe\s+a\s+case",               # case descriptions
    r"retrospectively\s+analyz",               # retrospective analysis
    r"was\s+initially\s+misdiagnos",            # misdiagnosis of an individual
    r"a\s+\d+[-\s]month\s+history",            # "a 6-month history"
    r"(led|lead)\s+to\s+partial\s+regression", # specific treatment outcomes
    r"at\s+\d+\s+weeks?\s+(,|and|post|after)", # time-specific outcomes
    r"prior\s+treatment\s+for\s+\w+\s+was\s+ineffective",  # specific treatment failure
    r"herein\s+we\s+report",                   # "herein we report"
    r"we\s+present\s+(a|the)\s+case",          # "we present a case"
    r"we\s+report\s+\d+\s+(case|patient)",     # "we report 3 cases"
    r"present\s+the\s+case\s+of\s+a",           # "present the case of a"
    r"had\s+been\s+treated\s+(for|with)",       # individual treatment history
    r"was\s+diagnosed\s+as\s+having",           # individual diagnosis
]

SYMPTOM_KEYWORDS = [
    # Neuromuscular symptoms
    "weakness", "muscle weakness", "progressive", "dystrophy", "myopathy",
    "contracture", "scoliosis", "cardiomyopathy", "arrhythmia",
    "respiratory", "pulmonary", "failure", "pneumonia",
    "cognitive", "intellectual", "seizure", "epilepsy",
    "gait", "ambulation", "motor", "atrophy", "hypertrophy",
    "skeletal", "pseudohypertrophy", "fasciculation", "spasticity",
    "paresthesia", "hyporeflexia", "hyperreflexia", "myasthenia",
    "ptosis", "ophthalmoplegia", "dysphagia", "dysarthria", "dysphonia",
    # Cancer/tumor symptoms
    "tumor", "cancer", "carcinoma", "sarcoma", "leukemia", "lymphoma",
    "neoplasm", "malignancy", "metastasis", "lesion", "oncogene",
    # Neurological symptoms
    "retinopathy", "neuropathy", "nephropathy", "chorea", "dystonia",
    "tremor", "ataxia", "parkinsonism", "dementia", "psychosis",
    "depression", "anxiety", "behavioral", "insomnia", "sleep",
    # Hematological/immune symptoms
    "immunodeficiency", "hemophilia", "anemia", "thalassemia",
    "coagulopathy", "bleeding", "thrombocytopenia", "pancytopenia",
    # Metabolic symptoms
    "fibrosis", "phenotype", "malformation", "deafness", "blindness",
    "developmental delay", "growth retardation", "short stature",
    "diabetes", "obstructive", "recurrent", "intractable",
    # Age of onset
    "infantile", "childhood", "adolescent", "neonatal", "congenital",
    "juvenile", "adult-onset", "late-onset", "early-onset",
    # Organ-specific symptoms
    "osteoporosis", "fracture", "joint", "limb", "proximal", "distal",
    "bulbar", "respiratory insufficiency", "cardiac", "hepatic", "renal",
    "pancreatic", "endocrine", "pulmonary", "gastrointestinal",
    # Additional clinical manifestations
    "fatigue", "exercise intolerance", "muscle pain", "myalgia",
    "cramps", "stiffness", "rigidity", "bradykinesia", "tremor",
    "coordination", "balance", "falling", " clumsiness",
    "speech", "language", "hearing", "vision", "swallowing",
    "feeding", "failure to thrive", "weight loss", "failure to thrive",
]

DIAGNOSTIC_KEYWORDS = [
    # Genetic testing methods
    "creatine kinase", "CK", "MLPA", "NGS", "next-generation sequencing",
    "whole exome", "WES", "whole genome", "WGS", "Sanger",
    "PCR", "multiplex ligation", "genetic testing", "sequencing",
    "deletion", "duplication", "carrier testing", "prenatal testing",
    "newborn screening", "genotyping", "karyotype", "FISH",
    "chromosomal microarray", "targeted gene panel", "diagnostic gene panel",
    # Tissue/biopsy methods
    "muscle biopsy", "immunohistochemistry", "western blot",
    "immunofluorescence", "histopathology", "electron microscopy",
    "liver biopsy", "skin biopsy", "nerve biopsy", "bone marrow biopsy",
    # Electrophysiology
    "EMG", "electromyography", "nerve conduction", "EEG", "electroencephalogram",
    "EKG", "ECG", "electrocardiogram", "echocardiography", "echocardiogram",
    "evoked potentials", "nerve conduction study",
    # Imaging
    "MRI", "magnetic resonance imaging", "CT scan", "computed tomography",
    "ultrasound", "X-ray", "radiography", "PET scan", "bone scan",
    "DEXA scan", "bone densitometry",
    # Laboratory tests
    "biomarker", "sweat test", "mass spectrometry", "liquid chromatography",
    "immunoassay", "ELISA", "enzyme assay", "metabolic screen",
    "amino acid analysis", "acylcarnitine profile", "organic acids",
    "lactate", "pyruvate", "ammonia", "uric acid",
    # Diagnosis terms
    "diagnosis", "diagnostic", "confirmed", "detect", "laboratory",
    "assay", "elevated", "increased", "plasma", "serum", "blood",
    "urine", "tissue", "genetic counseling", "prenatal diagnosis",
    "preimplantation diagnosis", "predictive testing", "presymptomatic",
    # Specific disease diagnostics
    "alpha-galactosidase", "hexosaminidase", "phenylalanine",
    "tyrosine", "galactose", "biotinidase", "17-hydroxyprogesterone",
    "immunoreactive trypsinogen", "TSH", "hemoglobin",
]

THERAPY_KEYWORDS = [
    # Pharmacological treatments
    "corticosteroid", "prednisone", "deflazacort", "steroid",
    "anti-inflammatory", "NSAID", "analgesic", "pain management",
    "anticonvulsant", "antiepileptic", "sedative", "anxiolytic",
    "antidepressant", "antipsychotic", "mood stabilizer",
    # Gene-based therapies
    "exon skipping", "eteplirsen", "golodirsen", "viltolarsen",
    "risdiplam", "nusinersen", "gene therapy", "AAV",
    "delandistrogene", "olmesartan", "losartan", "ACE inhibitor",
    "siRNA", "ASO", "antisense", "CRISPR", "gene editing",
    "base editing", "prime editing", "lentiviral", "retroviral",
    # Enzyme/protein replacement
    "enzyme replacement", "substrate reduction", "chaperone",
    "protein replacement", "metabolic cofactor", "vitamin",
    "coenzyme", "pharmacological chaperone",
    # Device/surgical interventions
    "pacemaker", "defibrillator", "ICD", "respiratory support",
    "ventilation", "BiPAP", "CPAP", "surgery", "orthosis",
    "bracing", "wheelchair", "assistive device", "prosthesis",
    "deep brain stimulation", "pallidotomy", "DBS",
    # Cell-based therapies
    "stem cell", "bone marrow transplant", "transplant",
    "hematopoietic stem cell", "mesenchymal stem cell",
    "autologous", "allogeneic", "gene corrected",
    # Immunotherapies
    "antibody", "monoclonal", "immunotherapy", "checkpoint inhibitor",
    "CAR-T", "adoptive cell transfer", "vaccine",
    # Supportive care
    "physiotherapy", "rehabilitation", "physical therapy",
    "occupational therapy", "speech therapy", "respiratory therapy",
    "nutritional support", "dietary management", "feeding supplement",
    "palliative care", "hospice", "pain management",
    # Clinical trial terms
    "treatment", "therapy", "therapeutic", "management",
    "clinical trial", "Phase", "FDA", "approved", "compassionate use",
    "expanded access", "investigational", "off-label",
    # Specific drug classes
    "VMAT2", "tetrabenazine", "deutetrabenazine", "valbenazine",
    "cholinesterase inhibitor", "dopamine", "serotonin",
    "GABA", "glutamate", "NMDA", "calcium channel",
    "potassium channel", "sodium channel",
    # Disease-specific treatments
    "chemotherapy", "radiation", "targeted therapy", "hormone therapy",
    "bisphosphonate", "growth hormone", "insulin", "metformin",
    "anticoagulant", "antiplatelet", "fibrinolytic",
]

SYMPTOM_DISQUALIFIERS = [
    r"\bclinical trial\b", r"\bphase\s+[ivx]+\b", r"\btrial(s)?\b",
    r"\btherapy\b", r"\btreatment\b", r"\bmanagement\b",
    r"\bmonitoring\b", r"\bdevice-based\b", r"\bgenetic testing\b",
    r"\bmolecular diagnosis\b", r"\bnext[- ]generation sequencing\b",
    r"\bNGS\b", r"\bMLPA\b", r"\bbiomarker\b", r"\bdiagnos(tic|is(es)?)\b",
    r"\bCRISPR\b", r"\bgene therapy\b", r"\bAAV\b", r"\bartificial intelligence\b",
    r"\bAI\b", r"\brehabilitation\b", r"\bassistive device\b",
    r"\bclinical trial\b", r"\bstudy\b", r"\bcohort\b",
    r"\bsurvey\b", r"\bquestionnaire\b", r"\brespondent(s)?\b",
    r"\bparticipants?\b", r"\bcaregiver(s)?\b", r"\bincidence\b",
    r"\bprevalence\b", r"\bawareness\b", r"\bimprove[d]?\s+medical\s+care\b",
]

DIAGNOSTIC_DISQUALIFIERS = [
    r"\bclinical trial\b", r"\bphase\s+[ivx]+\b",
    r"\bsurgery\b", r"\brehabilitation\b",
    r"\bpatient\s+(with|presenting|had|was|received)\b",
    r"\bcourse\s+of\s+treatment\b", r"\bfollow[- ]?up\b",
    r"\bin\s+this\s+study\b", r"\bthis\s+research\b",
    r"\bquestionnaire\b", r"\bparticipants?\b",
    r"\breported\s+(herein|previously)\b",
    r"\breported\s+elsewhere\b", r"\bof\s+\d+\s+patients\b",
    r"\bsix\s+to\s+eight\b", r"\boften\s+requires\b",
    r"\bhigh\s+index\s+of\s+suspicion\b",
    r"\bawareness\s+is\s+needed\b",
]

THERAPY_DISQUALIFIERS = [
    r"\bdifferential\s+diagnos", r"\bcase\s+report\b",
    r"\bmisdiagnos", r"\bclinical\s+distinction",
    r"\bearly\s+diagnosis\b", r"\bdiagnostic\s+accuracy\b",
    r"\bthree[- ]tier\b", r"\bdiagnostic\s+system\b",
    r"\breduce\s+long[- ]term\s+malignant\s+risk\b",
    r"\bretrospective\b", r"\bwas\s+ineffective\b",
    r"\bhad\s+no\s+effect\b",
    r"\bwe\s+present\b", r"\bwe\s+discuss\b",
    r"\bfragile\s+health\s+systems?\b",
    r"\bhigh\s+solar\s+gradients?\b",
    r"\bpreventive\s+measures\b",
    r"\bsub[- ]?Saharan\s+Africa\b",
    r"\bclose\s+communication\s+between\b",
    r"\boptimal\s+management\b",
    r"\biatrogenic\s+immunosuppressive\s+therapy\b",
    r"\bquality\s+of\s+life\b", r"\bpsychosocial\b",
]

SYMPTOM_CONTEXT_PATTERNS = [
    r"characterized by",
    r"present(?:s|ed|ing)? with",
    r"manifest(?:s|ed|ation)",
    r"symptom(?:s)? (?:include|are|is|such as)",
    r"signs? of",
    r"clinical feature(?:s)?",
    r"clinical presentation",
    r"reported (?:symptoms|features|manifestations)",
    r"patients? (?:with|presenting with)",
    r"characteristic(?:al)? of",
]

MECHANISM_KEYWORDS = [
    # Mutation types
    "mutation", "deletion", "duplication", "frameshift",
    "nonsense", "missense", "splicing", "reading frame",
    "point mutation", "insertion", "inversion", "translocation",
    "copy number variation", "CNV", "structural variant",
    # Functional consequences
    "loss of function", "gain of function", "pathogenic variant",
    "truncating", "premature termination", "premature stop",
    "protein deficiency", "enzyme deficiency", "protein misfolding",
    "haploinsufficiency", "dominant negative", "dominant-negative",
    "null allele", "amorphic", "hypomorphic", "hypermorphic",
    # Repeat disorders
    "repeat expansion", "trinucleotide repeat", "CAG repeat",
    "CTG repeat", "CGG repeat", "GCC repeat", "GAA repeat",
    "microsatellite instability", "dynamic mutation",
    # Expression defects
    "loss of expression", "reduced expression", "overexpression",
    "aberrant expression", "ectopic expression", "silencing",
    "promoter mutation", "enhancer mutation", "epigenetic",
    # Protein-level mechanisms
    "protein aggregation", "protein misfolding", "protein degradation",
    "protein trafficking", "protein localization", "protein interaction",
    "post-translational modification", "phosphorylation", "glycosylation",
    "acetylation", "methylation", "ubiquitination",
    # Cellular mechanisms
    "oxidative stress", "mitochondrial dysfunction", "calcium dysregulation",
    "apoptosis", "necrosis", "autophagy", "inflammation",
    "fibrosis", "fibrotic", "remodeling", "hypertrophy",
    "atrophy", "degeneration", "necrosis", "excitotoxicity",
    # Pathway disruption
    "signal transduction", "pathway disruption", "dysregulated",
    "aberrant signaling", "constitutive activation", "inhibited",
    "downregulated", "upregulated", "impaired", "defective",
    # Specific disease mechanisms
    "dystrophin deficiency", "sarcoglycan", "dystroglycan",
    "muscular dystrophy", "myotonic", "channelopathy",
    "storage disorder", "lysosomal", "peroxisomal", "mitochondrial",
    "metabolic", "neurodegenerative", "neurodevelopmental",
]


def _is_noise(sentence: str) -> bool:
    lower = sentence.lower()
    for pat in NOISE_PATTERNS:
        if re.search(pat, lower):
            return True
    return False


def _is_experimental_sentence(sentence: str) -> bool:
    lower = sentence.lower()
    return bool(re.search(
        r"\b(mice|mouse|murine|rat|rodent|zebrafish|xenograft|animal\s+model|in\s+vitro|in\s+vivo|cell\s+line|cell\s+culture|primary\s+cells|transgenic|knockout|irradiation|organoid|nonhuman|monkey|macaque|rabbit|dog|pig)\b",
        lower,
    ))


def _ncbi_search(query: str, db: str = "pubmed", retmax: int = 5) -> list:
    try:
        resp = requests.get(
            f"{NCBI_EUTILS}/esearch.fcgi",
            params={"db": db, "term": query, "retmax": retmax, "retmode": "json"},
            timeout=10,
        )
        resp.raise_for_status()
        return resp.json().get("esearchresult", {}).get("idlist", [])
    except Exception:
        return []


def _ncbi_fetch_abstracts(pmids: list) -> str:
    if not pmids:
        return ""
    try:
        resp = requests.get(
            f"{NCBI_EUTILS}/efetch.fcgi",
            params={
                "db": "pubmed",
                "id": ",".join(pmids),
                "rettype": "abstract",
                "retmode": "xml",
            },
            timeout=12,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        texts = []
        for article in root.findall(".//Article"):
            for abstract_text in article.findall(".//AbstractText"):
                # Abstract sections can contain nested XML tags (for example,
                # gene names or italicized terms).  ``.text`` stops at the
                # first nested tag and produces the incomplete fragments that
                # were appearing in the Disease Association card.
                texts.append("".join(abstract_text.itertext()).strip())
        return " ".join(texts)
    except Exception:
        return ""


def _extract_matching_terms(
    text: str,
    keywords: list,
    max_items: int = 5,
    disqualifiers: Optional[list] = None,
    require_context: bool = False,
    disease_terms: Optional[list[str]] = None,
) -> list:
    """Extract relevant sentences from text that contain keyword matches.

    Returns cleaned-up sentences or short phrases that are medically relevant.
    """
    matched = []
    seen = set()
    sentences = re.split(r'(?<=[.!?])\s+', text)

    for sentence in sentences:
        if len(matched) >= max_items:
            break
        cleaned = sentence.strip()
        if not cleaned or len(cleaned) < 15:
            continue
        if cleaned in seen:
            continue

        lower = cleaned.lower()

        # A PubMed search can still return a paper that merely mentions a
        # related gene. Do not present it as a disease manifestation unless
        # the sentence itself names the selected disease/phenotype.
        if disease_terms and not any(term.lower() in lower for term in disease_terms):
            continue

        # Check if any keyword matches
        has_match = any(kw.lower() in lower for kw in keywords)
        if not has_match:
            continue

        # Skip noisy sentences or sentences that are clearly experimental.
        if _is_noise(cleaned) or _is_experimental_sentence(cleaned):
            continue

        if disqualifiers and any(re.search(pat, lower) for pat in disqualifiers):
            continue

        if require_context and not any(re.search(pat, lower) for pat in SYMPTOM_CONTEXT_PATTERNS):
            continue
        # Skip pure statistics
        if re.match(r'^\d+[\s%(]', cleaned):
            continue
        # Skip very short fragments
        if len(cleaned) < 20:
            continue
        # Skip sentences that are mostly numbers
        if sum(c.isdigit() for c in cleaned) > len(cleaned) * 0.25:
            continue

        # Keep the complete sentence. The Disease Association card has its own
        # scroll area, so shortening evidence text here loses useful context.
        cleaned = cleaned.strip(".,;: ")

        seen.add(cleaned)
        matched.append(cleaned)

    return matched


def _extract_symptom_terms(
    text: str, max_items: int = 10, disease_terms: Optional[list[str]] = None
) -> list:
    """Extract specific clinical symptom terms from medical text.

    Looks for the SYMPTOM_KEYWORDS in text and returns concise descriptions.
    """
    return _extract_matching_terms(
        text,
        SYMPTOM_KEYWORDS,
        max_items,
        disqualifiers=SYMPTOM_DISQUALIFIERS,
        require_context=True,
        disease_terms=disease_terms,
    )


def _extract_diagnostic_terms(
    text: str, max_items: int = 8, disease_terms: Optional[list[str]] = None
) -> list:
    """Extract specific diagnostic test/biomarker names from medical text.

    Looks for the DIAGNOSTIC_KEYWORDS in text and returns concise descriptions.
    """
    return _extract_matching_terms(
        text, DIAGNOSTIC_KEYWORDS, max_items,
        disqualifiers=DIAGNOSTIC_DISQUALIFIERS,
        disease_terms=disease_terms,
    )


def _extract_therapy_terms(
    text: str, max_items: int = 8, disease_terms: Optional[list[str]] = None
) -> list:
    """Extract specific therapy/drug names from medical text.

    Looks for the THERAPY_KEYWORDS in text and returns concise descriptions.
    """
    return _extract_matching_terms(
        text, THERAPY_KEYWORDS, max_items,
        disqualifiers=THERAPY_DISQUALIFIERS,
        disease_terms=disease_terms,
    )


def _gene_linked_pmids(gene_id: Optional[str], retmax: int = 15) -> list:
    """PubMed ids NCBI has curated as being about this gene record.

    Uses elink from db=gene to db=pubmed. Unlike a free-text symbol search
    this cannot return a homonym: the link is to the gene record, not to the
    letters of its symbol.
    """
    if not gene_id:
        return []
    try:
        resp = requests.get(
            f"{NCBI_EUTILS}/elink.fcgi",
            params={
                "dbfrom": "gene",
                "db": "pubmed",
                "id": str(gene_id),
                "retmode": "json",
            },
            timeout=10,
        )
        resp.raise_for_status()
        linksets = resp.json().get("linksets", [])
        for ls in linksets:
            for db in ls.get("linksetdbs", []):
                if db.get("dbto") == "pubmed":
                    return [str(x) for x in db.get("links", [])][:retmax]
    except Exception as exc:
        logger.warning("Gene->PubMed elink failed for %s: %s", gene_id, exc)
    return []


def _get_gene_summary(gene_symbol: str) -> dict:
    result = {"summary": None, "omim_id": None, "description": None,
              "gene_id": None, "full_name": None}
    try:
        resp = requests.get(
            f"{NCBI_EUTILS}/esearch.fcgi",
            params={
                "db": "gene",
                "term": f"{gene_symbol}[Symbol] AND human[Organism]",
                "retmax": 5,
                "retmode": "json",
            },
            timeout=10,
        )
        resp.raise_for_status()
        gene_ids = resp.json().get("esearchresult", {}).get("idlist", [])
        if not gene_ids:
            return result

        # Fetch summaries for all candidates, pick exact symbol match
        for gid in gene_ids:
            resp = requests.get(
                f"{NCBI_EUTILS}/esummary.fcgi",
                params={"db": "gene", "id": gid, "retmode": "json"},
                timeout=10,
            )
            resp.raise_for_status()
            data = resp.json().get("result", {}).get(gid, {})
            name = (data.get("name") or "").upper().strip()
            if name != gene_symbol.upper():
                continue
            result["summary"] = data.get("summary")
            result["description"] = data.get("description")
            result["gene_id"] = str(gid)
            result["full_name"] = data.get("description")
            mim = data.get("mim")
            if mim and isinstance(mim, list) and mim:
                result["omim_id"] = str(mim[0])
            elif mim and isinstance(mim, str) and mim.isdigit():
                result["omim_id"] = mim
            break
    except Exception:
        pass
    return result


def _extract_mechanism_from_summary(summary: str) -> Optional[str]:
    """Extract disease mechanism from gene summary.

    Combines multiple relevant sentences to give a complete picture
    of the disease mechanism, without artificial truncation.

    Args:
        summary: NCBI gene summary text

    Returns:
        Extracted mechanism text or None if not found
    """
    if not summary:
        return None
    sentences = re.split(r'(?<=[.!?])\s+', summary)

    # Disease-related keywords
    disease_kw = [
        "disease", "disorder", "syndrome", "cancer", "tumor",
        "dystrophy", "myopathy", "deficiency", "insufficiency",
        "failure", "malformation", "abnormality", "pathology",
        "degeneration", "neurodegenerative", "disability", "condition",
        "lethal", "fatal", "severe", "progressive", "chronic",
    ]

    # Mechanism-related keywords (causal relationships)
    mechanism_kw = [
        "cause", "causes", "caused", "results in", "leads to",
        "responsible for", "mutations", "deletion", "duplication",
        "deficiency", "loss of function", "trinucleotide",
        "CAG repeat", "expansion", "frameshift", "nonsense",
        "missense", "splicing", "reading frame", "premature",
        "haploinsufficiency", "dominant negative", "pathogenic",
    ]

    # Collect ALL sentences that match disease+mechanism criteria
    mechanism_sentences = []
    for sentence in sentences:
        lower = sentence.lower()
        has_disease = any(dk in lower for dk in disease_kw)
        has_mechanism = any(mk in lower for mk in mechanism_kw)
        if has_disease and has_mechanism:
            cleaned = sentence.strip()
            if len(cleaned) > 20:
                mechanism_sentences.append(cleaned)

    # If we found disease+mechanism sentences, combine them
    if mechanism_sentences:
        combined = " ".join(mechanism_sentences)
        return combined if len(combined) <= 800 else combined[:797] + "..."

    # A generic gene-function sentence (for example, one describing
    # alternative splicing) is not a disease mechanism. Only return text when
    # the source explicitly ties a pathogenic mechanism to a disease.
    return None


def _fetch_omim_mechanism(omim_id: str) -> Optional[str]:
    """Fetch detailed disease mechanism from OMIM.

    OMIM provides the most authoritative disease mechanism descriptions.
    This fetches the OMIM entry page and extracts the mechanism text.
    """
    if not omim_id:
        return None

    # Clean up OMIM ID (remove # prefix)
    clean_id = omim_id.strip().lstrip("#")
    if not clean_id.isdigit():
        return None

    try:
        # Fetch OMIM entry via NCBI EFetch (OMIM is indexed in NCBI)
        resp = requests.get(
            f"{NCBI_EUTILS}/esearch.fcgi",
            params={
                "db": "omim",
                "term": f"{clean_id}[OMIM ID]",
                "retmode": "json",
            },
            timeout=10,
        )
        if resp.status_code != 200:
            return None

        ids = (resp.json().get("esearchresult") or {}).get("idlist", [])
        if not ids:
            return None

        # Fetch the OMIM entry summary
        resp = requests.get(
            f"{NCBI_EUTILS}/esummary.fcgi",
            params={"db": "omim", "id": ids[0], "retmode": "json"},
            timeout=10,
        )
        if resp.status_code != 200:
            return None

        data = resp.json().get("result", {}).get(ids[0], {})
        if not data:
            return None

        # Extract the clinical description / mechanism from OMIM
        # OMIM entries have a "titles" field and "text" field
        titles = data.get("titles", {})
        title = titles.get("preferredTitle", "") or titles.get("approvedTitle", "")

        # Also check for clinical features
        text = data.get("text", "")

        # Build mechanism from available OMIM data
        mechanism_parts = []
        if title:
            mechanism_parts.append(f"OMIM: {title}")
        if text:
            # Extract the first meaningful paragraph
            paragraphs = text.split("\n\n")
            for para in paragraphs:
                para = para.strip()
                if len(para) > 30 and any(kw in para.lower() for kw in [
                    "cause", "caused", "mutation", "deletion", "deficiency",
                    "loss of function", "pathogenic", "disease", "disorder",
                    "syndrome", "mechanism", "results in", "leads to",
                ]):
                    # Clean up
                    para = para.strip(".,;: ")
                    if len(para) > 500:
                        para = para[:497] + "..."
                    mechanism_parts.append(para)
                    break

        if mechanism_parts:
            return ". ".join(mechanism_parts)

    except Exception as e:
        logger.info(f"OMIM mechanism fetch failed for {omim_id}: {e}")

    return None


def _get_clinical_trials(gene_symbol: str, disease_name: Optional[str] = None) -> List[str]:
    """Fetch clinical trials information from ClinicalTrials.gov API.
    
    Args:
        gene_symbol: Official gene symbol
        disease_name: Associated disease name if known
        
    Returns:
        List of clinical trial-related therapeutic terms
    """
    trials = []
    try:
        search_term = gene_symbol
        if disease_name:
            search_term = f"{gene_symbol} AND {disease_name}"
        
        url = "https://clinicaltrials.gov/api/v2/studies"
        params = {
            "query.term": search_term,
            "pageSize": 15,
            "format": "json"
        }
        
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        
        studies = data.get("studies", [])
        for study in studies[:10]:
            try:
                protocol = study.get("protocolSection", {})
                ident = protocol.get("identificationModule", {})
                title = ident.get("briefTitle", "")
                
                # Extract intervention names
                interventions = protocol.get("armsInterventionsModule", {})
                for intervention in interventions.get("interventions", []):
                    name = intervention.get("name", "")
                    if name and len(name) > 5:
                        trials.append(f"Clinical trial: {name}")
                
                # Extract status/phase if relevant
                status_mod = protocol.get("statusModule", {})
                phase = status_mod.get("phases", [])
                if phase:
                    phase_str = ", ".join(phase)
                    trial_phase = f"Clinical trials in Phase {phase_str}"
                    if trial_phase not in trials:
                        trials.append(trial_phase)
            except Exception:
                continue
                    
    except Exception:
        pass
    
    return trials[:8]


def get_clinical_details(
    gene_symbol: str,
    disease_name: Optional[str] = None,
    omim_id: Optional[str] = None,
    phenotypes: Optional[List[str]] = None,
) -> dict:
    """Fetch comprehensive clinical details for a gene.

    This function retrieves disease mechanism, diagnostic tests, clinical symptoms,
    therapeutic options, and carrier manifestations for a given gene by searching
    PubMed clinical literature.

    Args:
        gene_symbol: Official gene symbol (e.g., "DMD", "TP53")
        disease_name: Associated disease name if known
        omim_id: OMIM identifier if available
        phenotypes: List of associated phenotype terms

    Returns:
        Dictionary containing:
        - diseaseMechanism: Full text describing the disease mechanism
        - diagnosticTests: List of specific diagnostic tests/biomarkers
        - clinicalSymptoms: List of specific clinical symptoms
        - carrierManifestations: List of carrier-related information
        - therapeuticOptions: List of specific treatment options
    """
    result = {
        "diseaseMechanism": None,
        "diagnosticTests": [],
        "clinicalSymptoms": [],
        "carrierManifestations": [],
        "therapeuticOptions": [],
    }

    gene_info = _get_gene_summary(gene_symbol)
    if not omim_id and gene_info["omim_id"]:
        omim_id = gene_info["omim_id"]

    if gene_info["summary"] and not result["diseaseMechanism"]:
        mech = _extract_mechanism_from_summary(gene_info["summary"])
        if mech:
            result["diseaseMechanism"] = mech

    # Also try the NCBI gene description (shorter but often more direct)
    if not result["diseaseMechanism"] and gene_info.get("description"):
        desc = gene_info["description"]
        if any(kw in desc.lower() for kw in ["cause", "associated", "linked", "disorder", "disease", "syndrome"]):
            result["diseaseMechanism"] = desc

    # Fetch detailed mechanism from OMIM if available
    if omim_id and not result["diseaseMechanism"]:
        omim_mech = _fetch_omim_mechanism(omim_id)
        if omim_mech:
            result["diseaseMechanism"] = omim_mech

    # Only use disease-specific literature. Searching by gene symbol alone
    # pulls in unrelated cell and animal experiments, which must not be shown
    # as human clinical symptoms or diagnostics.
    disease_terms = []
    if phenotypes:
        disease_terms.extend(phenotypes)
    if disease_name:
        disease_terms.extend(part.strip() for part in disease_name.split(";") if part.strip())
    disease_terms = list(dict.fromkeys(term for term in disease_terms if term))[:3]

    for disease_term in disease_terms:
        pmids = _ncbi_search(
            f'("{disease_term}"[Title/Abstract]) AND '
            "(diagnosis[Title/Abstract] OR clinical[Title/Abstract] OR "
            "treatment[Title/Abstract] OR management[Title/Abstract] OR phenotype[Title/Abstract]) AND "
            f"({HUMAN_FILTER})",
            retmax=20,
        )
        abstract_text = _ncbi_fetch_abstracts(pmids)
        if not abstract_text:
            continue
        if not result["clinicalSymptoms"]:
            result["clinicalSymptoms"] = _extract_symptom_terms(abstract_text, max_items=7)
        if not result["diagnosticTests"]:
            result["diagnosticTests"] = _extract_diagnostic_terms(abstract_text, max_items=8)
        if not result["therapeuticOptions"]:
            result["therapeuticOptions"] = _extract_therapy_terms(abstract_text, max_items=8)
        if result["clinicalSymptoms"] and result["diagnosticTests"] and result["therapeuticOptions"]:
            break

    # Fallback when the disease-term search came up short.
    #
    # This used to be `"{gene_symbol}"[Title/Abstract] AND (clinical terms)`,
    # with the comment "ensures every gene gets some clinical data". A bare
    # symbol is ambiguous in free text, and the search happily returned it:
    # HTT matched "hyalinizing trabecular tumours (HTTs) of the thyroid", so
    # huntingtin's clinical panel filled with thyroid-neoplasm diagnostics and
    # therapies. Filling a gap with a homonym is worse than leaving it empty.
    #
    # NCBI's curated Gene -> PubMed link is the disambiguated route: it returns
    # papers indexed against THIS gene record, so no amount of abbreviation
    # collision can leak in.
    if not result["clinicalSymptoms"] or not result["diagnosticTests"]:
        pmids = _gene_linked_pmids(gene_info.get("gene_id"), retmax=15)
        abstract_text = _ncbi_fetch_abstracts(pmids)
        if abstract_text:
            if not result["clinicalSymptoms"]:
                result["clinicalSymptoms"] = _extract_symptom_terms(abstract_text, max_items=7)
            if not result["diagnosticTests"]:
                result["diagnosticTests"] = _extract_diagnostic_terms(abstract_text, max_items=8)
            if not result["therapeuticOptions"]:
                result["therapeuticOptions"] = _extract_therapy_terms(abstract_text, max_items=8)

    # Second fallback: broader gene + disease association search
    if not result["clinicalSymptoms"] and not result["diagnosticTests"]:
        broad_query = (
            f'"{gene_symbol}"[Gene Name] AND '
            "(human[Organism]) AND "
            "(disease[Title/Abstract] OR disorder[Title/Abstract] OR syndrome[Title/Abstract] OR "
            "pathogenic[Title/Abstract] OR variant[Title/Abstract]) AND "
            f"({HUMAN_FILTER})"
        )
        pmids = _ncbi_search(broad_query, retmax=15)
        abstract_text = _ncbi_fetch_abstracts(pmids)
        if abstract_text:
            if not result["clinicalSymptoms"]:
                result["clinicalSymptoms"] = _extract_symptom_terms(abstract_text, max_items=7)
            if not result["diagnosticTests"]:
                result["diagnosticTests"] = _extract_diagnostic_terms(abstract_text, max_items=8)

    # Strategy 5: Carrier manifestations (X-linked genes)
    if gene_info["summary"] and any(w in (gene_info["summary"] or "").lower()
                                     for w in ["x-linked", "x linked", "xlink"]):
        result["carrierManifestations"] = [
            "X-linked inheritance: female carriers are typically asymptomatic or mildly affected due to skewed X-inactivation.",
            "Carrier testing available through molecular genetic testing of the causative gene.",
            "Routine cardiac and neurological screening recommended for female carriers.",
            "Variable expressivity in carriers depends on X-inactivation patterns.",
        ]

    # Strategy 6: Fetch clinical trials information for therapeutic options
    if not result["therapeuticOptions"] or len(result["therapeuticOptions"]) < 3:
        trial_therapies = _get_clinical_trials(gene_symbol, disease_terms[0] if disease_terms else disease_name)
        if trial_therapies:
            result["therapeuticOptions"].extend(trial_therapies)

    # Deduplicate results while preserving order
    result["clinicalSymptoms"] = list(dict.fromkeys(result["clinicalSymptoms"]))
    result["diagnosticTests"] = list(dict.fromkeys(result["diagnosticTests"]))
    result["therapeuticOptions"] = list(dict.fromkeys(result["therapeuticOptions"]))
    result["carrierManifestations"] = list(dict.fromkeys(result["carrierManifestations"]))

    # Final limits
    result["clinicalSymptoms"] = result["clinicalSymptoms"][:7]
    result["diagnosticTests"] = result["diagnosticTests"][:8]
    result["therapeuticOptions"] = result["therapeuticOptions"][:10]
    result["carrierManifestations"] = result["carrierManifestations"][:4]

    return result
