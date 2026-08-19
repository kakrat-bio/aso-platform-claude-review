import asyncio
import os
import sys
import time
import logging

# Add backend/ to sys.path so services resolve correctly
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dotenv import load_dotenv
load_dotenv()

import aiohttp
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional

try:  # ``uvicorn backend.main:app`` from the repository root
    from .services.notification_service import add_notification as _add_notification
    from .services.gene_service import EnsemblLookupUnavailable, clean_synonyms, get_gene_metadata, get_gene_phenotypes, ensembl_gene_url, build_gene_fallback_payload
    from .services.enrichment_service import get_gene_enrichment, get_aso_analysis
    from .services.constraint_service import get_human_constraint_metrics
    from .services.clinical_service import get_clinical_details
    from .services.protein_service import get_protein_db_ids
    from .services.variant_details_service import get_variant_details, get_clinvar_variants
    from .services.protein_properties_service import get_protein_properties
    from .services.single_cell_service import get_single_cell_expression
    from .services.tissue_expression_service import get_tissue_expression
    from .services.rna_halflife_service import get_rna_halflife
    from .services.dependency_service import get_gene_dependency
    from .services.fda_therapies_service import get_fda_therapies
    from .services.orphanet_service import get_orphanet_data
    from .services.mutation_breakdown_service import get_mutation_breakdown
    from .services.sequence_liability_service import get_sequence_liabilities
    from .api.mechanisms import router as mechanisms_router
    from .api.gene_silencing import router as gene_silencing_router
    from .api.gene_upregulation import router as gene_upregulation_router
    from .api.protein_replacement import router as protein_replacement_router
    from .api.isoform_engineering import router as isoform_engineering_router
    from .api.translational_regulation import router as translational_regulation_router
    from .api.rna_neutralization import router as rna_neutralization_router
    from .api.rna_editing import router as rna_editing_router
    from .api.upload import router as upload_router
    from .api.assistant import router as assistant_router
    from .api.notifications import router as notifications_router
    from .api.disease_search import router as disease_search_router
    from .api.auth import router as auth_router
    from .api.profile import router as profile_router
    from .api.reports import router as reports_router
    from .api.projects import router as projects_router
    from .api.bug_reports import router as bug_reports_router
    from .api.gene_search import router as gene_search_router
except ImportError:
    from services.notification_service import add_notification as _add_notification
    from services.gene_service import EnsemblLookupUnavailable, clean_synonyms, get_gene_metadata, get_gene_phenotypes, ensembl_gene_url, build_gene_fallback_payload
    from services.enrichment_service import get_gene_enrichment, get_aso_analysis
    from services.constraint_service import get_human_constraint_metrics
    from services.clinical_service import get_clinical_details
    from services.protein_service import get_protein_db_ids
    from services.variant_details_service import get_variant_details, get_clinvar_variants
    from services.protein_properties_service import get_protein_properties
    from services.single_cell_service import get_single_cell_expression
    from services.tissue_expression_service import get_tissue_expression
    from services.rna_halflife_service import get_rna_halflife
    from services.dependency_service import get_gene_dependency
    from services.fda_therapies_service import get_fda_therapies
    from services.orphanet_service import get_orphanet_data
    from services.mutation_breakdown_service import get_mutation_breakdown
    from services.sequence_liability_service import get_sequence_liabilities
    from api.mechanisms import router as mechanisms_router
    from api.gene_silencing import router as gene_silencing_router
    from api.gene_upregulation import router as gene_upregulation_router
    from api.protein_replacement import router as protein_replacement_router
    from api.isoform_engineering import router as isoform_engineering_router
    from api.translational_regulation import router as translational_regulation_router
    from api.rna_neutralization import router as rna_neutralization_router
    from api.rna_editing import router as rna_editing_router
    from api.upload import router as upload_router
    from api.assistant import router as assistant_router
    from api.notifications import router as notifications_router
    from api.disease_search import router as disease_search_router
    from api.auth import router as auth_router
    from api.profile import router as profile_router
    from api.reports import router as reports_router
    from api.projects import router as projects_router
    from api.bug_reports import router as bug_reports_router
    from api.gene_search import router as gene_search_router


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# NCBI E-utilities rate limit: 3 requests/sec without API key (lock created lazily for Python 3.9 compat)
_ncbi_req_lock = None
_last_ncbi_req = 0.0


def _get_ncbi_lock():
    global _ncbi_req_lock
    if _ncbi_req_lock is None:
        _ncbi_req_lock = asyncio.Lock()
    return _ncbi_req_lock

app = FastAPI(title="ASO Platform API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(mechanisms_router)
app.include_router(protein_replacement_router)
app.include_router(gene_silencing_router)
app.include_router(gene_upregulation_router)
app.include_router(isoform_engineering_router)
app.include_router(translational_regulation_router)
app.include_router(rna_neutralization_router)
app.include_router(rna_editing_router)
app.include_router(upload_router)
app.include_router(assistant_router)
app.include_router(notifications_router)
app.include_router(disease_search_router)
app.include_router(auth_router)
app.include_router(profile_router)
app.include_router(reports_router)
app.include_router(projects_router)
app.include_router(bug_reports_router)
app.include_router(gene_search_router)


@app.on_event("startup")
def _init_db():
    from database.db import init_db
    init_db()

SPECIES_TAXON_IDS = {
    "homo_sapiens": 9606,
    "mus_musculus": 10090,
    "rattus_norvegicus": 10116,
    "macaca_fascicularis": 9541,
    "macaca_mulatta": 9544,
    "danio_rerio": 7955,
    "drosophila_melanogaster": 7227,
    "caenorhabditis_elegans": 6239,
    "saccharomyces_cerevisiae": 4932,
    "schizosaccharomyces_pombe": 4896,
    "canis_lupus_familiaris": 9615,
    "felis_catus": 9685,
    "sus_scrofa": 9823,
    "bos_taurus": 9913,
    "equus_caballus": 9796,
    "ovis_aries": 9940,
    "capra_hircus": 9925,
    "gallus_gallus": 9031,
    # Tier 4 — Plants (Ensembl Plants)
    "arabidopsis_thaliana": 3702,
    "oryza_sativa": 39947,
    "zea_mays": 4577,
    "triticum_aestivum": 4565,
    "solanum_lycopersicum": 4081,
    # Tier 6 — Bacteria (NCBI Taxonomy)
    "escherichia_coli": 511145,
    "staphylococcus_aureus": 1280,
    "mycobacterium_tuberculosis": 83333,
    "pseudomonas_aeruginosa": 208964,
}

class TargetRequest(BaseModel):
    gene_symbol: str
    organism: str  
    disease_name: Optional[str] = None

def get_safe_ensembl_url(species: str, gene_id: str) -> str:
    """Safely resolves the Ensembl link, bypassing any type conversion conflicts."""
    try:
        if callable(ensembl_gene_url):
            return ensembl_gene_url(species, gene_id)
    except Exception:
        pass
    # Capitalize species name for Ensembl URL (e.g. "homo_sapiens" -> "Homo_sapiens")
    parts = species.split("_")
    formatted_species = "_".join(p.capitalize() for p in parts) if parts else "Homo_sapiens"
    return f"https://www.ensembl.org/{formatted_species}/Gene/Summary?g={gene_id}"

async def _ncbi_call(session: aiohttp.ClientSession, db: str, term: str) -> Optional[dict]:
    """Rate-limited NCBI ESearch (max ~2.5 req/sec)."""
    global _last_ncbi_req
    async with _get_ncbi_lock():
        now = time.monotonic()
        since_last = now - _last_ncbi_req
        if since_last < 0.4:
            await asyncio.sleep(0.4 - since_last)
        _last_ncbi_req = time.monotonic()
        async with session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": db, "term": term, "retmode": "json"},
            timeout=aiohttp.ClientTimeout(total=5),
        ) as response:
            return await response.json() if response.status == 200 else {}


async def _ncbi_fetch(session: aiohttp.ClientSession, url: str, params: dict) -> Optional[dict]:
    """Rate-limited NCBI EUtils GET (for esearch, esummary, etc.)."""
    global _last_ncbi_req
    async with _get_ncbi_lock():
        now = time.monotonic()
        since_last = now - _last_ncbi_req
        if since_last < 0.4:
            await asyncio.sleep(0.4 - since_last)
        _last_ncbi_req = time.monotonic()
        params.setdefault("retmode", "json")
        async with session.get(
            url, params=params, timeout=aiohttp.ClientTimeout(total=5),
        ) as response:
            return await response.json() if response.status == 200 else {}


async def get_pubmed_count(session: aiohttp.ClientSession, term: str) -> int:
    try:
        data = await _ncbi_call(session, "pubmed", term)
        return int((data.get("esearchresult") or {}).get("count", 0))
    except Exception:
        return 0


async def get_clinvar_count(session: aiohttp.ClientSession, symbol: str) -> Optional[int]:
    """Return the live ClinVar record count for a human gene."""
    try:
        data = await _ncbi_call(session, "clinvar", f"{symbol}[gene]")
        count = (data.get("esearchresult") or {}).get("count")
        return int(count) if count is not None else None
    except Exception:
        return None

async def get_dbsnp_count(session: aiohttp.ClientSession, symbol: str) -> Optional[int]:
    """Return the live dbSNP variant count for a gene."""
    try:
        data = await _ncbi_call(session, "snp", f"{symbol}[gene]")
        count = (data.get("esearchresult") or {}).get("count")
        return int(count) if count is not None else None
    except Exception:
        return None


async def fetch_ncbi_aliases(session: aiohttp.ClientSession, gene_symbol: str, taxon_id: int) -> list:
    """Fetch gene aliases from NCBI esummary by gene symbol + taxon ID."""
    try:
        search_term = f"{gene_symbol}[Gene Name] AND {taxon_id}[Taxonomy ID]"
        data = await _ncbi_fetch(
            session,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            {"db": "gene", "term": search_term},
        )
        if not data:
            return []
        ids = (data.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            return []

        summary = await _ncbi_fetch(
            session,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            {"db": "gene", "id": ids[0]},
        )
        if not summary:
            return []
        result = (summary.get("result") or {}).get(ids[0]) or {}
        aliases_str = result.get("otheraliases", "")
        if aliases_str:
            return [a.strip() for a in aliases_str.split(", ") if a.strip()]
    except Exception:
        pass
    return []

async def get_rxiv_count(session: aiohttp.ClientSession, symbol: str) -> int:
    """Fetch preprint count from bioRxiv/medRxiv via Europe PMC."""
    try:
        async with session.get(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query": f"{symbol}[title] AND SRC:PPR", "format": "json", "resultType": "lite"},
            timeout=aiohttp.ClientTimeout(total=8),
        ) as response:
            data = await response.json() if response.status == 200 else {}
            return int(data.get("hitCount", 0))
    except Exception:
        return 0


async def fetch_gene_from_ncbi(session: aiohttp.ClientSession, gene_symbol: str, species: str, taxon_id: int) -> Optional[dict]:
    """Fallback: look up a gene via NCBI Gene API when Ensembl doesn't cover the species."""
    try:
        search_term = f"{gene_symbol}[Gene Name] AND {taxon_id}[Taxonomy ID]"
        data = await _ncbi_fetch(
            session,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            {"db": "gene", "term": search_term},
        )
        if not data:
            return None
        ids = (data.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            return None
        gene_id = ids[0]

        summary = await _ncbi_fetch(
            session,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            {"db": "gene", "id": gene_id},
        )
        if not summary:
            return None
        result = (summary.get("result") or {}).get(gene_id) or {}
        if not result:
            return None

        official_symbol = result.get("name", gene_symbol)
        ncbi_gene_id = result.get("uid", gene_id)
        aliases = result.get("otheraliases", "").split(", ") if result.get("otheraliases") else []

        chromosome = result.get("chromosome")
        genomic_info = result.get("genomicinfo") or []
        if genomic_info and isinstance(genomic_info, list):
            genomic_start = genomic_info[0].get("chrstart")
            genomic_stop = genomic_info[0].get("chrstop")
        else:
            genomic_start = genomic_stop = None

        # NCBI reports minus-strand genes with chrstart > chrstop; infer strand
        # so the gene page no longer shows a blank Strand field for fallbacks.
        if genomic_start is not None and genomic_stop is not None and genomic_start != genomic_stop:
            inferred_strand = -1 if genomic_start > genomic_stop else 1
        else:
            inferred_strand = None

        return {
            "id": f"NCBI:{ncbi_gene_id}",
            "officialSymbol": official_symbol,
            "geneName": result.get("description", official_symbol),
            "seq_region_name": chromosome,
            "start": genomic_start,
            "end": genomic_stop,
            "cytoband": result.get("chromosomeLocation"),
            "genomeBuild": None,
            "strand": inferred_strand,
            "biotype": result.get("type_of_gene", "protein_coding"),
            "synonyms": aliases,
            "nomenclatureId": None,
            "canonicalTranscript": None,
            "otherTranscripts": [],
            "totalTranscripts": 0,
            "exonCount": 0,
            "proteinLength": 0,
            "proteinId": None,
            "entrezGeneId": str(ncbi_gene_id),
        }
    except Exception as e:
        logger.info(f"NCBI Gene fallback failed for {gene_symbol} in {species}: {e}")
        return None

async def fetch_expression_details(session: aiohttp.ClientSession, symbol: str, ensembl_gene_id: str, species: str) -> dict:
    """Fetch tissue expression data using the new multi-source tissue expression service."""
    try:
        tissue_data = await get_tissue_expression(
            symbol=symbol,
            ensembl_id=ensembl_gene_id,
            species=species,
        )

        # Map to the expected format
        expr_data = {
            "available": tissue_data.get("available", False),
            "top_tissue": tissue_data.get("top_tissue"),
            "tpm": tissue_data.get("tpm"),
            "gtex_level": f"{tissue_data['top_tissue']} ({tissue_data['tpm']} TPM)" if tissue_data.get("top_tissue") and tissue_data.get("tpm") else None,
            "hpa_level": None,
            "top_tissues": tissue_data.get("topTissues", []),
            "expression_cv": tissue_data.get("expressionStabilityCV"),
            "vital_organ_tpm": None,
            "vital_organ_tissues": tissue_data.get("vitalOrganTissues", []),
            "dominant_isoform_fraction": tissue_data.get("dominantIsoformFraction"),
            "dominant_isoform_id": tissue_data.get("dominantIsoformId"),
            "source": tissue_data.get("source"),
        }

        # Calculate vital_organ_tpm from top_tissues
        vital_keywords = ["Heart", "Kidney", "Lung", "Brain", "Liver"]
        vital_matches = [
            t for t in expr_data["top_tissues"]
            if any(kw in t.get("name", "") for kw in vital_keywords)
        ]
        expr_data["vital_organ_tpm"] = max((t.get("tpm", 0) or 0 for t in vital_matches), default=None)

        return expr_data
    except Exception as e:
        logger.info(f"Tissue expression lookup failed for {symbol}: {e}")
        return {
            "available": False,
            "top_tissue": None,
            "tpm": None,
            "gtex_level": None,
            "hpa_level": None,
            "top_tissues": [],
            "expression_cv": None,
            "vital_organ_tpm": None,
            "vital_organ_tissues": [],
            "dominant_isoform_fraction": None,
            "dominant_isoform_id": None,
            "source": None,
        }

async def fetch_disease_associations(session: aiohttp.ClientSession, symbol: str, ensembl_id: str, species: str) -> dict:
    result = {"diseases": [], "omim_id": None, "source": []}
    is_human = species == "homo_sapiens"

    def add_disease(name: Optional[str]) -> None:
        if not name:
            return
        cleaned = name.strip()
        if cleaned and cleaned.casefold() not in {item.casefold() for item in result["diseases"]}:
            result["diseases"].append(cleaned)

    def capture_omim(value: object) -> None:
        if result["omim_id"] or not isinstance(value, str):
            return
        normalized = value.strip()
        if normalized.upper().startswith("OMIM:"):
            result["omim_id"] = normalized.split(":", 1)[1].strip()
        elif normalized.isdigit() and len(normalized) == 6:
            result["omim_id"] = normalized

    if is_human and ensembl_id and ensembl_id.startswith("ENSG"):
        try:
            ot_url = "https://api.platform.opentargets.org/api/v4/graphql"
            query = """
            query targetInfo($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                associatedDiseases(page: {index: 0, size: 20}) {
                  rows {
                    disease { name dbXRefs }
                  }
                }
              }
            }
            """
            async with session.post(
                ot_url,
                json={"query": query, "variables": {"ensemblId": ensembl_id}},
                timeout=aiohttp.ClientTimeout(total=6),
            ) as ot_res:
                if ot_res.status == 200:
                    ot_data = await ot_res.json()
                    target_data = (ot_data.get("data") or {}).get("target")
                    if target_data:
                        rows = (target_data.get("associatedDiseases") or {}).get("rows", [])
                        for row in rows:
                            disease_node = row.get("disease") or {}
                            add_disease(disease_node.get("name"))
                            for xref in disease_node.get("dbXRefs", []) or []:
                                capture_omim(xref)
                        if rows:
                            result["source"].append("Open Targets Platform")
        except Exception as e:
            logger.warning(f"Open Targets lookup failed for {symbol}: {e}")

    # Open Targets returns ranked disease associations, while Ensembl carries
    # phenotype and model-organism annotations that often fill in gaps.  Merge
    # them rather than treating the first successful provider as exhaustive.
    phenotypes = await asyncio.to_thread(get_gene_phenotypes, symbol, species)
    if phenotypes:
        phenotype_sources = set()
        for phenotype in phenotypes:
            add_disease(phenotype.get("description"))
            source = phenotype.get("source")
            if source:
                phenotype_sources.add(source)
            for key in ("accession", "id", "ontology_accession", "external_reference"):
                capture_omim(phenotype.get(key))
        result["source"].extend(sorted(phenotype_sources) or ["Ensembl Phenotype"])

    result["diseases"] = result["diseases"][:20]
    result["source"] = list(dict.fromkeys(result["source"]))

    return result

@app.post("/api/pipeline/initialize-target")
async def initialize_target(payload: TargetRequest):
    try:
        symbol_upper = payload.gene_symbol.strip()
        species = payload.organism or "homo_sapiens"
        ORGANISM_ID_TO_ENSEMBL = {
            "human": "homo_sapiens", "mouse": "mus_musculus", "rat": "rattus_norvegicus",
            "cynomolgus": "macaca_fascicularis", "rhesus": "macaca_mulatta",
            "zebrafish": "danio_rerio", "fruitfly": "drosophila_melanogaster",
            "celegans": "caenorhabditis_elegans", "yeast": "saccharomyces_cerevisiae",
            "fissionyeast": "schizosaccharomyces_pombe", "dog": "canis_lupus_familiaris",
            "cat": "felis_catus", "pig": "sus_scrofa", "cow": "bos_taurus",
            "horse": "equus_caballus", "sheep": "ovis_aries", "goat": "capra_hircus",
            "chicken": "gallus_gallus",
            # Tier 4 — Plants
            "arabidopsis": "arabidopsis_thaliana", "rice": "oryza_sativa",
            "maize": "zea_mays", "wheat": "triticum_aestivum",
            "tomato": "solanum_lycopersicum",
            # Tier 6 — Bacteria
            "ecoli": "escherichia_coli", "saureus": "staphylococcus_aureus",
            "mtuberculosis": "mycobacterium_tuberculosis", "paeruginosa": "pseudomonas_aeruginosa",
        }
        species = ORGANISM_ID_TO_ENSEMBL.get(species.lower(), species)
        is_human = species == "homo_sapiens"

        try:
            meta = get_gene_metadata(symbol_upper, species)
        except EnsemblLookupUnavailable:
            meta = None
        if not meta or not meta.get("id"):
            # Try NCBI Gene fallback for species not in Ensembl
            async with aiohttp.ClientSession() as session:
                meta = await fetch_gene_from_ncbi(session, symbol_upper, species, SPECIES_TAXON_IDS.get(species, 0))
            if not meta or not meta.get("id"):
                raise HTTPException(
                    status_code=404,
                    detail=f'Gene "{symbol_upper}" was not found in {species.replace("_", " ")} (Ensembl or NCBI).',
                )

        gene_id = meta["id"]
        official_symbol = meta.get("officialSymbol") or symbol_upper
        taxon_id = SPECIES_TAXON_IDS.get(species, 0)

        async with aiohttp.ClientSession() as session:
            pubmed_task = get_pubmed_count(session, f"{official_symbol}[gene]")
            review_task = get_pubmed_count(session, f"{official_symbol}[gene] AND review[pt]")
            clinical_trial_task = get_pubmed_count(session, f"{official_symbol}[gene] AND clinical trial[pt]")
            case_report_task = get_pubmed_count(session, f"{official_symbol}[gene] AND case reports[pt]")
            biorxiv_task = get_rxiv_count(session, official_symbol)
            medrxiv_task = asyncio.sleep(0, result=0)
            disease_task = fetch_disease_associations(session, official_symbol, gene_id, species)
            expr_task = fetch_expression_details(session, official_symbol, gene_id, species)
            clinvar_task = get_clinvar_count(session, official_symbol) if is_human else None
            dbsnp_task = get_dbsnp_count(session, official_symbol)
            ncbi_aliases_task = fetch_ncbi_aliases(session, official_symbol, taxon_id) if not is_human and taxon_id else asyncio.sleep(0, result=[])

            pubmed_count, review_count, clinical_trial_count, case_report_count, biorxiv_count, medrxiv_count, disease_info, expr_details, clinvar_count, dbsnp_count, ncbi_aliases = await asyncio.gather(
                pubmed_task,
                review_task,
                clinical_trial_task,
                case_report_task,
                biorxiv_task,
                medrxiv_task,
                disease_task,
                expr_task,
                clinvar_task if clinvar_task else asyncio.sleep(0, result=None),
                dbsnp_task if dbsnp_task else asyncio.sleep(0, result=None),
                ncbi_aliases_task,
            )

        synonyms_list = clean_synonyms(
            [*(meta.get("synonyms") or []), *ncbi_aliases],
            official_symbol,
        )
        
        disease_resolved = "; ".join(disease_info["diseases"][:3]) if disease_info["diseases"] else None
        if not disease_resolved and payload.disease_name:
            disease_resolved = payload.disease_name.strip()
            
        omim_id = f"#{disease_info['omim_id']}" if disease_info.get("omim_id") else None

        async def _run_sync(fn, label, default=None, timeout_seconds: float = 18.0):
            try:
                return await asyncio.wait_for(asyncio.to_thread(fn), timeout=timeout_seconds)
            except asyncio.TimeoutError:
                logger.warning("%s timed out for %s", label, official_symbol)
                return default
            except Exception as e:
                logger.warning("%s failed for %s: %s", label, official_symbol, e)
                return default

        enrichment_data, constraint_data, aso_data, rna_halflife_data, dependency_data, variant_details, top_variants, single_cell, clinical_details, fda_therapies, orphanet_data, mutation_data = await asyncio.gather(
            _run_sync(lambda: get_gene_enrichment(gene_id, taxon_id, official_symbol), "Enrichment lookup", {}),
            _run_sync(lambda: get_human_constraint_metrics(official_symbol) if is_human else {}, "Constraint lookup", {}),
            _run_sync(lambda: get_aso_analysis(gene_id, taxon_id), "ASO analysis", {}),
            _run_sync(lambda: get_rna_halflife(official_symbol) if is_human else {}, "RNA half-life lookup", {}),
            _run_sync(lambda: get_gene_dependency(official_symbol), "Dependency lookup", {}),
            _run_sync(lambda: get_variant_details(gene_symbol=official_symbol, ensembl_gene_id=meta.get("id"), entrez_id=meta.get("entrezGeneId")) if is_human else {}, "Variant details lookup", {}),
            _run_sync(lambda: get_clinvar_variants(meta.get("id")) if is_human else [], "Top variants lookup", [], timeout_seconds=25.0),
            _run_sync(lambda: get_single_cell_expression(ensembl_id=gene_id, gene_symbol=official_symbol) if is_human else {}, "Single-cell lookup", {}),
            _run_sync(lambda: get_clinical_details(gene_symbol=official_symbol, disease_name=disease_resolved, omim_id=disease_info.get("omim_id"), phenotypes=disease_info.get("diseases")), "Clinical details lookup", {}),
            _run_sync(lambda: get_fda_therapies(official_symbol, disease_resolved) if is_human else {}, "FDA therapies lookup", {}, timeout_seconds=18.0),
            _run_sync(lambda: get_orphanet_data(official_symbol, ensembl_id=gene_id, disease_name=disease_resolved, phenotypes=disease_info.get("diseases")) if is_human else {}, "Orphanet lookup", {}),
            _run_sync(lambda: get_mutation_breakdown(official_symbol) if is_human else {}, "Mutation breakdown lookup", {}, timeout_seconds=60.0),
        )

        # This used to pass `aso_sequence=aso_data.get("cdsSequence")` — the
        # TARGET GENE's coding sequence, thousands of nucleotides long — into a
        # function whose descriptors (length, GC, Tm, CpG count) are meant for a
        # 16-25 nt oligonucleotide. Every CDS tripped "Length > 25 nt" as a
        # Lipinski violation and had a PBPK curve drawn for it. No candidate
        # oligo exists at this point in the pipeline, so no sequence liability
        # is computed here. What gene context alone genuinely supports is the
        # on-target pharmacology, and that is what remains.
        gene_pharmacology = await _run_sync(lambda: get_sequence_liabilities(
            aso_sequence=None,
            gene_context={
                "vitalOrganTpm": expr_details.get("vital_organ_tpm"),
                "vitalOrganTissues": expr_details.get("vital_organ_tissues", []),
                "essentialGene": dependency_data.get("essentialGene"),
                "loeufDecile": constraint_data.get("loeufDecile"),
            }
        ), "on-target pharmacology lookup", {})

        # Protein chain — depends on protein_db_ids → protein_properties
        protein_props = {}
        protein_db = {}
        try:
            protein_db = get_protein_db_ids(
                uniprot_id=meta.get("proteinId"),
                gene_symbol=official_symbol,
                entrez_id=meta.get("entrezGeneId"),
                taxon_id=taxon_id,
            )
            uniprot_acc = protein_db.get("uniprotAccession")
            protein_props = get_protein_properties(
                ensembl_protein_id=meta.get("proteinId"),
                gene_symbol=official_symbol,
                uniprot_accession=uniprot_acc,
                taxon_id=taxon_id,
            )
        except Exception as exc:
            logger.warning("Protein properties lookup failed for %s: %s", official_symbol, exc)
            protein_props = {}
            protein_db = {}

        _add_notification(
            "analysis",
            f"Gene lookup completed for {symbol_upper}",
            f"Retrieved data for {symbol_upper} in {species.replace('_', ' ')}.",
        )

        raw_strand = meta.get("strand")
        if str(raw_strand) in ["-1", "-"]:
            strand_display = "Reverse (−)"
        elif str(raw_strand) in ["1", "+1", "+"]:
            strand_display = "Forward (+)"
        else:
            strand_display = None

        hgnc_display = meta.get("nomenclatureId") or meta.get("hgncId")
        if not hgnc_display:
            hgnc_display = None

        gene_type_display = meta.get("biotype") or meta.get("geneType") or "protein_coding"

        start, end = meta.get("start"), meta.get("end")
        # abs() guards against providers (e.g. NCBI fallback) that report
        # minus-strand genes with start > end, which would otherwise render a
        # negative gene length.
        gene_length = (abs(end - start) + 1) if (start is not None and end is not None) else None
        exon_count = meta.get("exonCount")
        protein_length = meta.get("proteinLength")
        cds_length = (protein_length * 3 + 3) if protein_length else None

        ensembl_url = get_safe_ensembl_url(species, gene_id)

        deep_links = {
            "ensembl": ensembl_url,
            "ncbi": f"https://www.ncbi.nlm.nih.gov/gene/?term={official_symbol}",
            "gtex": f"https://gtexportal.org/home/gene/{official_symbol}" if is_human else None,
            "hpa": f"https://www.proteinatlas.org/search/{official_symbol}" if is_human else None,
            "uniprot": f"https://www.uniprot.org/uniprotkb?query={official_symbol}",
            "clinvar": f"https://www.ncbi.nlm.nih.gov/clinvar/?term={official_symbol}%5Bgene%5D",
            "kegg": f"https://www.genome.jp/dbget-bin/www_bget?q={official_symbol}",
            "reactome": f"https://reactome.org/content/query?q={official_symbol}",
            "pubmed": f"https://pubmed.ncbi.nlm.nih.gov/?term={official_symbol}%5Bgene%5D",
            "clinicaltrials": f"https://clinicaltrials.gov/search?cond={official_symbol}",
            "omim": f"https://www.omim.org/search?search={official_symbol}",
            "go": f"https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=ENSEMBL%3A{gene_id}",
            "string": f"https://string-db.org/cgi/network?identifiers={official_symbol}&species={taxon_id}",
        }

        tissue_level = None
        if expr_details["available"]:
            tpm_value = expr_details["tpm"] or 0
            tissue_level = "High" if tpm_value > 25 else ("Medium" if tpm_value > 5 else "Low")

        # Single-cell prevalence (% cell types with nCPM > 0)
        cell_type_all = single_cell.get("cellTypeAll", {})
        sc_total = len(cell_type_all)
        sc_positive = sum(1 for v in cell_type_all.values() if (v or 0) > 0)
        single_cell_prevalence = round(sc_positive / sc_total, 3) if sc_total > 0 else None

        # Developmental / age-dependent expression pattern
        gene_func = enrichment_data.get("geneFunction") or meta.get("geneName") or ""
        dev_keywords = [
            "development", "differentiation", "embryonic", "fetal", "morphogenesis",
            "organogenesis", "neurogenesis", "myogenesis", "angiogenesis",
        ]
        has_dev = any(kw in gene_func.lower() for kw in dev_keywords)
        developmental_expression = "Developmentally Regulated" if has_dev else "Ubiquitous (Age-Stable)"

        # Alternative polyadenylation — inferred from transcript isoform diversity
        total_transcripts = meta.get("totalTranscripts") or len(meta.get("otherTranscripts", [])) + 1
        if total_transcripts > 2:
            alternative_polyadenylation = f"Multiple 3' UTR Isoforms ({total_transcripts} transcripts)"
        elif total_transcripts > 1:
            alternative_polyadenylation = "Few 3' UTR Isoforms"
        else:
            alternative_polyadenylation = "Single 3' UTR"

        # Cytoplasmic vs nuclear retention index — heuristic from subcellular location + intron count
        subcell = protein_props.get("subcellularLocation") or ""
        is_nuclear = any(kw in subcell.lower() for kw in ["nucleus", "nucleolar", "nuclear", "nucleoplasm"])
        n_exons = meta.get("exonCount") or 0
        n_introns = max(int(n_exons) - 1, 0) if n_exons else max(int(meta.get("totalTranscripts", 5)) or 5, 1)
        if is_nuclear and n_introns > 5:
            nuclear_retention_index = round(min(0.4 + 0.4 * (1 - 1 / (1 + n_introns / 15)), 0.85), 2)
        elif is_nuclear:
            nuclear_retention_index = round(min(0.3 + 0.2 * (n_introns / 10), 0.6), 2)
        elif n_introns > 10:
            nuclear_retention_index = round(min(0.2 + 0.15 * (n_introns / 20), 0.45), 2)
        else:
            nuclear_retention_index = round(min(0.15 + 0.1 * (n_introns / 10), 0.3), 2)

        payload_dict = {
            "organism": species,
            "diseaseName": payload.disease_name.strip() if payload.disease_name else None,
            "geneSymbol": official_symbol,
            "geneName": meta.get("geneName"),
            "geneFunction": enrichment_data.get("geneFunction"),
            "geneId": gene_id,  
            "entrezGeneId": meta.get("entrezGeneId") or enrichment_data.get("entrezGeneId"),
            "hgncId": hgnc_display,
            "chromosome": meta.get("seq_region_name"),
            "location": f"{meta.get('seq_region_name')}:{start}-{end}" if start and end else None,
            "cytoband": meta.get("cytoband"),
            "genomeBuild": meta.get("genomeBuild"),
            "genomicStart": start,
            "genomicEnd": end,
            "strand": strand_display,
            "geneType": gene_type_display,
            "synonyms": synonyms_list,
            "source": ["Ensembl"] if str(gene_id).startswith("ENSG") else (["NCBI"] if str(gene_id).startswith("NCBI:") else ["Ensembl"]),
            "taxonId": str(taxon_id),

            "canonicalTranscript": meta.get("canonicalTranscript"),
            "canonicalTranscriptLabel": "Canonical (MANE Select)" if is_human else "Canonical",
            
            # Explicitly tie these to the parsed Ensembl values from service layer
            "otherTranscripts": meta.get("otherTranscripts", []),
            "totalTranscripts": meta.get("totalTranscripts") or len(meta.get("otherTranscripts", [])) + 1,

            "variantExamples": [],
            "totalKnownVariantsClinvar": None,
            "topVariants": top_variants[:10] if isinstance(top_variants, list) else [],

            "defaultTissue": expr_details["top_tissue"],
            "tissueExpressionLevel": tissue_level,
            "tissueTpm": expr_details["tpm"],
            "topTissues": expr_details["top_tissues"],

            "defaultCellType": single_cell.get("cellType"),
            "cellExpressionLevel": single_cell.get("cellType"),
            "cellTpm": single_cell.get("cellTpm"),
            "cellTypeAll": cell_type_all,

            "expressionStabilityCV": expr_details.get("expression_cv"),
            "vitalOrganTpm": expr_details.get("vital_organ_tpm"),
            "vitalOrganTissues": expr_details.get("vital_organ_tissues", []),
            "dominantIsoformFraction": expr_details.get("dominant_isoform_fraction"),
            "dominantIsoformId": expr_details.get("dominant_isoform_id"),
            "diseaseFoldChange": None,
            "singleCellPrevalence": single_cell_prevalence,
            "circadianAmplitude": expr_details.get("circadianAmplitude"),
            "intronRetentionRatio": None,
            "developmentalExpression": developmental_expression,
            "alternativePolyadenylation": alternative_polyadenylation,
            "nuclearRetentionIndex": nuclear_retention_index,

            "proteinId": meta.get("proteinId"),
            "proteinName": meta.get("geneName"),
            "proteinLength": protein_length,

            # Protein properties from UniProt
            "molecularWeight": protein_props.get("molecularWeight"),
            "isoelectricPoint": protein_props.get("isoelectricPoint"),
            "secondaryStructureDistribution": protein_props.get("secondaryStructureDistribution"),
            "criticalPhosphorylationSite": protein_props.get("criticalPhosphorylationSite"),
            "ubiquitinationTarget": protein_props.get("ubiquitinationTarget"),
            "quaternaryStructure": protein_props.get("quaternaryStructure"),
            "stabilityScore": protein_props.get("stabilityScore"),
            "subcellularLocation": protein_props.get("subcellularLocation"),
            "criticalFunctionalDomains": protein_props.get("criticalFunctionalDomains"),
            "disorderedContent": protein_props.get("disorderedContent"),
            "proteosomalTurnover": protein_props.get("proteosomalTurnover"),

            "alphafoldPlddt": protein_props.get("alphafoldPlddt"),
            "gravyIndex": protein_props.get("gravyIndex"),
            "proteinAbundance": protein_props.get("proteinAbundance"),
            "tractability": protein_props.get("tractability"),

            # Protein database IDs from UniProt + NCBI (available for all organisms)
            "interproId": protein_db.get("interproId"),
            "pfamId": protein_db.get("pfamId"),
            "pdbId": protein_db.get("pdbId"),
            "uniprotAccession": protein_db.get("uniprotAccession"),

            # Use gnomAD mutation rate if available (overrides protein service)
            **({"mutationRate": constraint_data["mutationRate"]} if constraint_data.get("mutationRate") else {}),

            # Top ClinVar variant details
            "topHgvsName": variant_details.get("topHgvsName"),
            "topRsId": variant_details.get("topRsId"),
            **({"populationFrequencyMaf": constraint_data["populationFrequencyMaf"]} if constraint_data.get("populationFrequencyMaf") else {}),

            "disease": disease_resolved,
            "diseaseAssociation": disease_resolved if disease_resolved else "None identified",
            "diseaseAssociationSource": disease_info.get("source", ["Ensembl Phenotype"] if disease_resolved else []),
            "phenotypes": disease_info.get("diseases", []),
            "associationStatus": "Reported" if disease_resolved else None,
            "omimId": omim_id,
            "diseaseMechanism": clinical_details.get("diseaseMechanism"),
            "diagnosticTests": clinical_details.get("diagnosticTests", []),
            "clinicalSymptoms": clinical_details.get("clinicalSymptoms", []),
            "carrierManifestations": clinical_details.get("carrierManifestations", []),
            "therapeuticOptions": clinical_details.get("therapeuticOptions", []),

            "exonCount": exon_count,
            "intronCount": (exon_count - 1) if exon_count else None,
            "cdsLength": cds_length,
            "geneLength": gene_length,

            # Provide actual counts from live API lookups
            "dbSnpCount": dbsnp_count,
            "gnomadAvailable": bool(constraint_data.get("intolerantToLossScore") or constraint_data.get("loeufScore")),
            "clinvarVariantCount": clinvar_count,

            "gtexAvailable": expr_details["available"],
            "humanProteinAtlasLevel": expr_details["hpa_level"] or ("Profile linked" if is_human else None),
            "gtexExpressionLevel": expr_details["gtex_level"],

            "haploinsufficiencyScore": constraint_data.get("haploinsufficiencyScore"),
            "intolerantToLossScore": constraint_data.get("intolerantToLossScore"),
            "recessiveConstraintZ": constraint_data.get("recessiveConstraintZ"),
            "hetExcessZ": constraint_data.get("hetExcessZ"),
            "compositeConstraintIndex": constraint_data.get("compositeConstraintIndex"),

            "loeufDecile": constraint_data.get("loeufDecile"),
            "triplosensitivity": constraint_data.get("triplosensitivity"),
            "activeIsoforms": aso_data.get("activeIsoforms"),
            "spliceSwitches": aso_data.get("spliceSwitches"),
            "structuralAccessibility": aso_data.get("structuralAccessibility"),
            "splicingMotifDensity": aso_data.get("splicingMotifDensity"),
            "preclinicalConservation": aso_data.get("preclinicalConservation"),
            "gQuadruplexes": aso_data.get("gQuadruplexes"),
            "cpgDensity": aso_data.get("cpgDensity"),
            "selfDimerRisk": aso_data.get("selfDimerRisk"),
            "polygTracts": aso_data.get("polygTracts"),
            "transcriptSpecificity": aso_data.get("transcriptSpecificity"),
            "codonUsageBias": aso_data.get("codonUsageBias"),

            # RNA half-life and dependency
            "rnaHalflife": rna_halflife_data.get("rnaHalflife"),
            "rnaHalflifeHours": rna_halflife_data.get("rnaHalflifeHours"),
            "rnaHalflifeSource": rna_halflife_data.get("rnaHalflifeSource"),
            "depmapDependency": dependency_data.get("depmapDependency"),
            "depmapDependencyScore": dependency_data.get("depmapDependencyScore"),
            "essentialGene": dependency_data.get("essentialGene"),
            "essentialGeneGeneTrap": dependency_data.get("essentialGeneGeneTrap"),
            "essentialGeneCrispr": dependency_data.get("essentialGeneCrispr"),
            "essentialGeneCrispr2": dependency_data.get("essentialGeneCrispr2"),
            "depmapSource": dependency_data.get("depmapSource"),

            "deepLinks": deep_links,

            "keggCount": enrichment_data.get("keggCount"),
            "reactomeCount": enrichment_data.get("reactomeCount"),
            "keggPathwayName": enrichment_data.get("keggPathwayName"),
            "reactomePathwayName": enrichment_data.get("reactomePathwayName"),
            "keggPathwayId": enrichment_data.get("keggPathwayId"),
            "reactomePathwayId": enrichment_data.get("reactomePathwayId"),
            "pathwayCommonsCount": enrichment_data.get("pathwayCommonsCount"),

            "goBiologicalProcess": enrichment_data.get("goBiologicalProcess"),
            "goMolecularFunction": enrichment_data.get("goMolecularFunction"),
            "goCellularComponent": enrichment_data.get("goCellularComponent"),
            "goBiologicalProcessAnnotations": enrichment_data.get("goBiologicalProcessAnnotations"),
            "goMolecularFunctionAnnotations": enrichment_data.get("goMolecularFunctionAnnotations"),
            "goCellularComponentAnnotations": enrichment_data.get("goCellularComponentAnnotations"),
            "pathwayHighlight": enrichment_data.get("pathwayHighlight"),
            "goBiologicalProcessHighlight": enrichment_data.get("goBiologicalProcessHighlight"),
            "goMolecularFunctionHighlight": enrichment_data.get("goMolecularFunctionHighlight"),
            "goCellularComponentHighlight": enrichment_data.get("goCellularComponentHighlight"),

            "stringHighConfidenceCount": enrichment_data.get("stringHighConfidenceCount"),
            "mediumConfidenceCount": enrichment_data.get("mediumConfidenceCount"),
            "totalInteractors": enrichment_data.get("totalInteractors"),
            "experimentalCount": enrichment_data.get("experimentalCount"),
            "databaseCount": enrichment_data.get("databaseCount"),
            "topInteractors": enrichment_data.get("topInteractors", []),
            "interactionNetworkDensity": enrichment_data.get("interactionNetworkDensity"),

            "pubmedArticleCount": pubmed_count,
            "reviewCount": review_count,
            "clinicalTrialsCount": clinical_trial_count,
            "caseReportsCount": case_report_count,
            "preprintCount": biorxiv_count,

            # New fields: genomic overview
            "genomicSize": gene_length,
            "mrnaLength": cds_length,
            "proteinMass": protein_props.get("molecularWeight"),
            "genomicOverviewDetails": {
                "canonicalTranscript": meta.get("canonicalTranscript"),
                "canonicalTranscriptLink": (ensembl_url.replace("/Gene/Summary?g=", "/Transcript/Summary?t=") if meta.get("canonicalTranscript") else None),
                "otherTranscripts": meta.get("otherTranscripts", [])[:5],
                "exonCount": exon_count,
                "proteinLength": protein_length,
                "proteinId": meta.get("proteinId"),
            },

            # FDA-approved ASO therapies
            "fdaApprovedTherapies": fda_therapies.get("fdaApprovedTherapies", []),
            "fdaMessage": fda_therapies.get("message"),
            "targetableExons": fda_therapies.get("targetableExons"),

            # Orphanet / ICD-11 / incidence
            "incidence": orphanet_data.get("incidence"),
            "orphanetCode": orphanet_data.get("orphanetCode"),
            "icd11Code": orphanet_data.get("icd11Code"),
            "orphanetDiseaseNames": orphanet_data.get("diseaseNames", []),

            # Known pathogenic variants and mutation breakdown
            "knownPathogenicVariants": mutation_data.get("knownPathogenicVariants"),
            "totalClinvarVariants": mutation_data.get("totalClinvarVariants"),
            "mutationBreakdown": mutation_data.get("mutationBreakdown", {
                "largeExonDeletions": None,
                "largeExonDuplications": None,
                "nonsensePointMutations": None,
                "frameshiftMutations": None,
                "spliceSiteMutations": None,
            }),
            "clinicalProfileAnnotations": [
                {
                    "description": p if isinstance(p, str) else (p.get("description") or p.get("phenotype") or None),
                    "source": None if isinstance(p, str) else (p.get("source") or None),
                    "id": None if isinstance(p, str) else (p.get("accession") or p.get("id") or p.get("ontology_accession") or p.get("external_reference"))
                }
                for p in (disease_info.get("diseases", []) or [])
            ],

            # On-target pharmacology for target assessment.
            #
            # This block used to flatten 30 "ADMET" fields — absorption,
            # distribution, metabolism, excretion, cell uptake, protein
            # binding, renal clearance, a PBPK time series, a charge/pH
            # profile, Lipinski violations and a 2-D chemical-space
            # projection — all derived from the target gene's CDS. Those
            # endpoints are properties of a finished oligonucleotide's
            # backbone chemistry and formulation, not of a coding sequence,
            # and no candidate oligo exists at this stage. See
            # services/sequence_liability_service.py for the full reasoning.
            #
            # What survives is the part gene context genuinely supports: the
            # consequence of modulating this gene, from gnomAD constraint,
            # essentiality and vital-organ expression.
            "onTargetPharmacology": gene_pharmacology.get("onTargetPharmacology"),
            "onTargetPharmacologyNotAssessed": gene_pharmacology.get("notAssessed"),
        }

        fallback_payload = build_gene_fallback_payload(
            meta=meta,
            official_symbol=official_symbol,
            gene_name=payload_dict.get("geneName") or meta.get("geneName") or official_symbol,
            gene_id=gene_id,
            is_human=is_human,
            enrichment_data=enrichment_data,
            protein_props=protein_props,
            protein_db=protein_db,
            clinical_details=clinical_details,
            disease_resolved=disease_resolved,
        )
        for key, value in fallback_payload.items():
            if value is not None and not payload_dict.get(key):
                payload_dict[key] = value

        return payload_dict
    except HTTPException:
        raise
    except EnsemblLookupUnavailable as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:
        logger.error(f"Error in initialize_target route: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

class AdmetRequest(BaseModel):
    """Request for /api/pipeline/sequence-liabilities.

    Name kept for the deprecated /api/pipeline/admet-prediction alias.
    `transcript_count` is retained for wire compatibility and is unused: it
    fed an off-target risk estimate that never aligned against a genome.
    """

    aso_sequence: str
    gene_symbol: Optional[str] = None
    chemistry: Optional[str] = None
    transcript_count: int = 5


@app.post("/api/pipeline/sequence-liabilities")
@app.post("/api/pipeline/admet-prediction")  # deprecated alias
async def sequence_liabilities(payload: AdmetRequest):
    """Sequence-determined liabilities for a candidate oligonucleotide.

    Was `/api/pipeline/admet-prediction`. The old path still resolves so an
    external caller is not silently broken, but the payload is the honest one:
    innate-immune and structural flags that the base sequence actually
    determines, plus a `notAssessed` map naming every ADMET endpoint that was
    withdrawn and the biological reason for it.
    """
    try:
        gene_context = {}
        if payload.gene_symbol:
            try:
                meta = get_gene_metadata(payload.gene_symbol.strip(), "homo_sapiens")
                if meta and meta.get("id"):
                    constraint_data = get_human_constraint_metrics(payload.gene_symbol.strip())
                    expr_data = get_tissue_expression(
                        symbol=payload.gene_symbol.strip(),
                        ensembl_id=meta.get("id"),
                        species="homo_sapiens",
                    )
                    dep_data = get_gene_dependency(payload.gene_symbol.strip())
                    gene_context = {
                        "vitalOrganTpm": expr_data.get("vital_organ_tpm"),
                        "vitalOrganTissues": expr_data.get("vital_organ_tissues", []),
                        "essentialGene": dep_data.get("essentialGene"),
                        "loeufDecile": constraint_data.get("loeufDecile"),
                    }
            except Exception:
                pass

        return get_sequence_liabilities(
            aso_sequence=payload.aso_sequence,
            chemistry=getattr(payload, "chemistry", None),
            gene_context=gene_context,
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in sequence_liabilities route: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/health")
def health():
    return {"status": "ok"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
