import asyncio
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as main_module
from services.enrichment_service import get_aso_analysis


@pytest.fixture
def client():
    with TestClient(main_module.app) as test_client:
        yield test_client


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Specifies a gene_metadata fast path that this codebase has never "
        "had -- no commit in the history adds the kwarg. Two of its four "
        "assertions cannot be satisfied honestly: structuralAccessibility is "
        "derived from the CDS sequence GC content "
        "(_compute_aso_metrics_from_sequence), and gene_metadata carries only "
        "counts and coordinates, so with the network forbidden there is no "
        "sequence to compute it from; and spliceSwitches would have to become "
        "totalTranscripts-1 here while the network path defines it as "
        "(distinct exon counts)-1, giving one field two meanings depending on "
        "which path ran. Implementing it means inventing both numbers. "
        "Left failing on purpose -- see docs/planning/model_training_results.md."
    ),
)
def test_get_aso_analysis_uses_gene_metadata_without_network(monkeypatch):
    def fail_ensembl_get(*args, **kwargs):
        raise AssertionError("Network lookup should not be attempted when gene metadata is present")

    monkeypatch.setattr("services.enrichment_service._ensembl_get", fail_ensembl_get)

    result = get_aso_analysis("ENSG00000198947", 9606, gene_metadata={
        "totalTranscripts": 95,
        "exonCount": 79,
        "start": 31097677,
        "end": 33339609,
    })

    assert result["activeIsoforms"] == 95
    assert result["spliceSwitches"] == 94
    assert result["structuralAccessibility"]
    assert result["transcriptSpecificity"] == "95 isoforms (Low)"


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Asks the endpoint to emit activeIsoforms=1, spliceSwitches=0 and a "
        "truthy structuralAccessibility when the ASO analysis TIMES OUT. "
        "Those would be fabricated values presented in the same fields as "
        "measured ones, with nothing marking them as a fallback -- the same "
        "defect class as the hash()-derived binding affinities and the "
        "always-True feature placeholders removed elsewhere in this review. "
        "On timeout _run_sync correctly returns {} and the fields stay None. "
        "Left failing on purpose -- see docs/planning/model_training_results.md."
    ),
)
def test_initialize_target_uses_generic_aso_fallback_for_any_gene_when_analysis_times_out(monkeypatch, client):
    def fake_get_gene_metadata(symbol, species):
        return {
            "id": "ENSG00000141510",
            "officialSymbol": "TP53",
            "geneName": "tumor protein p53",
            "biotype": "protein_coding",
            "strand": "+",
            "start": 7565097,
            "end": 7590856,
            "exonCount": 11,
            "proteinLength": 393,
            "canonicalTranscript": "ENST00000269305",
            "otherTranscripts": [],
            "proteinId": "ENSP00000269305",
            "synonyms": ["P53"],
            "entrezGeneId": "7157",
            "hgncId": "HGNC:11998",
            "nomenclatureId": "HGNC:11998",
        }

    async def fake_get_pubmed_count(*args, **kwargs):
        return 0

    async def fake_get_rxiv_count(*args, **kwargs):
        return 0

    async def fake_fetch_disease_associations(*args, **kwargs):
        return {"diseases": [], "omim_id": None, "source": []}

    async def fake_fetch_expression_details(*args, **kwargs):
        return {
            "available": False,
            "top_tissue": None,
            "tpm": None,
            "top_tissues": [],
            "cellType": None,
            "cellTpm": None,
            "expression_cv": None,
            "vital_organ_tpm": None,
            "vital_organ_tissues": [],
            "dominant_isoform_fraction": None,
            "dominant_isoform_id": None,
            "hpa_level": None,
            "gtex_level": None,
            "source": None,
            "expressionUnit": "TPM",
        }

    async def fake_get_clinvar_count(*args, **kwargs):
        return 0

    async def fake_get_dbsnp_count(*args, **kwargs):
        return 1212

    async def fake_fetch_ncbi_aliases(*args, **kwargs):
        return []

    def fake_get_gene_enrichment(*args, **kwargs):
        return {}

    def fake_get_human_constraint_metrics(*args, **kwargs):
        return {}

    def fake_get_aso_analysis(*args, **kwargs):
        raise TimeoutError("simulated timeout")

    def fake_get_rna_halflife(*args, **kwargs):
        return {}

    def fake_get_gene_dependency(*args, **kwargs):
        return {}

    def fake_get_variant_details(*args, **kwargs):
        return {}

    def fake_get_single_cell_expression(*args, **kwargs):
        return {}

    def fake_get_clinical_details(*args, **kwargs):
        return {}

    def fake_get_fda_therapies(*args, **kwargs):
        return {}

    def fake_get_orphanet_data(*args, **kwargs):
        return {}

    def fake_get_mutation_breakdown(*args, **kwargs):
        return {}

    def fake_get_protein_db_ids(*args, **kwargs):
        return {}

    def fake_get_protein_properties(*args, **kwargs):
        return {}

    monkeypatch.setattr(main_module, "get_gene_metadata", fake_get_gene_metadata)
    monkeypatch.setattr(main_module, "get_pubmed_count", fake_get_pubmed_count)
    monkeypatch.setattr(main_module, "get_rxiv_count", fake_get_rxiv_count)
    monkeypatch.setattr(main_module, "fetch_disease_associations", fake_fetch_disease_associations)
    monkeypatch.setattr(main_module, "fetch_expression_details", fake_fetch_expression_details)
    monkeypatch.setattr(main_module, "get_clinvar_count", fake_get_clinvar_count)
    monkeypatch.setattr(main_module, "get_dbsnp_count", fake_get_dbsnp_count)
    monkeypatch.setattr(main_module, "fetch_ncbi_aliases", fake_fetch_ncbi_aliases)
    monkeypatch.setattr(main_module, "get_gene_enrichment", fake_get_gene_enrichment)
    monkeypatch.setattr(main_module, "get_human_constraint_metrics", fake_get_human_constraint_metrics)
    monkeypatch.setattr(main_module, "get_aso_analysis", fake_get_aso_analysis)
    monkeypatch.setattr(main_module, "get_rna_halflife", fake_get_rna_halflife)
    monkeypatch.setattr(main_module, "get_gene_dependency", fake_get_gene_dependency)
    monkeypatch.setattr(main_module, "get_variant_details", fake_get_variant_details)
    monkeypatch.setattr(main_module, "get_single_cell_expression", fake_get_single_cell_expression)
    monkeypatch.setattr(main_module, "get_clinical_details", fake_get_clinical_details)
    monkeypatch.setattr(main_module, "get_fda_therapies", fake_get_fda_therapies)
    monkeypatch.setattr(main_module, "get_orphanet_data", fake_get_orphanet_data)
    monkeypatch.setattr(main_module, "get_mutation_breakdown", fake_get_mutation_breakdown)
    monkeypatch.setattr(main_module, "get_protein_db_ids", fake_get_protein_db_ids)
    monkeypatch.setattr(main_module, "get_protein_properties", fake_get_protein_properties)
    monkeypatch.setattr(main_module, "_add_notification", lambda *args, **kwargs: None)
    monkeypatch.setattr(main_module, "build_gene_fallback_payload", lambda **kwargs: {})
    monkeypatch.setattr(main_module, "get_safe_ensembl_url", lambda species, gene_id: "https://example.test/ensembl")

    response = client.post(
        "/api/pipeline/initialize-target",
        json={
            "gene_symbol": "TP53",
            "organism": "Human (Homo sapiens)",
            "disease_name": "Cancer",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["activeIsoforms"] == 1
    assert payload["spliceSwitches"] == 0
    assert payload["structuralAccessibility"]


def test_initialize_target_preserves_aso_metrics_when_analysis_is_slow(monkeypatch, client):
    def fake_get_gene_metadata(symbol, species):
        return {
            "id": "ENSG00000198947",
            "officialSymbol": "DMD",
            "geneName": "dystrophin",
            "biotype": "protein_coding",
            "strand": "+",
            "start": 31097677,
            "end": 33339609,
            "exonCount": 79,
            "proteinLength": 3685,
            "canonicalTranscript": "ENST00000357033",
            "otherTranscripts": [],
            "proteinId": "ENSP00000354923",
            "synonyms": ["BMD"],
            "entrezGeneId": "1756",
            "hgncId": "HGNC:2928",
            "nomenclatureId": "HGNC:2928",
        }

    async def fake_get_pubmed_count(*args, **kwargs):
        return 0

    async def fake_get_rxiv_count(*args, **kwargs):
        return 0

    async def fake_fetch_disease_associations(*args, **kwargs):
        return {"diseases": [], "omim_id": None, "source": []}

    async def fake_fetch_expression_details(*args, **kwargs):
        return {
            "available": False,
            "top_tissue": None,
            "tpm": None,
            "top_tissues": [],
            "cellType": None,
            "cellTpm": None,
            "expression_cv": None,
            "vital_organ_tpm": None,
            "vital_organ_tissues": [],
            "dominant_isoform_fraction": None,
            "dominant_isoform_id": None,
            "hpa_level": None,
            "gtex_level": None,
            "source": None,
            "expressionUnit": "TPM",
        }

    async def fake_get_clinvar_count(*args, **kwargs):
        return 0

    async def fake_get_dbsnp_count(*args, **kwargs):
        return 754712

    async def fake_fetch_ncbi_aliases(*args, **kwargs):
        return []

    def fake_get_gene_enrichment(*args, **kwargs):
        return {}

    def fake_get_human_constraint_metrics(*args, **kwargs):
        return {}

    def fake_get_aso_analysis(*args, **kwargs):
        return {
            "activeIsoforms": 40,
            "spliceSwitches": 30,
            "structuralAccessibility": "76% (Favorable)",
            "splicingMotifDensity": "65.3/kb (High)",
            "preclinicalConservation": "3/3 (Excellent)",
            "gQuadruplexes": "0 Blocks Found",
            "cpgDensity": "Medium Risk",
            "selfDimerRisk": "58/kb (High)",
            "polygTracts": "11 (High Risk)",
            "transcriptSpecificity": "40 isoforms (Low)",
            "codonUsageBias": "GC3=45.7% (Balanced)",
        }

    def fake_get_rna_halflife(*args, **kwargs):
        return {}

    def fake_get_gene_dependency(*args, **kwargs):
        return {}

    def fake_get_variant_details(*args, **kwargs):
        return {}

    def fake_get_single_cell_expression(*args, **kwargs):
        return {}

    def fake_get_clinical_details(*args, **kwargs):
        return {}

    def fake_get_fda_therapies(*args, **kwargs):
        return {}

    def fake_get_orphanet_data(*args, **kwargs):
        return {}

    def fake_get_mutation_breakdown(*args, **kwargs):
        return {}

    def fake_get_protein_db_ids(*args, **kwargs):
        return {}

    def fake_get_protein_properties(*args, **kwargs):
        return {}

    monkeypatch.setattr(main_module, "get_gene_metadata", fake_get_gene_metadata)
    monkeypatch.setattr(main_module, "get_pubmed_count", fake_get_pubmed_count)
    monkeypatch.setattr(main_module, "get_rxiv_count", fake_get_rxiv_count)
    monkeypatch.setattr(main_module, "fetch_disease_associations", fake_fetch_disease_associations)
    monkeypatch.setattr(main_module, "fetch_expression_details", fake_fetch_expression_details)
    monkeypatch.setattr(main_module, "get_clinvar_count", fake_get_clinvar_count)
    monkeypatch.setattr(main_module, "get_dbsnp_count", fake_get_dbsnp_count)
    monkeypatch.setattr(main_module, "fetch_ncbi_aliases", fake_fetch_ncbi_aliases)
    monkeypatch.setattr(main_module, "get_gene_enrichment", fake_get_gene_enrichment)
    monkeypatch.setattr(main_module, "get_human_constraint_metrics", fake_get_human_constraint_metrics)
    monkeypatch.setattr(main_module, "get_aso_analysis", fake_get_aso_analysis)
    monkeypatch.setattr(main_module, "get_rna_halflife", fake_get_rna_halflife)
    monkeypatch.setattr(main_module, "get_gene_dependency", fake_get_gene_dependency)
    monkeypatch.setattr(main_module, "get_variant_details", fake_get_variant_details)
    monkeypatch.setattr(main_module, "get_single_cell_expression", fake_get_single_cell_expression)
    monkeypatch.setattr(main_module, "get_clinical_details", fake_get_clinical_details)
    monkeypatch.setattr(main_module, "get_fda_therapies", fake_get_fda_therapies)
    monkeypatch.setattr(main_module, "get_orphanet_data", fake_get_orphanet_data)
    monkeypatch.setattr(main_module, "get_mutation_breakdown", fake_get_mutation_breakdown)
    monkeypatch.setattr(main_module, "get_protein_db_ids", fake_get_protein_db_ids)
    monkeypatch.setattr(main_module, "get_protein_properties", fake_get_protein_properties)
    monkeypatch.setattr(main_module, "_add_notification", lambda *args, **kwargs: None)
    monkeypatch.setattr(main_module, "build_gene_fallback_payload", lambda **kwargs: {})
    monkeypatch.setattr(main_module, "get_safe_ensembl_url", lambda species, gene_id: "https://example.test/ensembl")

    response = client.post(
        "/api/pipeline/initialize-target",
        json={
            "gene_symbol": "DMD",
            "organism": "Human (Homo sapiens)",
            "disease_name": "Duchenne Muscular Dystrophy",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["activeIsoforms"] == 40
    assert payload["structuralAccessibility"] == "76% (Favorable)"
    assert payload["dbSnpCount"] == 754712
