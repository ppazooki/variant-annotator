import pytest
import responses
from vep_client import (
    build_variant_region,
    build_hgvs_notation,
    create_error_response,
    parse_vep_response,
    parse_batch_vep_response,
    handle_vep_error,
    get_variant_effects,
    get_variant_effects_batch
)


def test_build_variant_region():
    assert build_variant_region('chr1', 100, 'A') == '1:100-100/A'
    assert build_variant_region('2', 200, 'G') == '2:200-200/G'


def test_build_hgvs_notation():
    assert build_hgvs_notation('chr1', 100, 'G', 'A') == '1:g.100G>A'
    assert build_hgvs_notation('2', 200, 'C', 'T') == '2:g.200C>T'


def test_create_error_response():
    result = create_error_response('API_ERROR')
    assert result['gene_id'] == 'API_ERROR'
    assert result['gene_symbol'] == 'API_ERROR'
    assert result['rsid'] == 'N/A'
    assert result['maf'] == 'N/A'


def test_parse_empty_response():
    result = parse_vep_response([])
    assert result['gene_id'] == 'N/A' and result['rsid'] == 'N/A'
    
def test_parse_valid_response():
    data = [{
        'transcript_consequences': [{
            'gene_id': 'ENSG00000139618',
            'gene_symbol': 'BRCA2',
            'biotype': 'protein_coding',
            'consequence_terms': ['missense_variant', 'splice_region_variant'],
            'impact': 'MODERATE',
            'strand': 1
        }],
        'colocated_variants': [{
            'id': 'rs123456',
            'frequencies': {
                'gnomad': 0.1234
            }
        }]
    }]
    result = parse_vep_response(data)
    assert result['gene_id'] == 'ENSG00000139618'
    assert result['gene_symbol'] == 'BRCA2'
    assert result['consequence_terms'] == 'missense_variant, splice_region_variant'
    assert result['rsid'] == 'rs123456'
    assert result['maf'] == 'N/A'  # MAF extraction removed, will be populated by Variation API
    
def test_parse_response_without_consequences():
    result = parse_vep_response([{'transcript_consequences': []}])
    assert result['gene_id'] == 'N/A' and result['rsid'] == 'N/A'


def test_handle_vep_error():
    assert handle_vep_error("Input reference allele matches reference", "1:100-100/A", "chr1", 100, "G", "A")['gene_id'] == 'REF_MISMATCH'
    assert handle_vep_error("Some other error", "1:100-100/A", "chr1", 100, "G", "A")['gene_id'] == 'API_ERROR'


@responses.activate
def test_successful_api_call():
    responses.add(
        responses.GET,
        'https://grch37.rest.ensembl.org/vep/human/region/1:100-100/A',
        json=[{
            'transcript_consequences': [{
                'gene_id': 'ENSG00000001',
                'gene_symbol': 'TEST',
                'biotype': 'protein_coding',
                'consequence_terms': ['missense_variant'],
                'impact': 'MODERATE',
                'strand': 1
            }]
        }],
        status=200
    )
    
    result = get_variant_effects('chr1', 100, 'G', 'A')
    assert result['gene_id'] == 'ENSG00000001'
    assert result['gene_symbol'] == 'TEST'
    
@responses.activate
def test_api_call_with_error_response():
    responses.add(
        responses.GET,
        'https://grch37.rest.ensembl.org/vep/human/region/1:100-100/A',
        json={'error': 'Input reference allele matches reference'},
        status=400
    )
    
    result = get_variant_effects('chr1', 100, 'G', 'A')
    assert result['gene_id'] == 'REF_MISMATCH'
    
def test_api_call_with_network_error(mocker):
    import requests
    
    mock_get = mocker.patch('requests.get')
    mock_get.side_effect = requests.exceptions.RequestException('Network error')
    
    result = get_variant_effects('chr1', 100, 'G', 'A')
    assert result['gene_id'] == 'API_ERROR'


def test_parse_batch_vep_response_valid():
    entry = {
        'input': '1:g.100G>A',
        'transcript_consequences': [{
            'gene_id': 'ENSG00000001',
            'gene_symbol': 'TEST',
            'biotype': 'protein_coding',
            'consequence_terms': ['missense_variant'],
            'impact': 'MODERATE',
            'strand': 1
        }]
    }
    result = parse_batch_vep_response(entry)
    assert result['gene_id'] == 'ENSG00000001'
    assert result['gene_symbol'] == 'TEST'
    assert result['consequence_terms'] == 'missense_variant'


def test_parse_batch_vep_response():
    assert parse_batch_vep_response({'input': '1:g.100G>A', 'error': 'Some error'})['gene_id'] == 'API_ERROR'
    result = parse_batch_vep_response({'input': '1:g.100G>A', 'transcript_consequences': []})
    assert result['gene_id'] == 'N/A' and result['gene_symbol'] == 'N/A'


@responses.activate
def test_get_variant_effects_batch_success():
    """Test batch API call with successful response"""
    responses.add(
        responses.POST,
        'https://grch37.rest.ensembl.org/vep/human/hgvs',
        json=[{
            'input': '1:g.100G>A',
            'transcript_consequences': [{
                'gene_id': 'ENSG00000001',
                'gene_symbol': 'TEST1',
                'biotype': 'protein_coding',
                'consequence_terms': ['missense_variant'],
                'impact': 'MODERATE',
                'strand': 1
            }]
        }, {
            'input': '2:g.200C>T',
            'transcript_consequences': [{
                'gene_id': 'ENSG00000002',
                'gene_symbol': 'TEST2',
                'biotype': 'protein_coding',
                'consequence_terms': ['synonymous_variant'],
                'impact': 'LOW',
                'strand': -1
            }]
        }],
        status=200
    )
    
    variants = [('chr1', 100, 'G', 'A'), ('chr2', 200, 'C', 'T')]
    results = get_variant_effects_batch(variants)
    
    assert len(results) == 2
    assert results[0]['gene_symbol'] == 'TEST1'
    assert results[1]['gene_symbol'] == 'TEST2'


@responses.activate
def test_get_variant_effects_batch_with_fallback(mocker):
    """Test batch API falls back to individual calls for failed variants"""
    # Batch response with one success and one error
    responses.add(
        responses.POST,
        'https://grch37.rest.ensembl.org/vep/human/hgvs',
        json=[{
            'input': '1:g.100G>A',
            'transcript_consequences': [{
                'gene_id': 'ENSG00000001',
                'gene_symbol': 'TEST1',
                'biotype': 'protein_coding',
                'consequence_terms': ['missense_variant'],
                'impact': 'MODERATE',
                'strand': 1
            }]
        }, {
            'input': '2:g.200ATTTT>GTTTC',
            'error': 'Unable to parse HGVS notation'
        }],
        status=200
    )
    
    # Mock individual call for the failed variant
    mock_get_individual = mocker.patch('vep_client.get_variant_effects')
    mock_get_individual.return_value = {
        'gene_id': 'ENSG00000002',
        'gene_symbol': 'TEST2',
        'biotype': 'protein_coding',
        'consequence_terms': 'frameshift_variant',
        'impact': 'HIGH',
        'strand': 1
    }
    
    variants = [('chr1', 100, 'G', 'A'), ('chr2', 200, 'ATTTT', 'GTTTC')]
    results = get_variant_effects_batch(variants)
    
    assert len(results) == 2
    assert results[0]['gene_symbol'] == 'TEST1'
    assert results[1]['gene_symbol'] == 'TEST2'  # From fallback
    mock_get_individual.assert_called_once()
