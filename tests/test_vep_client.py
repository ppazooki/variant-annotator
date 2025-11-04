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


def test_build_variant_region_with_chr_prefix():
    result = build_variant_region('chr1', 100, 'A')
    assert result == '1:100-100/A'
    
def test_build_variant_region_without_chr_prefix():
    result = build_variant_region('2', 200, 'G')
    assert result == '2:200-200/G'


def test_build_hgvs_notation_with_chr_prefix():
    result = build_hgvs_notation('chr1', 100, 'G', 'A')
    assert result == '1:g.100G>A'


def test_build_hgvs_notation_without_chr_prefix():
    result = build_hgvs_notation('2', 200, 'C', 'T')
    assert result == '2:g.200C>T'


def test_create_error_response():
    result = create_error_response('API_ERROR')
    assert result['gene_id'] == 'API_ERROR'
    assert result['gene_symbol'] == 'API_ERROR'


def test_parse_empty_response():
    result = parse_vep_response([])
    assert result['gene_id'] == 'N/A'
    assert result['gene_symbol'] == 'N/A'
    
def test_parse_valid_response():
    data = [{
        'transcript_consequences': [{
            'gene_id': 'ENSG00000139618',
            'gene_symbol': 'BRCA2',
            'biotype': 'protein_coding',
            'consequence_terms': ['missense_variant', 'splice_region_variant'],
            'impact': 'MODERATE',
            'strand': 1
        }]
    }]
    result = parse_vep_response(data)
    assert result['gene_id'] == 'ENSG00000139618'
    assert result['gene_symbol'] == 'BRCA2'
    assert result['consequence_terms'] == 'missense_variant, splice_region_variant'
    
def test_parse_response_without_consequences():
    data = [{'transcript_consequences': []}]
    result = parse_vep_response(data)
    assert result['gene_id'] == 'N/A'


def test_handle_reference_mismatch_error():
    result = handle_vep_error(
        "Input reference allele matches reference",
        "1:100-100/A",
        "chr1", 100, "G", "A"
    )
    assert result['gene_id'] == 'REF_MISMATCH'
    
def test_handle_generic_api_error():
    result = handle_vep_error(
        "Some other error",
        "1:100-100/A",
        "chr1", 100, "G", "A"
    )
    assert result['gene_id'] == 'API_ERROR'


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
    assert result['impact'] == 'MODERATE'


def test_parse_batch_vep_response_with_error():
    entry = {'input': '1:g.100G>A', 'error': 'Some error'}
    result = parse_batch_vep_response(entry)
    assert result['gene_id'] == 'API_ERROR'


def test_parse_batch_vep_response_empty_consequences():
    entry = {
        'input': '1:g.100G>A',
        'transcript_consequences': []
    }
    result = parse_batch_vep_response(entry)
    assert result['gene_id'] == 'N/A'
    assert result['gene_symbol'] == 'N/A'


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
