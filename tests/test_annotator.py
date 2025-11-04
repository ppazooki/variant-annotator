import csv
import os
from annotator import annotate_vcf, export_to_tsv, FIELDNAMES


def test_annotate_vcf_basic(tmp_path, mocker):
    """Test annotate_vcf uses batch API and combines results correctly"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs123\tA\tT\t30\tPASS\tDP=100;RO=60;AO=40;AF=0.4\n'
    )
    
    # Mock batch API response
    mock_batch_result = [{
        'gene_id': 'ENSG00000001',
        'gene_symbol': 'TEST',
        'biotype': 'protein_coding',
        'consequence_terms': 'missense_variant',
        'impact': 'MODERATE',
        'strand': 1,
        'rsid': 'N/A',
        'maf': 'N/A'
    }]
    
    mocker.patch('annotator.get_variant_effects_batch', return_value=mock_batch_result)
    mocker.patch('annotator.enrich_with_population_maf', side_effect=lambda x: x)
    
    annotations = annotate_vcf(str(vcf_file))
    
    assert len(annotations) == 1
    assert annotations[0]['chromosome'] == 'chr1'
    assert annotations[0]['position'] == 100
    assert annotations[0]['variant_id'] == 'rs123'
    assert annotations[0]['reference'] == 'A'
    assert annotations[0]['alternate'] == 'T'
    assert annotations[0]['quality'] == '30'
    assert annotations[0]['gene_symbol'] == 'TEST'


def test_annotate_vcf_with_limit(tmp_path, mocker):
    """Test annotate_vcf respects limit parameter"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs1\tA\tT\t30\tPASS\tDP=100\n'
        'chr1\t200\trs2\tG\tC\t30\tPASS\tDP=100\n'
        'chr1\t300\trs3\tT\tA\t30\tPASS\tDP=100\n'
    )
    
    mock_batch_results = [
        {'gene_id': 'ENSG1', 'gene_symbol': 'GENE1', 'biotype': 'pc', 'consequence_terms': 'synonymous', 'impact': 'LOW', 'strand': 1, 'rsid': 'N/A', 'maf': 'N/A'},
        {'gene_id': 'ENSG2', 'gene_symbol': 'GENE2', 'biotype': 'pc', 'consequence_terms': 'synonymous', 'impact': 'LOW', 'strand': 1, 'rsid': 'N/A', 'maf': 'N/A'}
    ]
    
    mocker.patch('annotator.get_variant_effects_batch', return_value=mock_batch_results)
    mocker.patch('annotator.enrich_with_population_maf', side_effect=lambda x: x)
    
    annotations = annotate_vcf(str(vcf_file), limit=2)
    
    assert len(annotations) == 2


def test_annotate_vcf_combines_variant_data(tmp_path, mocker):
    """Test that annotate_vcf correctly combines variant data with VEP annotations"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr2\t500\trs999\tG\tA\t25\tPASS\tDP=200;RO=150;AO=50;AF=0.25\n'
    )
    
    mock_batch_result = [{
        'gene_id': 'ENSG00000002',
        'gene_symbol': 'MYGENE',
        'biotype': 'protein_coding',
        'consequence_terms': 'stop_gained',
        'impact': 'HIGH',
        'strand': -1,
        'rsid': 'N/A',
        'maf': 'N/A'
    }]
    
    mocker.patch('annotator.get_variant_effects_batch', return_value=mock_batch_result)
    mocker.patch('annotator.enrich_with_population_maf', side_effect=lambda x: x)
    
    annotations = annotate_vcf(str(vcf_file))
    
    annotation = annotations[0]
    # Check variant data
    assert annotation['chromosome'] == 'chr2'
    assert annotation['position'] == 500
    assert annotation['variant_id'] == 'rs999'
    assert annotation['reference'] == 'G'
    assert annotation['alternate'] == 'A'
    assert annotation['quality'] == '25'
    # Check variant analysis
    assert annotation['depth'] == 200
    assert annotation['variant_reads'] == 50
    assert annotation['reference_reads'] == 150
    # Check VEP annotations
    assert annotation['gene_id'] == 'ENSG00000002'
    assert annotation['gene_symbol'] == 'MYGENE'
    assert annotation['consequence_terms'] == 'stop_gained'


def test_export_single_annotation(tmp_path):
    annotations = [{
        'chromosome': 'chr1',
        'position': 100,
        'variant_id': 'rs123',
        'reference': 'A',
        'alternate': 'T',
        'quality': '30',
        'variant_type': 'SNP (substitution)',
        'depth': 100,
        'variant_reads': 50,
        'reference_reads': 50,
        'variant_percentage': 50.0,
        'reference_percentage': 50.0,
        'allele_frequency': 0.5,
        'gene_id': 'ENSG00000001',
        'gene_symbol': 'TEST',
        'consequence_terms': 'missense_variant',
        'rsid': 'N/A',
        'maf': '0.0123'
    }]
    
    output_file = tmp_path / "output.tsv"
    export_to_tsv(annotations, str(output_file))
    
    with open(output_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    
    assert len(rows) == 1
    assert rows[0]['chromosome'] == 'chr1'
    assert rows[0]['gene_symbol'] == 'TEST'


def test_export_empty_annotations(tmp_path):
    output_file = tmp_path / "output.tsv"
    export_to_tsv([], str(output_file))
    
    # File should not be created or should be empty
    assert not output_file.exists() or output_file.stat().st_size == 0
