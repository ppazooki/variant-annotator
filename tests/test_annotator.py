import os
from annotator import annotate_variant, annotate_vcf


def test_annotate_basic_variant(mocker):
    variant = {
        'chrom': 'chr1',
        'pos': 100,
        'id': 'rs123',
        'ref': 'A',
        'alt': 'T',
        'qual': '30',
        'filter': 'PASS',
        'info': {
            'DP': 100,
            'RO': 60,
            'AO': 40,
            'AF': 0.4
        }
    }
    
    mocker.patch('annotator.get_variant_effects', return_value={
        'gene_id': 'ENSG00000001',
        'gene_symbol': 'TEST',
        'biotype': 'protein_coding',
        'consequence_terms': 'missense_variant',
        'impact': 'MODERATE',
        'strand': 1
    })
    
    result = annotate_variant(variant)
    
    assert result['chromosome'] == 'chr1'
    assert result['position'] == 100
    assert result['variant_id'] == 'rs123'
    assert result['reference'] == 'A'
    assert result['alternate'] == 'T'
    assert result['variant_type'] == 'SNP (substitution)'
    assert result['depth'] == 100
    assert result['variant_reads'] == 40
    assert result['reference_reads'] == 60
    assert result['gene_id'] == 'ENSG00000001'
    assert result['gene_symbol'] == 'TEST'


def test_annotate_insertion(mocker):
    variant = {
        'chrom': 'chr2',
        'pos': 200,
        'id': '.',
        'ref': 'A',
        'alt': 'AT',
        'qual': '.',
        'filter': '.',
        'info': {}
    }
    
    mocker.patch('annotator.get_variant_effects', return_value={
        'gene_id': 'ENSG00000002',
        'gene_symbol': 'TEST2',
        'biotype': 'protein_coding',
        'consequence_terms': 'frameshift_variant',
        'impact': 'HIGH',
        'strand': -1
    })
    
    result = annotate_variant(variant)
    
    assert result['variant_type'] == 'Insertion'
    assert result['impact'] == 'HIGH'


def test_annotate_deletion(mocker):
    variant = {
        'chrom': 'chr3',
        'pos': 300,
        'id': 'rs456',
        'ref': 'AT',
        'alt': 'A',
        'qual': '40',
        'filter': 'PASS',
        'info': {
            'DP': 200,
            'RO': 100,
            'AO': 100,
            'AF': 0.5
        }
    }
    
    mocker.patch('annotator.get_variant_effects', return_value={
        'gene_id': 'ENSG00000003',
        'gene_symbol': 'TEST3',
        'biotype': 'protein_coding',
        'consequence_terms': 'frameshift_variant',
        'impact': 'HIGH',
        'strand': 1
    })
    
    result = annotate_variant(variant)
    
    assert result['variant_type'] == 'Deletion'
    assert result['depth'] == 200


def test_annotate_vcf_basic(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs123\tA\tT\t30\tPASS\tDP=100;RO=60;AO=40;AF=0.4\n'
    )
    
    annotations = annotate_vcf(str(vcf_file))
    
    assert len(annotations) == 1


def test_annotate_vcf_with_limit(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs1\tA\tT\t30\tPASS\tDP=100\n'
        'chr1\t200\trs2\tG\tC\t30\tPASS\tDP=100\n'
        'chr1\t300\trs3\tT\tA\t30\tPASS\tDP=100\n'
    )
    
    annotations = annotate_vcf(str(vcf_file), limit=2)
    
    assert len(annotations) == 2
