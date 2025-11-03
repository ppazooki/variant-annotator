import csv
from exporter import export_to_tsv, FIELDNAMES


def test_export_single_annotation(tmp_path):
    annotations = [{
        'chromosome': 'chr1',
        'position': 100,
        'variant_id': 'rs123',
        'reference': 'A',
        'alternate': 'T',
        'quality': '30',
        'filter': 'PASS',
        'variant_type': 'SNP (substitution)',
        'depth': 100,
        'variant_reads': 50,
        'reference_reads': 50,
        'variant_percentage': 50.0,
        'reference_percentage': 50.0,
        'allele_frequency': 0.5,
        'gene_id': 'ENSG00000001',
        'gene_symbol': 'TEST',
        'biotype': 'protein_coding',
        'consequence_terms': 'missense_variant',
        'impact': 'MODERATE',
        'strand': 1
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
