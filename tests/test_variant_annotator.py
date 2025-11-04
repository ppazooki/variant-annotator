import sys
from unittest.mock import patch, MagicMock
from variant_annotator import main


def test_main_with_basic_arguments(tmp_path):
    """Test main function with basic VCF file"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs123\tA\tT\t30\tPASS\tDP=100;RO=60;AO=40;AF=0.4\n'
    )
    
    tsv_file = tmp_path / "output.tsv"
    
    test_args = ['variant_annotator.py', str(vcf_file), '--output', str(tsv_file), '--limit', '0']
    with patch('annotator.enrich_with_population_maf', side_effect=lambda x: x):
        with patch.object(sys, 'argv', test_args):
            main()
    
    # Check that output file exists (even if empty due to limit=0)
    assert tsv_file.exists()


def test_main_with_limit(tmp_path):
    """Test main function with limit parameter"""
    vcf_file = tmp_path / "test.vcf"
    lines = ['##fileformat=VCFv4.2\n', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n']
    for i in range(20):
        lines.append(f'chr1\t{100+i}\t.\tA\tT\t.\t.\tDP=100\n')
    vcf_file.write_text(''.join(lines))
    
    tsv_file = tmp_path / "output.tsv"
    
    test_args = ['variant_annotator.py', str(vcf_file), '--output', str(tsv_file), '--limit', '5']
    with patch('annotator.enrich_with_population_maf', side_effect=lambda x: x):
        with patch.object(sys, 'argv', test_args):
            main()
    
    # Verify output file was created
    assert tsv_file.exists()
