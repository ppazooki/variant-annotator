from vcf_parser import (
    parse_info,
    parse_genotype,
    parse_variant_line,
    parse_header,
    parse_variants,
    calculate_read_statistics,
    determine_variant_type
)


def test_parse_info():
    assert parse_info('DP=100')['DP'] == 100
    result = parse_info('DP=100;AF=0.5;AC=2')
    assert result['DP'] == 100 and result['AF'] == 0.5 and result['AC'] == 2
    assert parse_info('AF=0.5,0.3;DP=100')['AF'] == 0.5  # Takes first value
    assert parse_info('DB;DP=100')['DB'] is True and parse_info('DB;DP=100')['DP'] == 100
    assert parse_info('GENE=BRCA2;DP=100')['GENE'] == 'BRCA2'


def test_parse_genotype():
    result = parse_genotype('0/1:30:99', ['GT', 'DP', 'GQ'])
    assert result['GT'] == '0/1' and result['DP'] == '30' and result['GQ'] == '99'
    assert parse_genotype('1/1', ['GT'])['GT'] == '1/1'
    result = parse_genotype('./.:.', ['GT', 'DP'])
    assert result['GT'] == './.' and result['DP'] == '.'


def test_parse_variant_line_basic():
    line = 'chr1\t100\trs123\tA\tT\t30\tPASS\tDP=100;AF=0.5'
    samples = []
    result = parse_variant_line(line, samples)
    assert result['chrom'] == 'chr1'
    assert result['pos'] == 100
    assert result['id'] == 'rs123'
    assert result['ref'] == 'A'
    assert result['alt'] == 'T'
    assert result['qual'] == '30'
    assert result['filter'] == 'PASS'
    assert result['info']['DP'] == 100


def test_parse_variant_line_with_samples():
    line = 'chr1\t100\trs123\tA\tT\t30\tPASS\tDP=100\tGT:DP\t0/1:30\t1/1:40'
    samples = ['Sample1', 'Sample2']
    result = parse_variant_line(line, samples)
    assert result['genotypes']['Sample1']['GT'] == '0/1'
    assert result['genotypes']['Sample2']['GT'] == '1/1'


def test_parse_variant_line_no_id():
    line = 'chr1\t100\t.\tA\tT\t.\t.\tDP=100'
    samples = []
    result = parse_variant_line(line, samples)
    assert result['id'] == '.'
    assert result['qual'] == '.'
    assert result['filter'] == '.'


def test_parse_header_with_samples(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '##reference=hg19\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n'
    )
    
    header, samples = parse_header(str(vcf_file))
    
    assert '##fileformat=VCFv4.2' in header
    assert '##reference=hg19' in header
    assert samples == ['Sample1', 'Sample2']


def test_parse_header_without_samples(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )
    
    header, samples = parse_header(str(vcf_file))
    
    assert '##fileformat=VCFv4.2' in header
    assert samples == []


def test_parse_header_only_meta_lines(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '##reference=hg19\n'
    )
    
    header, samples = parse_header(str(vcf_file))
    
    assert len(header) == 2
    assert samples == []


def test_parse_multiple_variants(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\trs1\tA\tT\t30\tPASS\tDP=100\n'
        'chr1\t200\trs2\tG\tC\t40\tPASS\tDP=200\n'
        'chr2\t300\trs3\tT\tA\t50\tPASS\tDP=300\n'
    )
    
    _, samples = parse_header(str(vcf_file))
    variants = list(parse_variants(str(vcf_file), samples))
    
    assert len(variants) == 3
    assert variants[0]['id'] == 'rs1'
    assert variants[1]['info']['DP'] == 200
    assert variants[2]['chrom'] == 'chr2'


def test_parse_variants_with_samples(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n'
        'chr1\t100\trs1\tA\tT\t30\tPASS\tDP=100\tGT:DP\t0/1:30\n'
    )
    
    _, samples = parse_header(str(vcf_file))
    variants = list(parse_variants(str(vcf_file), samples))
    
    assert len(variants) == 1
    assert variants[0]['genotypes']['Sample1']['GT'] == '0/1'


def test_parse_no_variants(tmp_path):
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )
    
    _, samples = parse_header(str(vcf_file))
    variants = list(parse_variants(str(vcf_file), samples))
    
    assert len(variants) == 0


def test_calculate_with_valid_data():
    variant = {
        'info': {
            'DP': 100,
            'RO': 60,
            'AO': 40,
            'AF': 0.4
        }
    }
    
    result = calculate_read_statistics(variant)
    
    assert result['depth'] == 100
    assert result['variant_reads'] == 40
    assert result['reference_reads'] == 60
    assert result['variant_percentage'] == 40.0
    assert result['reference_percentage'] == 60.0
    assert result['allele_frequency'] == 0.4


def test_calculate_with_missing_fields():
    variant = {'info': {}}
    
    result = calculate_read_statistics(variant)
    
    assert result['depth'] == 0
    assert result['variant_reads'] == 0
    assert result['reference_reads'] == 0


def test_determine_variant_type():
    assert determine_variant_type('A', 'T') == "SNP (substitution)"
    assert determine_variant_type('A', 'AT') == "Insertion"
    assert determine_variant_type('AT', 'A') == "Deletion"
    assert determine_variant_type('AAA', '<A>') == "Structural variant"
    assert determine_variant_type('AAA', '[A]') == "CNV"
    assert determine_variant_type('AA', 'TT') == "Indel"
