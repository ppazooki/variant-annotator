from vcf_parser import (
    parse_info,
    parse_genotype,
    parse_variant_line,
    parse_header,
    parse_variants
)


def test_parse_info_single_value():
    info_str = 'DP=100'
    result = parse_info(info_str)
    assert result['DP'] == 100


def test_parse_info_multiple_values():
    info_str = 'DP=100;AF=0.5;AC=2'
    result = parse_info(info_str)
    assert result['DP'] == 100
    assert result['AF'] == 0.5
    assert result['AC'] == 2


def test_parse_info_with_comma_separated_values():
    info_str = 'AF=0.5,0.3;DP=100'
    result = parse_info(info_str)
    assert result['AF'] == 0.5  # Takes first value
    assert result['DP'] == 100


def test_parse_info_boolean_flag():
    info_str = 'DB;DP=100'
    result = parse_info(info_str)
    assert result['DB'] is True
    assert result['DP'] == 100


def test_parse_info_with_string_value():
    info_str = 'GENE=BRCA2;DP=100'
    result = parse_info(info_str)
    assert result['GENE'] == 'BRCA2'
    assert result['DP'] == 100


def test_parse_genotype_basic():
    format_fields = ['GT', 'DP', 'GQ']
    genotype_str = '0/1:30:99'
    result = parse_genotype(genotype_str, format_fields)
    assert result['GT'] == '0/1'
    assert result['DP'] == '30'
    assert result['GQ'] == '99'


def test_parse_genotype_homozygous():
    format_fields = ['GT']
    genotype_str = '1/1'
    result = parse_genotype(genotype_str, format_fields)
    assert result['GT'] == '1/1'


def test_parse_genotype_missing_value():
    format_fields = ['GT', 'DP']
    genotype_str = './.:.'
    result = parse_genotype(genotype_str, format_fields)
    assert result['GT'] == './.'
    assert result['DP'] == '.'


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
