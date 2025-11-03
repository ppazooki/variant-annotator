from variant_analysis import calculate_read_statistics, determine_variant_type


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


def test_snp_single_base_substitution():
    assert determine_variant_type('A', 'T') == "SNP (substitution)"


def test_insertion():
    assert determine_variant_type('A', 'AT') == "Insertion"


def test_deletion():
    assert determine_variant_type('AT', 'A') == "Deletion"


def test_structural_variant():
    assert determine_variant_type('AAA', '<A>') == "Structural variant"


def test_cnv_with_equal_length():
    assert determine_variant_type('AAA', '[A]') == "CNV"


def test_indel_multi_base():
    assert determine_variant_type('AA', 'TT') == "Indel"

