import sys
from typing import Dict, List, Optional

from vcf_parser import parse_header, parse_variants
from vep_client import get_variant_effects_batch
from variant_analysis import calculate_read_statistics, determine_variant_type


def annotate_vcf(vcf_file: str, limit: Optional[int] = None) -> List[Dict]:
    # Parse header to get samples
    _, samples = parse_header(vcf_file)
    
    # Collect all variants first
    variants = []
    for variant in parse_variants(vcf_file, samples):
        if limit and len(variants) >= limit:
            break
        variants.append(variant)
    
    if not variants:
        return []
    
    print(f"Processing {len(variants)} variants...", file=sys.stderr)
    
    variant_tuples = [(v['chrom'], v['pos'], v['ref'], v['alt']) for v in variants]
    vep_results = get_variant_effects_batch(variant_tuples)
    
    # Combine variant data with VEP annotations
    annotations = []
    for variant, vep_data in zip(variants, vep_results):
        stats = calculate_read_statistics(variant)
        variant_type = determine_variant_type(variant['ref'], variant['alt'])
        
        annotation = {
            'chromosome': variant['chrom'],
            'position': variant['pos'],
            'variant_id': variant['id'] if variant['id'] != '.' else f"{variant['chrom']}:{variant['pos']}",
            'reference': variant['ref'],
            'alternate': variant['alt'],
            'quality': variant['qual'],
            'filter': variant['filter'],
            'variant_type': variant_type,
            **stats,
            **vep_data
        }
        annotations.append(annotation)
    
    print(f"Total variants annotated: {len(annotations)}", file=sys.stderr)
    return annotations

