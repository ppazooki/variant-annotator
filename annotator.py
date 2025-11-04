import csv
import sys
from typing import Dict, List, Optional

from vcf_parser import parse_header, parse_variants, calculate_read_statistics, determine_variant_type
from vep_client import get_variant_effects_batch, enrich_with_population_maf


FIELDNAMES = [
    'chromosome', 'position', 'variant_id', 'reference', 'alternate',
    'quality', 'filter', 'variant_type',
    'depth', 'variant_reads', 'reference_reads', 
    'variant_percentage', 'reference_percentage', 'allele_frequency',
    'gene_id', 'gene_symbol', 'biotype', 'consequence_terms', 
    'impact', 'strand', 'rsid', 'maf'
]


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
        
        # Use VCF ID column as rsID if it starts with 'rs', otherwise use VEP rsID
        vcf_id = variant['id'] if variant['id'] != '.' else None
        rsid_from_vep = vep_data.get('rsid', 'N/A')
        final_rsid = rsid_from_vep
        if vcf_id and vcf_id.startswith('rs') and rsid_from_vep == 'N/A':
            final_rsid = vcf_id
        
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
            **vep_data,
            'rsid': final_rsid  # Override with VCF ID if available
        }
        annotations.append(annotation)
    
    print(f"Total variants annotated: {len(annotations)}", file=sys.stderr)
    
    # Enrich with MAF from Variation API for variants with rsIDs
    annotations = enrich_with_population_maf(annotations)
    
    return annotations


def export_to_tsv(annotations: List[Dict], output_file: str) -> None:
    if not annotations:
        print("No annotations to export", file=sys.stderr)
        return
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t')
        writer.writeheader()
        writer.writerows(annotations)
    
    print(f"Annotations exported to {output_file}", file=sys.stderr)

