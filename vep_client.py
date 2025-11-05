import sys
import requests
from typing import Dict, List, Tuple


BASE_URL = "https://grch37.rest.ensembl.org"


def build_variant_region(chrom: str, pos: int, alt: str) -> str:
    """Build variant region string for VEP region API."""
    formatted_chrom = chrom[3:] if chrom.startswith('chr') else chrom
    return f"{formatted_chrom}:{pos}-{pos}/{alt}"


def build_hgvs_notation(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build HGVS notation string for VEP batch API."""
    formatted_chrom = chrom[3:] if chrom.startswith('chr') else chrom
    return f"{formatted_chrom}:g.{pos}{ref}>{alt}"


def create_error_response(error_type: str = 'API_ERROR') -> Dict:
    """Create an error response dictionary."""
    return {
        'gene_id': error_type,
        'gene_symbol': error_type,
        'consequence_terms': error_type,
        'rsid': 'N/A',
        'maf': 'N/A'
    }


def parse_vep_response(data: list) -> Dict:
    """Parse VEP API response into annotation dictionary."""
    if not data:
        return {
            'gene_id': 'N/A',
            'gene_symbol': 'N/A',
            'consequence_terms': 'N/A',
            'rsid': 'N/A',
            'maf': 'N/A'
        }
    
    # Extract rsID from colocated_variants
    rsid = 'N/A'
    colocated = data[0].get('colocated_variants', [])
    if colocated:
        rsid = colocated[0].get('id', 'N/A')
    
    consequences = data[0].get('transcript_consequences', [])
    
    if consequences:
        consequence = consequences[0]
        return {
            'gene_id': consequence.get('gene_id', 'N/A'),
            'gene_symbol': consequence.get('gene_symbol', 'N/A'),
            'consequence_terms': ', '.join(consequence.get('consequence_terms', [])),
            'rsid': rsid,
            'maf': 'N/A'  # MAF will be populated by Variation API
        }
    
    return {
        'gene_id': 'N/A',
        'gene_symbol': 'N/A',
        'consequence_terms': 'N/A',
        'rsid': rsid,
        'maf': 'N/A'  # MAF will be populated by Variation API
    }


def handle_vep_error(error_msg: str, variant_region: str, chrom: str, pos: int, ref: str, alt: str) -> Dict:
    """Handle VEP API error messages and return appropriate error response."""
    if "matches reference" in error_msg:
        return create_error_response('REF_MISMATCH')
    else:
        return create_error_response('API_ERROR')


def get_variant_effects(
    chrom: str,
    pos: int,
    ref: str,
    alt: str
) -> Dict:
    """Get variant effects for a single variant using VEP region API."""
    variant_region = build_variant_region(chrom, pos, alt)
    endpoint = f"{BASE_URL}/vep/human/region/{variant_region}"
    params = {
        "content-type": "application/json",
        "assembly": "GRCh37"
    }
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(endpoint, headers=headers, params=params, timeout=10)
        data = response.json()
        
        if "error" in data:
            return handle_vep_error(data["error"], variant_region, chrom, pos, ref, alt)
        
        response.raise_for_status()
        return parse_vep_response(data)
        
    except requests.exceptions.RequestException as e:
        return create_error_response('API_ERROR')


def get_variant_effects_batch(
    variants: List[Tuple[str, int, str, str]],
    batch_size: int = 200
) -> List[Dict]:
    """Get variant effects for multiple variants using VEP batch API."""
    total = len(variants)
    if total == 0:
        return []
    
    print(f"Processing {total} variants (batches of {batch_size})...", file=sys.stderr)
    
    all_results = []
    
    # Build HGVS notation list and mapping
    hgvs_list = []
    hgvs_to_variant_idx = {}
    for idx, (chrom, pos, ref, alt) in enumerate(variants):
        hgvs = build_hgvs_notation(chrom, pos, ref, alt)
        hgvs_list.append(hgvs)
        hgvs_to_variant_idx[hgvs] = idx
    
    # Process in batches
    total_batches = (total + batch_size - 1) // batch_size
    for batch_num in range(0, total, batch_size):
        batch_hgvs = hgvs_list[batch_num:batch_num + batch_size]
        batch_idx = batch_num // batch_size + 1
        
        endpoint = f"{BASE_URL}/vep/human/hgvs"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        data = {"hgvs_notations": batch_hgvs}
        
        try:
            print(f"Batch {batch_idx}/{total_batches}: Processing {len(batch_hgvs)} variants...", file=sys.stderr)
            
            response = requests.post(
                endpoint,
                headers=headers,
                json=data,
                timeout=120
            )
            # Raise HTTPError for 4xx/5xx responses
            response.raise_for_status()
            
            batch_data = response.json()
            
            # Create mapping from HGVS to results
            hgvs_to_result = {}
            for entry in batch_data:
                input_hgvs = entry.get('input', '')
                hgvs_to_result[input_hgvs] = entry
            
            failed_variants = []
            for i, hgvs in enumerate(batch_hgvs):
                variant_idx = batch_num + i
                chrom, pos, ref, alt = variants[variant_idx]
                
                if hgvs in hgvs_to_result:
                    entry = hgvs_to_result[hgvs]
                    result = parse_batch_vep_response(entry)
                    all_results.append(result)
                    
                    # Check if this variant failed (API_ERROR response)
                    if result.get('gene_id') == 'API_ERROR':
                        failed_variants.append((variant_idx, chrom, pos, ref, alt))
                else:
                    all_results.append(create_error_response('API_ERROR'))
                    failed_variants.append((variant_idx, chrom, pos, ref, alt))
            
            # Fall back to single API calls for failed variants (complex variants)
            if failed_variants:
                print(f"  Falling back to individual calls for {len(failed_variants)} failed variants...", file=sys.stderr)
                for variant_idx, chrom, pos, ref, alt in failed_variants:
                    individual_result = get_variant_effects(chrom, pos, ref, alt)
                    all_results[variant_idx] = individual_result
            
            print(f"Batch {batch_idx}/{total_batches} completed", file=sys.stderr)
            
        except requests.exceptions.RequestException as e:
            print(f"Error: Batch API request failed: {e}", file=sys.stderr)
            # Return error responses for this batch
            for hgvs in batch_hgvs:
                all_results.append(create_error_response('API_ERROR'))
    
    print(f"Completed processing {len(all_results)}/{total} variants", file=sys.stderr)
    return all_results


def parse_batch_vep_response(entry: dict) -> Dict:
    """Parse a single entry from VEP batch API response."""
    if "error" in entry:
        return create_error_response('API_ERROR')
    
    # Extract rsID from colocated_variants
    rsid = 'N/A'
    colocated = entry.get('colocated_variants', [])
    if colocated:
        rsid = colocated[0].get('id', 'N/A')
    
    consequences = entry.get('transcript_consequences', [])
    
    if consequences:
        consequence = consequences[0]
        return {
            'gene_id': consequence.get('gene_id', 'N/A'),
            'gene_symbol': consequence.get('gene_symbol', 'N/A'),
            'consequence_terms': ', '.join(consequence.get('consequence_terms', [])),
            'rsid': rsid,
            'maf': 'N/A'  # MAF will be populated by Variation API
        }
    
    return {
        'gene_id': 'N/A',
        'gene_symbol': 'N/A',
        'consequence_terms': 'N/A',
        'rsid': rsid,
        'maf': 'N/A'  # MAF will be populated by Variation API
    }


def fetch_maf_from_variation_api(rsid: str) -> str:
    """Fetch Minor Allele Frequency (MAF) from Ensembl Variation API."""
    if rsid == 'N/A' or not rsid.startswith('rs'):
        return 'N/A'
    
    endpoint = f"{BASE_URL}/variation/human/{rsid}?pops=1"
    headers = {"Accept": "application/json"}
    
    try:
        response = requests.get(endpoint, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        if 'MAF' in data and data['MAF'] is not None:
            return f"{float(data['MAF']):.4f}"
        
        return 'N/A'
        
    except Exception as e:
        return 'N/A'


def enrich_with_population_maf(annotations: List[Dict]) -> List[Dict]:
    """Enrich annotations with MAF data from Ensembl Variation API."""
    variants_with_rsid = [(i, ann) for i, ann in enumerate(annotations) 
                          if ann.get('rsid') and ann.get('rsid') != 'N/A']
    
    if not variants_with_rsid:
        return annotations
    
    print(f"\nFetching MAF from Variation API for {len(variants_with_rsid)} variants with rsIDs...", file=sys.stderr)
    
    for idx, (ann_idx, ann) in enumerate(variants_with_rsid):
        # Skip if MAF already populated
        if ann.get('maf') != 'N/A':
            continue
        
        rsid = ann['rsid']
        maf = fetch_maf_from_variation_api(rsid)
        
        if maf != 'N/A':
            annotations[ann_idx]['maf'] = maf
        
    return annotations
