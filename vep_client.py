import sys
import requests
from typing import Dict


BASE_URL = "https://rest.ensembl.org"


def build_variant_region(chrom: str, pos: int, alt: str) -> str:
    formatted_chrom = chrom[3:] if chrom.startswith('chr') else chrom
    return f"{formatted_chrom}:{pos}-{pos}/{alt}"


def create_error_response(error_type: str = 'API_ERROR') -> Dict:
    return {
        'gene_id': error_type,
        'gene_symbol': error_type,
        'biotype': error_type,
        'consequence_terms': error_type,
        'impact': error_type,
        'strand': error_type
    }


def parse_vep_response(data: list) -> Dict:
    if not data:
        return {
            'gene_id': 'N/A',
            'gene_symbol': 'N/A',
            'biotype': 'N/A',
            'consequence_terms': 'N/A',
            'impact': 'N/A',
            'strand': 'N/A'
        }
    
    consequences = data[0].get('transcript_consequences', [])
    
    if consequences:
        consequence = consequences[0]
        return {
            'gene_id': consequence.get('gene_id', 'N/A'),
            'gene_symbol': consequence.get('gene_symbol', 'N/A'),
            'biotype': consequence.get('biotype', 'N/A'),
            'consequence_terms': ', '.join(consequence.get('consequence_terms', [])),
            'impact': consequence.get('impact', 'N/A'),
            'strand': consequence.get('strand', 'N/A')
        }
    
    return {
        'gene_id': 'N/A',
        'gene_symbol': 'N/A',
        'biotype': 'N/A',
        'consequence_terms': 'N/A',
        'impact': 'N/A',
        'strand': 'N/A'
    }


def handle_vep_error(error_msg: str, variant_region: str, chrom: str, pos: int, ref: str, alt: str) -> Dict:
    if "matches reference" in error_msg:
        print(f"Warning: Reference mismatch for {chrom}:{pos} ref={ref} alt={alt}: {error_msg}", file=sys.stderr)
        return create_error_response('REF_MISMATCH')
    else:
        print(f"Warning: VEP API error for {variant_region}: {error_msg}", file=sys.stderr)
        return create_error_response('API_ERROR')


def get_variant_effects(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    species: str = "human"
) -> Dict:
    variant_region = build_variant_region(chrom, pos, alt)
    endpoint = f"{BASE_URL}/vep/{species}/region/{variant_region}"
    params = {
        "content-type": "application/json",
        "assembly": "GRCh37"
    }
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(endpoint, headers=headers, params=params, timeout=10)
        data = response.json()
        
        # Handle API error responses
        if "error" in data:
            return handle_vep_error(data["error"], variant_region, chrom, pos, ref, alt)
        
        response.raise_for_status()
        return parse_vep_response(data)
        
    except requests.exceptions.RequestException as e:
        print(f"Warning: VEP API request failed for {variant_region}: {e}", file=sys.stderr)
        return create_error_response('API_ERROR')

