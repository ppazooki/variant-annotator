# Variant Annotator

Annotate variants from VCF files using the Ensembl REST APIs.

## Features

- Parses VCF files and extracts variant information
- Annotates variants with gene information and consequence terms
- Extracts rsIDs and Minor Allele Frequency (MAF) data
- Uses batch processing for efficient annotation of large variant sets
- Exports results to TSV format

## Ensembl APIs Used

This tool integrates with **three Ensembl REST API endpoints**:

1. **VEP Batch API** (`POST /vep/human/hgvs`) - Primary method for batch variant effect prediction using HGVS notation
2. **VEP Region API** (`GET /vep/human/region/{region}`) - Fallback method for individual complex variants that cannot be processed by the batch endpoint
3. **Variation API** (`GET /variation/human/{rsid}`) - Fetches Minor Allele Frequency (MAF) data from population databases (1000 Genomes, gnomAD)

## Installation

### Using pip (for local development):
```bash
make install
# or
pip install -e ".[test]"
```

### Using Docker:
```bash
make build
```

## Usage

### Command Line (Local):
```bash
python variant_annotator.py input.vcf [--output output.tsv] [--limit N]
```

**Arguments:**
- `input.vcf` - Input VCF file (required)
- `--output` - Output TSV file path (default: `output.tsv`)
- `--limit` - Limit number of variants to process (optional, for testing)

**Example:**
```bash
python variant_annotator.py data/sample_annotations.vcf --output data/output.tsv
```

### Using Make (Docker):
```bash
# Build the Docker image
make build

# Run annotation (specify full path relative to project root)
make run VCF=data/sample_annotations.vcf OUTPUT=data/output.tsv

# Run with limit
make run VCF=data/sample_annotations.vcf OUTPUT=data/output.tsv LIMIT=100
```

**Note:** When using `make run`, specify the path to your VCF file relative to the project root. The output file path is also relative to the project root. Complex variants (indels, multi-allelic variants) that cannot be processed by the batch endpoint are automatically retried using individual VEP API calls.

## Testing

```bash
make test
```

## Available Make Commands

- `make install` - Install the package and test dependencies
- `make test` - Run the test suite
- `make build` - Build the Docker image
- `make run VCF=path/to/input.vcf [OUTPUT=path/to/output.tsv] [LIMIT=N]` - Run annotation using Docker

## Output

The tool generates a TSV file with the following columns:
- Variant info: chromosome, position, variant_id, reference, alternate, quality, variant_type
- Read statistics: depth, variant_reads, reference_reads, variant_percentage, reference_percentage, allele_frequency
- Annotations: gene_id, gene_symbol, consequence_terms, rsid, maf
