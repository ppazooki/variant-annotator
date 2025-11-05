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
python variant_annotator.py data/input.vcf --output data/output.tsv
```

### Using Make (Docker):
```bash
# Build the Docker image
make build

# Run annotation
make run VCF=data/input.vcf OUTPUT=data/output.tsv

# Run with limit
make run VCF=data/input.vcf OUTPUT=data/output.tsv LIMIT=100

# Use default output filename
make run VCF=data/input.vcf
```

**Parameters:**
- `VCF` - Path to input VCF file (relative to project root, required)
- `OUTPUT` - Path to output TSV file (relative to project root, optional, defaults to `output.tsv`)
- `LIMIT` - Number of variants to process (optional, for testing)

**Troubleshooting:**
- **Error: VCF file not specified** - Make sure to include `VCF=path/to/file.vcf` in your command
- **File not found** - Ensure the VCF file path is relative to the project root directory
- **Docker not running** - Start Docker Desktop or your Docker daemon before running
- **Build errors** - Run `make build` first to create the Docker image
- **Permission errors** - Ensure Docker has access to the project directory

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
  1. `depth` - Depth of sequence coverage at the site of variation
  2. `variant_reads` - Number of reads supporting the variant
  3. `variant_percentage`, `reference_percentage` - Percentage of reads supporting variant vs reference
  4. `gene_id`, `gene_symbol`, `variant_type`, `consequence_terms` - Gene information and variant effects from Ensembl VEP
  5. `maf` - Minor allele frequency (if available)
- **Additional annotations:**
  - `chromosome`, `position`, `variant_id`, `reference`, `alternate`, `quality`, `reference_reads`, `allele_frequency`, `rsid`
