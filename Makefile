.PHONY: install test build run

install:
	pip install -e ".[test]"

test:
	pytest

build:
	docker build -t variant-annotator .

run:
	@if [ -z "$(VCF)" ]; then \
		echo "Error: VCF file not specified. Usage: make run VCF=path/to/input.vcf [OUTPUT=path/to/output.tsv] [LIMIT=N]"; \
		exit 1; \
	fi
	@OUTPUT_FILE=$(if $(OUTPUT),$(OUTPUT),output.tsv); \
	LIMIT_ARG=$(if $(LIMIT),--limit $(LIMIT),); \
	docker run --rm -v $(PWD):/workspace -v $(PWD):/app variant-annotator /workspace/$(VCF) --output /workspace/$$OUTPUT_FILE $$LIMIT_ARG

