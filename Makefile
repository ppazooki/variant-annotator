.PHONY: install test build run

install:
	pip install -e ".[test]"

test:
	pytest

build:
	docker build -t variant-annotator .

run:
	@test -n "$(VCF)" || (echo "Error: VCF parameter required"; exit 1)
	docker run --rm -v $(PWD):/workspace -v $(PWD):/app variant-annotator \
		/workspace/$(VCF) \
		--output /workspace/$(or $(OUTPUT),output.tsv) \
		$(if $(LIMIT),--limit $(LIMIT))

