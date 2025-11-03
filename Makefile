.PHONY: install test clean

install:
	pip install -e ".[test]"

test:
	pytest

clean:
	rm -rf __pycache__ tests/__pycache__ htmlcov .pytest_cache .coverage

