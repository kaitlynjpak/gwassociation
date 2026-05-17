# Makefile for GWAssociation Framework

.PHONY: help install test clean docs

help:
	@echo "GWAssociation Framework Makefile"
	@echo "============================"
	@echo "Available targets:"
	@echo "  install      Install the package in development mode"
	@echo "  test         Run all tests"
	@echo "  test-quick   Run quick tests only"
	@echo "  test-minimal Run minimal tests (no optional deps)"
	@echo "  clean        Remove build artifacts and cache files"
	@echo "  docs         Generate documentation"
	@echo "  format       Format code with black"
	@echo "  lint         Check code with flake8"

install:
	pip install -e .
	pip install -r requirements.txt

install-full:
	pip install -e ".[full]"

install-dev:
	pip install -e ".[dev]"

test:
	python test_gwassociation.py

test-quick:
	python test_gwassociation.py --quick

test-minimal:
	python test_gwassociation.py --minimal

test-coverage:
	pytest --cov=gwassociation tests/

clean:
	rm -rf build dist *.egg-info src/*.egg-info
	rm -rf __pycache__ */__pycache__ */*/__pycache__
	rm -rf .pytest_cache .coverage htmlcov
	rm -rf test_data test_plots test_output
	find . -name "*.pyc" -delete
	find . -name "*.pyo" -delete

format:
	black src/gwassociation/

lint:
	flake8 src/gwassociation/ --max-line-length=100

docs:
	@echo "Generating documentation..."
	PYTHONPATH=src python -m pydoc -w gwassociation
	@echo "Documentation generated as HTML files"