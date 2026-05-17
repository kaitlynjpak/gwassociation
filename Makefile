# Makefile for GW-Assoc Framework

.PHONY: help install test clean docs

help:
	@echo "GW-Assoc Framework Makefile"
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
	python test_gw_assoc.py

test-quick:
	python test_gw_assoc.py --quick

test-minimal:
	python test_gw_assoc.py --minimal

test-coverage:
	pytest --cov=gw_assoc tests/

clean:
	rm -rf build dist *.egg-info
	rm -rf __pycache__ */__pycache__ */*/__pycache__
	rm -rf .pytest_cache .coverage htmlcov
	rm -rf test_data test_plots test_output
	find . -name "*.pyc" -delete
	find . -name "*.pyo" -delete

format:
	black gw_assoc/

lint:
	flake8 gw_assoc/ --max-line-length=100

docs:
	@echo "Generating documentation..."
	python -m pydoc -w gw_assoc
	@echo "Documentation generated as HTML files"