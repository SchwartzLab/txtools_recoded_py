# txtools developer tasks
PYTHON ?= python

.PHONY: install dev test binary clean

## Editable install into the current (conda) environment
install:
	pip install -e .

## Install with dev + binary extras
dev:
	pip install -e ".[dev,binary]"

## Run the test suite
test:
	$(PYTHON) -m pytest -q

## Build a standalone single-file executable (on-demand, platform-specific).
## Requires the `binary` extra (pyinstaller). The htslib shared libs that pysam
## relies on are collected explicitly so the frozen binary is self-contained.
binary:
	pyinstaller --onefile --name txtools \
		--collect-all pysam \
		--copy-metadata txtools \
		src/txtools/__main__.py
	@echo "Standalone binary written to dist/txtools"

clean:
	rm -rf build dist *.spec src/txtools/__pycache__ src/txtools.egg-info src/*.egg-info
