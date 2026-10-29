SHELL := /usr/bin/env bash

CONFIG ?= config/iso_slam_as.yaml
PYTHON_BIN ?= python3
SNAKEMAKE ?= snakemake
RUNTIME_TMP ?= tmp/snakemake_runtime
SNAKEMAKE_CACHE ?= tmp/snakemake_cache
MPL_CACHE ?= tmp/mplconfig
SNAKEMAKE_ENV = TMPDIR=$(RUNTIME_TMP) TMP=$(RUNTIME_TMP) TEMP=$(RUNTIME_TMP) XDG_CACHE_HOME=$(SNAKEMAKE_CACHE) MPLCONFIGDIR=$(MPL_CACHE)

.PHONY: help lint validate layout-audit dry-run run smoke tree

help:
	@printf "%s\n" \
	"make validate     validate the configured AS sample sheet" \
	"make lint         compile active Python helpers" \
	"make layout-audit audit the seven-directory contract" \
	"make dry-run      resolve the Snakemake DAG without running" \
	"make run          execute the configured AS workflow" \
	"make smoke        run lightweight syntax and layout checks"

lint:
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/as_core.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/check_samples.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/stage_input_audit.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/stage_as_tables.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/stage_as_figures.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/stage_as_report.py
	"$(PYTHON_BIN)" -m py_compile pipeline/scripts/stage_customer_delivery.py

validate:
	"$(PYTHON_BIN)" pipeline/scripts/check_samples.py \
		--samples config/samples.iso_slam_as.tsv --analysis-type lr_as_rebuild

layout-audit:
	echo 'optional layout checker is not part of this repository'

dry-run:
	mkdir -p "$(RUNTIME_TMP)" "$(SNAKEMAKE_CACHE)" "$(MPL_CACHE)"
	$(SNAKEMAKE_ENV) "$(SNAKEMAKE)" --snakefile pipeline/Snakefile --configfile "$(CONFIG)" --cores 1 --dry-run

run:
	mkdir -p "$(RUNTIME_TMP)" "$(SNAKEMAKE_CACHE)" "$(MPL_CACHE)"
	$(SNAKEMAKE_ENV) "$(SNAKEMAKE)" --snakefile pipeline/Snakefile --configfile "$(CONFIG)" --cores 1 --rerun-incomplete

smoke: lint layout-audit

tree:
	find . -maxdepth 3 -not -path './.git/*' | sort
