# Build the reversibility paper: figures from data, then the PDF.
#
#   make figures    regenerate every data-derived figure
#   make papers     compile main.pdf into manuscript/build/
#   make check      report undefined citations / references / missing figures
#   make clean      remove LaTeX aux files
#
# Mirrors the Makefile of PubCode-2021-IDPP: same conda environment, same
# TEXINPUTS trick so the .tex can name a figure that lives under analysis/,
# same vendored texmf so no root TeX install is needed.

SHELL := /bin/bash
ROOT  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
MS    := $(ROOT)/manuscript
ANA   := $(ROOT)/analysis
BUILD := $(MS)/build

CONDA_ENV := env-2021-LLDP
PY        := conda run --no-capture-output -n $(CONDA_ENV) python

export TEXMFHOME := $(MS)/texmf
export TEXINPUTS := .:$(MS)/figures:$(ANA)/SAXS_12id:$(ANA)/SAXPCS_8id:$(ANA)/ramps:
export BIBINPUTS := .:$(MS):$(BUILD):
export MPLBACKEND := Agg

F_12ID  := $(ANA)/SAXS_12id/Figure1_SAXS_WAXS.pdf
F_CYCLE := $(ANA)/SAXPCS_8id/Figure2_Thermal_Cycle.pdf
XPCS_DATA := $(wildcard $(ANA)/SAXPCS_8id/data/*.hdf) \
             $(ANA)/SAXPCS_8id/data/thermal_cycle_temperature.csv
F_RAMPS := $(ANA)/ramps/Figure3_Chain_Length_Ramps.pdf
RAMP_DATA := $(wildcard $(ANA)/ramps/data/*.hdf) $(ANA)/ramps/data/trace_2021-2.csv
FIGURES := $(F_12ID) $(F_CYCLE) $(F_RAMPS)

.PHONY: all figures papers check clean

all: figures papers

figures: $(FIGURES)
	@echo "all figures up to date"

# Plot_12ID.py writes Figure2_SAXS_WAXS.pdf, the name it had in the combined
# Letter; this paper calls it Figure 1.  Renaming the output is a Monday job --
# see MONDAY_TODO.md -- so for now the copy is made here and the provenance is
# explicit rather than silent.
$(F_12ID): $(ANA)/SAXS_12id/Plot_12ID.py $(wildcard $(ANA)/SAXS_12id/reduced_data/*.csv) $(ANA)/common/acs_style.py
	cd $(ANA)/SAXS_12id && $(PY) Plot_12ID.py && cp Figure2_SAXS_WAXS.pdf Figure1_SAXS_WAXS.pdf

# thermal_cycle.py writes FigureS4_Thermal_Cycle.pdf, the name it had in the
# combined Letter; this paper calls it Figure 2.  Renaming the output in the
# script is a Monday job -- see MONDAY_TODO.md.
$(F_CYCLE): $(ANA)/SAXPCS_8id/thermal_cycle.py $(ANA)/SAXPCS_8id/abs_xsec.py \
            $(XPCS_DATA) $(ANA)/common/acs_style.py $(ANA)/SAXPCS_8id/xpcs_fit.py \
            $(ANA)/SAXPCS_8id/nexus_read.py
	cd $(ANA)/SAXPCS_8id && $(PY) thermal_cycle.py && mv FigureS4_Thermal_Cycle.pdf Figure2_Thermal_Cycle.pdf

# Reads only ramps/data/, written by reduce_ramps.py on a /gdata host and
# committed, so this runs from a clone.
$(F_RAMPS): $(ANA)/ramps/plot_ramps.py $(ANA)/SAXPCS_8id/abs_xsec.py \
            $(RAMP_DATA) $(ANA)/common/acs_style.py $(ANA)/SAXPCS_8id/nexus_read.py
	cd $(ANA)/ramps && $(PY) plot_ramps.py

define compile
	@mkdir -p $(BUILD)
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@cd $(MS) && if grep -q '\\bibdata' build/$(1).aux 2>/dev/null; then \
	    bibtex build/$(1) >/dev/null; \
	    pdflatex -interaction=nonstopmode -halt-on-error \
	        -file-line-error -output-directory=build $(1).tex >/dev/null; \
	fi
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@echo "  -> manuscript/build/$(1).pdf ($$(pdfinfo $(BUILD)/$(1).pdf | awk '/^Pages/{print $$2}') pages)"
endef

papers: $(FIGURES) $(MS)/main.tex $(MS)/reference.bib
	$(call compile,main)

check:
	@echo "== undefined citations / references =="
	@grep -hoE "(Citation|Reference) \`[^']+'.*undefined" $(BUILD)/main.log 2>/dev/null | sort -u || true
	@echo "== bibtex problems =="
	@grep -hE "couldn't open|I didn't find a database entry" $(BUILD)/*.blg 2>/dev/null | sort -u | head -20 || true
	@echo "== missing figures =="
	@grep -hE "No file|not found" $(BUILD)/*.log 2>/dev/null | grep -iE '\.pdf|\.png' || true
	@echo "== done =="

clean:
	rm -f $(BUILD)/*.aux $(BUILD)/*.log $(BUILD)/*.out $(BUILD)/*.bbl $(BUILD)/*.blg $(BUILD)/*.toc
