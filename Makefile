# Build the reversibility paper: figures from data, then the PDFs.
#
#   make figures           regenerate every data-derived figure of main.tex and si.tex
#   make papers            compile main.pdf and si.pdf into manuscript/build/
#   make check             report undefined citations / references / missing figures
#   make clean             remove LaTeX aux files
#   make previous-figures  the figures of the 28 Oct version that are not in this one
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

GEL  := $(ANA)/gelation
S12  := $(ANA)/SAXS_12id
XPCS := $(ANA)/SAXPCS_8id
RAD  := $(ANA)/Rad_Dam_Check

export TEXMFHOME := $(MS)/texmf
export TEXINPUTS := .:$(MS)/figures:$(GEL):$(S12):$(XPCS):$(RAD):$(ANA)/ramps:
export BIBINPUTS := .:$(MS):$(BUILD):
export MPLBACKEND := Agg

STYLE    := $(ANA)/common/acs_style.py
GEL_LIB  := $(GEL)/gelation_data.py $(XPCS)/abs_xsec.py $(STYLE)
GUINIER  := $(S12)/guinier.py $(S12)/reduced_data/Merged_Reference_10C.csv
data      = $(addprefix $(GEL)/data/,$(addsuffix .h5,$(1)))

# --- figures of main.tex and si.tex -----------------------------------------
F_SAXS    := $(GEL)/Figure2_SAXS_Scaling.pdf
F_ARREST  := $(GEL)/Figure3_Dynamic_Arrest.pdf
F_TWOSTEP := $(GEL)/Figure4_Two_Step_Decay.pdf
T_TWOSTEP := $(GEL)/Table_Two_Step_Fit.tex
F_RAMPS   := $(GEL)/FigureS2_Ramp_VPAVG10.pdf $(GEL)/FigureS3_Ramp_VPAVG20.pdf \
             $(GEL)/FigureS4_Ramp_VPAVG30.pdf
F_12ID    := $(S12)/FigureS5_12ID_Reversibility.pdf
F_CYCLE   := $(XPCS)/FigureS6_Thermal_Cycle.pdf
F_IC      := $(XPCS)/FigureS10_IC_Calibration.pdf
F_STAB    := $(GEL)/FigureS11_SAXS_Stability.pdf $(GEL)/FigureS12_XPCS_Stability.pdf
F_FLUX    := $(RAD)/FigureS13_Flux_Control.pdf
FIGURES   := $(F_SAXS) $(F_ARREST) $(F_TWOSTEP) $(T_TWOSTEP) $(F_RAMPS) $(F_12ID) \
             $(F_CYCLE) $(F_IC) $(F_STAB) $(F_FLUX)
ARTWORK   := $(addprefix $(MS)/figures/,Setup.pdf VideoS1_Stills.jpg VideoS2_Stills.jpg \
             FigureS1_Sample_Cells.jpg FigureS7_Irreversibility.jpg FigureS8_MD.png)
XPCS_DATA := $(wildcard $(XPCS)/data/*.hdf) $(XPCS)/data/thermal_cycle_temperature.csv

.PHONY: all figures papers check clean previous-figures

all: figures papers

figures: $(FIGURES)
	@echo "all figures up to date"

# Every rule below reads only what is committed and never rebuilds it.  The
# data were reduced on the /gdata hosts (CLAUDE.md) from /gdata/s8id-dmdtn alone:
# analysis/gelation/data/ by reduce_gelation.py, analysis/SAXPCS_8id/data/ by
# average_ranges.py, SAXS_12id/reduced_data/ by Read_12ID_SAWAXS.py, and
# Rad_Dam_Check/data/ by reduce_flux_control.py.
$(F_SAXS): $(GEL)/saxs_scaling.py $(GEL_LIB) $(GUINIER) $(call data,E0110_n10 D0138)
	cd $(GEL) && $(PY) saxs_scaling.py

$(F_ARREST): $(GEL)/g2_arrest.py $(GEL_LIB) $(call data,B039_n10)
	cd $(GEL) && $(PY) g2_arrest.py

$(F_TWOSTEP) $(T_TWOSTEP) &: $(GEL)/two_step_fit.py $(GEL_LIB) \
                             $(call data,B039_n20_skip70 T214_n20 E0110_n20)
	cd $(GEL) && $(PY) two_step_fit.py

$(F_RAMPS) &: $(GEL)/ramp_si.py $(GEL_LIB) \
              $(call data,B039_n10 B040_n10 T214_n20 T215_n20 E0110_n10 E0111_n10 D029 D0138)
	cd $(GEL) && $(PY) ramp_si.py

$(F_12ID): $(S12)/reversibility_12id.py $(S12)/saxs_waxs_fits.py $(GUINIER) \
           $(wildcard $(S12)/reduced_data/*.csv) $(STYLE)
	cd $(S12) && $(PY) reversibility_12id.py

$(F_CYCLE): $(XPCS)/thermal_cycle.py $(XPCS)/abs_xsec.py $(XPCS_DATA) $(STYLE) \
            $(XPCS)/xpcs_fit.py $(XPCS)/nexus_read.py
	cd $(XPCS) && $(PY) thermal_cycle.py

$(F_IC): $(XPCS)/ic_calibration.py $(XPCS)/abs_xsec.py $(STYLE)
	cd $(XPCS) && $(PY) ic_calibration.py

$(F_STAB) &: $(GEL)/stability_si.py $(GEL_LIB) \
             $(call data,B039_segments B039_n10 B039_part1_n10 B039_part2_n10)
	cd $(GEL) && $(PY) stability_si.py

$(F_FLUX): $(RAD)/flux_control.py $(RAD)/Pind_calibration.csv \
           $(RAD)/data/flux_control.h5 $(STYLE)
	cd $(RAD) && $(PY) flux_control.py

# --- the 28 Oct 2026 version of main.tex ----------------------------------------
# Two of its three figures are not in the present manuscript: the 12-ID-B
# SAXS/WAXS figure as a single panel (its fits are Figure S5a here) and the
# chain-length ramps in 2 C windows.  Its third, the seven-cycle thermal
# cycling, is Figure S6.  The scripts and data are kept and still build.
P_12ID  := $(S12)/Figure1_SAXS_WAXS.pdf
P_RAMPS := $(ANA)/ramps/Figure3_Chain_Length_Ramps.pdf
RAMP_DATA := $(wildcard $(ANA)/ramps/data/*.hdf) $(ANA)/ramps/data/trace_2021-2.csv

previous-figures: $(P_12ID) $(P_RAMPS)

# Plot_12ID.py writes Figure2_SAXS_WAXS.pdf, its name in the combined Letter.
$(P_12ID): $(S12)/Plot_12ID.py $(S12)/saxs_waxs_fits.py $(GUINIER) \
           $(wildcard $(S12)/reduced_data/*.csv) $(STYLE)
	cd $(S12) && $(PY) Plot_12ID.py && cp Figure2_SAXS_WAXS.pdf Figure1_SAXS_WAXS.pdf

$(P_RAMPS): $(ANA)/ramps/plot_ramps.py $(XPCS)/abs_xsec.py $(RAMP_DATA) $(STYLE) \
            $(XPCS)/nexus_read.py
	cd $(ANA)/ramps && $(PY) plot_ramps.py

# --- LaTeX ------------------------------------------------------------------------
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

# si -> main -> si, so the cross-references resolve in both directions.
papers: $(FIGURES) $(ARTWORK) $(MS)/main.tex $(MS)/si.tex $(MS)/reference.bib
	$(call compile,si)
	$(call compile,main)
	$(call compile,si)

check:
	@echo "== undefined citations / references =="
	@grep -hoE "(Citation|Reference) \`[^']+'.*undefined" \
	    $(BUILD)/main.log $(BUILD)/si.log 2>/dev/null | sort -u || true
	@echo "== bibtex problems =="
	@grep -hE "couldn't open|I didn't find a database entry" $(BUILD)/*.blg 2>/dev/null | sort -u | head -20 || true
	@echo "== missing figures =="
	@grep -hE "No file|not found" $(BUILD)/*.log 2>/dev/null | grep -iE '\.pdf|\.png|\.jpg' || true
	@echo "== done (no lines under a heading means clean) =="

clean:
	rm -f $(BUILD)/*.aux $(BUILD)/*.log $(BUILD)/*.out $(BUILD)/*.bbl $(BUILD)/*.blg $(BUILD)/*.toc
