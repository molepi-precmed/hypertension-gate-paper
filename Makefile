# Makefile for hypertension-gate-paper
# Usage:
#   make          → build everything (deps + data check + manuscript)
#   make deps     → restore R package environment via renv
#   make data     → fetch data from GitHub Release (skipped if local data present)
#   make build    → render hypertension_paper.Rmd to output/hypertension_paper.pdf
#   make clean    → remove compiled output
#   make clean-data → remove fetched data (not local data dir)
#
# Use system R libraries only (no renv restore, no installs):
#   export USE_SYSTEM_R_LIBS=1
# Then: make deps  → no-op (skips renv); make build  → uses preinstalled packages only.
# Other hosts: leave unset for normal make deps (renv::restore) and install-if-missing.

RSCRIPT    := Rscript --vanilla
MANUSCRIPT := hypertension_paper.Rmd
OUTPUT_PDF := output/hypertension_paper.pdf

.PHONY: all deps data build clean clean-data

all: deps data build

## Restore R package environment from renv.lock (skipped when USE_SYSTEM_R_LIBS=1)
deps:
	@if [ -n "$$USE_SYSTEM_R_LIBS" ]; then \
	  echo ">>> USE_SYSTEM_R_LIBS is set: using system R libraries only (skipping renv restore)."; \
	else \
	  echo ">>> Restoring R package environment..."; \
	  $(RSCRIPT) -e "if (!requireNamespace('renv', quietly = TRUE)) \
	    install.packages('renv', repos = 'https://cloud.r-project.org'); \
	    renv::restore(prompt = FALSE)"; \
	fi

## Fetch data from GitHub Release if local data dir is absent
data:
	@echo ">>> Checking / fetching data..."
	$(RSCRIPT) R/fetch_data.R

## Run the analysis pipeline, which renders the manuscript
## (Always runs; data/ is not in the dependency list so we do not use file timestamps.)
## If ~/.TinyTeX exists, prepend it to PATH so the LaTeX step uses TinyTeX (e.g. tabu.sty).
build:
	@echo ">>> Running analysis pipeline..."
	@mkdir -p output
	@([ -d "$$HOME/.TinyTeX/bin" ] && export PATH="$$HOME/.TinyTeX/bin/$$(ls "$$HOME/.TinyTeX/bin" 2>/dev/null | head -1):$$PATH"; $(RSCRIPT) pipeline.R) && echo ">>> Output: $(OUTPUT_PDF)"

## Remove compiled outputs
clean:
	@echo ">>> Cleaning output..."
	rm -f output/hypertension_paper.pdf output/hypertension_paper.tex

## Remove fetched data (leaves local sibling dir untouched)
clean-data:
	@echo ">>> Removing fetched data files from data/..."
	find data/ -type f ! -name '.gitkeep' -delete

