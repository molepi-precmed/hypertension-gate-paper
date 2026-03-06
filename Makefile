# Makefile for hypertension-gate-paper
# Usage:
#   make          → build everything (deps + data check + manuscript)
#   make deps     → restore R package environment via renv
#   make data     → fetch data from GitHub Release (skipped if local data present)
#   make build    → render hypertension_paper.Rmd to output/hypertension_paper.pdf
#   make clean    → remove compiled output
#   make clean-data → remove fetched data (not local data dir)

RSCRIPT    := Rscript --vanilla
MANUSCRIPT := hypertension_paper.Rmd
OUTPUT_PDF := output/hypertension_paper.pdf

.PHONY: all deps data build clean clean-data

all: deps data build

## Restore R package environment from renv.lock
deps:
	@echo ">>> Restoring R package environment..."
	$(RSCRIPT) -e "if (!requireNamespace('renv', quietly = TRUE)) \
	  install.packages('renv', repos = 'https://cloud.r-project.org'); \
	  renv::restore(prompt = FALSE)"

## Fetch data from GitHub Release if local data dir is absent
data:
	@echo ">>> Checking / fetching data..."
	$(RSCRIPT) R/fetch_data.R

## Run the analysis pipeline, which renders the manuscript
build: $(OUTPUT_PDF)

$(OUTPUT_PDF): pipeline.R $(MANUSCRIPT) R/config.R R/utils.R
	@echo ">>> Running analysis pipeline..."
	@mkdir -p output
	$(RSCRIPT) pipeline.R
	@echo ">>> Output: $(OUTPUT_PDF)"

## Remove compiled outputs
clean:
	@echo ">>> Cleaning output..."
	rm -f output/hypertension_paper.pdf output/hypertension_paper.tex

## Remove fetched data (leaves local sibling dir untouched)
clean-data:
	@echo ">>> Removing fetched data files from data/..."
	find data/ -type f ! -name '.gitkeep' -delete

