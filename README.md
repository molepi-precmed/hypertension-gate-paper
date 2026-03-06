# hypertension-gate-paper

Reproducible manuscript build system for the GATE hypertension analysis.

## Directory layout

This repository (`hypertension-gate-paper/`) expects the large data directory
to exist as a **sibling** directory on your machine:

```text
<parent>/
├── hypertension-gate-paper/
└── hypertension-gate-data/
```

- If `../hypertension-gate-data/` exists, the manuscript will read data from
  there (local development).
- Otherwise, `make data` will fetch a copy into `hypertension-gate-paper/data/`
  from the GitHub Release tag `data-release` using `piggyback`.

## Build the manuscript

From the repo root:

```bash
make
```

Useful targets:

- `make deps`: restore the R environment from `renv.lock`
- `make data`: download data into `data/` (skipped if local sibling data exists)
- `make build`: render `hypertension_paper.Rmd` → `output/hypertension_paper.pdf`
- `make clean`: remove compiled outputs
- `make clean-data`: remove fetched data files under `data/` (keeps `.gitkeep`)

## Data uploads (authors only)

To upload data from the local sibling directory to GitHub Releases:

```bash
Rscript R/upload_data.R
```

This requires `GITHUB_PAT` with write access to the repository.
