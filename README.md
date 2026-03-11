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
- Otherwise, `make data` will download a single archive from the GitHub Release
  tag `data-release` and extract it into `hypertension-gate-paper/data/`
  (preserving `eqtl/`, `pqtl/`, `mrresults/`, etc.). One asset keeps the layout
  intact and avoids API rate limits; **set `GITHUB_PAT`** if you see 403 errors.

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

## Setting GITHUB_PAT (for data download)

If you run `make data` without a local `../hypertension-gate-data/` directory, the
script downloads one archive from the release. If you see **403** errors (e.g. on
private repos or under rate limits), set a GitHub Personal Access Token (PAT).

### Create a token

1. On GitHub: **Settings** → **Developer settings** → **Personal access tokens** →
   **Tokens (classic)** → **Generate new token (classic)**.
2. Give it a name (e.g. `hypertension-gate-paper data`), choose an expiry, and
   leave all scopes **unchecked** (no scopes are required for public release
   downloads).
3. Generate the token and copy it (it starts with `ghp_`). You won’t see it again.

### Use the token when running the build

**Option A — one-off (current terminal only):**

```bash
GITHUB_PAT="ghp_xxxxxxxxxxxx" make data
```

**Option B — in your shell profile** (e.g. `~/.zshrc` or `~/.bashrc`), so it’s
set whenever you open a terminal:

```bash
export GITHUB_PAT="ghp_xxxxxxxxxxxx"
```

Then run `make data` (or `make`) as usual. With the token set, 403 rate-limit
errors from the download step should stop.

## Data uploads (authors only)

To upload data from the local sibling directory to GitHub Releases (as a single
`hypertension-gate-data.tar.gz` asset):

```bash
Rscript R/upload_data.R
```

This requires `GITHUB_PAT` with write access to the repository. The script
archives the contents of `../hypertension-gate-data/` and uploads that one file;
`make data` then downloads and extracts it so `data/eqtl/`, `data/pqtl/`, etc.
match the original layout.

## Using system R libraries (skip installing packages from scratch)

On a host where R and the required packages are already installed (e.g. in
`/usr/local/lib/R/site-library`, `/usr/lib/R/site-library`, `/usr/lib/R/library`),
you can rely on those and avoid running `renv::restore()` or installing anything:

```bash
export USE_SYSTEM_R_LIBS=1
make
```

With `USE_SYSTEM_R_LIBS=1`:

- `make deps` does nothing (no renv restore, no downloads).
- `make build` loads packages from the system library only; it does not install
  missing packages. If a required package is missing, the pipeline stops with an
  error.

Leave `USE_SYSTEM_R_LIBS` unset on other machines if you want the usual
behaviour: `make deps` restores from `renv.lock`, and the pipeline may install
missing packages when needed.
