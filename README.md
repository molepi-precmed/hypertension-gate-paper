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
  from the GitHub Release tag `data-release` using `piggyback`. For many files,
  GitHub’s API rate limit (60 req/hr unauthenticated) can cause 403 errors;
  **set `GITHUB_PAT`** (a personal access token with no extra scopes) for higher
  limits when running `make data`.

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
script downloads release assets from GitHub. Unauthenticated requests are limited
to 60 API calls per hour, which can trigger **403** errors when many files are
downloaded. Setting a GitHub Personal Access Token (PAT) raises the limit and
avoids this.

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

To upload data from the local sibling directory to GitHub Releases:

```bash
Rscript R/upload_data.R
```

This requires `GITHUB_PAT` with write access to the repository.
