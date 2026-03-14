# Project guidelines for Claude

## Building the manuscript from Rmarkdown
- Use `make build` from the command line. 
- Ignore ghostscript errors


## Scientific writing conventions
- HUGO gene symbols should always be in italic (e.g. *ABCA1*), unless the text refers to the encoded protein and the protein has the same acronym as the gene — in that case use roman type
- Regression coefficients in table headers should be labelled "Slope" rather than the Greek letter β, as "Slope" is more interpretable
- Legend titles in ggplot figures should appear above the plot symbols (`guide_legend(title.position = "top")`)
- Table footnote text should be left-aligned (`\raggedright`); kableExtra tables default to centred footnotes

## R Markdown / pandoc / knitr
- To place raw LaTeX (e.g. `\begin{adjustwidth}`) around a knitr code chunk without suppressing pandoc's processing of the chunk output, use pandoc raw blocks (` ```{=latex} `) rather than bare LaTeX commands — bare LaTeX causes pandoc to treat the entire surrounding block as raw LaTeX and leave markdown image syntax unconverted
- To make a figure wider than the text column, use `\begin{adjustwidth}{-Xcm}{-Xcm}` (from the `changepage` package) as pandoc raw blocks around the chunk, `fig.pos='H'` (requires `\usepackage{float}`), and `out.width='100%'` — the figure then fills the wider `\linewidth`
- `fig.width` controls the R graphics device size (aspect ratio and resolution); `out.width` controls the printed width in the document — they are independent
- To suppress a table or figure chunk without deleting it, add `eval=FALSE` to the chunk header
- When `eqtl.start`/`eqtl.end` row indices for `pack_rows()` or `row_spec()` are computed before a filter is applied to the table data, they must be recomputed from the filtered data

## kableExtra
- kableExtra's fixed-width column padding can introduce trailing spaces, causing `},` sequences to break LaTeX; fix with `gsub("}[ \t]+,", "},", ...)` before `knitr::asis_output()`
- `footnote()` text alignment is centred by default inside `table` environments; add `\raggedright` before `\fontsize{...}` in the injected footnote LaTeX to left-align it
- Footnotes inserted by `sub("\\end{table}", ...)` or `sub("\\end{longtable}", ...)` should use `\par\noindent{\footnotesize ...}` for proper paragraph-level formatting

## ggraph / bipartite graphs
- For bipartite graph layouts, fix layer y-coordinates (clumps at y=1, targets at y=0) after `layout_with_stress()` to enforce a two-row layout; then evenly re-space x-coordinates within each layer
- When a GWAS clump has no attributed gene symbol, label it `<chr:pos>` (chromosome and start position in Mb to one decimal place) to provide location context

## LaTeX / document formatting
- To remove the fancyhdr footer rule, use `\renewcommand{\footrulewidth}{0pt}` (not just removing the `\footrule` definition)
- Page date belongs on the title page (e.g. as a bullet item computed via `format(Sys.Date(), "%d %B %Y")`), not in the running footer
- `\usepackage{float}` is required for `fig.pos='H'` (forced "here" figure placement) to work
