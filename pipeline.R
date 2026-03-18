source("R/config.R")

# Ensure required packages are available (load from .libPaths(), install if missing).
for (pkg in REQUIRED_PACKAGES)
  ensure_library(pkg)

# Data and results live in the sibling hypertension-gate-data directory.
working.dir <- DATA_DIR
reports.dir <- getwd()
basegatedir <- DATA_DIR

library(flextable)
library(data.table)
library(ggplot2)
library(writexl)

source("report.helpers.R")
source("shared.helpers.R")

phenoname <- "meanarterialpressure_imputed"
name.string.eqtl <- "map_score-excl_strict-thresh"
name.string.pqtl <- "map_score-excl_strict-thresh"
coeffs.wide.eqtl.file <- sprintf("coeffs.wide.%s.%s.%s.csv",
                                  "eqtl", phenoname, name.string.eqtl)
coeffs.wide.pqtl.file <- sprintf("coeffs.wide.%s.%s.%s.csv",
                                 "pqtl", phenoname, name.string.pqtl)

annotate.gwashits <- TRUE
compute.correlations <- FALSE
run.create.clumps <- FALSE
test.proteins <- FALSE
newrun <- FALSE
generate.summaries <- FALSE

## filters for associations
min.eql.diversity <- 5
min.pqtl.diversity <- 5
min.pvalue <- 1e-5
min.pvalue.pqtl <- 1e-5
select.genes <- NULL
overlap.flank <- 2e5

proteins.associations.file <- "coeffs.protein.all.csv"
mr.results.file <- "mrres.RDS"

#' Helper function to load a given Rdata object.
load.rdata <- function(filename) {
    load(filename)
    get(ls()[ls() != "filename"])
}

## Pipeline to create automated report from summary stats
## -----------------------------------------------------------------------------
coeffs.eqtl <- process.coeffs(study="eqtl", gate.dir=basegatedir, working.dir,
                              coeffs.wide.file=coeffs.wide.eqtl.file,
                              min.diversity=min.eql.diversity,
                              min.pvalue=min.pvalue, select.genes=select.genes)
coeffs.pqtl <- process.coeffs(study="pqtl", gate.dir=basegatedir, working.dir,
                              coeffs.wide.file=coeffs.wide.pqtl.file,
                              min.diversity=min.pqtl.diversity,
                              min.pvalue=min.pvalue.pqtl, select.genes=select.genes)

## Q-Q plots
save.qqplot(study="eqtl", coeffs.wide=coeffs.eqtl$coeffs.wide, working.dir)
save.qqplot(study="pqtl", coeffs.wide=coeffs.pqtl$coeffs.wide, working.dir)

## -----------------------------------------------------------------------------
## Read GWAS hits
## -----------------------------------------------------------------------------
if (annotate.gwashits) {
    gwas.dt.file <- file.path(working.dir, "gwas-association-downloaded_2025-08-14-EFO_0000537-withChildTraits.tsv")
    if (file.exists(gwas.dt.file)) {
        gwas.dt <- fread(gwas.dt.file)
    } else {
        assoc <- get_associations(efo_id=efo.id, warnings=FALSE, interactive=FALSE)
        variants <- get_variants(efo_id=efo.id, warnings=FALSE, interactive=FALSE)
        traits <- get_traits(efo_id=efo.id, warnings=FALSE)
        studies <- get_studies(efo_id=efo.id, warnings=FALSE)

        assoc_trait_map <- association_to_trait(assoc@associations$association_id, warnings = FALSE)
        assoc_trait_dt <- setDT(as.data.table(assoc_trait_map))[, .(association_id, efo_id)]

        variants_dt <- setDT(as.data.table(variants@variants))[, .(variant_id, chromosome_name, chromosome_position)]
        studies_dt <- setDT(as.data.table(studies@studies))[, .(study_id, reported_trait)]

        assoc_study_map <- association_to_study(assoc@associations$association_id, warnings = FALSE)
        assoc_study_dt <- setDT(as.data.table(assoc_study_map))[, .(association_id, study_id)]

        ## create main join with chromosomal and position info
        gwas.dt <- setDT(as.data.table(assoc@associations))[
            setDT(as.data.table(assoc@loci)), on = "association_id"][
            setDT(as.data.table(assoc@risk_alleles)), on = c("association_id", "locus_id")][
            setDT(as.data.table(assoc@genes)), on = c("association_id", "locus_id")][
            variants_dt, on = "variant_id"][
            assoc_trait_dt, on = "association_id"][
            assoc_study_dt, on = "association_id"][
            studies_dt, on = "study_id"]

        ## select and rename columns to match website download format
        gwas.dt <- gwas.dt[, .(
            SNPS = variant_id,
            CHR_ID = chromosome_name,
            CHR_POS = chromosome_position,
            `P-VALUE` = pvalue,
            `REPORTED GENE(S)` = gene_name,
            `MAPPED_GENE` = gene_name,
            `DISEASE/TRAIT` = reported_trait,
            `MAPPED_TRAIT` = reported_trait,
            STUDY=study_id,
            association_id,
            locus_id
        )]
        fwrite(gwas.dt, file=gwas.dt.file)
    }

    gwas.hits <- clump.GWAShits(gwas.dt, clump.gap=overlap.flank,
                                phenotypes=NULL)
    gwas.hits.genes <- unique(unlist(strsplit(gwas.hits[, reported.genes], split=", ")))
    gwas.hits.genes <- trimws(gwas.hits.genes, whitespace="[\\h\\v]")

    gwas.hits <- gwas.hits[!grepl("HLA-", reported.genes)]
    gwas.hits <- na.omit(gwas.hits,
                         cols=c("chrom", "x.clumpStart", "x.clumpEnd"))

    gwas.hits <- gwas.hits[,
        .(chrom, startpos, endpos, x.clumpStart, x.clumpEnd, reported.genes),
        by=clump]
    gwas.hits <- unique(gwas.hits,
                        by=c("chrom", "x.clumpStart", "x.clumpEnd"))

    coeffs.eqtl$coregenes <- overlap.coeffs.hits(coeff.dt=coeffs.eqtl$coregenes,
                                                 hits.dt=gwas.hits,
                                                 gene.column="reported.genes",
                                                 overlap.flank=overlap.flank)

    coeffs.pqtl$coregenes <- overlap.coeffs.hits(coeff.dt=coeffs.pqtl$coregenes,
                                                 hits.dt=gwas.hits,
                                                 gene.column="reported.genes",
                                                 overlap.flank=overlap.flank)
} else {
    gwas.dt <- NULL
    gwas.hits <- NULL
    coeffs.eqtl$coregenes[, reported.genes := NA]
    coeffs.pqtl$coregenes[, reported.genes := NA]
}

## -----------------------------------------------------------------------------
## QTL clumps
## -----------------------------------------------------------------------------
if (run.create.clumps) {
    create.clumps(study="eqtl", coregenes=coeffs.eqtl$coregenes,
              gate.meta=coeffs.eqtl$gate.meta,
              locus.meta=coeffs.eqtl$locus.meta, gwas.hits=gwas.hits,
                working.dir)
    clumpseqtl.print <- fread(file.path(working.dir, study="eqtl", "clumps.csv"))
    create.clumps(study="pqtl", coregenes=coeffs.pqtl$coregenes,
                gate.meta=coeffs.pqtl$gate.meta,
                locus.meta=coeffs.pqtl$locus.meta, gwas.hits=gwas.hits,
                working.dir)
    clumpspqtl.print <- fread(file.path(working.dir, study="pqtl", "clumps.csv"))
} else {
    clumpseqtl.print <- fread(file.path(working.dir, study="eqtl", "clumps.csv"))
    clumpspqtl.print <- fread(file.path(working.dir, study="pqtl", "clumps.csv"))
}

## -----------------------------------------------------------------------------
## Core genes
## -----------------------------------------------------------------------------
coregenes.eqtl <- unique(coeffs.eqtl$coregenes, by="scoreid_trans")
coregenes.pqtl <- unique(coeffs.pqtl$coregenes, by=c("scoreid_trans", "studyid"))

## -----------------------------------------------------------------------------
## Score correlations
## -----------------------------------------------------------------------------
if (compute.correlations) {
    compute.score.correlations(study="eqtl", phenoname, pheno.dt,
                               gate.dir=basegatedir,
                               working.dir=working.dir,
                               coregenes.file="coregenes.csv")

    compute.score.correlations(study="pqtl", phenoname, pheno.dt,
                               gate.dir=basegatedir,
                               working.dir=working.dir,
                               coregenes.file="coregenes.csv")
}

## -----------------------------------------------------------------------------
## Test for protein associations
## -----------------------------------------------------------------------------
proteins.results.file <- file.path(working.dir, proteins.associations.file)
if (file.exists(proteins.results.file)) {
    coeffs.protein.all <- fread(proteins.results.file)
}

if (test.proteins) {
    proteins.local.file <- file.path(output.base, basename(proteins.file))
    system(sprintf("dx download %s --output %s --overwrite",
                proteins.file, proteins.local.file))

    if (grepl("\\.rds$", proteins.local.file, ignore.case=TRUE)) {
        pheno.proteins <- readRDS(proteins.local.file)
    } else {
        pheno.proteins <- fread(proteins.local.file)
        pheno.proteins[, V1 := NULL]
    }
    setDT(pheno.proteins)

    genenames.proteins <- setdiff(names(pheno.proteins), "eid")

    pheno.proteins <- pheno.proteins[pheno.dt, on=.(eid)]

    y.isbinary <- setequal(length(unique(na.omit(pheno.dt[[phenoname]]))), 2)

    coeffs.protein.all <- get.coeffs.protein(genenames.proteins, pheno.proteins,
                                            formula=proteins.formula,
                                            y.isbinary, cores=cores)
    coeffs.protein.all[, Gene := gsub("_", ":", Gene)]
    fwrite(coeffs.protein.all, file=proteins.results.file)
}

## -----------------------------------------------------------------------------
## Run MRhevo for core genes
## -----------------------------------------------------------------------------
mrres.file <- file.path(working.dir, mr.results.file)

if (file.exists(mrres.file)) {
    all.mrres <- readRDS(mrres.file)
}

if (newrun) {
    if (!"plink.id" %in% names(pheno.dt))
        pheno.dt[, plink.id := sprintf("%s_%s", eid, eid)]

    mrres <- report.mr(core.genes=coregenes.eqtl$gene_symbol,
            scoresinfo=coregenes.eqtl,
            pheno.dt=pheno.dt,
            phenoname=phenoname,
            covariates=mr.covariates,
            study.name="eQTLGen",
            working.dir=working.dir,
            gwas.dt=gwas.dt,
            basegatedir=gate.dir.eqtl,
            ld.refplinkfile=ld.refplinkfile,
            cores=cores)

    mrres <- report.mr(core.genes=coregenes.pqtl$gene_symbol,
            scoresinfo=coregenes.pqtl,
            pheno.dt=pheno.dt,
            phenoname=phenoname,
            covariates=mr.covariates,
            study.name="deCODE",
            working.dir=working.dir,
            gwas.dt=gwas.dt,
            basegatedir=gate.dir.pqtl,
            ld.refplinkfile=ld.refplinkfile,
            cores=cores)

    mrres <- report.mr(core.genes=coregenes.pqtl$gene_symbol,
            scoresinfo=coregenes.pqtl,
            pheno.dt=pheno.dt,
            phenoname=phenoname,
            covariates=mr.covariates,
            study.name="UKBPPP",
            working.dir=working.dir,
            gwas.dt=gwas.dt,
            basegatedir=gate.dir.pqtl,
            ld.refplinkfile=ld.refplinkfile,
            cores=cores)
}

mr.studies <- c("eQTLGen", "deCODE", "UKBPPP")
mr.exposures <- c(coregenes.eqtl$gene_symbol, coregenes.pqtl$gene_symbol)
all.mrres <- list()
for (this.study in mr.studies) {
    for (this.exposure in mr.exposures) {
        this.mr.dir <- file.path(working.dir, "mrresults", this.study)
        mr.file <- list.files(this.mr.dir,
                                # pattern=sprintf("estimators\\.(%s)+\\.\\d+(\\.divergent_transitions_\\d+)?\\.csv",
                                pattern=sprintf("^estimators\\.(%s)\\.\\d+(\\.divergent_transitions_\\d+)?\\.csv",
                                                this.exposure),
                                full.names=TRUE)
        if (length(mr.file)==0) next
        dt <- fread(mr.file)
        dt[, study := this.study]
        dt[, exposure := this.exposure]
        all.mrres <- append(all.mrres, list(dt))
    }
}
all.mrres <- rbindlist(all.mrres, use.names=TRUE, fill=TRUE)
saveRDS(all.mrres, file=mrres.file)

summaries.file <- file.path(working.dir, "summaries.Rdata.gz")
if (generate.summaries) {
    all.cases <- pheno.dt[, table(get(phenoname))]
    included.cases <- pheno.dt[exclude==FALSE, table(get(phenoname))]
    n.variants <- count.variants(basedir=basename,
                                 genotype.file.prefix="ukb22828_c",
                                 genotype.file.suffix="_b0_v3")
    n.variants.forscores <- count.variants(basedir=basename,
                                    genotype.file.prefix=genotype.file.prefix,
                                    genotype.file.suffix=genotype.file.suffix)
    save(all.cases, included.cases, n.variants, n.variants.forscores,
         file=summaries.file)
} else {
    load(summaries.file)
}

 coeffs.protein.all[, Gene := toupper(Gene)]
genenames.proteins.core <- c(coeffs.eqtl$coregenes$gene_symbol,
                              coeffs.pqtl$coregenes$gene_symbol)
coeffs.protein.core <- coeffs.protein.all[Gene %in% genenames.proteins.core]
 
if ("cases" %in% names(coeffs.protein.all)) {
    binary <- TRUE
    numobs <- c("noncases", "cases")
    keep.names <- c("Gene", numobs, "coeff", "pvalue.formatted")
} else {
    binary <- FALSE
    keep.names <- c("Gene", "coeff", "pvalue.formatted")
}
coeffs.protein.core <- coeffs.protein.core[, ..keep.names]
coeffs.protein.core[, coeff := round(coeff, 3)]

## Join coregens and protein measurements
coregenes.eqtl[coeffs.protein.all, on=.(gene_symbol=Gene), `:=`(protein.coeff = i.coeff, protein.pvalue = i.pvalue)]
coregenes.pqtl[coeffs.protein.all, on=.(gene_symbol=Gene), `:=`(protein.coeff = i.coeff, protein.pvalue = i.pvalue)]
# coregenes.pqtl[abs(protein.coeff)>2*abs(estimate_trans), uniqueN(gene_symbol)]

coregenes.eqtl.strict <- coregenes.eqtl[(locus.diversity>=8&pvalue_trans<=1e-5), ]
coregenes.pqtl.strict <- coregenes.pqtl[(locus.diversity>=20&pvalue_trans<=min.pvalue.pqtl), ]

setorder(coregenes.eqtl.strict, pvalue_trans)

setorder(coregenes.pqtl.strict, pvalue_trans)

coregenes.eqtl.rest <- coregenes.eqtl[(locus.diversity<8&locus.diversity>=5&pvalue_trans<=1e-5), ]
coregenes.pqtl.rest <- coregenes.pqtl[(locus.diversity<20&locus.diversity>=5&pvalue_trans<=min.pvalue.pqtl), ]

phenoname <- "mean arterial pressure"

## Prepare table with validations

## prepare coregenes table
## -----------------------------------------------------------------------------

coregenes.tbl <- rbind(coregenes.eqtl.strict, coregenes.pqtl.strict)

## Join with MR
coregenes.tbl[all.mrres, on=.(gene_symbol=exposure), mr.pvalue:=pvalue]

cols.select <- c("gene_symbol", "reported.genes", "estimate_trans", 
                 "protein.coeff", "protein.pvalue", "mr.pvalue", "pvalue_cis", "locus.diversity")

all.evidence.tbl <- coregenes.tbl[, ..cols.select]
all.evidence.tbl[, cis.validated := ifelse(pvalue_cis < 0.001, TRUE, FALSE)]
all.evidence.tbl[, gwas.hit := fifelse(is.na(reported.genes), FALSE, TRUE)]
all.evidence.tbl[, gwas.validated := fifelse((cis.validated | gwas.hit), "+", "-")]
all.evidence.tbl[is.na(gwas.validated), gwas.validated := "-"]
## 4. validated by the association with the gene product
all.evidence.tbl[is.na(protein.coeff), protein.validated := "."]
all.evidence.tbl[protein.pvalue > 1e-4, protein.validated := "0"]
all.evidence.tbl[is.na(protein.validated), protein.validated := fifelse(sign(protein.coeff) != sign(estimate_trans), "-", "+")]
all.evidence.tbl[, mr.validated := fifelse(mr.pvalue<0.01&locus.diversity>20, "+", ".")]
all.evidence.tbl[, mr.validated := fifelse(mr.pvalue>=0.01&locus.diversity>20, "-", mr.validated)]

## Read other evidence from the files
add.evidence.eqtl <- fread(file.path(working.dir, "gene_hypertension_eQTL.csv"))
setnames(add.evidence.eqtl,
         old = tail(names(add.evidence.eqtl), 3),
         new = c("model.validated", "monogenic.cause", "drug.validated"))
add.evidence.eqtl[,source:="eQTL"]

add.evidence.pqtl <- fread(file.path(working.dir, "gene_hypertension_pQTL.csv"))
setnames(add.evidence.pqtl,
         old = tail(names(add.evidence.pqtl), 3),
         new = c("model.validated", "monogenic.cause", "drug.validated"))
add.evidence.pqtl[,source:="pQTL"]
add.evidence <- rbind(add.evidence.eqtl, add.evidence.pqtl)

all.evidence.tbl <- all.evidence.tbl[add.evidence, on=.(gene_symbol=Gene), .(gene_symbol, gwas.validated, protein.validated, mr.validated, model.validated, monogenic.cause, drug.validated, source)]
# all.evidence.tbl[, model.validated := fifelse(model.validated=="+", "+", "-")]
all.evidence.tbl[, monogenic.cause := fifelse(monogenic.cause=="+", "+", "-")]
all.evidence.tbl[, drug.validated := fifelse(drug.validated=="+", "+", "-")]

cols.select <- c(cols.select, "z_trans", "pvalue_trans", "studyid", "gw.variance")
coeffs.protein.core.joined <- coeffs.protein.core[coregenes.tbl[, ..cols.select], on="Gene==gene_symbol", nomatch=NA] 
coeffs.protein.core.joined[,study:=car::recode(studyid,
                                     "21='FIN'; 31='IMP'; 26='INT'; 1008='DeC'; 50='UKB'; 36='eQTLGen'")]
coeffs.protein.core.joined[,studyid:=NULL]                                    
coeffs.protein.core.joined[, coeff := as.character(coeff)]
coeffs.protein.core.joined[is.na(coeff), coeff := "."]
coeffs.protein.core.joined[is.na(pvalue.formatted), pvalue.formatted := "."]
setorder(coeffs.protein.core.joined, pvalue_trans)

rmarkdown::render(
    'hypertension_paper.Rmd',
    output_dir = 'output',
    output_file = 'hypertension_paper.pdf',
    params = list(output_format = 'pdf')
)
