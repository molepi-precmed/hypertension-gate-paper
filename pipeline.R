source("R/config.R")

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
                              min.pvalue=min.pvalue, select.genes=select.genes)

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
coregenes.pqtl.strict <- coregenes.pqtl[(locus.diversity>=20&pvalue_trans<=1e-5), ]

## Apply correct ordering 

# gene_order_eqtl <- c("IRAG1", "RRP12", "EMP1", "MS4A6A", "DPYSL2", "GAS7", "MEF2C", 
#                 "BCL11A", "COBLL1", "BANK1", "CXCR5", "TNFRSF13B", "TREM1", 
#                 "TP53INP2", "SLC31A2", "SFTPD", "CCN3", "TESC", "UBB", "CD36", 
#                 "RTN1", "KCNJ15", "ALOX5AP", "GPBAR1", "ABCA1", "IL2RB", "COLQ", "DHRS9")

# coregenes.eqtl.strict <- coregenes.eqtl.strict[match(gene_order_eqtl, coregenes.eqtl.strict$gene_symbol), ]
setorder(coregenes.eqtl.strict, pvalue_trans)

# gene_order_pqtl <- c("ELANE", "COCH", "APOD", "ADIPOQ", "CKB", "IGFBP2", "CD300LG", 
#                 "DSG2", "BTN2A1", "VASN", "PTPRF", "IL17RB", "CEACAM16", "CLEC4C", 
#                 "TCL1A", "GCNT1", "CD209", "CPM", "CNTN3", "CTRB1", "CPA1", "REG3G", 
#                 "MYOM3", "EDA2R")

# coregenes.pqtl.strict <- coregenes.pqtl.strict[match(gene_order_pqtl, coregenes.pqtl.strict$gene_symbol), ]

setorder(coregenes.pqtl.strict, pvalue_trans)

## Genes exctarcted from Guzik's review
immune_list <- c(
  "SH2B3",      # also known as LNK
  "CD247",      # T cell receptor component
  "CD86",
  "CD80",
  "CD209",
  "ADGRE5",     # also known as CD97
  "FGFBP2",     # related to cytotoxic T cells
  "PRF1",       # perforin 1, cytotoxic/apoptotic pathway
  "GNLY",       # granulysin, NK/cytotoxic pathways
  "NKG7",       # NK/granzyme related
  "IL2",        # Interleukin-2
  "IL2RA",      # Interleukin-2 receptor alpha
  "IL2RB",      # Interleukin-2 receptor beta
  "IL15",       # Interleukin-15
  "LGALS9",     # Galectin-9, T cell regulation
  "HAVCR2",     # TIM-3
  "IRF5",
  "IRAK1BP1",
  "TRAF1",
  "TXNDC17",
  "PSMA3", "PSMA4", "PSMC3", "PSMC4", "PSMD3", "PSMD5", "SEC31A", # proteasome and antigen presentation genes
  "NLRP3",
  "IL1R2",
  "IL1RAP",
  "IL10RA",
  "IL6", "IL7", "IL9", "IL10",  "IL17A", "IL17RB", "IL18", "IL21",           # key cytokines with evidence for causal or regulatory roles
  "IFNG",       # Interferon gamma
  "TBX21",      # T-bet, TH1 cell regulator
  "TNF",        # Tumor necrosis factor
  "RORC",       # RORgammaT (TH17 differentiation)
  "SGK1",       # salt-sensing kinase
  "CD70", "CD83", "CD28", "CD25", "CTLA4",                        # T cell costimulation and regulatory molecules
  "TLR2", "TLR4", "TLR9",                                         # Toll-like receptors
  "VEGFC", "VEGFD",                                               # lymphangiogenesis
  "CXCR2",                                                        # chemokine receptor
  "NOX1", "NOX2", "NOX4", "NOX5",                                 # NADPH oxidases (oxidative stress)
  "ENaC", "NCC", "NKCC", "NHE3", "Kir4.1", "ClC-K"                # renal sodium transporters (sometimes family/complex, not always a single gene)
)

## Add genes that Paul identified as immune
immune_genes <- fread(file.path(working.dir, "htens_genes.csv"))

immune_list <- unique(c(immune_list, intersect(immune_genes[system=="immune"]$V2, coregenes.eqtl$gene_symbol)))
immune_list <- unique(c(immune_list, intersect(immune_genes[system=="immune"]$V2, coregenes.pqtl$gene_symbol)))

## Adding mannually from the discussion
immune_list <- unique(c(immune_list, "CKB", "BCL11A", "BANK1", "MEF2C", "TREM1", "ALOX5AP", "TNFRSF13B", "IL2RB", "CPA1", "RTN1", "KCNJ15"))



coregenes.guzik <- coregenes.eqtl[gene_symbol %in% immune_list]
coregenes.guzik.pqtl <- coregenes.pqtl[gene_symbol %in% immune_list]
coregenes.guzik <- rbind(coregenes.guzik, coregenes.guzik.pqtl, fill=TRUE)

setorder(coregenes.pqtl.strict, pvalue_trans)

coregenes.eqtl.rest <- coregenes.eqtl[(locus.diversity<8&locus.diversity>=5&pvalue_trans<=1e-5), ]
coregenes.pqtl.rest <- coregenes.pqtl[(locus.diversity<20&locus.diversity>=5&pvalue_trans<=1e-5), ]

phenoname <- "mean arterial pressure"

## Prepare table with validations

## prepare coregenes table
## -----------------------------------------------------------------------------

coregenes.tbl <- rbind(coregenes.eqtl.strict, coregenes.pqtl.strict)

## Join with MR
coregenes.tbl[all.mrres, on=.(gene_symbol=exposure), mr.pvalue:=pvalue]

cols.select <- c("gene_symbol", "reported.genes", "estimate_trans", 
                 "protein.coeff", "protein.pvalue", "mr.pvalue", "pvalue_cis", "locus.diversity")

# gene_symbol, gwas.hit, product.validated,
#                                   mr.validated, model.validated,
#                                   drug.validated

all.evidence.tbl <- coregenes.tbl[, ..cols.select]
all.evidence.tbl[, cis.validated := ifelse(pvalue_cis < 0.001, TRUE, FALSE)]
all.evidence.tbl[, gwas.hit := fifelse(is.na(reported.genes), FALSE, TRUE)]
all.evidence.tbl[, gwas.validated := fifelse((cis.validated | gwas.hit), "+", "-")]
all.evidence.tbl[is.na(gwas.validated), gwas.validated := "-"]
## 4. validated by the association with the gene product
all.evidence.tbl[is.na(protein.coeff), protein.validated := "."]
all.evidence.tbl[protein.pvalue > 1e-4, protein.validated := "0"]
all.evidence.tbl[is.na(protein.validated), protein.validated := fifelse(sign(protein.coeff) != sign(estimate_trans), "-", "+")]
# all.evidence.tbl[protein.coeff>0, protein.validated :=
#             ifelse(protein.coeff>=2*estimate_trans, "+", ".")]
# all.evidence.tbl[protein.coeff<0, protein.validated :=
#             ifelse(protein.coeff<=2*estimate_trans, "+", ".")]

# all.evidence.tbl[, protein.validated := fifelse(abs(protein.coeff)>2*abs(estimate_trans), "+", ".")]
# all.evidence.tbl[is.na(protein.validated), protein.validated := "."]
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
all.evidence.tbl[, model.validated := fifelse(model.validated=="+", "+", "-")]
all.evidence.tbl[, monogenic.cause := fifelse(monogenic.cause=="+", "+", "-")]
all.evidence.tbl[, drug.validated := fifelse(drug.validated=="+", "+", "-")]

## Evidence for immune genes
## Join with MR
coregenes.guzik[all.mrres, on=.(gene_symbol=exposure), mr.pvalue:=pvalue]
all.evidence.immune <- coregenes.guzik[, ..cols.select]
all.evidence.immune[, cis.validated := ifelse(pvalue_cis < 0.001, TRUE, FALSE)]
all.evidence.immune[, gwas.hit := fifelse(is.na(reported.genes), FALSE, TRUE)]
all.evidence.immune[, gwas.validated := fifelse((cis.validated | gwas.hit), "+", "-")]
all.evidence.immune[is.na(gwas.validated), gwas.validated := "-"]
## 4. validated by the association with the gene product
# all.evidence.immune[, protein.validated := "."]
# all.evidence.immune[protein.coeff>0, protein.validated :=
#             ifelse(protein.coeff>=2*estimate_trans, "+", ".")]
# all.evidence.immune[protein.coeff<0, protein.validated :=
#             ifelse(protein.coeff<=2*estimate_trans, "+", ".")]
# all.evidence.immune[, protein.validated :=
#             ifelse(protein.coeff>=2*estimate_trans, "+", ".")]

all.evidence.immune[is.na(protein.coeff), protein.validated := "."]
all.evidence.immune[protein.pvalue > 1e-4, protein.validated := "0"]
all.evidence.immune[is.na(protein.validated), protein.validated := fifelse(sign(protein.coeff) != sign(estimate_trans), "-", "+")]

all.evidence.immune[, mr.validated := fifelse(mr.pvalue<0.01&locus.diversity>20, "+", ".")]
all.evidence.immune[, mr.validated := fifelse(mr.pvalue>=0.01&locus.diversity>20, "-", mr.validated)]

xls.tbl.immune <- all.evidence.immune[, .(gene_symbol, gwas.validated, protein.validated,
                                  mr.validated)]
write_xlsx(xls.tbl.immune, "Core_genes_immune_evidence.xlsx")

cols.select <- c(cols.select, "z_trans", "pvalue_trans", "studyid", "gw.variance")
coeffs.protein.core.joined <- coeffs.protein.core[coregenes.tbl[, ..cols.select], on="Gene==gene_symbol", nomatch=NA] 
coeffs.protein.core.joined[,study:=car::recode(studyid,
                                     "21='FIN'; 31='IMP'; 26='INT'; 1008='DeC'; 50='UKB'; 36='eQTLGen'")]
coeffs.protein.core.joined[,studyid:=NULL]                                    
coeffs.protein.core.joined[, coeff := as.character(coeff)]
coeffs.protein.core.joined[is.na(coeff), coeff := "."]
coeffs.protein.core.joined[is.na(pvalue.formatted), pvalue.formatted := "."]
setorder(coeffs.protein.core.joined, pvalue_trans)

## Add references
## eQTL
# all.evidence.tbl[gene_symbol == "TREM1", `:=`(
#   drug.validated = "+\\cite{francoisProspectiveEvaluationEfficacy2023}"
# )]
# # all.evidence.tbl[gene_symbol == "GPBAR1", `:=`(
# #   drug.validated = "+\\cite{heTakedaProteinCoupled2025}"
# # )]
# all.evidence.tbl[gene_symbol == "ABCA1", `:=`(
#   monogenic.cause = "+\\cite{rhyneMultipleSpliceDefects2009}"
# #   ,
# #   drug.validated = "+\\cite{choiBiomedicalAdvancesABCA12023}"
# )]

# ## pQTL
# all.evidence.tbl[gene_symbol == "IL17RB", `:=`(
#   drug.validated = "+\\cite{davisInterleukin17AKey2021}"
# )]
# all.evidence.tbl[gene_symbol == "ADIPOQ", `:=`(
#   monogenic.cause = "+\\cite{simeoneDominantNegativeADIPOQ2022}"
# #   ,
# #   drug.validated = "+\\cite{raisulabedin303ORNotchSignaling2024}"
# )]
# # all.evidence.tbl[gene_symbol == "CKB", `:=`(
# #   drug.validated = "+\\cite{brewsterCreatineKinaseEnergy2018}"
# # )]
# all.evidence.tbl[gene_symbol == "DSG2", `:=`(
#   monogenic.cause = "+\\cite{sumidaFourCardiomyopathyPatients2024}"
# )]
# all.evidence.tbl[gene_symbol == "MYOM3", `:=`(
#   monogenic.cause = "+\\cite{rouillonSerumProteomicProfiling2015}"
# )]
# all.evidence.tbl[gene_symbol == "ELANE", `:=`(
#   monogenic.cause = "+\\cite{tidwellNeutropeniaassociatedELANEMutations2014}"
# #   ,
# #   drug.validated = "+\\cite{makaryanElastaseInhibitorsPotential2017}"
# )]

# all.evidence.tbl[gene_symbol == "ELANE", `:=`(
#   monogenic.cause = "+\\cite{tidwellNeutropeniaassociatedELANEMutations2014}"
# #   ,
# #   drug.validated = "+\\cite{makaryanElastaseInhibitorsPotential2017}"
# )]

# all.evidence.tbl[gene_symbol == "EDA2R", `:=`(
#   monogenic.cause = "+\\cite{farringtonRolesEDA2RAgeing2025}"
# )]

# # all.evidence.tbl[gene_symbol == "COCH", `:=`(
# #   model.validated = "+\\cite{carreonInteractionCochlinMechanosensitive2017}"
# # #   ,
# # #   drug.validated = "+\\cite{vriezeAllelespecificAntisenseOligonucleotide2020}"
# #  )]

rmarkdown::render(
    'hypertension_paper.Rmd',
    output_dir = 'output',
    output_file = 'hypertension_paper.pdf',
    params = list(output_format = 'pdf')
)
