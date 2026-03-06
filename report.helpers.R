#' Convert vector of chromosome symbols to numbers
#'
#' @param x Vector of chromosome symbols (1:22, X, Y, XY, MT).
#' @return
#' A vector of chromosome numbers (1:26).
#' @keywords internal
#' @export
chrom.to.integer <- function(x) {
    ## from https://www.cog-genomics.org/plink/2.0/input:
    ## when there are n autosome pairs, the X chromosome is assigned numeric code n+1,
    ## Y is n+2, XY (pseudo-autosomal region of X) is n+3, and MT (mitochondria) is n+4.
    symbols <- c(1:22, "X", "Y", "XY", "MT")
    return(match(x, symbols))
}

#' Convert vector of chromosome numbers to symbols
#'
#' @param x Vector of chromosome numbers (1:26).
#' @return
#' A vector of chromosome symbols (1:22, X, Y, XY, MT).
#' @keywords internal
#' @export
chrom.to.symbol <- function(x) {
    symbols <- c(1:22, "X", "Y", "XY", "MT")
    return(symbols[x])
}

#' Create clumps of SNPs or QTLs. This can be used in one of two instances.
#' 1. To look through all regional QTL scores and create QTL clumps. Each clump
#' can be a QTL for 1 GWAS (i.e., 1 trait) or consist of overlapping QTLs for
#' multiple traits.
#' 2. To create GWAS hit clumps.
#'
#' @param qtls Data table containing the necessary genomic information such as
#'        chromosome number ('chrom_int'), start and end positions
#'        ('startpos', 'endpos') for each locus score or GWAS hit.
#' @param maxgap Integer value that specifies the maximum allowable gap
#'        (in base pairs, BP) between loci for them to be considered in the same
#'        clump. Default is 10,000 BP (10 kb).
#'
#' @details
#' The function checks for the presence of the required columns
#' ('chrom_int', 'startpos', 'endpos'). It converts the 'chrom_int' values into
#' a symbolic representation (e.g., 1-22, X, Y) and orders the data by chromosome
#' and start position. Absolute genomic positions are calculated followed by
#' clump assignment where any loci separated by less than or equal to the
#' `maxgap` are grouped together. Changes in chromosome or excessive gaps
#' between loci initiate a new clump. The function also prints out the count of
#' regions that overlap and those that do not.
#'
#' @examples
#' # Assuming 'QTLdata' is a pre-loaded data table with the necessary columns:
#' clumped_QTLs <- clump(QTLdata, maxgap=10000)
#'
#' @importFrom data.table data.table
#' @importFrom genoscores chrom.to.symbol
#' @import data.table
clump <- function(qtls, maxgap=1E4) {
    require(data.table)
    required.cols <- c("chrom_int", "startpos", "endpos")
    check.cols(required.cols, qtls)
    if (class(qtls$chrom_int) != "integer")
        stop("Column chrom_int must be integer.")

    qtls[, chrom := factor(chrom.to.symbol(chrom_int),
                           levels=c(1:22, "X", "Y"))]
    setorder(qtls, chrom, startpos)
    qtls[, startpos := as.integer(startpos)]
    qtls[, endpos := as.integer(endpos)]
    qtls[, x.start := position.absolute(chrom, startpos)]
    qtls[, x.end := position.absolute(chrom, endpos)]
    qtls[, x.end.previous.max := data.table::shift(cummax(x.end))]
    qtls[1, x.end.previous.max := 0]
    ## allow gaps of maxgap (in kb) for overlap
    qtls[, overlap := x.start < x.end.previous.max + maxgap]
    qtls[1, overlap := 0]
    qtls[, newclump := 1 - overlap]
    ## create a new clump if the chromosome changes
    qtls[, previous.chrom := c(0, qtls[-.N, chrom])]
    qtls[previous.chrom != chrom, newclump := 1]
    qtls[, chrom := as.character(chrom)]

    cat("The following ranges overlap\n")
    print(qtls[, .N, by=overlap])

    qtls[, clump := cumsum(newclump)]
    return(qtls)
}

#' Clump GWAS Hits by Proximity
#'
#' This function clumps trait-associated SNPs from GWAS data based on spatial
#' proximity. It creates a 20kb window around each GWAS hit and then clumps
#' these windows together if they appear within the specified distance from each
#' other. Such clumps will be used to annotate core genes and trans-
#' QTLs contributing to these core genes.
#'
#' The function assumes GWAS hits were downloaded from the GWAS catalog.
#'
#' @param GWAShits.dt Data table containing GWAS hits. At least the following
#'        columns are required:  `P-VALUE`, `CHR_ID`, `CHR_POS`,
#'        `REPORTED GENE(S)` `MAPPED_GENE`, `SNPS`, `DISEASE/TRAIT`, `STUDY`.
#' @param clump.gap A numeric value defining the max distance (in base pairs)
#'        between GWAS hits for clumping them togather.
#' @param phenotypes Character vector with the names of the phenotypes in
#'        `DISEASE/TRAIT` column. Only GWAS hits for specified phenotypes will
#'        be considered for clumping. If is NULL, the function will assume only
#'        one phenotype and consider all GWAS hits in GWAShits.dt for clumping.
#'
#' @return A data table with clumped GWAS hits.
#'
#' @details
#' The function starts by counting all SNP associations and studies in the
#' dataset. It then filters the dataset based on the provided phenotypes, if
#' any. All hits are then filtered to retain only those with a p-value below a
#' genome-wide significance threshold of 5E-8. The function ensures valid
#' genomic positions for these SNPs and defines a clumping region around each
#' SNP based on the provided gap size. It then clumps together SNPs that are
#' within this gap. For each clump, it summarizes the data and provides a
#' detailed annotation of involved genes and traits. Finally, the function
#' returns the clumped SNP data as a data table.
#'
#' @examples
#' # Assume `GWAShitsData` is a data.table already loaded with GWAS hits data.
#' clumped_data <- clump.GWAShits(GWAShitsData, clump.gap=20,000, phenotypes=NULL)
#'
#' @importFrom data.table data.table
#' @importFrom genoscores chrom.to.integer
#' @importFrom clump clump
clump.GWAShits <- function(GWAShits.dt, clump.gap, phenotypes) {
    required.cols <- c("P-VALUE", "CHR_ID", "CHR_POS", "REPORTED GENE(S)",
                       "MAPPED_GENE", "SNPS", "DISEASE/TRAIT", "STUDY")
    check.cols(required.cols, GWAShits.dt)

    cat("Processing and clumping GWAS hits.\n")
    msg <- sprintf("Started from %s SNP associations in %s studies.\n",
                   GWAShits.dt[, uniqueN(SNPS)], GWAShits.dt[, uniqueN(STUDY)])
    cat(msg)

    ## subset relevant phenotypes if specified
    if (is.null(phenotypes)) {
        msg <- paste("The list of phenotypes of interest is not provided.",
            "Assuming all hits are for the relevant phenotypes.\n")
        cat(msg)
    }

    if (!is.null(phenotypes)) {
        msg <- sprintf("Proceeding for specified phenotypes: \n %s\n",
                       paste0(phenotypes, collapse="; "))
        cat(msg)

        ## ensure we match traits in GWAS Catalog with phenotypes
        GWAShits.dt <- GWAShits.dt[tolower(`DISEASE/TRAIT`) %in%
                                   tolower(phenotypes)]

        if (nrow(GWAShits.dt)==0) {
            stop("Specified phenotypes were not found in `DISEASE/TRAIT` column.
            Please specify correct phenotype names or set it to NULL.")
        }

        msg <- sprintf("%s SNPs remain in %s studies of relevant phenotypes.\n",
                       GWAShits.dt[, uniqueN(SNPS)],
                       GWAShits.dt[, uniqueN(STUDY)])
        cat(msg)
    }

    GWAShits.dt <- GWAShits.dt[as.numeric(`P-VALUE`) < 5E-8]
    msg <- sprintf("%s SNPs are genome-wide significant.\n",
                   GWAShits.dt[, uniqueN(SNPS)])
    cat(msg)

    ## ensure chrom_int is an integer
    GWAShits.dt[, `:=`(chrom_int=chrom.to.integer(CHR_ID),
                       pos=CHR_POS)]

    ## take only those GWAS hits that have both chromosome and position
    GWAShits.dt <- GWAShits.dt[!is.na(chrom_int) & !is.na(pos)]
    msg <- sprintf("%s SNPs have valid genomic positions.\n",
                   GWAShits.dt[, uniqueN(SNPS)])
    cat(msg)

    ## define a 20kb window around top SNP
    GWAShits.dt[, `:=`(startpos=as.integer(pos)-2E4, endpos=as.integer(pos)+2E4)]

    ## clump the SNP hit windows
    GWAS.clumps <- clump(GWAShits.dt, maxgap=clump.gap)
    if (any(is.na(GWAS.clumps$clump))) stop("GWAS clump cannot be NA.")

    msg <- sprintf("Defined %s trait associated regions from %s SNPs.\n",
                   GWAS.clumps[, uniqueN(clump)], GWAS.clumps[, uniqueN(SNPS)])
    cat(msg)

    ## reported.genes are reported by the study, mapped.genes are assigned by
    ## position by GWAS catalog
    setnames(GWAS.clumps, old=c("REPORTED GENE(S)", "MAPPED_GENE"),
             new=c("reported.genes", "mapped.genes"))

    ## if reported genes are NULL or empty, we replace them with 'NR'
    GWAS.clumps[is.na(reported.genes) | is.null(reported.genes) |
        reported.genes == "", reported.genes := "NR"]

    ## large clumps can map to multiple SNP hits, we collapse the information
    ## per clump
    GWAS.clumps <- GWAS.clumps[,
        .(
            chrom = chrom[[1]],
            startpos = min(startpos),
            endpos = max(endpos),
            x.clumpStart = min(x.start),
            x.clumpEnd = max(x.end),
            reported.genes = paste0(unique(
                reported.genes[!is.na(reported.genes)]), collapse=", "),
            mapped.genes = paste0(unique(mapped.genes[!is.na(mapped.genes)]),
                                  collapse=", "),
            pvalue = min(`P-VALUE`),
            snps = paste0(SNPS, collapse=", "),
            trait = paste0(unique(MAPPED_TRAIT), collapse=", ")
        ),
        by = clump
    ]
    cat("Finished clumping GWAS hits.\n")
    return(GWAS.clumps)
}

#' For each gene, whose GATE score was tested for association with disease and
#' whose transcription site overlaps a clump of GWAS hit SNPs or a monogenic
#' cause of disease, this function reports GWAS gene names or monogenic genes
#' near tested gene transcription site.
#'
#' @param coeff.dt Data table with summary statystics from the tests of
#'        association between phenotype and GATE scores.
#'        The required columns are:
#'        - `gene_chrom` integer specifying the chromosome of the transcription
#'          site of the gene corresponding to the GATE score.
#'        - `gene_startpos` numeric value specifying location of the start of
#'          the gene transcription site.
#'        - `gene_endpos` numeric value specifying location of the end of
#'          the gene transcription site.
#' @param hits.dt Data table with GWAS hits clumped by position or
#'        monogenic causes of a diseases. The required columns are:
#'        - `x.clumpStart` numeric specifying the absolute position of the
#'           beginning of the clump
#'        - `x.clumpEnd` numeric specifying the absolute position of the
#'           end of the clump
#'        - `reported.genes` or `monogenic` character specifying gene symbols
#'          reported as GWAS hits or as monogenic causes of disease.
#' @param gene.column Character specifying the name of the column containing
#'        gene name for the hit. If overlapping with GWAS hits, this column
#'        should be called 'reported.genes'. If overlapping with monogenic
#'        disease causes, this column should be called 'monogenic'.
#' @param overlap.flank Numeric specifying the window to search for hits or
#'        monogenic causes that lie just outside the clump. This window extends
#'        to both sides of a clump.
#'
#' @return Data table with GWAS hit genes (if `gene.column="reported.genes"`) or
#' monogenic causes of disease (if `gene.column="monogenic"`) which overlap
#' with the transcription site of a gene corresponding to each GATE score.
overlap.coeffs.hits <- function(coeff.dt, hits.dt,
                                gene.column=c("reported.genes", "monogenic"),
                                overlap.flank) {
    cat("Overlapping gene annotations with GATE associations.\n")
    ## check coeffs.dt table
    if (class(coeff.dt)[[1]] != "data.table")
        stop("A data.table containing associations with scores is required.")

    cols.required <- c("gene_chrom", "gene_startpos", "gene_endpos")
    check.cols(cols.required, coeff.dt)

    ## check provided hits table
    if (class(hits.dt)[[1]] != "data.table")
        stop("A data.table containing hits is required.")

    cols.required <- c("x.clumpStart", "x.clumpEnd")
    check.cols(cols.required, hits.dt)

    ## overlap coeffs.wide with hits.dt
    if (!coeff.dt[, all(is.numeric(gene_startpos))] ||
        !coeff.dt[, all(is.numeric(gene_endpos))]) {
        stop("Columns gene_startpos and gene_endpos should be numeric.")
    }
    coeff.dt[, x.gene_startpos := position.absolute(gene_chrom, gene_startpos)]
    coeff.dt[, x.gene_endpos := position.absolute(gene_chrom, gene_endpos)]

    ## add new columns for foverlaps with the flank
    coeff.dt[, x.gene.startpos.with.flank := x.gene_startpos - overlap.flank]
    coeff.dt[, x.gene.endpos.with.flank := x.gene_endpos + overlap.flank]
    coeff.dt[x.gene.startpos.with.flank < 0, x.gene.startpos.with.flank := 0]

    coeff.dt[x.gene_startpos<0, x.gene_startpos := 0]

    setkey(hits.dt, x.clumpStart, x.clumpEnd)
    setkey(coeff.dt, x.gene.startpos.with.flank, x.gene.endpos.with.flank)
    coeff.hit.overlaps <- foverlaps(coeff.dt,
                                    hits.dt[, .(x.clumpStart, x.clumpEnd,
                                                get(gene.column))],
                                    by.x=c("x.gene.startpos.with.flank",
                                           "x.gene.endpos.with.flank"),
                                    by.y=c("x.clumpStart", "x.clumpEnd"),
                                    type="any")
    coeff.hit.overlaps[, `:=`(x.clumpStart=NULL, x.clumpEnd=NULL)]
    setnames(coeff.hit.overlaps, old="V3", new=eval(gene.column))

    ## remove columns added for foverlaps
    coeff.hit.overlaps[, `:=`(x.gene.startpos.with.flank=NULL,
                              x.gene.endpos.with.flank=NULL)]

    return(coeff.hit.overlaps)
}

#' Overlap clumps table with either a table of GWAS hits or a monogenic table.
#'
#' @param trans.clumps Data table with  QTL clumps for a given study (e.g.
#'        protein, transcripts). Minimum columns expected - x.clumpStart,
#'        x.clumpEnd.
#' @param hits.dt Data table with clumped GWAS hits. Minimum columns expected -
#'        x.clumpStart, x.clumpEnd
#' @param monogenic Data table with monogenic causes of disease.
#'        Minimum columns expected - x.clumpStart, x.clumpEnd
#' @param overlap.flank Numeric specifying the window to search for hits or
#'        monogenic causes that lie just outside the clump. This window extends
#'        to both sides of a clump.
#' @param monogenic.available Logical to indicate if the disease being analysed
#'        has monogenic causes.
#'
#' @return updated clumps data.table with reported GWAS hits and monogenic
#' causes of disease.
overlap.qtls.hits <- function(trans.clumps, hits.dt, monogenic, overlap.flank,
                              monogenic.available = TRUE) {
    ## trans.clumps and hits must have required columns in correct format
    required.cols <- c("chrom", "startpos", "endpos", "x.clumpStart",
                       "x.clumpEnd")
    check.cols(required.cols, trans.clumps)
    check.cols(required.cols, hits.dt)
    if (monogenic.available) {
        check.cols(required.cols, monogenic)
    }

    ## hits must be unique, if they are not, throw an error and stop
    if (any(duplicated(hits.dt, by=c("chrom", "x.clumpStart", "x.clumpEnd"))))
        stop("Hits table has duplicated clumps. Check annotation.")

    ## add columns with foverlaps with the flank
    trans.clumps[, x.clumpStart.with.flank := x.clumpStart - overlap.flank]
    trans.clumps[, x.clumpEnd.with.flank := x.clumpEnd + overlap.flank]
    trans.clumps[x.clumpStart.with.flank < 0, x.clumpStart.with.flank := 0]

    ## find the overlaps between QTL clumps and GWAS hit clumps keeping all
    ## trans.clumps irrespective of whether they overlap hits
    setkey(trans.clumps, x.clumpStart.with.flank, x.clumpEnd.with.flank)
    setkey(hits.dt, x.clumpStart, x.clumpEnd)

    clump.hit.overlaps <- foverlaps(trans.clumps, hits.dt,
        by.x = c("x.clumpStart.with.flank", "x.clumpEnd.with.flank"),
        by.y = c("x.clumpStart", "x.clumpEnd"),
        type = "any"
    )
    clump.hit.overlaps <- clump.hit.overlaps[,
        .(clump=i.clump, chrom=i.chrom, x.clumpStart=i.x.clumpStart,
          x.clumpEnd=i.x.clumpEnd, startpos=i.startpos, endpos=i.endpos,
          target.genes, reported.genes, nearby_genes,
          x.clumpStart.with.flank, x.clumpEnd.with.flank)]

    if (monogenic.available) {
        setkey(monogenic, x.clumpStart, x.clumpEnd)

        clump.hit.overlaps <- foverlaps(
            clump.hit.overlaps,
            monogenic[, .(x.clumpStart, x.clumpEnd, monogenic)],
            by.x = c("x.clumpStart.with.flank", "x.clumpEnd.with.flank"),
            by.y = c("x.clumpStart", "x.clumpEnd"),
            type="any")
        clump.hit.overlaps <- clump.hit.overlaps[,
            .(clump, chrom, x.clumpStart=i.x.clumpStart,
            x.clumpEnd=i.x.clumpEnd, startpos, endpos,
            target.genes, reported.genes, monogenic, nearby_genes)]
    }

    ## overlap produces duplicated clumps were several GWAS hit regions map
    ## to the same clump, so we should combine these into 1 record per clump
    clumps.with.hits <- clump.hit.overlaps[,
        .(
            chrom = chrom[[1]],
            x.clumpStart = x.clumpStart[[1]],
            x.clumpEnd = x.clumpEnd[[1]],
            startpos = min(startpos),
            endpos = max(endpos),
            target.genes = paste0(unique(target.genes), collapse=", "),
            reported.genes = paste0(unique(
                reported.genes[!is.na(reported.genes)]), collapse=", "),
            monogenic = paste0(
                unique(monogenic[!is.na(monogenic)]),
                collapse=", "),
            nearby.genes = paste0(unique(nearby_genes), collapse=", ")
        ),
        by = clump
    ]
    if(!monogenic.available) {
        clumps.with.hits[, monogenic := NA]
    }
    return(clumps.with.hits)
}

get.metadata.transqtls.gwasids <- function(all.scoresinfo, ans, gwasids) {
  scoresinfo.bylocus <- all.scoresinfo ## locus-specific scores
  scoresinfo.aggregated <- ans # cis-, cis-x- and genomewide trans- scores only

  # get metadata only where minpvalue < threshold
  transqtls.metadata <-
    scoresinfo.bylocus[!is.na(startpos) &
                       !(scoreid %in% ans[substr(qtl_type, 1, 3)=="cis", score.name])]
  transqtls.metadata.gwasids <- transqtls.metadata[gwasid %in% gwasids,
                                                     .(gwasid, chrom_int, chrom_min, startpos, endpos,
                                                       nearby_genes)]
  transqtls.metadata.gwasids[, nearby_genes := gsub(",havana_tagene", "", nearby_genes)]
  transqtls.metadata.gwasids[, nearby_genes := gsub(",havana", "", nearby_genes)]

  return(transqtls.metadata.gwasids)
}

process.coeffs <- function (study, gate.dir, working.dir,
                            coeffs.wide.file,
                            min.diversity=5, min.pvalue=1e-6,
                            select.genes=NULL) {
    study.dir <- file.path(gate.dir, study)

    coeffs.wide <- fread(file.path(gate.dir, coeffs.wide.file))
    locus.meta <- load.rdata(file.path(study.dir,
                                    "all.scoresinfo.annotated.Rdata.gz"))
    gate.meta <- load.rdata(file.path(study.dir,
                                    "trans.scoresinfo.annotated.Rdata.gz"))

    coregenes <- coeffs.wide[pvalue_trans<min.pvalue &
                             locus.diversity>min.diversity]
    if (!is.null(select.genes)) {
        coregenes <- rbind(coregenes,
                           coeffs.wide[gene_symbol %in% select.genes])
    }
    ## select unique genes
    # coregenes <- coregenes[, .SD[which.min(pvalue_trans)], by=gene_symbol]

    if (nrow(coregenes)<5) {
        coregenes.top <-
            coeffs.wide[order(pvalue_trans)][locus.diversity>min.diversity][1:5]
        coregenes.top <- coregenes.top[!gene_symbol %in% coregenes]
        coregenes <- rbind(coregenes, coregenes.top)
    }
    dir.create(file.path(working.dir, study), showWarnings=FALSE,
               recursive=TRUE)
    fwrite(coregenes, file=file.path(working.dir, study, "coregenes.csv"))

    index.name <- switch(study,
                         eqtl="trans.genotypicscore.1e-5.col.idx",
                         pqtl="trans.genotypicscore.1e-6.col.idx")
    col.idx <- fread(file.path(study.dir, index.name))
    col.idx.coregenes <- col.idx[score.name %in% coregenes$scoreid_trans]
    col.idx.coregenes[coregenes, on=.(score.name=scoreid_trans),
                      coregene.colnames:=gene_symbol]
    fwrite(col.idx.coregenes, file=file.path(working.dir, study, "col.idx.coregenes.csv"))
    return(list(coregenes=coregenes, coeffs.wide=coeffs.wide,
                locus.meta=locus.meta, gate.meta=gate.meta))
}

annotate.genelists <- function(genelists, hits) {
  ## genelists is a character vector in which each element is a comma-separated list of genes
  x <- strsplit(genelists, ",") # split into a list of character vectors
  for(i in 1:length(x)) {
    y <- unique(x[[i]])
    y <- y[grep("\\\\$|ensemble|havana|tagene|Y\\_RNA", y, invert=TRUE)]
    hits.indices <- which(y %in% hits)
    nohits.indices <- which(!y %in% hits)
    truncate <- length(nohits.indices) > 5
    y[hits.indices] <- paste0("\\textbf{", y[hits.indices], "}")
    indices.keep <- unique(c(hits.indices, nohits.indices[1:5])) # keep up to 5 other
    indices.keep <- sort(indices.keep)
    y <- y[indices.keep]
    y <- gsub("\\_", "\\\\_", y)
    x[[i]] <- paste(y, collapse=", ")
    if(truncate) {
      x[[i]] <- paste0(x[[i]], ", ...")
    }
  }
  x <- unlist(x)
  x[x=="" | x=="NA" | is.na(x)] <- "."
  return(x)
}

create.clumps <- function(study, coregenes, gate.meta, locus.meta, gwas.hits,
                          working.dir) {
    metadata <- process.gate.metadata(gate.meta, locus.meta)
    metadata <- metadata$collapsed.scores.info
    metadata.transqtls.coregenes <-
        get.metadata.transqtls.gwasids(all.scoresinfo=locus.meta, ans=gate.meta,
                                       gwasids=coregenes$gwasid)

    metadata.transqtls.coregenes <- clump(metadata.transqtls.coregenes)
    cat("trans-QTLs for core genes fall into",
        max(metadata.transqtls.coregenes$clump), "clumps\n")

    metadata.transqtls.coregenes <-
        metadata[, .(gwasid, gene_symbol)][metadata.transqtls.coregenes,
                                           on="gwasid"]

    metadata.transqtls.coregenes[, chrom_min := as.integer(chrom_min)]
    clumps.transqtls.coregenes <-
    metadata.transqtls.coregenes[, .(clump,
                                    chrom_int=min(chrom_int), chrom=min(chrom_min),
                                    clumpStart=min(startpos),
                                    clumpEnd=max(endpos),
                                    target.genes=paste(unique(gene_symbol), collapse=","),
                                    nearby_genes=paste(nearby_genes[1], collapse=",")),
                                by=clump]
                                            # clumps.transqtls.coregenes[, nearby_genes := deduplicate.split(nearby_genes)]
    clumps.transqtls.coregenes[, .(clump, chrom, clumpStart, clumpEnd, target.genes, nearby_genes)]
    clumps.print <-
    clumps.transqtls.coregenes[, .(clump, chrom,
                                    start=round(clumpStart * 1E-6, 2),
                                    stop=round(clumpEnd * 1E-6, 2),
                                    targetgenes=gsub(",", ", ", target.genes),
                                    nearby_genes=annotate.genelists(nearby_genes,
                                                                    hits=gwas.hits))]
    clumps.file <- file.path(working.dir, study, "clumps.csv")
    fwrite(clumps.print, file=clumps.file)
}

dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriate(d_corr, method = method, ...)
  i = get_order(s)
  return(i)
}

#' Helper function to compute correlations between GATE scores for core genes
#' discovered in a specific study (e.g., eQTLGen, deCODE, UKB-PPP)
#'
#' @param study String specifying the name of the study in which core genes were
#'        discovered. This assumes that GATE scores for each study are saved in
#'        the respective subdirectory: <parent_dir>/<study>.
#'        For instance, if study is eQTLGen, the scores should be in
#'        <parent_dir>/eQTLGen.
#' @param phenoname String specifying the name of the phenotype in question. It
#'        should exists in column names of pheno.dt.
#' @param pheno.dt Data table with the phenotype to be used to subset scores to
#'        compute correlations in controls only.
#' @param gate.dir Full path to the parent directory containing GATE results for
#'        each study.
#' @param working.dir Full path to the working directory which should contian
#'        study subdirectory containing a file with the associations between
#'        GATE scores and the phenotype. Score correlations will be saved into
#'        this directory as well.
#' @param coregenes.file Name of the file containing associations between
#'        GATE scores and the phenotype.
compute.score.correlations <- function(study, phenoname, pheno.dt, gate.dir,
                                       working.dir, coregenes.file) {
    require(data.table)
    require(bigstatsr)
    require(ggcorrplot)

    ## read coeffs wide for identified core genes from a file
    coregenes <- fread(file.path(working.dir, study, coregenes.file))
    if (nrow(coregenes)<2) {
        cat("Less than 2 core genes. Correlations won't be computed.\n")
        return(NULL)
    }
    coregenes <- coregenes[,
        .SD[locus.diversity==max(locus.diversity) |
            pvalue_trans == min(pvalue_trans)],
        by=gene_symbol]

    ## attach big matrix of GATE scores
    gate.dir <- file.path(gate.dir, study)

    file1 <- file.path(gate.dir, 'trans.genotypicscore.1e-6.rds')
    file2 <- file.path(gate.dir, 'trans.genotypicscore.1e-5.rds')

    if (file.exists(file1)) {
        selected_file <- file1
    } else if (file.exists(file2)) {
        selected_file <- file2
    }
    genome.wide.scores <- big_attach(selected_file)

    ## read column idex of GATE matrix and subset the columns to keep only those
    ## corresponding to GATE scores
    file1 <- file.path(gate.dir, "trans.genotypicscore.1e-6.col.idx")
    file2 <- file.path(gate.dir, "trans.genotypicscore.1e-5.col.idx")

    if (file.exists(file1)) {
        selected_file <- file1
    } else if (file.exists(file2)) {
        selected_file <- file2
    }
    col.idx.coregenes <- fread(selected_file)
    col.idx.coregenes <- col.idx.coregenes[score.name %in% coregenes$scoreid_trans]
    col.idx.coregenes[coregenes, on=.(score.name=scoreid_trans),
                      gene_symbol:=gene_symbol]

    ## read row index of GATE matrix and subset the rows to keep only those
    ## corresponding to controls
    file1 <- file.path(gate.dir, "trans.genotypicscore.1e-6.row.idx")
    file2 <- file.path(gate.dir, "trans.genotypicscore.1e-5.row.idx")
    if (file.exists(file1)) {
        selected_file <- file1
    } else if (file.exists(file2)) {
        selected_file <- file2
    }
    row.idx <- fread(selected_file)
    y.isbinary <- setequal(length(unique(na.omit(pheno.dt[[phenoname]]))), 2)
    if (y.isbinary) {
        row.idx.control <- row.idx[row.name %in%
                               pheno.dt[get(phenoname)==0 & exclude==FALSE,
                               row.name]]
    } else {
        row.idx.control <- row.idx
    }

    ## subset the GATE matrix and convert to data.tableassigning gene names
    ## to columns which will be used as labels on correlation plots
    coregene.scores <- as.data.table(genome.wide.scores[row.idx.control$row.idx,
                             col.idx.coregenes$gw.container.id])
    colnames(coregene.scores) <- col.idx.coregenes$gene_symbol

    ## compute and plot correlations
    allqtl.corr.core <- cor(coregene.scores)
    allqtl.corr.core <- allqtl.corr.core - diag(nrow(allqtl.corr.core)) # set diag elements to 0
    maxcorr <- apply(abs(allqtl.corr.core), 1, max)
    maxcorr <- data.table(gene_symbol=names(maxcorr), maxcorr=as.numeric(maxcorr))
    coregenes <- maxcorr[coregenes, on="gene_symbol"]
    d_corr <- as.dist(1 - allqtl.corr.core)  # 1) distance matrix
    hc <- hclust(d_corr, method = "complete")      # 2) hierarchical clustering
    reorder.i <- hc$order                           # 3) extract order of leaves
    # reorder.i <- dist2order(abs(allqtl.corr.core), method="HC")
    allqtl.corr.core <- allqtl.corr.core[rev(reorder.i), reorder.i] # reverse rows for ggcorrplot

    coregenes <- coregenes[match(rownames(allqtl.corr.core), gene_symbol)]
    coregenes[maxcorr > 0.8, gene_symbol := paste0("(", gene_symbol, ")")]
    options(warn=-1)
    ## FIXME: https://stackoverflow.com/questions/67530405/change-orientation-of-diagonal-of-correlation-plot-using-ggcorrplot-package-if

    p.allqtl.corr.core <- ggcorrplot::ggcorrplot(allqtl.corr.core, hc.order=FALSE,
                                                 method="circle", lab=TRUE,
                                                 lab_size=2,
                                                 sig.level=0.3, digits=1,
                                                 type="full",
                                                 show.diag=TRUE) +
    theme(axis.text.y = element_text(face="italic", size=8)) +
    theme(axis.text.x = element_text(face="italic", size=8)) +
    theme(legend.position="bottom")
    options(warn=2)

    ## save the plot
    plot.file <- file.path(working.dir, study,
                              sprintf("%s_hits_correlation.png", study))
    ggsave(filename=plot.file,
            p.allqtl.corr.core, device="png",
            dpi=300, width=7, height=7, units="in")
}

#' pheno.test is name of outcome variable
#' pheno.proteins is dataset
get.coeff.protein <- function(gene, pheno.proteins,
                              formula=formula, y.isbinary) {
    options(warn=1) # tryCatch fails to catch a warning if it is converted to an error
    gene.formula <- as.formula(sprintf("%s + %s", formula, gene))
    coeff.dt <- tryCatch( {
        gene.model <- glm(pheno.proteins,
                          formula=gene.formula,
                          family=ifelse(y.isbinary, "binomial", "gaussian"))
        coeff.dt <- as.data.table(summary(gene.model)$coefficients,
                                  keep.rownames="variable")
        if(y.isbinary) {
            y.freqs <- table(gene.model$y)
            coeff.dt[, noncases := y.freqs[1]]
            coeff.dt[, cases := y.freqs[2]]
        } else {
            coeff.dt[, N := length(gene.model$y)]
        }
        coeff.dt[variable==gene]
    },
    error=function(cond) {
        message(paste("model fitting failed for protein", gene))
        message(paste(cond, "\n"))
        return(data.table(variable=gene))
    },
    warning=function(cond) {
        message(paste("model fitting for protein", gene, "caused a warning"))
        message(paste(cond, "\n"))
        return(data.table(variable=gene))
    },
    finally={}
    )
    return(coeff.dt)
}

get.coeffs.protein <- function(genenames.proteins, pheno.proteins,
                               formula, y.isbinary, cores=10) {
    require(foreach)
    require(doParallel)
    registerDoParallel(cores=cores)
    coeffs.protein <- foreach(i = 1:length(genenames.proteins),
                              .combine=function(...) rbind(..., fill = TRUE),
                              .multicombine = TRUE) %do% {
                                  get.coeff.protein(gene=genenames.proteins[i],
                                                    pheno.proteins, formula,
                                                    y.isbinary)
                              }
    if (y.isbinary) {
        colnames(coeffs.protein) <- c("Gene", "coeff", "SE", "z", "pvalue",
                                      "noncases", "cases")
    } else {
        colnames(coeffs.protein) <- c("Gene", "coeff", "SE", "z", "pvalue", "N")
    }
    proteins.sd <- apply(X=pheno.proteins[, ..genenames.proteins],
                         MARGIN=2, FUN=sd, na.rm=TRUE)
    coeffs.protein[, sd.protein := proteins.sd]
    coeffs.protein[, pvalue.formatted := format.z.aspvalue(z)]
    coeffs.protein[, coeff := coeff * sd.protein]
    coeffs.protein[, SE := SE * sd.protein]
    return(coeffs.protein)
}

report.mr <- function(core.genes, scoresinfo, pheno.dt, phenoname, covariates,
                      study.name, working.dir, gwas.dt, basegatedir=basegatedir,
                      ld.refplinkfile=ld.refplinkfile, cores=20) {
    source("mr.R")
    n_workers <- cores
    ## Use socket clusters (better memory management than fork) and output
    ## messages directly to console
    cl <- makeCluster(n_workers, type="SOCK", outfile="")
    registerDoParallel(cl)

    output.dir <- file.path(working.dir, "mrresults")
    dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
    study.output.dir <- file.path(output.dir, study.name)
    dir.create(study.output.dir, showWarnings=FALSE, recursive=TRUE)

    if (!is.null(gwas.dt)) {
        gwas.hits.genes <- unique(unlist(strsplit(gwas.dt[, `REPORTED GENE(S)`],
                                                split=", ")))
        gwas.hits.genes <- trimws(gwas.hits.genes, whitespace="[\\h\\v]")
    } else {
        gwas.hits.genes <- NULL
    }

    pheno.name <- phenoname
    mr.exposures <- core.genes

    scoresinfo <- scoresinfo[!is.na(scoreid_trans) & locus.diversity > 5]
    scoresinfo[, pvalue := pvalue_trans]

    scoresbychr.dir <- switch(study.name,
                              eQTLGen = file.path(basegatedir, "eqtl"),
                              deCODE = file.path(basegatedir, "pqtl"),
                              UKBPPP = file.path(basegatedir, "pqtl"))
    scoresall.dir <- switch(study.name,
                            eQTLGen = file.path(basegatedir, "eqtl"),
                            deCODE = file.path(basegatedir, "pqtl"),
                            UKBPPP = file.path(basegatedir, "pqtl"))
    this.studyid <- switch(study.name,
                           eQTLGen = 36,
                           deCODE = 1008,
                           UKBPPP = 50)
    scores.file <- file.path(scoresbychr.dir, "all.genotypicscore")
    meta.file <- file.path(scoresbychr.dir, "all.scoresinfo.annotated.Rdata.gz")
    kg.file <- sprintf("%s.stats", ld.refplinkfile)

    ## Fixed parameters for analysis
    ## -----------------------------------------------------------------------------
    ## exclude locus-specific scores overlapping with the HLA region from the MR
    ## analysis and the MR plots
    exclude.hla <- TRUE

    ## scale each instrument to have zero mean and unit variance (recommended).
    scale.instruments <- TRUE

    ## Set MR-Hevo parameters to control MCMC sampling and inference
    ## -----------------------------------------------------------------------------
    ## initial guess on the proportion of pleiotropic instruments
    fraction.pleiotropic <- 0.5
    ## scale of the slab portion of spike-and-lab prior on the proportion of
    ## pleiotropic instruments
    slab.scale <- 0.05
    ## slab degrees of freedom
    slab.df <- 4
    ## number of MCMC iterations
    num.iter.mcmc <- 8000
    ## number of iterations for sampler warmup
    num.warmup.mcmc <- 4000

    ## Perform MR for all putative core genes (eQTL/pQTL detected) in sequence,
    ## using summary statistics from any of the three large QTL studies (if the
    ## gene product is quantified and has at least 5 locus-specific trans scores)
    ## --------------------------------------------------------------------------

    ## initialise list to store the results from MR
    if (!exists("mrres")) {
        ## only initialise if it does not exist
        ## this allows to restart the job from an existing session
        mrres <- vector("list", length(mr.exposures))
        names(mrres) <- mr.exposures
    }

    ## start MR

    ## specify threshold for filtering out weak instruments.
    ## Recommended thresholds:
    ## -- 1E-6 when the genetic instruments are based on the pQTL studies
    ##    (tested all SNPs genome-wide)
    ## -- 1E-5 when the genetic instruments are based on the eQTL study
    ##    (pre-selected trait-associated SNPs)
    instrument.minpvalue <- switch(study.name,
                                    eQTLGen = 1E-5,
                                    deCODE = 1E-6,
                                    UKBPPP = 1E-6)

    ## loop over exposures
    newrun <- TRUE
    results <- foreach(exposure.name=mr.exposures,
            .combine='c',
            .packages=(.packages()),
            .errorhandling='pass') %dopar% {
        source("mr.R")
        source("helper-functions/report.helpers.R")
        source("helper-functions/aggregation.helpers.R")
        source("helper-functions/shared.helpers.R")

        tryCatch({
            if (exposure.name %in% scoresinfo[studyid==this.studyid, gene_symbol]) {

                ## get gwasid for the exposure
                exposure.scoresinfo <- scoresinfo[gene_symbol == exposure.name &
                                                  studyid==this.studyid, ]
                exposure.gwasid <- unique(exposure.scoresinfo[pvalue==min(pvalue),
                                                              gwasid])

                ## check if this MR analysis has already been run
                output.file <- file.path(study.output.dir, paste0(
                    "mrplot.", exposure.name, ".", exposure.gwasid, ".png"))

                re <- sprintf(".*%s.*.csv", exposure.name)
                result_file <- list.files(path=study.output.dir, pattern=re,
                                          full.names=TRUE)[1]

                ## Try to read existing result first
                if (!is.na(result_file) & file.exists(result_file)) {
                    cat("Reading existing MR results for", exposure.name, "\n")
                    result <- fread(result_file)

                } else if (!newrun & file.exists(output.file)) {
                    cat("MR plot exists but no result file. Skipping", exposure.name, "\n")
                    return(NULL)

                } else {
                    cat("Running MR analysis for", exposure.name, "\n")
                    result <- run.mrhevo(scoresbychr.dir, scoresall.dir,
                                        scores.file, meta.file, kg.file,
                                        exposure.name, exposure.gwasid,
                                        pheno.dt, pheno.name, covariates,
                                        study.output.dir, gwas.hits.genes,
                                        instrument.minpvalue, exclude.hla,
                                        scale.instruments, fraction.pleiotropic,
                                        slab.scale, slab.df, num.iter.mcmc,
                                        num.warmup.mcmc, cores)
                    gc()
                    if (!is.data.table(result)) {
                        cat("Results were not generated. Check other errors.\n")
                        result <- data.table()
                    }
                }
                # Add metadata columns
                result[, exposure := exposure.name]
                result[, study := study.name]
                return(list(result))
            } else {
                return(NULL)
            }

        }, error = function(e) {
            cat("Error processing", exposure.name, ":",
                conditionMessage(e), "\n")
            return(NULL)
        })
    }
    stopCluster(cl)
    return(mrres)
}

count.variants <- function(basedir, genotype.file.prefix,
                           genotype.file.suffix) {
    snps.total <- data.table(chr=paste0("chr", 1:22), n_variants=NA_integer_)
    for (id in 1:22) {
        bim <- file.path(basename, id, paste0(genotype.file.prefix, id,
                                              genotype.file.suffix, ".bim"))
        nvar <- as.integer(system2("wc", args=c("-l", bim), stdout=TRUE) |>
                           sub(pattern=" .*", replacement=""))
        snps.total[id, n_variants := nvar]
    }
    return(snps.total[, sum(n_variants)])
}
