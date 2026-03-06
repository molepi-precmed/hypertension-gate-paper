#' This script contains shared helper functions that are used by multiple helper
#' and main scripts.

#! convert chr bp to absolute position x, using hg38 to set length of chr
position.absolute <- function(CHR, BP, padding=0) {
    positions <- data.table(CHR, BP)
    positions[, rownum := .I]
    setkey(positions, CHR)
    setkey(chr.lengths, CHR)
    positions <- chr.lengths[, .(CHR, cumbp)][positions]
    positions[, x := BP + cumbp + padding]
    setorder(positions, rownum)
    return(positions[, x])
}

pnorm.extreme <- function(z, upper=TRUE) {
    ## https://www.johndcook.com/blog/norm-dist-bounds/
    if(upper) { # upper bound on log10 tail prob
        c <- 8 / pi
    } else { # lower bound
        c = 4
    }
    x <-  2 * sqrt(2 / pi) / (z + sqrt(z^2 + c))
    ln_p = -0.5 * z^2 + log(x)
    log10p <- ln_p / log(10)
    exponent <- floor(log10p)
    coeff <- 10^(log10p - exponent)
    string <- paste0(round(coeff), "E", exponent)
    return(string)
}

#' reformat for LaTeX p-values that are in scientific notation
format.scinot.pvalue <- function(x, nexp=1) {
    x <- toupper(x)
    x.split <- as.numeric(unlist(strsplit(as.character(x), "E"))) # split x.split at E
    x.split <- signif(as.numeric(x.split, 1))
    x.split <- t(matrix(x.split, nrow=2))
    roundedto10 <- x.split[, 1]==10
    x.split[, 1][roundedto10] <- 1
    x.split[, 2][roundedto10] <- x.split[, 2][roundedto10] - 1
    p.latex <- sprintf("\\ensuremath{%.*f \\times 10^{%0*d}}", 0, x.split[, 1], nexp, round(x.split[, 2]))
    return(p.latex)
}

#' generate formatted p-values from z value
#' global variables: sigfig, neglogp.threshold.scinot, neglogp.threshold
format.z.aspvalue <- function(z) {
    if(!exists("sigfig")) {
        sigfig <- 1
    } else {
        if(sigfig==0) {
            sigfig <- 1
        }
    }


    ## neglogp.threshold.scinot is threshold for using scientific notation
    if(!exists("neglogp.threshold.scinot")) {
        neglogp.threshold.scinot <- 3
    } else {
        if(neglogp.threshold.scinot==0) {
            neglogp.threshold.scinot <- 3
        }
    }

    p <- signif(2 * pnorm(-abs(z)), sigfig)
    p.char <- toupper(as.character(p))
    ## pnorm.extreme returns a character string of form "NE-NNN"
    p.char[!is.na(p.char) & p.char=="0"] <- pnorm.extreme(z[!is.na(p.char) & p.char=="0"])    # where R outputs 0
    sci.revert <- grepl("E", p.char) & p > 10^-neglogp.threshold.scinot
    p.char[sci.revert] <-  format(p[sci.revert], scientific=FALSE)

    if(exists("neglogp.threshold")) {
        if(neglogp.threshold > 0) { # thresholding of p values
            p.char[p < 10^-neglogp.threshold] <- paste0("<",
                                                        format(10^-neglogp.threshold,
                                                               scientific=FALSE))
        }
    } else {
        }
    p.char[grep("E", p.char)] <- format.scinot.pvalue(p.char[grep("E", p.char)])
    return(p.char)
}

format.pvalue <- function(z, pvalue=NULL) {
    format.z.aspvalue(z)
}

## format a vector of pvalues in LaTeX and return a vector of mode character
pvalue.latex <- function(x, n=1, nexp=1, s.threshold=s.threshold) {
    ## this function has to be able to handle x whether numeric or character
    pvalue <- sapply(x, function(z) { # sapply returns a vector applying FUN to each element of x

        if (is.na(z) | is.nan(z)) {
            return(NA)
        } else if(as.numeric(z) >= 10^s.threshold) {
            ## return character string to one sig fig, not in scientific notation
            return(as.character(signif(as.numeric(z), 1)))
        } else {
            if(is.numeric(z)) {
                ## rounds to 1 sig fig and convert to character string
                z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
            } else {
                z <- toupper(z)
            }
            z <- as.numeric(unlist(strsplit(as.character(z), "E"))) # split z at E
            sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
        }
    }
    )
    pvalue <- as.character(pvalue)
    if (s.threshold==-4) {
        pvalue[grep("\\\\times", pvalue)] <- "<0.0001" # fix for thresholding at 0.0001
    }
    return(pvalue)
}

dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriation::seriate(d_corr, method = method, ...)
  i = seriation::get_order(s)
  return(i)
}

chr.lengths <- fread(text = "CHR,length
1, 248956422
2, 242193529
3,198295559
4,190214555
5,181538259
6,170805979
7,159345973
8,145138636
9,138394717
10,133797422
11,135086622
12,133275309
13,114364328
14,107043718
15,101991189
16,90338345
17,83257441
18,80373285
19,58617616
20,64444167
21,46709983
22,50818468
X,156040895
Y,57227415")
chr.lengths[, CHR := factor(CHR, levels = c(1:22, "X", "Y"))]
chr.lengths[, cumbp := c(0, cumsum(as.numeric(length)[-.N]))]
chr.lengths[, midpoint.chr := 0.5 * length]
chrom.labels <- as.character(1:22)
chrom.labels[c(19, 21)] <- " "
chrom.midpoints <- position.absolute(chr.lengths$CHR, chr.lengths$midpoint.chr)[1:22]
chrom.breaks <- chr.lengths[1:22, cumbp]
chrom.ends <- position.absolute(chr.lengths$CHR, chr.lengths$length)

unique.csv <- function(x) {
    unique.elementwise <- function(x) {
        x.vector <- unlist(strsplit(x, split="[ ,]"))
        x.vector <- x.vector[nchar(x.vector) > 0]
        x.vector <- unique(x.vector)
        x.string <- paste(x.vector, collapse=", ")
        return(x.string)
    }
    return(sapply(x, unique.elementwise))
}

#' Helper function to load a given Rdata object.
load.rdata <- function(filename) {
    load(filename)
    get(ls()[ls() != "filename"])
}

load.matrix <- function(matrix.file) {
    if (!file.exists(matrix.file))
        stop("Specified file:", matrix.file, " does not exist.")
    m <- load.rdata(matrix.file)
    return(m)
}

check.cols <- function(required.cols, object) {
    if (!all(required.cols %in% names(object))) {
        cols.missing <- required.cols[!required.cols %in% names(object)]
        stop("Required columns: ", paste0(cols.missing, collapse=", "),
             " are missing from required columns: ",
             paste(paste0(required.cols, collapse=", ")), ".")
    }
}

#' Check if a data.table contains columns with NA values or columns where the
#' number of unique valus is less than 2. If a column has NA values, print a
#' warning and continue. If a column has only one value for all rows, print
#' error message and stop (used to check score dt before fitting GLM for
#' association calculation)
#'
#' @param dt data.table object with values to check.
check.cols.values <- function(dt) {
    ## check for NA values
    na.check <- lapply(dt, function(col) any(is.na(col)))
    if (any(unlist(na.check))) {
        warning(
            "Warning: NA values found in column(s): ",
            paste(names(dt)[unlist(na.check)], collapse = ", ")
        )
    }

    ## check for number of unique values
    single.val.check <- lapply(dt, function(col) length(unique(col)) < 2)
    if (any(unlist(single.val.check))) {
        stop(
            "Error: Column(s) with less than 2 unique values found: ",
            paste(names(dt)[unlist(single.val.check)], collapse = ", ")
        )
    }
}

#' Report a message to the terminal.
#'
#' @param mode One of \sQuote{info} (normal messages), \sQuote{note} (messages
#'        that require some highlighting), \sQuote{warn} (important information
#'        the user should definitely notice).
#' @param ... Strings to be reported.
#' @param LF Whether a newline character should be added at the end of the
#'        message (\code{TRUE} by default).
#'
#' @import crayon
#' @export
msg <- function(mode, ..., LF=TRUE) {
    message(mode(...), appendLF=LF)
}
info <- crayon::reset
note <- crayon::green
warn <- crayon::yellow
bold <- crayon::bold
qqplot <- function(z) {
  qq <- data.table(z= sort(z))
  qq[, p := ppoints(.N)]
  qq[, z.exp := qnorm(p)]

  ticks <- c(0, 2, 5, 10, 20, 50, 100)
  breaks <- c(qnorm(0.5 * 10^(-rev(ticks))), -qnorm(0.5 * 10^(-ticks)))
  p.qq <- ggplot(data=qq, aes(z.exp, z)) +
    geom_point() +
    geom_abline(slope=(quantile(qq$z, 0.75) - quantile(qq$z, 0.25)) /
                  (quantile(qq$z.exp, 0.75) - quantile(qq$z.exp, 0.25)),
                intercept=0) +
    geom_abline(slope=1, intercept=0, linetype="dotted", linewidth=0.5) +
    scale_y_continuous(breaks=breaks, labels=c(rev(ticks), ticks)) +
    xlab(paste0("Quantile of test statistic under the null, based on ", length(z), " tests")) +
    ylab("Minus log10 p-value")
  return(p.qq)
}

save.qqplot <- function(study, coeffs.wide, working.dir, diversity=5) {
    z_trans <- coeffs.wide[!is.na(z_trans) & locus.diversity>diversity, z_trans]
    p.qq <- qqplot(z_trans)
    plot.file <- file.path(working.dir, study, "qqplot.png")
    ggsave(filename=plot.file,
           p.qq, device="png",
           dpi=300, width=7, height=7, units="in")
}
