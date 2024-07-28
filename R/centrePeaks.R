#' Re-estimate peak centres from coverage
#'
#' Use coverage to estimate peak centres
#'
#' @details
#' Use coverage to estimate the centre of a set of peaks or GenomicRanges.
#'
#' If using the mean or median, the point of maximum coverage for each sample
#' will be found within each peak and these positions will be averaged to return
#' a position representing an estimated peak centre.
#'
#' If using weighted.cov, positions are weighted by the combined coverage across
#' all samples to return the weighted mean position. In this case coverage will
#' be scaled by total alignments within each bam file before summing across files
#'
#' @return A GRanges object with all widths set to one
#'
#' @param x A set of GRanges representing peaks
#' @param y A suitable set of files with methods defined
#' @param f The function to use when estimating a combined peak centre
#' @param BPPARAM An object of class BPPARAM
#' @param ... Used to pass arguments between methods
#'
#' @examples
#' ## Define some peaks
#' f <- system.file("extdata/peaks.bed.gz", package = "extraChIPs")
#' peaks <- importPeaks(f, type = "bed")[[1]]
#' peaks
#'
#' ## Use a bam file to re-centre the regions using highest coverage
#' bf <- system.file("extdata/bam/ex1.bam", package = "extraChIPs")
#' centres <- centrePeaks(peaks, bf, BPPARAM = SerialParam())
#' centres
#'
#' @name centrePeaks
#' @rdname centrePeaks-methods
#' @export
#'
setGeneric("centrePeaks", function(x, y, ...) standardGeneric("centrePeaks"))
#' @importFrom GenomicAlignments coverage
#' @importFrom Rsamtools ScanBamParam idxstatsBam
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "BamFileList"),
    function(
        x, y, f = c("weighted.cov", "mean", "median"), BPPARAM = bpparam(), ...
    ) {
        ## Set the function to use when calculating positions
        f <- match.arg(f)

        ## Check BiocParallel is ready to go
        if (!bpisup(BPPARAM)) {
            bpstart(BPPARAM)
            on.exit(bpstop(BPPARAM))
        }

        ## Get the coverage & scale by total reads
        message(
            sprintf("Getting coverage across all %i peaks...", length(x)),
            appendLF = FALSE
        )
        t1 <- Sys.time()
        sbp <- ScanBamParam(which = x)
        cov <- bplapply(y, coverage, param = sbp, BPPARAM = BPPARAM)
        tot_aln <- vapply(y, \(i) sum(idxstatsBam(i)$mapped), numeric(1))
        t2 <- Sys.time()
        message(sprintf("done (%.2fs)", t2 - t1))

        .centreFromCov(cov, x, f, tot_aln, BPPARAM)

    }
)
#' @importFrom Rsamtools BamFileList
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "BamFile"),
    function(x, y, ...) {
        bfl <- BamFileList(y)
        names(bfl) <- basename(path(y))
        centrePeaks(x, bfl, ...)
    }
)
#' @importClassesFrom IRanges RleList
#' @importClassesFrom rtracklayer BigWigFileList
#' @importFrom rtracklayer import.bw
#' @importFrom S4Vectors splitAsList Rle
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "BigWigFileList"),
    function(
        x, y, f = c("weighted.cov", "mean", "median"), BPPARAM = bpparam(), ...
    ) {
        ## Set the function to use when calculating positions
        f <- match.arg(f)

        ## Check BiocParallel is ready to go
        if (!bpisup(BPPARAM)) {
            bpstart(BPPARAM)
            on.exit(bpstop(BPPARAM))
        }

        ## Get the coverage without any scaling
        message(
            sprintf("Getting coverage across all %i peaks...", length(x)),
            appendLF = FALSE
        )
        t1 <- Sys.time()
        grl <- bplapply(y, import.bw, which = x, BPPARAM = BPPARAM)
        ## Add zero scores to all missing positions
        grl_exp <- lapply(
            grl,
            \(gr) {
                gnm <- setdiff(GRanges(seqinfo(gr)), gr)
                gnm$score <- 0
                sort(c(gr, gnm))
            }
        )
        cov <- lapply(
            grl_exp,
            \(gr) {
                grl <- splitAsList(gr, seqnames(gr))
                rl <- lapply(grl, \(x) Rle(values = x$score, lengths = width(x)))
                RleList(rl, compress = FALSE)
            }
        )
        t2 <- Sys.time()
        message(sprintf("done (%.2fs)", t2 - t1))

        .centreFromCov(cov, x, f, 1, BPPARAM)

    }
)
#' @importClassesFrom rtracklayer BigWigFile
#' @importFrom rtracklayer BigWigFileList
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "BigWigFile"),
    function(x, y, ...) {
        bwfl <- BigWigFileList(y@resource)
        centrePeaks(x, bwfl, ...)
    }
)
#' @importClassesFrom rtracklayer BigWigFile
#' @importClassesFrom Rsamtools BamFileList
#' @importFrom rtracklayer BigWigFileList
#' @importFrom Rsamtools BamFileList
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "character"),
    function(x, y, ...) {
        suff <- tolower(unique(gsub(".+\\.([A-Za-z]+$)", "\\1", y)))
        if (length(suff) > 1) stop("paths must be to files of the same type")
        f <- NULL
        if (suff == "bam") f <- BamFileList(y)
        if (suff %in% c("bw", "bigwig")) f <- BigWigFileList(y)
        if (is.null(f)) stop("Couldn't detect file type from suffix")
        centrePeaks(x, f, ...)
    }
)

#' @importClassesFrom IRanges RleViewsList RleList
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges RleList viewWhichMaxs viewMaxs RleViewsList Views
#' @importFrom matrixStats rowMedians rowMeans2 weightedMean
#' @keywords internal
.centreFromCov <- function(cov, peaks, f, totals = 1, BPPARAM) {
    ## Split everything by chromosome
    sq <- seqinfo(peaks)
    grl <- splitAsList(peaks, seqnames(peaks))
    grl <- grl[vapply(grl, length, integer(1)) > 0]
    n <- length(cov)
    totals <- rep_len(totals, n)

    message("Finding centres by chromosome...", appendLF = FALSE)
    t1 <- Sys.time()
    new <- bplapply(
        grl,
        \(gr) {
            ## Just work with the coverage from each chromosome
            chr <- as.character(unique(seqnames(gr)))
            if (!length(chr)) return(GRanges(seqinfo = sq))
            rl <- RleList(lapply(cov, \(i) i[[chr]]))
            if (f == "weighted.cov") {
                scaled <- lapply(seq_along(rl), \(i) rl[[i]] / totals[[i]])
                rvl <- RleViewsList(lapply(scaled, Views, start(gr), end(gr)))
                ## Calculate combined weights by adding scaled coverage
                ## across positions
                w <- lapply(
                    seq_along(gr),
                    \(j) rowSums(
                        do.call(
                            "cbind", lapply(rvl, \(i) as.numeric(i[[j]]))
                        )
                    )
                )
                ## Define positions within each peak
                pos <- mapply(
                    \(i, j) seq(i, length.out = j),
                    i = start(gr), j = width(gr),
                    SIMPLIFY = FALSE
                )
                ## Return centres as the mean position weighted by coverage
                centres <- mapply(weightedMean, x = pos, w = w)
            } else{
                rvl <- RleViewsList(lapply(rl, Views, start(gr), end(gr)))
                pos <- viewWhichMaxs(rvl)
                posm <- do.call("cbind", as.list(pos))
                if (f == "mean") centres <- rowMeans2(posm)
                if (f == "median") centres <- rowMedians(posm)
            }
            GRanges(paste0(chr, ":", as.integer(centres)), seqinfo = sq)
        }, BPPARAM = BPPARAM
    )
    t2 <- Sys.time()
    message(sprintf("done (%.2fs)", t2 - t1))
    new <- unlist(GRangesList(new))
    mcols(new) <- mcols(peaks)
    new
}


