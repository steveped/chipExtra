#' Re-estimate peak centres from coverage
#'
#' Use coverage to estimate peak centres
#'
#' @details
#' Use coverage to estimate the centre of a set of peaks or GenomicRanges.
#' The point of maximum coverage for each sample will be found and these
#' positions will be averaged to return a position representing an estimated
#' peak centre.
#' If using weighted.mean, positions are weighted by the maximum coverage
#' within each sample.
#'
#' @return A GRanges object with all widths set to one
#'
#' @param x A set of GRanges representing peaks
#' @param y A BamFileList
#' @param f The function to use when estimating a combined peak centre
#' @param BPPARAM An object of class BPPARAM
#' @param ... Not used
#'
#' @name centrePeaks
#' @rdname centrePeaks-methods
#' @export
#'
setGeneric("centrePeaks", function(x, y, ...) standardGeneric("centrePeaks"))
#' @importClassesFrom IRanges RleViewsList RleList
#' @importFrom GenomicAlignments coverage
#' @importFrom Rsamtools ScanBamParam
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges RleList viewWhichMaxs viewMaxs RleViewsList Views
#' @importFrom stats weighted.mean
#' @importFrom matrixStats rowMedians rowMeans2
#' @rdname centrePeaks-methods
#' @export
setMethod(
    "centrePeaks",
    signature = signature(x = "GRanges", y = "BamFileList"),
    function(
        x, y, f = c("mean", "median", "weighted.mean"), BPPARAM = bpparam(), ...
    ) {
        ## Set the function to use when calculating positions
        f <- match.arg(f)

        ## Check BiocParallel is ready to go
        if (!bpisup(BPPARAM)) {
            bpstart(BPPARAM)
            on.exit(bpstop(BPPARAM))
        }

        ## Get the coverage
        sbp <- ScanBamParam(which = x)
        cov <- bplapply(y, coverage, param = sbp, BPPARAM = BPPARAM)

        ## Split everything by chromosome
        sq <- seqinfo(x)
        grl <- splitAsList(x, seqnames(x))
        grl <- grl[vapply(grl, length, integer(1)) > 0]
        new <- bplapply(
            grl,
            \(x) {
                chr <- as.character(unique(seqnames(x)))
                if (!length(chr)) return(GRanges())
                rl <- RleList(lapply(cov, \(i) i[[chr]]))
                rvl <- RleViewsList(lapply(rl, Views, start(x), end(x)))
                pos <- viewWhichMaxs(rvl)
                posm <- do.call("cbind", pos)
                if (f == "weighted.mean") {
                    maxs <- viewMaxs(rvl)
                    wm <- do.call("cbind", maxs)
                    wm <- wm / max(wm) # Avoids integer overflows
                    i <- seq_len(nrow(wm))
                    centres <- vapply(
                        i, \(i) weighted.mean(posm[i,], wm[i,]), numeric(1)
                    )
                }
                if (f == "mean") centres <- rowMeans2(posm)
                if (f == "median") centres <- rowMedians(posm)
                GRanges(paste0(chr, ":", as.integer(centres)), seqinfo = sq)
            }, BPPARAM = BPPARAM
        )
        new <- unlist(GRangesList(new))
        mcols(new) <- mcols(x)
        new
    }
)


