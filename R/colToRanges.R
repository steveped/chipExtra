#' @title Coerce a column to a GRanges object
#'
#' @description Coerce a column to a GRanges object from a rectangular object
#'
#' @details
#' Take a data.frame-like object and coerce one column to a GRanges object,
#' setting the remainder as the `mcols`.
#' A particularly useful application of this is when you have a GRanges object
#' with one mcol being a secondary GRanges object.
#'
#' Alternatively, if you have a data.frame with GRanges represented as a
#' character column, this provides a simple method of coercion.
#' In this case, no Seqinfo element will be applied to the GRanges element.
#'
#' @return
#' A GenomicRanges object
#'
#' @examples
#' set.seed(73)
#' x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
#' seqinfo(x) <- Seqinfo("chr1", 60, FALSE, "Example")
#' df <- data.frame(logFC = rnorm(3), logCPM = rnorm(3,8), p = 10^-rexp(3))
#' mcols(x) <- df
#' gr <- mergeByCol(x, col = "logCPM", pval = "p")
#' colToRanges(gr, "keyval_range")
#'
#' @param x A data-frame or GRanges object containing the column to coerce
#' @param var The name of the column to coerce
#' @param seqinfo A seqinfo object to be applied to the new GRanges object.
#' Ignored if the column is already a GRanges object
#' @param ... Used to pass arguments to lower-level functions
#'
#' @importFrom methods as
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo<- seqlevels
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
    "colToRanges", signature = signature(x = "DataFrame"),
    function(x, var, seqinfo = NULL, ...) {
        stopifnot(var %in% colnames(x))
        gr <- GRanges(x[[var]])
        if (!is.null(seqinfo) & !is(x[[var]], "GRanges")) {
            map <- match(seqlevels(seqinfo), seqlevels(gr))
            seqinfo(gr, new2old = map) <- seqinfo
        }
        keep <- setdiff(colnames(x), var)
        mcols(gr) <- x[keep]
        gr
    }
)
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
    "colToRanges", signature = signature(x = "GRanges"),
    function(x, var, ...) {
        df <- mcols(x)
        colToRanges(df, var, ...)
    }
)
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo<- seqlevels
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom IRanges CompressedList
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
    "colToRanges", signature = signature(x = "data.frame"),
    function(x, var, seqinfo = NULL, ...) {
        stopifnot(var %in% colnames(x))
        gr <- GRanges(x[[var]])
        if (!is.null(seqinfo) & !is(x[[var]], "GRanges")) {
            map <- match(seqlevels(seqinfo), seqlevels(gr))
            seqinfo(gr, new2old = map) <- seqinfo
        }
        keep <- setdiff(colnames(x), var)
        DF <- x[keep]
        list_cols <- vapply(DF, is.list, logical(1))
        if (any(list_cols)) {
            DF <- as.list(DF)
            DF[list_cols] <- lapply(DF[list_cols], as, Class = "CompressedList")
            DF <- DataFrame(DF)
        }
        mcols(gr) <- DF
        gr
    }
)
