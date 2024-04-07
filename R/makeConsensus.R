#' @title Make a set of consensus peaks
#'
#' @description Make a set of consensus peaks based on the number of replicates
#'
#' @details
#' This takes a list of GRanges objects and forms a set of consensus peaks.
#'
#' When using method = "union" the union ranges of all overlapping peaks will
#' be returned, using the minimum proportion of replicates specified.
#' When using method = "coverage", only the regions within each overlapping
#' range which are 'covered' by the minimum proportion of replicates specified
#' are returned. This will return narrower peaks in general, although some
#' artefactual very small ranges may be included (e.g. 10bp). Careful setting
#' of the min_width and merge_within parameters may be very helpful for these
#' instances. It is also expected that setting method = "coverage" should return
#' the region within each range which is more likely to contain the true binding
#'  site for the relevant ChIP targets
#'
#'
#' @param x A GRangesList
#' @param p The minimum proportion of samples (i.e. elements of `x`) required
#' for a peak to be retained in the output. By default all merged peaks will
#' be returned
#' @param var Additional columns in the `mcols` element to retain
#' @param method Either return the union of all overlapping ranges, or the
#' regions within the overlapping ranges which are covered by the specified
#' proportion of replicates. When using p = 0, both methods will return identical
#' results
#' @param ignore.strand,simplify,... Passed to \link{reduceMC} or
#' \link{intersectMC} internally
#' @param min_width Discard any regions below this width
#' @param merge_within Passed to \link[GenomicRanges]{reduce} as `min.gapwidth`
#'
#' @return
#' `GRanges` object with mcols containing a logical vector for every element of
#' x, along with the column `n` which adds all logical columns. These columns
#' denote which replicates contain an overlapping peak for each range
#'
#' If any additional columns have been requested using `var`, these will be
#' returned as CompressedList objects as produced by `reduceMC()` or
#' `intersectMC()`.
#'
#' @seealso \link{reduceMC} \link{intersectMC}
#'
#' @examples
#' data("peaks")
#' ## The first three replicates are from the same treatment group
#' grl <- peaks[1:3]
#' names(grl) <- gsub("_peaks.+", "", names(grl))
#' makeConsensus(grl)
#' makeConsensus(grl, p = 2/3, var = "score")
#'
#' ## Using method = 'coverage' finds ranges based on the intersection
#' makeConsensus(grl, p = 2/3, var = "score", method = "coverage")
#'
#'
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps overlapsAny
#' @importFrom S4Vectors mcols<- subset mcols endoapply
#' @importFrom methods is
#' @export
makeConsensus <- function(
        x, p = 0, var = NULL, method = c("union", "coverage"),
        ignore.strand = TRUE, simplify = FALSE, min_width = 0,
        merge_within = 1L, ...
) {

    ## Starting with a GRList
    if (!is(x, "GRangesList")) stop("Input must be a GRangesList")
    nm <- names(x)
    if (length(nm) != length(x)) stop("Each element of 'x' must be named")
    method <- match.arg(method)
    if (!is.null(var)) {
        mc_names <- .mcolnames(x[[1]])
        if (length(setdiff(var, mc_names))) {
            d <- paste(setdiff(var, mc_names), collapse = ", ")
            stop("Couldn't find column ", d)
        }
        x <- endoapply(
            x,
            function(gr) {
                mcols(gr) <- mcols(gr)[var]
                gr
            }
        )
    } else {
        x <- endoapply(x, granges)
    }
    unlisted <- unlist(x)

    if (method == "coverage") {
        ## Find how many elements each range is covered by
        cov_gr <- GRanges(coverage(x))
        cov_gr <- subsetByOverlaps(cov_gr, unlisted)
        ## Coverage will be returned as the column 'score'
        keep <- mcols(cov_gr)[["score"]] >= p * length(x)
        cov_gr <- cov_gr[keep]
        ## If merging within a certain range, merge here. This ensures that
        ## returned regions still have coverage but deals with small holes
        ## that appear due to drops in coverage, without messing up the *MC
        cov_gr <- GenomicRanges::reduce(
            cov_gr, min.gapwidth = merge_within, ignore.strand = ignore.strand
        )
        red_ranges <- intersectMC(
            unlisted, cov_gr, ignore.strand = ignore.strand,
            simplify = simplify, ...
        )

    }
    if (method == "union") {
        ## For now, remove all mcols, however reduceMC may be useful if
        ## wishing to retain these for use with plotOverlaps()
        red_ranges <- reduceMC(
            unlisted, ignore.strand = ignore.strand, simplify = simplify,
            min.gapwidth = merge_within, ...
        )
    }
    ## Return the columns for overlaps
    ol <- lapply(x, function(x) overlapsAny(red_ranges, x))
    ol <- as.data.frame(ol)

    ## Ensure names are strictly retained
    names(ol) <- names(x)
    ol$n <- rowSums(ol)
    mcols(red_ranges) <- cbind(mcols(red_ranges), ol)
    subset(red_ranges, n >= p * length(x) & width >= min_width)

}
