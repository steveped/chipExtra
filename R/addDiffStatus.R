#' @title Add a status column
#' @description
#' Add a status column based on significance and estimated change
#'
#' @details
#' This takes a simple object and adds a new column classifying entries into
#' one of three categories, as specified using `up`, `down` or `other`.
#' Results in the new column will always be returned as a factor with levels in
#' order of the values provided in the arguments `other`, `down` and `up`
#'
#' @param x Object to be classified
#' @param fc_col Name of the fold-change column
#' @param sig_col Name of the column with significance values
#' @param alpha significance threshold
#' @param cutoff minimum estimated change to be considered in either of the up
#' or down categories
#' @param up,down,other factor levels to annotate regions based on the above
#' criteria
#' @param missing Value to add when either fc_col or sig_col has NA values
#' @param new_col name of the new column to be added
#' @param drop logical(1) Drop unused factor levels from the status column
#' @param ... Used to pass arguments between methods
#'
#' @return An object of the same type as provided
#'
#' @examples
#' ## Working with a data.frame
#' set.seed(101)
#' df <- data.frame(logFC = rnorm(20), p = rbeta(20, shape1 = 1, shape2 = 20))
#' df$FDR <- p.adjust(df$p, "fdr")
#' addDiffStatus(df)
#'
#' ## This works identically with a GRanges object, amongst others
#' gr <- GRanges(paste0("chr1:", seq_len(20)))
#' mcols(gr) <- df
#' addDiffStatus(gr)
#'
#' @name addDiffStatus
#' @rdname addDiffStatus-methods
#' @export
#'
setGeneric("addDiffStatus", function(x, ...) standardGeneric("addDiffStatus"))
#' @importFrom dplyr case_when
#' @rdname addDiffStatus-methods
#' @export
setMethod(
  "addDiffStatus",
  signature = signature(x = "data.frame"),
  function(
    x, fc_col = "logFC", sig_col = c("FDR", "hmp_fdr", "p_fdr", "adj.P.Value"),
    alpha = 0.05, cutoff = 0, up = "Increased", down = "Decreased",
    other = "Unchanged", missing = "Undetected", new_col = "status",
    drop = FALSE, ...
  ) {

    # Start with a df
    stopifnot(is(x, "data.frame"))
    nm <- colnames(x)
    fc_col <- match.arg(fc_col, nm)
    sig_col <- intersect(sig_col, nm)[[1]]
    stopifnot(length(sig_col) == 1)
    fc <- x[[fc_col]]
    stopifnot(is.numeric(fc))
    stopifnot(is.numeric(x[[sig_col]]))
    sig <- x[[sig_col]] < alpha
    status <- case_when(
      is.na(sig) | is.na(fc) ~ missing[[1]],
      !sig ~ other[[1]],
      fc > abs(cutoff) ~ up[[1]],
      fc < -abs(cutoff) ~ down[[1]]
    )
    ## Do we need to add an explicit NA value here?
    lv <- unique(c(other, down, up, missing))
    if (drop) lv <- intersect(lv, status)
    x[[new_col]] <- factor(status, levels = lv)
    x
  }
)
#' @rdname addDiffStatus-methods
#' @export
setMethod(
  "addDiffStatus", signature = signature(x = "DataFrame"),
  function(x, new_col = "status", ...) {
    df <- as.data.frame(x)
    orig_names <- colnames(df)
    df <- addDiffStatus(df, new_col = new_col, ...)
    x[[new_col]] <- df[[new_col]]
    x
  }
)
#' @rdname addDiffStatus-methods
#' @export
setMethod(
  "addDiffStatus", signature = signature(x = "GRanges"), function(x, ...) {
    df <- mcols(x)
    df <- addDiffStatus(df, ...)
    mcols(x) <- df
    x
  }
)
#' @rdname addDiffStatus-methods
#' @export
setMethod(
  "addDiffStatus", signature = signature(x = "GRangesList"), function(x, ...) {
    endoapply(x, addDiffStatus, ...)
  }
)
#' @rdname addDiffStatus-methods
#' @export
setMethod(
  "addDiffStatus", signature = signature(x = "SummarizedExperiment"),
  function(x, ...) {
    df <- rowData(x)
    df <- addDiffStatus(df, ...)
    rowData(x) <- df
    x
  }
)
