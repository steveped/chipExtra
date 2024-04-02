#' @title Merge Sliding Windows using the Harmonic Mean P
#'
#' @description Merge overlapping windows using harmonic mean p-values from significance
#' testing
#'
#' @details
#' When using sliding windows to test for differential signal, overlapping
#' windows can be merged based on the significance of results.
#' `mergeByHMP()` merges overlapping windows using the asymptotically exact
#' harmonic mean p-value \link[harmonicmeanp]{p.hmp} from the individual,
#' window-level tests. This tests the Null Hypothesis that there is no
#' significance amongst the initial set of p-values, and returns a summarised
#' value which controls the FDR within a set of tests (Wilson, PNAS, 2019).
#' Multilevel testing across the set of results is currently implemented using
#' `p_adj_method = "fwer"`
#'
#' Given that the harmonic mean p-value is calculated from the inverse p-values,
#' these are used to provide a *weighted average* of expression and logFC values
#' in the returned object. Any weights provided in `w` are ignored for these
#' values as they are simple representative estimates.
#' The representative range returned in `keyval_range` corresponds to the window
#' with the lowest p-value.
#'
#' The total number of windows is also returned in the final object, with the
#' summarised values n_up and n_down indicating the number of windows with raw
#' p-values below the calculated harmonic mean p-value, and with the
#' corresponding direction of change.
#'
#' The column containing the harmonic mean p-values is returned as 'hmp'.
#' An additional column with adjusted hmp-values is returned with the suffix
#'  '_*' added where the p-value adjustment method is added after the underscore.
#'
#' @param x GenomicRanges object
#' @param df data.frame with results of differential binding analysis performed
#' using a sliding window strategy. If not provided, the columns in the
#' `mcols()` element of `x` will be used
#' @param w vector of weights to applied when calculating harmonic mean p-values
#' @param logfc,pval,cpm Column names for the values holding window specific
#' estimates of change in binding (logfc), overall signal intensity (cpm) and
#' the significance from statistical testing (pval).
#' @param inc_cols (Optional) Character vector of any additional columns in
#' `df` to return. Values will correspond to the range in the `keyval_range`
#' column
#' @param p_adj_method One of `p.adjust.methods` or "fwer". If "fwer" is
#' specified the adjusted harmonic-mean p-value will be returned in a form
#' which strictly controls the experiment-wide FWER. Please see
#' vignette("harmonicmeanp") for more details
#' @param merge_within Merge any non-overlapping windows within this distance
#' @param ignore_strand Passed internally to \link[GenomicRanges]{reduce} and
#' \link[GenomicRanges]{findOverlaps}
#' @param min_win Only keep merged windows derived from at least this number
#' @param keyval Return the key-value range as the window associated with the
#' minimum p-value, or by merging the ranges from all windows with raw p-values
#' below the merged harmonic-mean p-value
#' @param hm_pre Prefix to add to the beginning of all HMP-derived columns
#' @param ... Not used
#'
#' @return
#' A GenomicRanges object with merged ranges from the original object along with
#' summarised or representative values from the relevant columns. The range
#' corresponding to a representative values is also returned as described above
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
#' set.seed(1001)
#' df <- DataFrame(logFC = rnorm(3), logCPM = rnorm(3,8), p = rexp(3, 10))
#' mergeByHMP(x, df, pval = "p")
#' mcols(x) <- df
#' x
#' mergeByHMP(x, pval = "p", p_adj_method = "fwer")
#'
#' @name mergeByHMP
#' @rdname mergeByHMP-methods
#' @export
setGeneric("mergeByHMP", function(x, ...) standardGeneric("mergeByHMP"))
#' @importFrom S4Vectors DataFrame mcols mcols<- subjectHits
#' @importFrom dplyr group_by summarise arrange distinct left_join bind_rows
#' @importFrom dplyr across
#' @importFrom tidyselect all_of
#' @importFrom rlang := !! sym
#' @importFrom stats weighted.mean
#' @import GenomicRanges
#' @rdname mergeByHMP-methods
#' @export
setMethod(
    "mergeByHMP",
    signature = signature(x = "GenomicRanges"),
    function(
        x, df = NULL, w = NULL,
        logfc = "logFC", pval = "P", cpm = "logCPM", inc_cols = NULL,
        p_adj_method = "fdr", merge_within = 1L, ignore_strand = TRUE,
        min_win = 1, keyval = c("min", "merged"), hm_pre = "hm", ...
    ){

        ## Checks & defining the key columns
        if (is.null(df)) df <- mcols(x)
        stopifnot(nrow(df) == length(x))
        df <- as.data.frame(df)
        df_cols <- colnames(df)
        logfc <- match.arg(logfc, df_cols)
        pcol <- match.arg(pval, df_cols, several.ok = TRUE)
        cpm <- match.arg(cpm, df_cols)
        p_adj_method <- match.arg(p_adj_method, c(p.adjust.methods, "fwer"))
        keyval <- match.arg(keyval)
        if (is.null(w)) {
            ## Need to sum \leq one using the HMP algorithms
            df[["weights"]] <- 1 / nrow(df)
        } else {
            df[["weights"]] <- w / sum(w)
        }

        ## Define the columns to return
        ret_cols <- c("keyval_range", cpm, logfc, pcol)
        if (!is.null(inc_cols)) {
            inc_cols <- vapply(
                inc_cols, match.arg, character(1), choices = df_cols
            )
            ret_cols <- unique(c(inc_cols, ret_cols))
        }
        ## Always return the columns counting the total, up & down windows
        n_cols <- c("n_windows", "n_up", "n_down")
        if ("index" %in% inc_cols)
            warning("The column named 'index' will not be returned")

        ## Merge the ranges and get the map back to the original windows
        ranges_out <- GenomicRanges::reduce(
            x, min.gapwidth = merge_within, ignore.strand = ignore_strand
        )
        ol <- findOverlaps(x, ranges_out, ignore.strand = ignore_strand)
        stopifnot(length(ol) == length(x))

        ## Summarise the key values
        df[["subjectHits"]] <- subjectHits(ol)
        grp_df <- group_by(df, subjectHits)
        ref_p <- pcol[[1]]
        ref_hmp <- paste0(hm_pre, ref_p)
        ret_df <- summarise(
            grp_df,
            n_windows = dplyr::n(),
            across(
                all_of(c(cpm, logfc)), \(x) weighted.mean(x, 1 / !!sym(ref_p))
            ),
            across(
                all_of(pcol), \(x) .ec_HMP(x, !!sym("weights")),
                .names = "{hm_pre}{.col}"
            ),
            n_up = sum(!!sym(logfc) > 0 & !!sym(ref_p) < !!sym(ref_hmp)),
            n_down = sum(!!sym(logfc) < 0 & !!sym(ref_p) < !!sym(ref_hmp))
        )
        ## Replace the 'pval' column in the return columns with 'hmp'
        new_cols <- paste0(hm_pre, ret_cols[ret_cols %in% pcol])
        ret_cols[ret_cols %in% pcol] <- new_cols

        ## Remaining columns, including the keyval range are based on the
        ## minimal p-value. Form a separate df with these columns, then merge
        inc_df <- as.data.frame(df[c(pcol, inc_cols)])
        inc_df[["keyval_range"]] <- as.character(x)
        inc_df[["subjectHits"]] <- subjectHits(ol)
        inc_df <- arrange(inc_df, !!sym("subjectHits"), !!sym(ref_p))
        inc_df <- distinct(inc_df, !!sym("subjectHits"), .keep_all = TRUE)
        inc_df <- inc_df[!names(inc_df) %in% pcol]

        ## If keyval = "merged" (as opposed to "min")
        if (keyval == "merged") { ## Reform the keyvalue range by merging
            df$range <- as.character(x)
            df <- left_join(
                df, ret_df[c("subjectHits", new_cols)], by = "subjectHits"
            )
            df <- dplyr::filter(
                df,
                (
                    !!sym(ref_p) <= !!sym(ref_hmp) |
                     !!sym(ref_p) == min(!!sym(ref_p))
                ),
                .by = !!sym("subjectHits")
            )
            kv_gr <- colToRanges(df, "range")
            kv_grl <- splitAsList(kv_gr, kv_gr$subjectHits)
            kv_grl <- range(kv_grl)
            kv_gr <- unlist(kv_grl)
            mcols(kv_gr)[["subjectHits"]] <- as.integer(names(kv_gr))
            kv_df <- as_tibble(kv_gr, name = "keyval_range")
            inc_df <- inc_df[names(inc_df) != "keyval_range"]
            inc_df <- left_join(
                inc_df, kv_df, by = "subjectHits", relationship = "one-to-one"
            )
        }

        ## Merge the two, select the final columns, adjust p & return
        ret_df <- left_join(
            ret_df, inc_df, by = "subjectHits", relationship = "one-to-one"
        )
        ret_df <- ret_df[c(n_cols, ret_cols)]
        if (p_adj_method == "fwer") {
            ## Apply the FWER adjusted version if requested
            L <- length(x)
            adj_df <- summarise(
                grp_df,
                across(
                    all_of(pcol), \(x) .ec_HMP_adj(x, !!sym("weights"), L),
                    .names = "{hm_pre}{.col}_fwer"
                ), .groups = "drop"
            )
            ## This could be revisited for better integration with filtering
            stopifnot(nrow(ret_df) == nrow(adj_df))
            ret_df <- bind_cols(ret_df, adj_df)
        }

        ## Apply the filter based on window size
        stopifnot(is.numeric(min_win))
        keep_win <- ret_df$n_windows >= min_win
        ranges_out <- ranges_out[keep_win]
        ret_df <- ret_df[keep_win,]

        ## Adjust p-values using p.adj.methods (after window size filtering)
        if (p_adj_method != "fwer") {
            ## This can now create multiple p-adjust columns.
            ret_df <- mutate(
                ret_df,
                across(
                    all_of(new_cols), \(x) p.adjust(x, method = p_adj_method),
                    .names = "{.col}_{p_adj_method}"
                )
            )
        }

        ## Add to mcols to the ranges
        mcols(ranges_out) <- as.data.frame(ret_df)
        ranges_out$keyval_range <- GRanges(
            ranges_out$keyval_range, seqinfo = seqinfo(ranges_out)
        )

        ranges_out
    }
)
#' @rdname mergeByHMP-methods
#' @export
setMethod(
    "mergeByHMP",
    signature = signature(x = "RangedSummarizedExperiment"),
    function(
        x, df = NULL, w = NULL,
        logfc = "logFC", pval = "P", cpm = "logCPM", inc_cols = NULL,
        p_adj_method = "fdr", merge_within = 1L, ignore_strand = FALSE,
        hm_pre = "hm", ...
    ) {

        gr <- rowRanges(x)
        mergeByHMP(
            gr, df, w, logfc, pval, cpm, inc_cols, p_adj_method, merge_within,
            ignore_strand, hm_pre = hm_pre, ...
        )

    }

)


# This is a modified version of harmonicmeanp::p.hmp developed by Prof Daniel
# Wilson, and hardwired to simply return a combined asymptotically exact HMP.
# Hardwiring like this gives a 10-fold speed-up. Further modifications may be
# possible, but this seems enough for now
# @param p vector of p-values
# @param w vector of weights
#' @useDynLib extraChIPs, .registration = TRUE
#' @keywords internal
.ec_HMP <- function(p, w) {
    n <- length(p)
    hmp <- sum(w) / sum(w / p)
    loc <- log(n) + 1 + digamma(1) - log(2/pi)
    dbl <- double(1)
    # cout <- .C(
    #     "RtailsMSS", 1, 0, 1, loc, log(pi/2), 1.0, 1/hmp,
    #     out1 = dbl, out2 = dbl, out3 = dbl, out4 = dbl,
    #     out5 = dbl, out6 = dbl, COPY = rep(c(FALSE, TRUE), c(7, 6)),
    #     PACKAGE = "FMStable"
    # )
    # Still need to figure out if all of those outputs above can be dropped
    Rcout <- .C(
        "my_RtailsMSS", loc, 1/hmp, d = dbl, logd = dbl, `F` = dbl, logF = dbl,
        cF = dbl, logcF = dbl,
        COPY = rep(c(FALSE, TRUE), c(2, 6)), PACKAGE = "extraChIPs"
    )
    Rcout$cF
}

# Similar to the above, this produces the FWER-controlled version in a
# streamlined way
# @param p vector of p-values
# @param w vector of weights
# @param L Number of global tests
#' @useDynLib extraChIPs, .registration = TRUE
#' @keywords internal
.ec_HMP_adj <- function(p, w, L) {
    hmp <- sum(w) / sum(w / p)
    w.sum <- sum(w)
    loc <- log(L[[1]]) + 1 + digamma(1) - log(2/pi)
    dbl <- double(1)
    Rcout <- .C(
        "my_RtailsMSS", loc, w.sum/hmp, d = dbl, logd = dbl, `F` = dbl,
        logF = dbl, cF = dbl, logcF = dbl,
        COPY = rep(c(FALSE, TRUE), c(2, 6)), PACKAGE = "extraChIPs"
    )
    Rcout$cF
}


