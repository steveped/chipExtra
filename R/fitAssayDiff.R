#' @title Detect Differential ChIP Signal
#'
#' @description Detect differential ChIP signal using one of many approaches
#'
#' @details
#' Starting with a SummarizedExperiment object this function fits either a
#' \link[edgeR]{glmQLFit} model to count data, performs the
#' \link[DESeq2]{nbinomWaldTest} on count data, or applies the
#' \link[limma:eBayes]{limma-trend} model to logCPM data.
#'
#' If fitting Generalised Linear Models via glmQLFit, options for normalisation
#' are "none", which normalises to library size. Existing library sizes are
#' commonly found in the "totals" column of the colData element and this is
#' attempted by default. All methods provided in \link[edgeR]{normLibSizes}
#' are also implemented, with the added possibility of normalising within groups
#' instead of across the entire dataset. To enable this, the column with the
#' grouping factor is expected to be in the colData element and is simply
#' called by column name.
#' No normalisation is applied when using the limma-trend model, as this allows
#' for previous normalisation strategies to be performed on the data.
#'
#' If testing with \link[DESeq2]{nbinomWaldTest}, applying RLE normalisation
#' without groups, and using colSums for library sizes (instead of total
#' alignments), the standard normalisation factors from
#' \link[DESeq2]{estimateSizeFactors} will be used.
#' In all other scenarios, normalisation factors as returned by
#' \link[edgeR]{normLibSizes} will be used.
#' The fitType is set to 'local' when estimating dispersions, and this can be
#' easily modified by passing fitType via the dot arguments.
#' Results are additionally returned after applying \link[DESeq2]{lfcShrink},
#' including the svalue returned by this approach.
#'
#' Normalising to ChIP Input samples is not yet implemented.
#' Similarly, the use of offsets when applying the Wald test is not yet
#' implemented.
#'
#' Range-based hypothesis testing is implemented using
#' \link[edgeR]{glmTreat} or \link[limma]{treat}. Setting fc to 1 (or lfc to 0)
#' will default to a point-based null hypothesis, equivalent to either
#' \link[edgeR]{glmQLFTest} (method = "qlf") or
#' \link[limma]{eBayes} (method = "lt").
#' When applying \link[DESeq2]{nbinomWaldTest}, \link[DESeq2]{lfcShrink} will
#' be applied.
#'
#' It should also be noted that this is primarily a convenience function and
#' if requiring intermediate output from any steps, then these can be run
#' individually as conventionally specified.
#'
#' @return
#' A SummarizedExperiment object with results set as the `rowData` element.
#' Any existing columns not contained in the differential ChIP results will be
#' retained.
#' Results from testing will contain logCPM, logFC, PValue and any t/F
#' statistic as appropriate, along with an FDR-adjusted p-value.
#'
#' If specifying a range-based H0 by setting lfc != 0, an additional column
#' p_mu0 will be included which is the p-value for the point H0: logFC = 0.
#' These are not used for FDR-adjusted p-values but can be helpful when
#' integrating multiple ChIP targets due to the increase in false-negatives when
#' using a range-based H0, and when requiring more accurate identification of
#' truly unchanged sites, as opposed to those which simply fail to achieve
#' significance using a range-based H0 where arbitrary cutoff values are used.
#'
#' @param x a SummarizedExperiment object
#' @param assay The assay to use for analysis
#' @param design The design matrix to use for analysis
#' @param coef The required column from the design matrix
#' @param lib.size The column within the colData element which contains the
#' library size information. If set to NULL, column summaries will be used.
#' @param method the analytic method to be used. Can be 'qlf' which will fit
#' counts using the \link[edgeR]{glmQLFit} strategy , or 'lt' which fits the
#' \link[limma:eBayes]{limma-trend} model on logCPM, or pre-processed logCPM
#' values. Setting method = 'wald' will call \link[DESeq2]{nbinomWaldTest}
#' @param norm The normalisation strategy to use when running the
#' glmQLF model or the Wald test. The value 'none' relies solely on
#' library-size normalisation, and is the default. All methods available in
#' \link[edgeR]{normLibSizes} are implemented. Ignored when using method = "lt"
#' @param groups character(1) If a column name is supplied here, group-based
#' normalisation will be applied to GLM models treating data in this column
#' as a grouping factor. Ignored when using method = "lt"
#' @param fc,lfc Thresholds passed to \link[limma]{treat},
#' \link[edgeR]{glmTreat} or \link[DESeq2]{lfcShrink}
#' @param asRanges logical(1). By default, the returned object will be a
#' `SummarizedExperiment` object with the results added to the `rowData`
#' element. Setting `asRanges = TRUE` will only return the GRanges object from
#' this element
#' @param offset If provided will be used as the offset when a DGEList object
#' is created during model fitting for method = 'qlf'
#' @param weighted logical(1) Passed to  \link[edgeR]{normLibSizes}. Only used
#' when applying a TMM-type normalisation strategy
#' @param ... Passed to \link[edgeR]{normLibSizes} and
#' @param null Passed to \link[edgeR]{glmTreat}
#' \link[edgeR]{glmQLFit} when method = "qlf".
#' If method = "lt", instead passed to \link[limma]{lmFit}
#' @param robust Passed to \link[limma]{treat} and \link[limma]{eBayes}
#' @param type Passed to \link[DESeq2]{lfcShrink}
#'
#' @examples
#' nrows <- 200; ncols <- 6
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' colnames(counts) <- paste0("Sample_", seq_len(ncols))
#' df <- DataFrame(treat = c("A", "A", "A", "B", "B", "B"))
#' df$treat <- as.factor(df$treat)
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts), colData = df
#' )
#' X <- model.matrix(~treat, colData(se))
#' se <- fitAssayDiff(se, design = X, lib.size = NULL)
#' rowData(se)
#'
#'
#' @name fitAssayDiff
#' @rdname fitAssayDiff-methods
#' @export
#'
setGeneric("fitAssayDiff", function(x, ...) standardGeneric("fitAssayDiff"))
#' @import SummarizedExperiment
#' @importFrom edgeR glmTreat topTags glmQLFTest
#' @importFrom stats model.matrix
#' @rdname fitAssayDiff-methods
#' @export
setMethod(
    "fitAssayDiff",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", design = NULL, coef = NULL,
        lib.size = "totals", method = c("qlf", "lt", "wald"),
        norm = c("none", "TMM", "RLE", "TMMwsp", "upperquartile"),
        groups = NULL, fc = 1, lfc = log2(fc), asRanges = FALSE,
        offset = NULL, weighted = FALSE, ...,
        null = c("interval", "worst.case"), robust = FALSE,
        type = c("apeglm", "ashr", "normal")
    ) {
        method <- match.arg(method)
        norm <- match.arg(norm)
        args <- colnames(colData(x))
        n <- nrow(x)
        if (is.null(design)) stop("A design matrix must be specified")
        stopifnot(nrow(design) == ncol(x))
        if (is.null(coef)) coef <- colnames(design)[ncol(design)]
        if (is.numeric(coef)) stopifnot(coef <= ncol(design))
        if (is.character(coef)) stopifnot(coef %in% colnames(design))
        if (!is.null(groups)) groups <- match.arg(groups, args)
        if (!is.null(lib.size)) lib.size <- match.arg(lib.size, args)

        if (method == "qlf") {
            ## Only required for GLM fits
            fit <- .se2DGEGLM(
                x, assay, design, lib.size, norm, groups, offset, weighted, ...
            )
            fit0 <- glmQLFTest(fit, coef = coef) # fits mu0
            res0 <- topTags(
                fit0, n = n, adjust.method = "none", sort.by = "none"
            )$table
            res <- res0
            p_mu0 <- res0$PValue
            if (lfc != 0) {
                null <- match.arg(null)
                fit <- glmTreat(fit, coef, lfc = lfc, null = null)
                res <- topTags(
                    fit, n = n, adjust.method = "none", sort.by = "none"
                )$table
                res <- res[!colnames(res) %in% "unshrunk.logFC"]
            }
        }
        if (method == "lt") {
            if (!requireNamespace('limma', quietly = TRUE))
                stop("Please install 'limma' to use this function.")
            fit <- .se2LT(x, assay, design, ...)
            fit0 <- limma::eBayes(fit, trend = TRUE, robust = robust)
            res0 <- limma::topTable(
                fit0, coef = coef, number = n, sort.by = "none",
                adjust.method = "none"
            )
            res <- res0
            p_mu0 <- res0$P.Value
            if (lfc != 0) {
                fit <- limma::treat(
                    fit, lfc = lfc, trend = TRUE, robust = robust
                )
                res <- limma::topTreat(
                    fit, coef = coef, number = n, sort.by = "none",
                    adjust.method = "none"
                )
            }
            res <- res[!colnames(res) %in% c("adj.P.Val", "B")]
            colnames(res) <- gsub("AveExpr", "logCPM", colnames(res))
            colnames(res) <- gsub("P.Value", "PValue", colnames(res))
        }
        if (method == "wald") {
            if (!requireNamespace('DESeq2', quietly = TRUE))
                stop("Please install 'DESeq2' to use this function.")
            type <- match.arg(type)
            dds <- .se2Wald(
                x, assay, design, lib.size, norm, groups, offset, weighted, ...
            )
            res <- DESeq2::results(dds, name = coef, lfcThreshold = 0)
            if (abs(lfc) > 0) {
                p_mu0 <- res$pvalue
                res <- DESeq2::results(dds, name = coef, lfcThreshold = abs(lfc))
            }
            ## Always return shrunken estimates. Maybe make optional later
            shrink <- DESeq2::lfcShrink(
                dds, coef = coef, res = res, type = type,
                lfcThreshold = abs(lfc), svalue = TRUE
            )
            ## Collate into consistent output
            res <- data.frame(
                logFC = shrink$log2FoldChange, logCPM = log2(shrink$baseMean),
                svalue = shrink$svalue, PValue = res$pvalue,
                row.names = rownames(shrink)
            )

        }
        res[["FDR"]] <- p.adjust(res[["PValue"]], "fdr")
        ## Add the p-value for H0: fc = 0
        if (lfc != 0) res$p_mu0 <- p_mu0

        keep_cols <- setdiff(colnames(rowData(x)), colnames(res))
        any_gr <- vapply(
            keep_cols, function(i) is(rowData(x)[[i]], "GRanges"), logical(1)
        )
        gr_cols <- keep_cols[any_gr]
        orig <- as_tibble(rowData(x))[,keep_cols]
        rowData(x) <- cbind(orig, res)
        list_cols <- vapply(rowData(x), function(x) is(x, "list"), logical(1))
        rowData(x)[list_cols] <- lapply(
            rowData(x)[list_cols], function(x) as(x, "CompressedList")
        )
        if (any(any_gr) & is(x, "RangedSummarizedExperiment")) {
            sq <- seqinfo(x)
            rowData(x)[gr_cols] <- lapply(
                rowData(x)[gr_cols], GRanges, seqinfo = sq
            )
        }
        if (asRanges) {
            if (is.null(rowRanges(x))) {
                warning(
                    "No ranges found. Results will be returned in the rowData",
                    "element of the original object"
                )
            } else {
                return(rowRanges(x))
            }
        }
        x
    }
)

#' @import SummarizedExperiment
.se2LT <- function(x, assay, design, ...){
    ## 1. Create an MArrayM/Elist object
    ## 2. Don't normalise
    mat <- assay(x, assay)
    diff <- mat - as.integer(mat)
    if (all(diff == 0)) stop("Expected non-integer values")
    limma::lmFit(mat, design = design, ...)
}

#' @import SummarizedExperiment
#' @importFrom edgeR DGEList normLibSizes estimateDisp glmQLFit
.se2DGEGLM <- function(
        x, assay, design, lib.size, norm, groups, offset, weighted, ...
) {

    ## 1. Create a DGE list
    ## 2. Normalise as requested
    ## 3. Apply the model returning an object of class DGEGLM
    mat <- assay(x, assay)
    if (any(mat < 0)) stop("Counts cannot contain negative values")
    if (anyNA(mat)) stop("Missing values detected")

    col_df <- as_tibble(colData(x), rownames = "colnames")
    ls <- colSums(mat)
    if (!is.null(lib.size)) ls <- col_df[[lib.size]]
    message("Creating DGE list...")
    dge <- DGEList(counts = mat, lib.size = ls, samples = col_df)
    if (!is.null(offset)) dge$offset <- edgeR::scaleOffset(ls, offset)
    if (!is.null(groups)) {
        grp_fac <- as.factor(col_df[[groups]])
        dge$samples$group <- grp_fac
        split_df <- split(dge$samples, dge$samples[[groups]])
        list_cols <- lapply(split_df, rownames)
        message("Calculating group-wise normalisation factors...")
        nf <- lapply(
            list_cols,
            function(i) normLibSizes(dge$counts[,i], method = norm, ...)
        )
        names(nf) <- NULL
        nf <- unlist(nf)
        dge$samples$norm.factors <- nf[colnames(dge)]
    } else {
        message("Calculating experiment-wide normalisation factors...")
        ## Both DiffBind and csaw set doWeighting = FALSE. They're smart people
        dge <- normLibSizes(dge, method = norm, doWeighting = weighted, ...)
    }
    message("Estimating dispersions...")
    # dge <- estimateDisp(dge, design = design, ...)
    ## Switch to the method used by DiffBind
    dge <- edgeR::estimateGLMTrendedDisp(dge, design = design)
    dge <- edgeR::estimateGLMTagwiseDisp(dge, design = design)
    message("Running glmQLFit...")
    glmQLFit(dge, design = design, ...)

}

#' @import SummarizedExperiment
.se2Wald <- function(
        x, assay, design, lib.size, norm, groups, offset, weighted, ...
) {

    mat <- assay(x, assay)
    if (any(mat < 0)) stop("Counts cannot contain negative values")
    if (anyNA(mat)) stop("Missing values detected")

    col_df <- as_tibble(colData(x), rownames = "colnames")
    ls <- colSums(mat)
    if (!is.null(lib.size)) ls <- col_df[[lib.size]]

    dotArgs <- list(...)

    message("Creating DESeqDataSet...")
    dds <- DESeq2::DESeqDataSetFromMatrix(mat, colData = col_df, design = design)
    ## According to DiffBind, now we:
    ## 1. Tidy up the offsets for DESeq2
    if (!is.null(offset))
        warning("The use of offsets with DESeq2 is not currently implemented")

    ## 2. set Normalisaztion
    ## The only time to use DESeq2 is when lib.size is NULL as this method only
    ## allows for colSums as library sizes
    nf <- rep_len(NA, ncol(mat)) # Should error out if something fails below
    if (is.null(lib.size) & norm == "RLE" & is.null(groups)) {
        message("Calculating default DESeq2 normalisation factors...")
        dds <- DESeq2::estimateSizeFactors(dds)
    } else {
        if (is.null(groups)) {
            message("Calculating ", norm, " normalisation factors...")
            nf <- edgeR::normLibSizes(
                mat, method = norm, lib.size = ls, doWeighting = weighted
            )
        } else {
            message("Calculating ", norm, " normalisation factors by group...")
            grp_fac <- as.factor(col_df[[groups]])
            grp_lev <- levels(grp_fac)
            for (i in grp_lev) {
                idx <- grp_fac == i
                nf[idx] <- edgeR::normLibSizes(
                    mat[,idx], method = norm, lib.size = ls[idx],
                    doWeighting = weighted
                )
            }
        }
        DESeq2::sizeFactors(dds) <- nf
    }

    ## 3. Estimate dispersions
    ft <- "local"
    if ("fitType" %in% names(dotArgs)) ft <- dotArgs$fitType
    message("Estimating dispersions using fitType = ", ft, "...")
    dds <- DESeq2::estimateDispersions(dds, fitType = ft)

    ## 4. Fit (using nbinomWaldTest as default)
    dds <- DESeq2::nbinomWaldTest(dds, modelMatrix = design)

    ## Apply the fc threshold in the parent function
    dds
}
