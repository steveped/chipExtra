#' @title Import peaks
#'
#' @description Import peaks in narrowPeak, broadPeak or bed format
#'
#' @details
#' Peaks are imported from narrowPeak, broadPeak or bed format as
#' GenomicRanges objects.
#'
#' If importing bed files, only the default 3-6 columns will imported.
#'
#' @param x One or more files to be imported. All files must be of the same
#' type, i.e. narrow or broad
#' @param type The type of peaks to be imported
#' @param blacklist A set of ranges to be excluded
#' @param seqinfo A seqinfo object to be applied to the GRanges objects
#' @param pruning.mode How to handle conflicts if supplying a seqinfo object.
#' Defaults to `pruning.mode = "coarse"`. Only "coarse" and "error" are
#' implemented. See \link[GenomeInfoDb]{seqinfo}.
#' @param sort logical. Should the ranges be sorted during import
#' @param setNames logical Set basename(x) as the name for each element of the
#' GRangesList
#' @param glueNames \link[glue]{glue} syntax for naming list elements
#' @param centre Add the estimated peak centre. Ignored unless type = "narrow"
#' @param nameRanges Place any values in the name column as range names within
#' each file.
#' @param ... passed to `sort`
#'
#' @return
#' A GRangesList
#'
#' @examples
#' fl <- system.file(
#'     c("extdata/ER_1.narrowPeak", "extdata/ER_2.narrowPeak"),
#'     package = "extraChIPs"
#' )
#' peaks <- importPeaks(fl)
#' peaks
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom glue glue
#' @import GenomicRanges
#' @export
importPeaks <- function(
        x, type = c("narrow", "broad", "bed"), blacklist, seqinfo,
        pruning.mode = c("coarse", "error"), sort = TRUE, setNames = TRUE,
        glueNames = "{basename(x)}", centre = FALSE, nameRanges = TRUE, ...
) {

    ## Argument checks
    stopifnot(file.exists(x))
    type <- match.arg(type)
    pruning.mode <- match.arg(pruning.mode)
    stopifnot(is.logical(sort) & is.logical(setNames))
    n <- length(x)
    if (type == "bed") {
        out <- lapply(
            x,
            .importBedFile,
            seqinfo = seqinfo, blacklist = blacklist,
            pruning.mode = pruning.mode, nameRanges = nameRanges
        )
    } else {
        out <- lapply(
            x,
            .importPeakFile,
            type = type, seqinfo = seqinfo, blacklist = blacklist,
            pruning.mode = pruning.mode, centre = centre, nameRanges = nameRanges
        )
    }

    out <- GRangesList(out)
    if (sort) sort(out, ...)
    if (setNames) {
        nm <- as.character(glue(glueNames))
        stopifnot(length(nm) == length(out))
        names(out) <- nm
    }
    out

}

#' @import GenomicRanges
#' @importFrom IRanges overlapsAny
#' @importFrom methods is
#' @importFrom utils read.table
.importPeakFile <- function(
        x, type, seqinfo, blacklist, pruning.mode, centre, nameRanges
) {

    stopifnot(length(x) == 1)
    if (!missing(seqinfo)) {
        stopifnot(is(seqinfo, "Seqinfo"))
    } else {
        seqinfo <- NULL
    }

    ## Handle empty files separately first
    if (.isEmpty(x)) return(GRanges(NULL, seqinfo = seqinfo))

    ## Define the parameters for the two types
    colNames <- c(
        "seqnames", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "peak"
    )
    if (type == "narrow") stopifnot(.isValidNarrow(x))
    if (type == "broad") {
        stopifnot(.isValidBroad(x))
        colNames <- setdiff(colNames, "peak")
    }
    classes <- c("character", "numeric", "numeric", "character", "numeric",
                 "character", rep("numeric", 4))[seq_along(colNames)]

    ## Parse
    df <- read.table(x, sep = "\t", col.names = colNames, colClasses = classes)
    if (nameRanges) {
        if (length(unique(df[[4]])) == nrow(df)) rownames(df) <- df[["name"]]
        df <- df[-4] # The name column has been set as rownames
    }

    ## Deal with the Seqinfo object
    if (!is.null(seqinfo)) {

        ## Fail if there is a mismatch between any seqnames and seqinfo
        if (pruning.mode == "error")
            stopifnot(all(df[["seqnames"]] %in% seqnames(seqinfo)))
        ## Subset to remove any ranges which are not included in seqinfo
        if (pruning.mode == "coarse")
            df <- subset(df, seqnames %in% seqnames(seqinfo))
        if (nrow(df) == 0)
            message("No ranges match the supplied seqinfo object")
    }

    ## Form the object
    gr <- makeGRangesFromDataFrame(
        df, TRUE, seqinfo = seqinfo, starts.in.df.are.0based = TRUE
    )
    if (type == "narrow" & centre)
        gr$centre <- start(gr) + gr$peak

    ## Apply the blacklist if supplied
    if (missing(blacklist)) blacklist <- GRanges(seqinfo = seqinfo)
    stopifnot(is(blacklist, "GRanges"))
    gr[!overlapsAny(gr, blacklist)]

}

#' @import GenomicRanges
#' @importFrom IRanges overlapsAny
#' @importFrom methods is
#' @importFrom utils read.table
.importBedFile <- function(x, seqinfo, blacklist, pruning.mode, nameRanges) {

    stopifnot(length(x) == 1)
    if (!missing(seqinfo)) {
        stopifnot(is(seqinfo, "Seqinfo"))
    } else {
        seqinfo <- NULL
    }

    ## Handle empty files separately first
    if (.isEmpty(x)) return(GRanges(NULL, seqinfo = seqinfo))

    ## Define the colnames
    stopifnot(.isValidBed(x))
    nCol <- min(ncol(read.table(x, sep = "\t", nrows = 1)), 6)
    colNames <- c("seqnames", "start", "end", "name", "score", "strand")
    classes <- c(
        "character", "numeric", "numeric", "character", "numeric", "character"
    )

    ## Parse
    df <- read.table(x, sep = "\t", header = FALSE)[seq_len(nCol)]
    df <- lapply(
        seq_len(nCol), function(i) df[[i]] <- as(df[[i]], classes[[i]])
    )
    df <- as.data.frame(df)
    colnames(df) <- colNames[seq_len(nCol)]

    if (nCol >= 4 & nameRanges) {
        if (length(unique(df[[4]])) == nrow(df)) rownames(df) <- df[[4]]
        df <- df[-4] # The name column has been set as rownames
    }

    ## Deal with the Seqinfo object
    if (!is.null(seqinfo)) {

        ## Fail if there is a mismatch between any seqnames and seqinfo
        if (pruning.mode == "error")
            stopifnot(all(df[["seqnames"]] %in% seqnames(seqinfo)))
        ## Subset to remove any ranges which are not included in seqinfo
        if (pruning.mode == "coarse")
            df <- subset(df, seqnames %in% seqnames(seqinfo))
        if (nrow(df) == 0)
            message("No ranges match the supplied seqinfo object")
    }

    ## Form the object
    gr <- makeGRangesFromDataFrame(
        df, TRUE, seqinfo = seqinfo, starts.in.df.are.0based = TRUE
    )

    ## Apply the blacklist if supplied
    if (missing(blacklist)) blacklist <- GRanges(seqinfo = seqinfo)
    stopifnot(is(blacklist, "GRanges"))
    gr[!overlapsAny(gr, blacklist)]

}


#' @importFrom utils read.table
.isValidNarrow <- function(x) {
    r1 <- read.table(x, sep = "\t", nrows = 1)
    nCols <- ncol(r1) == 10
    if (!nCols) return(FALSE)
    charCols <- c(1, 4)
    allChars <- all(vapply(r1[charCols], is.character, logical(1)))
    numericCols <- c(2, 3, 5, seq(7, 10))
    allNumerics <- all(vapply(r1[numericCols], is.numeric, logical(1)))
    strandOK <- r1[[6]] %in% c("+", "-", ".")
    all(allChars, allNumerics, strandOK)
}

#' @importFrom utils read.table
.isValidBroad <- function(x) {
    r1 <- read.table(x, sep = "\t", nrows = 1)
    nCols <- ncol(r1) == 9
    if (!nCols) return(FALSE)
    charCols <- c(1, 4)
    allChars <- all(vapply(r1[charCols], is.character, logical(1)))
    numericCols <- c(2, 3, 5, seq(7, 9))
    allNumerics <- all(vapply(r1[numericCols], is.numeric, logical(1)))
    strandOK <- r1[[6]] %in% c("+", "-", ".")
    all(allChars, allNumerics, strandOK)
}

#' @importFrom utils read.table
.isValidBed <- function(x){
    r1 <- read.table(x, sep = "\t", nrows = 1)
    nCols <- ncol(r1)
    if (nCols < 3) return(FALSE)
    numericCols <- c(2, 3)
    allNumerics <- all(vapply(r1[numericCols], is.numeric, logical(1)))
    strandOK <- TRUE
    if (nCols >= 6) {
        strandOK <- r1[[6]] %in% c("+", "-", ".")
    }
    all(allNumerics, strandOK)
}

.isEmpty <- function(x) {
    ## This is a better approach then checking file.size, given that a
    ## compressed empty file will have size ~20B
    !length(readLines(x, n = 1))
}
