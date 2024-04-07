#' @title Plot Overlaps Between List Elements
#'
#' @description Plot Overlaps between list elements as an upset or Venn diagram
#'
#' @details
#' This function should give the capability to show overlaps for any number of
#' replicates or groups, or a list of items such as gene names.
#' For n = 2, a scaled Venn Diagram will be produced, however no scaling is
#' implemented for n = 3
#'
#' UpSet plots are possible for any lists with length > 1, and are the only
#' implemented possibility for lists > 3.
#'
#' If the input is a `GRangesList` an additional boxplot can be requested
#' using any numeric column within the existing `mcols()` element.
#' Values will be summarised across all elements using the requested function
#' and the boxplot will be included as an upper panel above the intersections
#'
#' @return
#' Either a VennDiagram (i.e. grid) object, or a ComplexUpset plot
#'
#' @param x GRangesList of S3 list able to be coerced to character vectors
#' @param type The type of plot to be produced
#' @param var Column to summarised as a boxplot in an upper panel
#' (UpSet plot only)
#' @param f Summarisation function. Must return a single value from any
#' numeric vector
#' @param set_col Colours to be assigned to each set
#' @param ... Passed to \link[VennDiagram]{draw.pairwise.venn} (or
#' `draw.single/triple.venn`) for Venn Diagrams, and to
#' \link[ComplexUpset]{upset} for UpSet plots
#' @param .sort_sets passed to `sort_sets` in \link[ComplexUpset]{upset}
#' @param ignore.strand Passed to \link[GenomicRanges]{reduce}
#' @param merge_within Passed to \link{makeConsensus}
#' @param sz_sets Text size for set size labels. Passed internally to
#' `geom_text(size = sz_sets)`
#' @param hj_sets Horizontal adjustment of set size labels
#' @param exp_sets X-axis expansion for set size panel
#'
#' @examples
#' ## Examples using a list of character vectors
#' ex <- list(
#'   x = letters[1:5], y = letters[c(6:15, 26)], z = letters[c(2, 10:25)]
#' )
#' plotOverlaps(ex, type = "upset")
#' plotOverlaps(ex, type = "venn", set_col = 1:3, alpha = 0.3)
#' plotOverlaps(ex, type = "upset", set_col = 1:3, labeller = stringr::str_to_title)
#' plotOverlaps(ex[1:2])
#'
#' ## GRangesList object will produce a boxplot of summarised values in the
#' ## upper panel
#' data("peaks")
#' grl <- peaks[1:3]
#' names(grl) <- gsub("_peaks.+", "", names(grl))
#' plotOverlaps(grl, type = 'upset', var = 'score', f = 'max')
#'
#' ## If only two samples are present, a VennDiagram will be produced
#' plotOverlaps(grl[1:2], set_col = c("green", "blue"))
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply mcols
#' @importFrom dplyr bind_cols
#' @importFrom rlang list2 := !! sym
#' @importFrom ComplexUpset upset
#' @import ggplot2
#' @rdname plotOverlaps-methods
#' @aliases plotOverlaps
#' @export
setMethod(
  "plotOverlaps",
  signature = "GRangesList",
  function(
    x, type = c("auto", "venn", "upset"), var = NULL,
    f = c("mean", "median", "max", "min", "sd"),
    set_col = NULL, ..., .sort_sets = "ascending", hj_sets = 1.15,
    sz_sets = 3.5, exp_sets = 0.25, merge_within = 1L, ignore.strand = TRUE
  ) {

    stopifnot(is(x, "GRangesList"))
    nm <- names(x)
    n <- length(x)
    stopifnot(length(nm) == n)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")
    if (n == 1 & type == "upset")
      stop("UpSet plots can only be drawn using more than one group")
    if (!is.null(var)) var <- match.arg(var[[1]], .mcolnames(x[[1]]))

    # Collapse as required
    gr <- makeConsensus(
      x, var = var, merge_within = merge_within, ignore.strand = ignore.strand
    )

    if (is.null(var) | type == "venn") {

      ## Form a character list & plot
      l <- lapply(nm, function(x) as.character(gr)[mcols(gr)[[x]]])
      names(l) <- nm
      plotOverlaps(
        l, type, set_col = set_col, .sort_sets = .sort_sets,
        hj_sets = hj_sets, sz_sets = sz_sets, exp_sets = exp_sets, ...
      )

    } else {

      if (!is.numeric(mcols(x[[1]])[[var]]))
        stop(var, " must contain numeric values")

      ## Setup the df
      tbl <- as_tibble(gr)
      f <- match.arg(f)
      f <- match.fun(f)
      if (is(tbl[[var]], "list"))
        tbl[[var]] <- vapply(tbl[[var]], f, numeric(1))

      ## Setup the boxplot & key inputs
      ann <- list2(
        "{var}" := list(
          aes = aes(x = !!sym("intersection"), y = !!sym(var)),
          geom = geom_boxplot(na.rm = TRUE)
        )
      )
      ip <- list(data = tbl, intersect = nm, annotations = ann)

      ## Add default arguments, respecting any supplied
      dotArgs <- .parseDotArgs(set_col, n, nm, ...)
      if (!"set_sizes" %in% names(dotArgs))
        dotArgs$set_sizes <- .makeSetSizes(hj_sets, sz_sets, exp_sets)
      dotArgs$sort_sets <- .sort_sets
      ip <- ip[!names(ip) %in% names(dotArgs)]
      p <- do.call("upset", c(ip, dotArgs))
      return(p)
    }

  }
)
#'
#' @importFrom methods is
#' @importFrom ComplexUpset upset
#' @importFrom grid grid.newpage
#' @rdname plotOverlaps-methods
#' @aliases plotOverlaps
#' @export
setMethod(
  "plotOverlaps",
  signature = "list",
  function(
    x, type = c("auto", "venn", "upset"), set_col = NULL, ...,
    .sort_sets = 'ascending', hj_sets = 1.15, sz_sets = 3.5, exp_sets = 0.25
  ) {

    stopifnot(length(names(x)) == length(x))
    n <- length(x)
    nm <- names(x)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")
    x <- lapply(x, as.character)
    x <- lapply(x, unique)

    if (type == "upset") {
      count <- c()
      if (n == 1)
        stop("UpSet plots can only be drawn using more than one group")
      ## Setup the df
      all_vals <- unique(unlist(x))
      df <- lapply(x, function(i) as.integer(all_vals %in% i))

      ## Ensure colnames are respected
      col_names <- names(df)
      df <- as.data.frame(df, row.names = all_vals)
      colnames(df) <- col_names
      ip <- list(data = df, intersect = names(df))

      ## Add default arguments, respecting any supplied
      dotArgs <- .parseDotArgs(set_col, n, nm, ...)
      if (!"set_sizes" %in% names(dotArgs))
        dotArgs$set_sizes <- .makeSetSizes(hj_sets, sz_sets, exp_sets)
      dotArgs$sort_sets <- .sort_sets
      p <- do.call("upset", c(ip, dotArgs))
      return(p)
    }

    if (type == "venn") {
      grid.newpage()
      if (n == 1) p <- .plotSingleVenn(x, fill = set_col, ...)
      if (n == 2) p <- .plotDoubleVenn(x, fill = set_col, ...)
      if (n == 3) p <- .plotTripleVenn(x, fill = set_col, ...)
    }
    invisible(p)

  }
)

#' @importFrom VennDiagram draw.single.venn
.plotSingleVenn <- function(x, ...) {
  stopifnot(length(x) == 1)
  draw.single.venn(area = length(x[[1]]), category = names(x)[[1]], ...)
}

#' @importFrom VennDiagram draw.pairwise.venn
.plotDoubleVenn <- function(x, ...) {
  stopifnot(length(x) == 2)
  plotArgs <- setNames(lapply(x, length), c("area1", "area2"))
  plotArgs$cross.area <- sum(duplicated(unlist(x)))
  plotArgs$category <- names(x)
  allowed <- c("gList1", "margin", names(formals(draw.pairwise.venn)))
  dotArgs <- list(...)
  dotArgs <- dotArgs[names(dotArgs) %in% allowed]
  do.call("draw.pairwise.venn", c(plotArgs, dotArgs))

}

#' @importFrom VennDiagram draw.triple.venn
.plotTripleVenn <- function(x, ...) {
  stopifnot(length(x) == 3)
  plotArgs <- setNames(lapply(x, length), paste0("area", seq_len(3)))
  plotArgs$n12 <- sum(duplicated(unlist(x[c(1, 2)])))
  plotArgs$n13 <- sum(duplicated(unlist(x[c(1, 3)])))
  plotArgs$n23 <- sum(duplicated(unlist(x[c(2, 3)])))
  plotArgs$n123 <- sum(table(unlist(x)) == 3)
  plotArgs$category <- names(x)
  plotArgs$overrideTriple <- TRUE
  allowed <- c("gList1", "margin", names(formals(draw.triple.venn)))
  dotArgs <- list(...)
  dotArgs <- dotArgs[names(dotArgs) %in% allowed]
  do.call("draw.triple.venn", c(plotArgs, dotArgs))
}

#' @importFrom ComplexUpset upset_set_size
#' @importFrom scales comma
#' @importFrom rlang sym !!
.makeSetSizes <- function(hj, sz, exp) {
  upset_set_size() +
    geom_text(
      aes(label = comma(after_stat(!!sym("count")))),  stat = 'count',
      hjust = hj, size = sz
    ) +
    scale_y_reverse(expand = expansion(c(exp, 0)))
}

#' @importFrom ComplexUpset upset upset_default_themes upset_data upset_query
.parseDotArgs <- function(set_col, n, nm, ...) {
  dotArgs <- list(...)
  allowed <- unique(names(c(formals(upset), formals(upset_data))))
  dotArgs <- dotArgs[names(dotArgs) %in% allowed]
  if (!'themes' %in% names(dotArgs)) {
    dotArgs$themes <- upset_default_themes(panel.grid = element_blank())
  }
  ## There is currently an issue with ComplexUpset. This places theme arguments
  ## which are not supported beyond ggplot2 3.5.0. This will remove them
  valid_theme_args <- names(formals(theme))
  dotArgs$themes <- lapply(
      dotArgs$themes,
      \(x) lapply(x, \(y) do.call("theme", y[names(y) %in% valid_theme_args]))
  )

  if (!is.null(set_col)) {
    ## Respect any existing set queries
    existing_sets <- unlist(
      lapply(dotArgs$queries, function(x) x$set)
    )
    set_col <- rep(set_col, n)
    names(set_col)[seq_len(n)] <- nm
    ql <- lapply(
      setdiff(nm, existing_sets),
      function(i) upset_query(set = i, fill = set_col[[i]])
    )
    dotArgs$queries <- c(dotArgs$queries, ql)
  }
  dotArgs
}

#' @importFrom S4Vectors mcols
#' @keywords internal
.mcolnames <- function(x) {colnames(mcols(x))}
