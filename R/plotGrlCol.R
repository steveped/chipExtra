#' @title Draw a plot from a GRangesList column
#'
#' @description
#' Draw a plot from a GRangesList column using ggplot2
#'
#' @details
#' Using a common column or the width of the ranges, produces a boxplot or
#' violinplot from each element of the provided GRangesList.
#' The names of the GRangesList will be passed to the x-axis using the `.id`
#' argument.
#' A data frame containing annotations corresponding to each element can be
#' supplied, ensuring that the column associated with each elements is the name
#' passed to the `.id` argument.
#'
#' If q is > 0, a horizontal line will be draw corresponding to this percentile
#' across the complete dataset, with parameters for this line able to be set
#' using the qline_* arguments.
#' The digits argument controls how many decimal points will be shown for the
#' associated label.
#'
#' The total length of each element will be added by default as a total, and is
#' able to be placed across the median values, or at the top and bottom
#' extremes of the plot.
#'
#' @return A ggplot object
#'
#' @param x A GRangesList
#' @param var The variable to plot. Either a column in the mcols element or
#' width. Can be quoted or unquoted
#' @param geom Choose between different geoms, or even provide a geom_*()
#' function
#' @param .id The column name to place the element names. Passed internally to
#' the same argument in \link[dplyr]{bind_rows}
#' @param df Optional data.frame with columns to be passed to the colour or
#' fill parameters. Must contain a column with the same name as the
#' value passed to the `.id` argument.
#' @param fill,colour Optional column names found in the df. Can be quoted or
#' unquoted
#' @param q The overall percentile to be drawn as a labelled, horizontal line.
#' Set q = 0 to hide this line
#' @param q_size Text size of percentile label
#' @param qline_type,qline_col Linetype and colour arguments for the horizontal
#' line showing the specified percentile(s)
#' @param total Glue syntax for totals, representing the length of each
#' GRangesList element
#' @param total_geom Passed to \link[ggplot2]{annotate}. Set to `none` to hide
#' totals
#' @param total_size,total_alpha Size and transparency of totals
#' @param total_pos Position for placing totals
#' @param total_adj Adjustment for labels
#' @param ... Passed to the geom if selecting via character string. Ignored
#' otherwise
#' @param digits Number of decimal places for the horizontal line label
#'
#' @examples
#' ## Load some peaks
#' data('peaks')
#' names(peaks) <- gsub("_peaks.+", "", names(peaks))
#'
#' ## The default boxplot
#' plotGrlCol(peaks)
#'
#' ## A customised violin plot
#' df <- data.frame(sample = names(peaks), treat = rep(c("A", "B"), each = 3))
#' plotGrlCol(
#'   peaks, geom = "violin", total_pos = "bottom", total_adj = 0.05,
#'   df = df, fill = "treat",
#'   draw_quantiles = 0.5, trim = FALSE, width = 0.7, alpha = 0.7
#' ) +
#' scale_y_log10()
#'
#' plotGrlCol(
#'   peaks, var = score, geom = "jitter", total_pos = "bottom", total_adj = 0.05,
#'   df = df, colour = treat, width = 0.2, height = 0
#' )
#'
#' plotGrlCol(
#'   peaks, geom = geom_boxplot(colour = "grey70"), df = df, fill = treat,
#'   total_pos = "bottom", total_adj = 0.05,
#' ) +
#' scale_y_log10()
#'
#' @importFrom tibble tibble
#' @importFrom S4Vectors mcols
#' @importFrom dplyr bind_rows left_join summarise
#' @importFrom rlang := !! sym enquo
#' @importFrom glue glue
#' @importFrom scales comma
#' @export
plotGrlCol <- function(
        x, var = "width", geom = c("boxplot", "violin", "point", "jitter"),
        .id = "sample", df, fill, colour,
        q = 0.1, q_size = 3.5, qline_type = 2, qline_col = "blue",
        total = "{comma(n)}", total_geom = c("label", "text", "none"),
        total_pos = c("median", "top", "bottom"),
        total_size = 3.5, total_alpha = 1, total_adj = 0.025, ..., digits = 0
) {

    stopifnot(is(x, "GRangesList"))

    ## Create basic data.frame
    mcnames <- .mcolnames(x[[1]])
    var <- ensym(var)
    var <- match.arg(as.character(var), c("width", mcnames))
    var <- sym(var)
    if (var == sym("width")) {
        df_list <- lapply(x, function(x) tibble(width = width(x)))
    } else {
        df_list <- lapply(x, function(gr) as_tibble(mcols(gr)))
        df_list <- lapply(df_list, dplyr::select, !!var)
    }
    plot_df <- bind_rows(df_list, .id = .id)
    if (!missing(fill)) fill <- ensym(fill)
    if (!missing(colour))  colour <- ensym(colour)

    ## Add any annotations
    if (!missing(df)) {
        stopifnot(is(df, "data.frame"))
        ann_cols <- colnames(df)
        if (!.id %in% ann_cols)
            stop(.id, " must be in the supplied annotations")
        plot_df <- left_join(plot_df, df, by = .id)
    }

    ## Summary values
    rng <- range(plot_df[[var]])
    n <- vapply(x, length, integer(1))
    totals_df <- tibble("{.id}" := names(x), label = glue(total))

    ## The plotting geometry
    if (is.character(geom)) {
        geom <- match.arg(geom)
        geom <- match.fun(paste0("geom_", geom))(...)
    }

    ## The basic plot
    p <- ggplot(
        plot_df,
        aes(!!sym(.id), !!sym(var), fill = {{ fill }}, colour = {{ colour }})
    ) +
        geom

    ## Totals
    total_geom <- match.arg(total_geom)
    if (total_geom != "none") {
        total_pos <- match.arg(total_pos)
        yval <- list(
            top = max(rng) + total_adj * diff(rng),
            bottom = min(rng) - total_adj * diff(rng),
            median = summarise(
                plot_df, "{var}" := median(!!sym(var)), .by = !!sym(.id)
            )[[var]]
        )[[total_pos]]
        p <- p + annotate(
            total_geom, x = names(x), y = yval,
            label = totals_df$label, size = total_size, alpha = total_alpha
        )
    }

    ## Percentiles as a horizontal line
    if (all(q > 0)) {
        qval <- quantile(plot_df[[var]], probs = q)
        p <- p + geom_hline(
            yintercept = qval, linetype = qline_type, colour = qline_col
        )
        p <- p +
            annotate(
                "label", x = 0.55, y = qval, label = round(qval, digits),
                colour = qline_col, size = q_size
            )
    }

    p
}

