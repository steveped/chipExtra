#' @title Plot a Genomic Region showing HiC, Features, Genes and Coverage
#'
#' @description Plot a region with showing HiC, Features, Genes and Coverage
#'
#' @details
#' Convenience function for plotting a common set of tracks. All tracks are
#' optional. For more fine control, users are advised to simply use Gviz
#' directly.
#'
#' The primary tracks defined in this function are H (HiC), F (features), G
#' (genes), and C (coverage). Axis and Ideogram tracks are an additional part of
#' this visualisation, with the Ideogram also being optional
#'
#' Use all tracks specific to this dataset to generate a simple visualisation.
#' In descending order the tracks displayed will be:
#' \enumerate{
#'   \item HiC Interactions (if supplied)
#'   \item Regulatory features
#'   \item Genes/genes
#'   \item Coverage tracks as supplied
#' }
#'
#' All tracks are optional and will simply be omitted if no data is supplied.
#' See individual sections below for a more detailed explanation of each track
#'
#' If wanting a single track of genes, simply pass a GRanges object in the
#' format specified for a \link[Gviz]{GeneRegionTrack}. Passing a GRangesList
#' with the same format will yield an individual track for each list element,
#' with each track shown by default as a separate colour. This can be used for
#' showing Up/Down-regulated genes, or Detected/Undetected genes.
#'
#' If passing a BigWigFileList for the coverage track, each file within the
#' object will be drawn on a separate track. If specified, the same y-limits
#' will be applied to each track
#' If passing a list of BigWigFileList
#' objects, each list element will be drawn as a single track with the
#' individual files within each BigWigFileList overlaid within each track.
#'
#' Cytogenetic band information must be in the structure required by
#' \link[Gviz]{IdeogramTrack}, with data for both GRCh37 and GRCh38 provided in
#' this package (\link{grch37.cytobands}, \link{grch38.cytobands}).
#'
#' A highlight overlay over the GRanges provided as the `gr` argument will be
#' added if a colour is provided. If set to NULL, no highlight will be added.
#'
#' @section Displaying HiC Interactions:
#'
#' The available arguments for displaying HiC Interactions are defined below.
#' If `hic` is supplied, a single \link[GenomicInteractions]{InteractionTrack}
#' will be added displaying
#' all interactions with an anchor within the range specified by `gr`.
#' Only interactions with an anchor explicitly overlapping `gr` will be shown.
#' If no interactions are found within `gr`, the track will not be displayed.
#' The **plotting range will expand to incorporate these interactions**, with
#' the paramater `max` providing an upper limit on the displayed range.
#'
#' \describe{
#'   \item{hic}{This is the `GInteractions` object required for inclusion of
#'   a HiC track in the final output. Will be ignored if not supplied}
#'   \item{hiccol}{Determines the colours used for display of anchors and
#'   interactions}
#'   \item{hicsize}{Relative size of the track compared to others}
#'   \item{hicname}{The name to display on the LHS panel}
#'   \item{max}{The maximum width of the plotted region. If multiple long-range
#'   interactions are identified, this provides an upper limit for the display.
#'   This defaults to `10Mb`.}
#' }
#'
#'
#' @section Displaying Features:
#'
#' If wanting to add an \link[Gviz]{AnnotationTrack} with regions defined as
#' 'features', the following arguments are highly relevant.
#' All are ignored if `features` is not provided.
#'
#' \describe{
#'   \item{features}{A named `GRangesList`. Each element will be considered as
#'   a separate feature and drawn as a block in a distinct colour. Any `mcols`
#'   data will be ignored.}
#'   \item{featcol}{A **named** vector (or list) providing a colour for each
#'   element of `features`}
#'   \item{featname}{The name to display on the LHS panel}
#'   \item{featstack}{Stacking to be applied to all supplied features}
#'   \item{featsize}{Relative size of the track compared to others}
#' }
#'
#' @section Displaying Genes And Transcripts:
#'
#' To display genes or transcripts, simply provide a single `GRanges` object if
#' you wish to display all genes on a single track.
#' The `mcols` element of this object should contain the columns `feature`,
#' `gene`, `exon`, `transcript` and `symbol` as seen on the
#' \link[Gviz]{GeneRegionTrack} help page.
#'
#' Alternatively, a `GRangesList` can be provided to display genes on separate
#' tracks based on their category.
#' This can be useful for separating and colouring Up/Down regulated genes in a
#' precise way.
#' All elements should be as described above.
#' Again, all parameters associated with this track-set will be ignored of no
#' object is supplied to this argument.
#'
#' \describe{
#'   \item{genes}{A `GRanges` or `GRangesList` object as described above}
#'   \item{genecol}{A single colour if supplying a `GRanges` object, or a
#'   **named** vector/list of colours matching the `GRangesList`}
#'   \item{genesize}{Relative size of the track compared to others}
#'   \item{collapseTranscripts}{Passed to all tracks. See the GeneRegionTrack
#'   section in \link[Gviz]{settings} for detail regarding possible arguments.
#'   If genes is a `GRangesList`, can be a **named** vector/list with names
#'   matching the names of the `genes` object.
#'   }
#' }
#'
#' @section Displaying Coverage Tracks:
#'
#' This section contains the most flexibility and can take two types of input.
#' The first option is a `BigWigFileList`, which will lead to each BigWig file
#' being plotted on it's own track.
#' An alternative is a list of `BigWigFileList` objects.
#' In this case, each list element will be plotted as a separate track,
#' with all individual `BigWig` files within each list element
#' overlaid within the relevant track.
#'
#' In addition to the coverage tracks, annotations can be added to each
#' `BigWigFileList` in the form of coloured ranges, indicating anything of the
#' users choice. Common usage may be to indicate regions with binding of a
#' ChIP target is found to be detected, unchanged, gained or lost.
#'
#' \describe{
#'   \item{coverage}{A `BigWigFileList` or `list` of `BigWigFileList` objects.
#'   A single `BigWigFileList` will be displayed with each individual file on a
#'   separate track with independent y-axes. Each element of the
#'   `BigWigFileList` **must be named** and these names will be displayed on the
#'   LHS panels
#'   A list of `BigWigFileList` objects will be displayed with each list element
#'   as a separate track, with any `BigWig` files overlaid using the same
#'   y-axis. The list **must be named** with these names displayed on the LHS
#'   panel. Each internal `BigWig` within a `BigWigFileList` must also be named.
#'   }
#'   \item{covtype}{Currently only lines (`covtype = "l"`) and
#'   heatmaps (`covtype = "heatmap"`) are supported. Colours can be
#'   specified using the arguments below}
#'   \item{linecol}{Can be a single colour applied to all tracks, or a *named*
#'   vector (or list) of colours. If `coverage` is a single `BigWigFileList`,
#'   these names should match the names of this object exactly.
#'   If `coverage` is a list of `BigWigFileList` objects, `linecol` should be
#'   a list with matching names. Each element of this list should also be a
#'   **named** vector with names that exactly match those of each corresponding
#'   `BigWigFileList`.}
#'   \item{gradient}{A colour gradient applied to all heatmap tracks. No
#'   specific structure is required beyond a vector of colours.}
#'   \item{covsize}{Relative size of the tracks compared to others}
#'   \item{ylim}{Can be a vector of length 2 applied to all coverage tracks.
#'   Alternatively, if passing a list of `BigWigFlieList` objects to the
#'   `coverage` argument, this can be a **named** list of numeric vectors with
#'   names matching `coverage`}
#'   \item{annotation}{Each `BigWigFileList` needs annotations to be passed to
#'   this argument as a **named** `GRangesList`, with names being used to
#'   associate unique colours with that set of ranges. If `coverage` is a
#'   `BigWigFileList` a simple `GRangesList` would be supplied and a single
#'   'annotation' track will appear at the top of the set of coverage tracks.
#'   If `coverage` is a `list`, then a **named** list of `GRangesList` objects
#'   should be supplied, with each being displayed above the corresponding track
#'   from the `coverage` object.}
#'   \item{annotcol}{A vector of colours corresponding to all names within all
#'   `GRangesList` elements supplied as `annotation`. It is assumed that the
#'   same colour scheme will be applied to all annotation tracks and, as such,
#'   the colours should **not** be provided as a list which matches the
#'   coverage tracks. Instead, every named element anywhere in the annotation
#'   GRanges, across all of the tracks must be included as a colour}
#'   \item{annotsize}{Relative size of the tracks compared to others}
#' }
#'
#'
#' @examples
#' \donttest{
#' library(rtracklayer)
#' ## Make sure we have the cytobands active
#' data(grch37.cytobands)
#'
#' ## Prepare the HiC, promoter & transcript information
#' data(ex_hic, ex_trans, ex_prom)
#' ex_features <- GRangesList(Promoter = ex_prom)
#' featcol <- c(Promoter = "red")
#'
#' ## Prepare the coverage
#' fl <- system.file(
#' "extdata", "bigwig", c("ex1.bw", "ex2.bw"), package = "extraChIPs"
#' )
#' bwfl <- BigWigFileList(fl)
#' names(bwfl)  <- c("ex1", "ex2")
#' bw_col <- c(ex1 = "#4B0055", ex2 = "#007094")
#'
#' ## Define the plotting range
#' gr <- GRanges("chr10:103862000-103900000")
#'
#' ## Now create the basic plot
#' plotHFGC(
#'   gr,
#'   hic = ex_hic, features = ex_features, genes = ex_trans, coverage = bwfl,
#'   featcol = featcol, linecol = bw_col, cytobands = grch37.cytobands
#' )
#'
#' plotHFGC(
#'   gr,
#'   hic = ex_hic, features = ex_features, genes = ex_trans, coverage = bwfl,
#'   featcol = featcol, linecol = bw_col, cytobands = grch37.cytobands,
#'   maxTrans = 1
#' )
#' }
#'
#' @return
#' A Gviz object
#'
#' @param gr The range(s) of interest. Must be on a single chromosome
#' @param hic Any HiC interactions to be included as a GenomicInteractions
#' object. If not supplied, no HiC track will be drawn.
#' @param features A named GRangesList or list of GRangesList objects. Each
#' GRangesList should contain features in each element which will drawn on the
#' same track. If providing a list, each GRangesList within the list will drawn
#' on a separate track. If this argument is not specified, no feature track will
#' be drawn. Features will be drawn with colours provided in  `featcol`.
#' @param genes A GRanges object with exon structure for each transcript/gene.
#' If not included, no track will be drawn for gene/transcript structure.
#' The expected mcols in this object are `type`, `gene`, `exon` `transcript`
#' and `symbol`. See data(ex_trans) for an example.
#' @param coverage A named list of BigWigFileList objects containing the
#' primary tracks to show coverage for. Each list element will be drawn on a
#' separate track, with elements within each BigWigFileList shown on the same
#' track. List names will become track names. Alternatively, a single
#' BigWigFileList will plot all individual files on separate tracks.
#' If not included, no coverage tracks will be drawn.
#' @param annotation Annotations for the coverage track(s).
#' A single GRangesList if coverage is a BigWigListList.
#' If coverage is supplied as a list of BigWigFileLists, a named list of
#' GRangesList objects for each coverage track being annotatated. Names must
#' match those given for coverage.
#' @param zoom Multiplicative factor for zooming in and out
#' @param shift Shift the plot. Applied after zooming
#' @param axistrack logical. Add an AxisTrack()
#' @param cytobands Cytogenetic bands to be displayed on each chromosome.
#' See data('grch37.cytobands') for the correct format. Only drawn if a
#' cytobands data.frame is provided.
#' @param max The maximum width of the plotting region. Given that the width of
#' the final plotting window will be determined by any HiC interactions, this
#' argument excludes any interactions beyond this distance. Plotting can be
#' somewhat slow if any long range interactions are included. Ignored if no HiC
#' interactions are supplied.
#' @param fontsize Applied across all tracks
#' @param hiccol list with names `"anchors"` and `"interactions"`. Colours
#' are passed to these elements
#' @param featcol Named vector (or list) of colours for each feature. Must be
#' provided if drawing features
#' @param genecol Named vector (or list) of colours for each gene category
#' @param annotcol Colours matching the coverage annotations
#' @param covtype The plot type for coverage. Currently only lines ("l")
#' and heatmaps ("heatmap") are supported
#' @param linecol If passing a BigWigFileList to coverage, a vector of colours.
#' If passing a list of BigWigFileList objects to coverage, a list of colours
#' with structure that matches the object being passed to coverage, i.e. a
#' named list of the same length, with elements who's length matches each
#' BigWigFileList. Only used if covtype = "l".
#' @param gradient Colour gradient for heatmaps
#' @param highlight Outline colour for the highlight track. Setting this to
#' `NULL` will remove the highlight
#' @param hicsize,featsize,genesize,covsize,annotsize
#' Relative sizes for each track (hic, features, genes, coverage & annotation)
#' @param hicname,featname Names displayed in the LHS panel
#' @param featstack Stacking for the fature track
#' @param ylim If a numeric vector, this will be passed to all coverage tracks.
#' Alternatively, a named list of y-limits for each coverage track with names
#' that match those in each element of the coverage list.
#' @param ... Passed to \link[Gviz]{DataTrack} for the **coverage tracks** only.
#' Useful arguments may be things like `legend`
#' @param cex.title Passed to all tracks
#' @param rotation.title Passed to all tracks
#' @param col.title Passed to all tracks
#' @param background.title Passed to all tracks
#' @param collapseTranscripts Passed to \link[Gviz]{GeneRegionTrack} for the
#' genes track. Defaults to `"auto"` for automatic setting. If the number of
#' transcripts to be plotted is > `maxtrans`, the argument will be
#' automatically set to `"meta"`, otherwise this will be passed as `FALSE` which
#' will show all transcripts.
#' @param maxTrans Only used if `collapseTranscripts` is set to "auto".
#' @param title.width Expansion factor passed to \link[Gviz]{plotTracks}, and
#' used to widen the panels on the LHS of all tracks.
#' Can have unpredictable effects on the font
#' size of y-axis limits due to the algorithm applied by `plotTracks`
#'
#' @importFrom methods slot
#' @importFrom InteractionSet anchors
#' @importFrom GenomicInteractions calculateDistances
#' @importFrom grDevices hcl.colors
#' @import GenomicRanges
#'
#' @export
plotHFGC <- function(
        gr, hic, features, genes, coverage, annotation,
        zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, cytobands,
        covtype = c("l", "heatmap"),
        linecol = c(), gradient = hcl.colors(101, "viridis"),
        hiccol = list(anchors = "lightblue", interactions = "red"),
        featcol, genecol, annotcol, highlight = "blue",
        hicsize = 1, featsize = 1, genesize = 1, covsize = 4, annotsize = 0.5,
        hicname = "HiC", featname = "Features",
        featstack = c("full", "hide", "dense", "squish", "pack"),
        collapseTranscripts = "auto", maxTrans = 12,
        ylim = NULL, ...,
        fontsize = 12, cex.title = 0.8, rotation.title = 0, col.title = "white",
        background.title = "lightgray",
        title.width = 1.5
) {

    if (!requireNamespace('Gviz', quietly = TRUE))
        stop("Please install 'Gviz' to use this function.")

    ## Argument checks
    ## Add checks for standard chromosomes!!!
    covtype <- match.arg(covtype)
    checkArgs <- .checkHFGCArgs(
        gr = gr, zoom = zoom, shift = shift, hic = hic, features = features,
        genes = genes, coverage = coverage, annotation = annotation,
        axistrack = axistrack, cytobands = cytobands, max = max,
        hiccol = hiccol, linecol = linecol, genecol = genecol,
        featcol = featcol, annotcol = annotcol, type = covtype, ylim = ylim,
        collapseTranscripts = collapseTranscripts, maxTrans = maxTrans
    )
    stopifnot(checkArgs)

    ## Add a step for restricting to standard chromosomes!!!

    ## Form the HiC track, including all interactions beyond the max
    hic_track <- .makeHiCTrack(
        hic, gr, fontsize, hicsize, cex.title, rotation.title, hiccol, hicname,
        col.title, background.title
    )

    ## If interactions were found, reset the plot range. This should be the
    ## maximum of all interactions < max or the initial range
    plot_range <- range(gr)
    if (width(plot_range) > max)
        warning("Provided range is wider than the max permitted. Highlights may be unpredictable")
    if (length(hic_track)) {
        hic <- slot(hic_track, "giobject")
        anchors <- anchors(hic[calculateDistances(hic) < max])
        anchors <- unlist(GRangesList(anchors))
        plot_range <- range(c(anchors, gr), ignore.strand = TRUE)
    }
    ## Now resize/shift as required
    max_width <- min(max, zoom*width(plot_range))
    plot_range <- resize(plot_range, max_width, fix = "center")
    plot_range <- shift(plot_range, shift)

    ## Form the IdeogramTrack
    ideo_track <- .makeIdeoTrack(plot_range, cytobands, fontsize)

    ## Form the features track
    featstack <- match.arg(featstack)
    if (missing(features)) {
        feature_track <- NULL
    } else {
        if (is(features, "GRangesList")) {
            feature_track <- .makeFeatureTrack(
                features, plot_range, fontsize, featcol, featsize, cex.title,
                rotation.title, featname, featstack, col.title, background.title
            )
        }
        if (is(features, "list")) {
            ## Only update using featnames if names are provided for every element
            if (length(featname) != length(features))
                featname <- names(features)
            feature_track <- lapply(
                seq_along(features),
                function(i) {
                    x <- names(features)[[i]]
                    .makeFeatureTrack(
                        features[[x]], plot_range, fontsize, featcol[[x]],
                        featsize, cex.title, rotation.title, featname[[i]],
                        featstack, col.title, background.title
                    )
                }
            )
        }
    }

    ## Form the genes tracks. NB: This will be a list of tracks
    gene_tracks <- .makeGeneTracks(
        genes, plot_range, collapseTranscripts, fontsize, genecol, genesize,
        cex.title, rotation.title, col.title, background.title, maxTrans
    )

    ## The coverage tracks NB: This will be list of tracks
    cov_tracks <- .makeCoverageTracks(
        coverage, plot_range, fontsize, covtype, linecol, gradient, covsize,
        cex.title, rotation.title, ylim, col.title, background.title, ...
    )
    cov_tracks <- .addAnnotations(
        annotation, plot_range, cov_tracks, coverage, annotcol, annotsize,
        col.title, background.title
    )

    ## Add the highlight track if wanted. Include features, genes & coverage
    hl_track <- c(hic_track, feature_track, gene_tracks, cov_tracks)
    hl_track <- hl_track[!vapply(hl_track, is.null, logical(1))]
    if (!is.null(highlight)) {
        gr <- subsetByOverlaps(gr, plot_range)
        if (length(gr)) {
            hl_track <- Gviz::HighlightTrack(
                trackList = hl_track, range = gr, col = highlight,
                fill = "#FFFFFF00", inBackground = FALSE
            )
        } else {
            message("No ranges able to be highlighted")
        }
    }

    plot_list <- list()
    if (!is.null(ideo_track)) plot_list <- list(ideo_track)
    if (axistrack)
        plot_list <- c(
            plot_list, Gviz::GenomeAxisTrack(plot_range, fontsize = fontsize)
        )

    plot_list <- c(plot_list, hl_track)
    Gviz::plotTracks(
        plot_list,
        from = start(plot_range), end(plot_range), title.width = title.width
    )

}


#' @importFrom GenomeInfoDb genome seqnames
#' @importFrom rtracklayer ucscGenomes
.makeIdeoTrack <- function(.gr, .bands, .fontsize) {
    if (missing(.bands)) return(NULL)
    ## Checks have been done in .checkHFGCArgs
    gen <- as.character(unique(genome(.gr)))
    chr <- as.character(unique(seqnames(.gr)))
    Gviz::IdeogramTrack(
        chromosome = chr, genome = gen, name = chr, bands = .bands,
        fontsize = .fontsize
    )
}

#' @importFrom GenomicInteractions anchorOne anchorTwo
#' @importFrom GenomicInteractions GenomicInteractions InteractionTrack
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames
.makeHiCTrack <- function(
        .hic, .gr, .fontsize, .tracksize, .cex, .rot, .col, .name, .col.title,
        .bg.title
) {

    if (missing(.hic)) return(NULL)

    ## Checks have been performed previously on any provided object & params

    ## Just keep cis interactions on the chromosome of interest
    .hic <- subsetByOverlaps(.hic, .gr)
    if (length(.hic) == 0) return(NULL)
    chr <- as.character(seqnames(.gr))[[1]]
    cis <- seqnames(anchorOne(.hic)) == chr & seqnames(anchorTwo(.hic)) == chr
    if (sum(cis) == 0) return(NULL)
    .hic <- .hic[cis]

    track <- InteractionTrack(
        x = GenomicInteractions(.hic),
        chromosome = chr,
        name = .name
    )
    track@dp@pars$size <- .tracksize
    track@dp@pars$fontsize <- .fontsize
    track@dp@pars$fontcolor.title  <- .col.title
    track@dp@pars$background.title  <- .bg.title
    track@dp@pars$cex.title <- .cex
    track@dp@pars$rotation.title <- .rot
    track@dp@pars$col.anchors.line <- .col[["anchors"]]
    track@dp@pars$col.anchors.fill <- .col[["anchors"]]
    track@dp@pars$col.interactions <- .col[["interactions"]]
    track
}

#' @importFrom methods is
#' @importFrom stringr str_trim
.addAnnotations <- function(
        .annotation, .gr, .cov_tracks, .coverage, .fill, .size, .col.title,
        .bg.title
) {
    if (missing(.annotation) | is.null(.cov_tracks)) return(.cov_tracks)
    ## Set everything grey if no colour is specified
    if (missing(.fill)) .fill <- "grey"

    if (is(.coverage, "BigWigFileList")) {
        ## Here we just add another feature track. Input will be a single GRList
        fontsize <- .cov_tracks[[1]]@dp@pars$fontsize
        cex <- .cov_tracks[[1]]@dp@pars$cex.title
        annot_track <- .makeFeatureTrack(
            .annotation, .gr, fontsize, .fill, .size, cex, 0, .name = c(),
            .stacking = "full", .col.title, .bg.title
        )
        annot_track@name <- ""
        tracks <- c(list(annot_track), .cov_tracks)
    }

    if (is(.coverage, "list")) {
        ## Here we just add a feature track before any coverage tracks
        ## Find the common names
        cov2add <- intersect(names(.annotation), names(.coverage))
        if (length(cov2add) == 0) return(.cov_tracks)
        tracks <- lapply(
            .cov_tracks,
            function(x) {
                nm <- str_trim(x@name)
                if (!nm %in% cov2add) return(x)
                fontsize <- x@dp@pars$fontsize
                cex <- x@dp@pars$cex.title
                annot_track <- .makeFeatureTrack(
                    .annotation[[nm]], .gr, fontsize, .fill, .size, cex, 0, "",
                    "full",.col.title, .bg.title
                )
                # annot_track@name <- nm
                c(list(annot_track), x)
            }
        )
    }

    unlist(tracks)

}

#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom stats setNames
#' @importFrom S4Vectors endoapply
.makeFeatureTrack <- function(
        .features, .gr, .fontsize, .fill, .tracksize, .cex, .rot, .name,
        .stacking, .col.title, .bg.title
) {

    if (missing(.features)) return(NULL)

    if (missing(.fill)) {
        pos <- ((seq_along(.features) - 1) %% 8) + 1 # Use modulo for recursion
        .fill <- brewer.pal(8, "Set2")[pos]
        names(.fill) <- names(.features)
    }
    .fill <- unlist(.fill)[names(.features)]

    .features <- endoapply(.features, setNames, NULL) ## Makes the unlist safer
    .features <- granges(unlist(.features))
    .features$feature <- names(.features)
    .features <- subsetByOverlaps(.features, .gr)
    if (length(.features) == 0) return(NULL)
    Gviz::AnnotationTrack(
        ## Change the name later
        range = .features, name = .name, stacking = .stacking,
        col = "transparent", fill = .fill[.features$feature],
        feature = .features$feature,
        ## Tidy setting this up later & also tidy the colour setting
        ## These can likely go in the @dp@pars
        fontsize = .fontsize, cex.title = .cex, col.title = .col.title,
        background.title = .bg.title, rotation.title = .rot,
        size = .tracksize
    )

}

#' @import GenomicRanges
#' @importFrom stringr str_to_title str_pad str_count
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods is
#' @importFrom IRanges subsetByOverlaps
.makeGeneTracks <- function(
        .genes, .gr, .collapse, .fontsize, .col, .tracksize, .cex, .rot,
        .col.title, .bg.title, .max_trans
) {

    if (missing(.genes)) return(list(NULL))
    gene <- c() # Avoid an R CMD check error

    if (is(.genes, "GRanges")) {
        ## If a single GRanges object, just return a single track
        if (missing(.col)) .col <- "#FFD58A"
                ids <- subsetByOverlaps(.genes, .gr)$gene
                if (length(ids) == 0) return(list(NULL))
                .genes <- subset(.genes, gene %in% ids)

                if (.collapse == "auto") {
                    ## Check for the maximum at every genomic position
                    n_trans <- max(max(coverage(.genes)))
                    if (n_trans > .max_trans) {
                        .collapse <- "meta"
                    } else {
                        .collapse <- FALSE
                    }
                }

                trackList <- Gviz::GeneRegionTrack(
                    .genes, name = "Genes", transcriptAnnotation = "symbol",
                    collapseTranscripts = .collapse, size = .tracksize,
                    fontsize = .fontsize, col = "transparent", fill = .col[[1]],
                    cex.title = .cex, rotation.title = .rot,
                    col.title = .col.title, background.title = .bg.title
                )
                return(trackList)
    }

    ## Below will only execute if .genes is a GRangesList
    ids <- lapply(.genes, subsetByOverlaps, .gr)
    ids <- unlist( lapply(ids, function(x) x$gene) )
    if (length(ids) == 0) return(list(NULL))
    .genes <- lapply(.genes, subset, gene %in% ids)
    .genes <- .genes[vapply(.genes, length, integer(1)) > 0]

    ## Make a default set of colours for the gene tracks
    if (missing(.col)) {
        n <- 9L # Set1 contains 9 colours
        pos <- ((seq_along(.genes) - 1) %% n) + 1 # Use modulo for recursion
        .col <- structure(brewer.pal(n, "Set1")[pos], names = names(.genes))
    }
    .col <- .col[names(.genes)]

    ## Setup the collapseTranscripts argument if we just have a vector
    if (length(.collapse) == 1) {
        .collapse <- rep(.collapse[[1]], length(.genes))
        names(.collapse) <- names(.genes)
    }

    trackList <- lapply(
        names(.genes),
        function(x) {
            nm <- x
            ids <- subsetByOverlaps(.genes[[x]], .gr)$gene
            gt <- subset(.genes[[x]], gene %in% ids)
            cl <- .collapse[[x]]
            if (cl == "auto") {
                n_trans <- max(max(coverage(gt)))
                if (n_trans > .max_trans) {
                    cl <- "meta"
                } else {
                    cl <- FALSE
                }
            }
            ## This just helps the weird alignment algorithm
            if (.rot == 0) nm <- str_pad(x, min(str_count(x) + 4, 12))
            Gviz::GeneRegionTrack(
                gt, name = str_to_title(nm),
                transcriptAnnotation = "symbol",
                collapseTranscripts = cl,
                size = .tracksize, fontsize = .fontsize, col = .col[[x]],
                fill = .col[[x]], cex.title = .cex, rotation.title = .rot,
                col.title = .col.title, background.title = .bg.title
            )
        }
    )

    return(trackList)

}

#' @importFrom rtracklayer import.bw
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_count str_pad
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @import GenomicRanges
.makeCoverageTracks <- function(
        .coverage, .gr, .fontsize, .type = c("l", "heatmap"), .linecol,
        .gradient, .tracksize, .cex, .rot, .ylim, .col.title, .bg.title, ...
) {

    ## Always returns a list
    if (missing(.coverage)) return(list(NULL))
    if (is.null(.coverage)) return(list(NULL))
    .type <- match.arg(.type)

    ## Split a single BigWigFileList. This will give a list of BigWigFileLists
    ## (which is a little unexpected). After this step, we will always(!) have a
    ## named S3 list of BigWigFileList objects
    is_single_list <- is(.coverage, "BigWigFileList")
    if (is_single_list) .coverage <- split(.coverage, f = names(.coverage))
    # Make sure all seqinfo objects are compatible with the ranges
    seq_levels <- lapply(.coverage, function(x) lapply(x, seqlevels))
    seq_levels <- unique(unlist(seq_levels))
    gr_levels <- as.character(seqnames(.gr))
    .gr <- .gr[gr_levels %in% seq_levels]
    stopifnot(length(.gr) > 0)
    seq_levels <- intersect(seq_levels, gr_levels)
    .gr <- keepSeqlevels(.gr, seq_levels)
    ## Ensure .linecol is a list of the same structure as .coverage
    .linecol <- .assignColours(.coverage, .linecol, .type)
    ## Ensure the .ylim is a list of the same structure as .coverage
    .ylim <- .assignLimits(.coverage, .ylim)

    tracks <- lapply(
        names(.coverage),
        function(x) {
            ## groups need to be unset if drawing a heatmap
            nm <- names(.coverage[[x]])
            grp <- NULL
            if (.type == "l") grp <- factor(nm, levels = nm)
            ## Now import each file within each list element
            cov <- lapply(nm, function(f) {
                data <- import.bw(.coverage[[x]][[f]], which = .gr)
                colnames(mcols(data)) <- f
                data
            })
            cov <- cov[vapply(cov, length, integer(1)) > 0]
            if (length(cov) == 0) return(NULL)
            nm <- x
            ## This just helps the weird alignment algorithm
            if (.rot == 0) nm <- str_pad(x, min(str_count(x) + 4, 12))
            Gviz::DataTrack(
                range = unlist(GRangesList(cov)),
                name = nm, groups = grp, type = .type,
                fontsize = .fontsize, col = .linecol[[x]], gradient = .gradient,
                showSampleNames = .type == "heatmap", # Always set for heatmaps
                size = .tracksize, cex.title = .cex, rotation.title = .rot,
                .cex.axis = 0.95 * .cex, ylim = .ylim[[x]],
                col.title = .col.title, col.axis = .col.title,
                background.title = .bg.title, ...
            )
        }
    )
    tracks
}

.assignLimits <- function(.coverage, .ylim) {

    lim <- list() # An empty list
    if (is.list(.ylim)) {
        ## If .ylim is already a list, just make sure the names are present
        lim <- lapply(.ylim, range)
        if (is.null(names(.ylim))) names(lim) <- names(.coverage)
    }

    ## If .ylim is a simple vector, return a list with this as each element
    if (is.numeric(.ylim)) lim <- lapply(.coverage, function(x) range(.ylim))

    ## If no limits are provided, return a list which matches .coverage
    ## Each element will be a NULL which is then passed to ylim for auto-limits
    if (is.null(.ylim)) {
        lim <- vector("list", length(.coverage))
        names(lim) <- names(.coverage)
    }
    lim
}

.assignColours <- function(.coverage, .linecol, .type) {

    ## Coverage will always be a list of BigWigFileLists
    ## If the initial object was a single BWFL, then all will internally
    ## have a length of 1. If this is the case, the required structure is a list
    ## of length == length(coverage), with all elements as a single value
    ## Alternatively, if the original object was a list of BWFL objects, we need
    ## a list of colours with the same names as coverage, and with matching
    ## lengths. This gives the structure in essence:
    ## We need a list with the same names as coverage, and elements which have
    ## matching lengths. Even in the case of a heatmap. However, in this case
    ## all elements will be null
    if (.type == "heatmap") {
        cols <- vector("list", length(.coverage))
        names(cols) <- names(.coverage)
        return(cols)
    }
    if (is.null(.linecol)) {
        ## Assign default colours from GViz using a recursion strategy
        defCols <- c(
            "#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00",
            "brown"
        )
        cols <- lapply(
            .coverage,
            function(x) {
                n <- length(x)
                # Use modulo for recursion if n > 6
                ind <- ((seq_len(n) - 1) %% 6) + 1
                defCols[ind]
            })
        return(cols)
    }
    ## If the function has progressed this far, we are using the supplied colours
    ## If the initial input was a BWFL, a vector will have been provided. Split
    ## this into a named list. If a named list was provided, return the same
    all_single <- all(vapply(.coverage, length, integer(1)) == 1)
    if (all_single) {
        cols <- unlist(.linecol) # In case we have a list
        ## Names will have been checked if they are present
        if (is.null(names(cols))) names(cols) <- names(.coverage)
        return(as.list(cols))
    }
    ## If the coverage was a list of BWFL objects, .linecol must have been a
    ## list which matches the structure exactly. This will have been checked in
    ## the first steps of the function and can be returned as is
    .linecol

}



#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
.checkHFGCArgs <- function(
        gr, zoom, shift, hic, features, genes, coverage, annotation, axistrack,
        cytobands, max, hiccol, linecol, genecol, featcol, annotcol, type, ylim,
        collapseTranscripts, maxTrans
) {

    msg <- c()
    if (!is(gr, "GRanges")) {
        msg <- c(msg, "'gr' must be a GRanges object.\n")
    } else {
        if (length(gr) == 0) msg <- c(msg, "Cannot be an empty range\n")
        seq_names <- as.character(unique(seqnames(gr)))
        if (length(seq_names) > 1)
            msg <- c(msg, "All ranges must be on the same chromosome\n")
    }

    nums <- as.numeric(c(zoom, shift, max))
    if (any(is.na(nums)))
        msg <- c(
            msg, "zoom/shift/max must be numeric or coercible to numeric\n"
        )

    if (!missing(hic)) {
        if (!is(hic, "GInteractions"))
            msg <- c(msg, "'hic' must be provided as a GInteractions object\n")
        colargs <- c("anchors", "interactions")
        if (!is(hiccol, "list") | !all(colargs %in% names(hiccol)))
            msg <- c(
                msg, paste("'hiccols' must be a list with name:", colargs, "\n")
            )
    }

    msg <- .checkFeatures(msg, features, featcol)
    msg <- .checkGenes(msg, genes, genecol, collapseTranscripts, maxTrans)
    msg <- .checkCoverage(
        msg, coverage, linecol, type, annotation, annotcol, ylim
    )
    bools <- c(axistrack)
    if (!is.logical(bools)) msg <- c(msg, "axistrack must be logical values")

    if (!missing(cytobands)) {
        band_cols <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
        if (!is(cytobands, "data.frame"))
            msg <- c(msg, "cytobands must be supplied as a data.frame\n")
        if (!all(band_cols %in% colnames(cytobands)))
            msg <- c(
                msg,
                paste(
                    "Columns in cytobands must be exactly:",
                    paste(band_cols, collapse = ", "), "\n"
                )
            )
        if (!seq_names %in% cytobands$chrom)
            msg <- c(msg, "seqnames not found in cytobands\n")
    }

    if (is.null(msg)) return(TRUE)
    message(msg)
    FALSE
}

.checkCoverage <- function(
        msg, coverage, linecol, type, annotation, annotcol, ylim
) {

    if (missing(coverage)) return(msg)

    if (!any(is(coverage, "list"), is(coverage, "BigWigFileList")))
        return(c(msg, "'coverage' should be a BigWigFileList or a list\n"))

    ## Add checks for coverage. This can be a list of BigWigFileLists, or a
    ## single BigWigFileList
    if (is(coverage, "list")) {
        chk <- vapply(coverage, is, logical(1), class2 = "BigWigFileList")
        if (!all(chk))
            msg <- c(
                msg, "All elements of 'coverage' must be a BigWigFileList\n"
            )
        if ("" %in% names(coverage))
            msg <- c(msg, "'coverage' must a be a named list\n")
        nm <- lapply(coverage, names)
        if (any(vapply(nm, is.null, logical(1))))
            msg <- c(msg, "All individual bigwig files must be named\n")
    }

    if (!is.null(linecol) & type == "l") {
        ## If we have BigWigFileList, we need a vector or list of colours
        ## the same length. If named, names must match
        if (is(coverage, "BigWigFileList")) {
            linecol <- unlist(linecol)
            covnames <- names(coverage)
            if (length(linecol) != length(coverage))
                msg <- c(
                    msg,
                    "linecol must be the same length as the coverage tracks\n"
                )

            if (!is.null(names(linecol)) & !all(covnames %in% names(linecol)))
                msg <- c(
                    msg, "All elements of coverage must be named in linecol\n"
                )
        }
        ## If passing a list of BigWigFileList objects colours must be a list of
        ## the same length. Each element must also be NULL, or match
        if (is(coverage, "list")) {
            nm <- names(coverage)
            if (!is.list(linecol)) {
                msg <- c(msg, "linecol must be a named list\n")
            } else {
                if (!all(nm %in% names(linecol))) {
                    msg <- c(
                        msg,
                        "All elements of coverage must be named in linecol\n"
                    )
                } else {
                    cov_len <- vapply(coverage, length, integer(1))
                    col_len <- vapply(linecol, length, integer(1))
                    cov_len <- cov_len[names(which(col_len > 0))]
                    if (!all(cov_len == col_len[names(cov_len)]))
                        msg <- c(
                            msg,
                            paste(
                                "All elements of linecol must be NULL,",
                                "or match coverage\n"
                            )
                        )
                }
                match_names <- mapply(
                    function(x, y) all(names(x) %in% names(y)),
                    x = coverage, y = linecol
                )
                if (!all(match_names)) msg <- c(
                    msg,
                    "All elements within coverage must have a matching colour\n"
                )
            }
        }
    }

    if (!is.null(ylim)) {
        if (is(coverage, "list")) {
            # If coverage is a list, ylim can be a vector passed to all elements
            # or a named list of numerics
            if (!is(ylim, "list")) {
                if (length(ylim) < 2 | !is.numeric(ylim))
                    msg <- c(
                        msg,
                        paste(
                            "ylim can be a named list or numeric vector of",
                            "length >= 2\n"
                        )
                    )
            } else {
                if (!all(names(coverage) %in% names(ylim)))
                    msg <- c(
                        msg,
                        "All elements of coverage must also be named in ylim\n"
                    )
                if (!all(vapply(ylim, length, integer(1)) >= 2))
                    msg <- c(
                        msg,
                        "All supplied elements of ylim must be of length >= 2\n"
                    )
            }
        }
        if (is(coverage, "BigWigFileList")) {
            if (length(ylim) < 2 | !is.numeric(ylim))
                msg <- c(
                    msg,
                    "ylim should be passed as a numeric vector of length >= 2"
                )
        }
    }

    if (missing(annotation)) return(msg)

    all_annot <- c()
    ## We need to check that it's a GRangesList if coverage is a BigWigFileList
    ## This would then be plotted as a single track above all coverage tracks
    if (is(coverage, "BigWigFileList")) {
        if (!is(annotation, "GRangesList")) {
            msg <- c(msg, "annotation must be a GRangesList\n")
        } else {
            all_annot <- names(annotation)
        }
    }

    ## Or a list of GRangesLists if coverage is a list of BigWigFileLists
    if (is(coverage, "list")) {
        if (!is(annotation, "list")) {
            msg <- c(msg, "annotation must be a list of GRangesList objects\n")
        } else {
            is_grl <- vapply(annotation, is, logical(1), class2 = "GRangesList")
            if (!all(is_grl)) {
                msg <- c(
                    msg, "annotation must be a list of GRangesList objects\n"
                )
            } else {
                all_annot <- unlist(lapply(annotation, names))
            }
        }
    }

    if (missing(annotcol)) return(msg)

    ## Also we need to check that colours match if specified. For a single BWFL
    ## we need a simple vector of colours. For a list of BWFLs we need the same!
    if (!missing(annotcol)) {
        all_cols <- names(annotcol)
        miss <- setdiff(all_annot, all_cols)
        if (length(miss) > 0)
            msg <- c(msg, "Colours not specified for ", miss, "\n")
    }

    msg

}


.checkGenes <- function(msg, genes, genecol, collapseTranscripts, maxTrans) {

    if (missing(genes)) return(msg)

    if (!is(genes, "GenomicRanges_OR_GRangesList"))
        msg <- c(msg, "genes must be a 'GRanges' or 'GRangesList' object\n")

    reqd_mcols <- c("gene", "exon", "transcript", "symbol")
    trans_vals <- c(
        "TRUE", "T", "FALSE", "F", "gene", "longest", "shortest", "meta", "auto"
    )

    if (is(genes, "GRangesList")) {
        nm <- names(genes)
        if (length(nm) != length(genes)) {
            msg <- c(
                msg, "All elements of the 'genes' GRangesList must be named \n"
            )
            return(msg)
        }

        chk_mcols <- vapply(
            genes,
            function(x) all(reqd_mcols %in% .mcolnames(x)),
            logical(1)
        )
        if (!all(chk_mcols))
            msg <- c(
                msg,
                paste(
                    "The 'genes' GRangesList must have the columns",
                    collapseGenes(
                        c("gene", "exon", "transcript", "symbol"), format = ""
                    ),
                    "\n"
                )
            )

        if (!missing(genecol)) {

            chkColNames <- (
                all(names(genes) %in% names(genecol)) | length(genecol) == 1
            )
            if (!chkColNames) msg <- c(
                msg,
                paste(
                    "All elements of the 'genes' GRangesList should have a",
                    "named colour in 'genecol'\n"
                )
            )

        }

        ## collapseTranscripts can be a logical vector or list
        ## An logical vector is fine & only alternatives need to be checked
        if (!is.logical(collapseTranscripts)) {

            if (length(collapseTranscripts) > 1) {
                ## Check the names
                if (!all(nm %in% names(collapseTranscripts)))
                    msg <- c(
                        msg,
                        paste(
                            "All elements of the 'genes' GRangesList must be",
                            "named in collapseTranscripts\n"
                        )
                    )
            }

            ## Check the values which aren't logical
            is_log <- vapply(collapseTranscripts, is.logical, logical(1))
            ct <- unlist(collapseTranscripts[!is_log])
            if (!all(ct %in% trans_vals)) {
                msg <- c(
                    msg,
                    paste(
                        "collapseTranscripts can only be logical or one of",
                        "gene, longest, shortest or meta\n"
                    )
                )
            }
            if (any(ct == "auto") & !is.numeric(maxTrans)) {
                msg <- c(msg, "'maxTrans' must be numeric")
            }
        }

    }

    if (is(genes, "GRanges")) {
        chk_cols <- all(reqd_mcols %in% .mcolnames(genes))
        if (!chk_cols)
            msg <- c(
                msg,
                paste(
                    "The 'genes' GRanges object must have the columns",
                    collapseGenes(
                        c("gene", "exon", "transcript", "symbol"), format = ""
                    ),
                    "\n"
                )
            )
        ## Check the values
        if (!is.logical(collapseTranscripts)) {
            if (!all(collapseTranscripts %in% trans_vals))
                msg <- c(
                    msg,
                    paste(
                        "collapseTranscripts can only be logical or one of",
                        "gene, longest, shortest or meta\n"
                    )
                )
        }
    }

    msg

}

.checkFeatures <- function(msg, features, featcol) {

    if (missing(features)) return(msg)

    if (!any(is(features, "GRangesList"), is(features, "list"))) {
        msg <- c(msg, "'features' must be a GRangesList or list")
        return(msg)
    }

    ## These two checks apply if features is either a GRangesList or list
    if ("" %in% names(features) | length(names(features)) != length(features))
        msg <- c(
            msg, "All elements of 'features' must be explicitly named\n"
        )

    if (!missing(featcol)) {
        if (!all(names(features) %in% names(featcol)))
            msg <- c(
                msg, "All elements of 'features' must be in 'featcol'\n"
            )
    }

    if (is(features, "list")) {
        all_grl <- all(vapply(features, is, logical(1), class2 = "GRangesList"))
        if (!all_grl)
            msg <- c(msg, "All elements of features must be a GRangesList")

        # All elements must be a named GRangesList
        all_named <- all(
            c(
                vapply(features, function(x) !"" %in% names(x), logical(1)),
                vapply(
                    features,
                    function(x) length(names(x)) == length(x),
                    logical(1)
                )
            )
        )
        if (!all_named)
            msg <- c(msg, "All GRangesList elements of features must be named")

        # featcol must be a named list with identical names in each element
        named_cols <- mapply(
            function(x, y) all(names(x) %in% names(y)),
            x = features, y = featcol
        )
        if (!all(named_cols))
            msg <- c(msg, "All features must have a named colour")
    }

    msg

}
