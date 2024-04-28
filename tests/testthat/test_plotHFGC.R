gr1 <- GRanges("chr1:1-10")
gr2 <- GRanges("chr1:21-30")
hic <- InteractionSet::GInteractions(gr1, gr2)
hiccol <- list(anchors = "lightblue", interactions = "red")
feat <- GRangesList(a = GRanges("chr1:8"))
genes <- gr1
mcols(genes) <- DataFrame(
    feature = "exon", gene = "ENSG", exon = "ENSE", transcript = "ENST",
    symbol = "ID"
)
cyto_df <- data.frame(
    chrom = "chr1", chromStart = 1, chromEnd = 100, name = "p1", gieStain = "gneg"
)

test_path <- system.file("tests", "test.bw", package = "rtracklayer")
test_bw <- rtracklayer::BigWigFileList(test_path)
names(test_bw) <- "a"

test_that(".checkHFGCArgs catches GRanges issues", {
    expect_error(suppressMessages(plotHFGC(gr = NULL)))
    expect_error(
        .checkHFGCArgs(zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l")
    )
    expect_message(
        .checkHFGCArgs(
            gr = "", zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l"
        ),
        "must be a GRanges object"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges(), zoom = 1, shift = 0, max = 1e7, axistrack = TRUE,
            type = "l"
        ),
        "Cannot be an empty range"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges(c("chr1:1", "chr2:1")), zoom = 1, shift = 0, max = 1e7,
            axistrack = TRUE, type = "l"
        ),
        "All ranges must be on the same chromosome"
    )
    expect_message(
        .checkHFGCArgs(
            gr = gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = c()
        ),
        "genes must be a 'GRanges' or 'GRangesList' object"
    )
    expect_true(
        .checkHFGCArgs(
            gr = gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l"
        )
    )

})

test_that("Malformed numerics/bools are caught", {
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = "", shift = 0, max = 1e7, axistrack = TRUE, type = "l"
        ),
        "zoom/shift/max must be numeric or coercible to numeric"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = "", type = "l"
        ),
        "axistrack must be logical"
    )
})

test_that("Malformed HiC Interactions are caught", {
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            hic = hic, hiccol = hiccol
        )
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            hic = c(), hiccol = hiccol
        ),
        "'hic' must be provided as a GInteractions object"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            hic = hic, hiccol = hiccol[1]
        ),
        "'hiccols' must be a list with name: interactions"
    )

})

test_that("Malformed features are caught", {
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            features = GRangesList()
        )
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            features = c()
        ),
        "'features' must be a GRangesList"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            features = GRangesList(a = GRanges(), GRanges())
        ),
        "All elements of 'features' must be explicitly named"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            features = GRangesList(a = GRanges()), featcol = c(b = "blue")
        ),
        "All elements of 'features' must be in 'featcol'"
    )
})

test_that("Malformed genes are caught", {
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = genes, collapseTranscripts = TRUE
        )
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = GRanges(), collapseTranscripts = TRUE
        ),
        "The 'genes' GRanges object must have the columns"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = GRangesList(a = gr1), genecol = "blue",
            collapseTranscripts = TRUE
        ),
        "The 'genes' GRangesList must have the columns"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = GRangesList(a = genes, b = genes), genecol = "blue",
            collapseTranscripts = list(a = TRUE, c = TRUE), maxTrans = 12
        ),
        "All elements of the 'genes' GRangesList must be named in collapseTranscr"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            genes = genes, genecol = "blue", collapseTranscripts = ""
        ),
        "collapseTranscripts can only be logical or one of gene, longest, shortest"
    )
})

test_that("Malformed coverage parameters are caught", {
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            coverage = test_bw, linecol = c(a = "blue"), ylim = c(0, 1)
        )
    )
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            coverage = list(a = test_bw), linecol = list(a = c(a = "blue")),
            ylim = list(a = c(0, 1))
        )
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            coverage = list(a = c()), ylim = list(a = c(0, 1)), linecol = c()
        ),
        "All elements of 'coverage' must be a BigWigFileList"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            coverage = c(), ylim = c(0, 1), linecol = c()
        ),
        "'coverage' should be a BigWigFileList or a list"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            coverage = list(a = test_bw), linecol =  "blue", ylim = list(a = c(0, 1))
        ),
        "linecol must be a named list"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = test_bw, annotation = GRanges(),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = c(0, 1), linecol = c()
        ),
        "annotation must be a GRangesList"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = list(a = test_bw),
            annotation = GRangesList(),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = list(a = c(0, 1)), linecol = c()
        ),
        "annotation must be a list of GRangesList objects"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = list(a = test_bw),
            annotation = list(NULL),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = list(a = c(0, 1)), linecol = c()
        ),
        "annotation must be a list of GRangesList objects"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = test_bw,
            annotation = GRangesList(a = GRanges()), annotcol = c(b = "blue"),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = c(0, 1), linecol = c()
        ),
        "Colours not specified for a"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = list(a = test_bw, test_bw),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = list(a = c(0, 1)), linecol = c()
        ),
        "'coverage' must a be a named list"
    )
    expect_message(
        .checkHFGCArgs(
            gr = GRanges("chr2:1-1000"), coverage = list(a = test_bw),
            zoom = 1, shift = 0, max = Inf, type = "l", axistrack = TRUE,
            ylim = 1, linecol = c()
        ),
        "ylim can be a named list or numeric vector of length >= 2"
    )
})

test_that(".checkCoverage checks correctly", {
    expect_equal(
        .checkCoverage(
            msg = c(),
            coverage = list(a = setNames(test_bw, c())),
            linecol = list(a = c(a = "blue")), type = "l",
            annotation = list(a = GRangesList()), annotcol = NULL,
            ylim = c(0, 1)
        ),
        "All individual bigwig files must be named\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = test_bw, linecol = rep("blue", 2), type = "l",
            annotation = GRangesList(), annotcol = NULL, ylim = c(0, 1)
        ),
        "linecol must be the same length as the coverage tracks\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = test_bw, linecol = setNames("blue", "b"),
            type = "l", annotation = GRangesList(), annotcol = NULL,
            ylim = c(0, 1)
        ),
        "All elements of coverage must be named in linecol\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = list(a = test_bw),
            linecol = setNames("blue", "b"), type = "l",
            annotation = list(a = GRangesList()), annotcol = NULL,
            ylim = c(0, 1)
        ),
        "linecol must be a named list\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = list(a = test_bw),
            linecol = list("b" = "blue"), type = "l",
            annotation = list(a = GRangesList()), annotcol = NULL,
            ylim = c(0, 1)
        ),
        c(
            "All elements of coverage must be named in linecol\n",
            "All elements within coverage must have a matching colour\n"
        )
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = list(a = test_bw),
            linecol = list(a = c(a = "blue")), type = "l",
            annotation = list(a = GRangesList()), annotcol = NULL, ylim = 0
        ),
        "ylim can be a named list or numeric vector of length >= 2\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = list(a = test_bw),
            linecol = list(a = c(a = "blue")), type = "l",
            annotation = list(a = GRangesList()), annotcol = NULL,
            ylim = list(a = 0)
        ),
        "All supplied elements of ylim must be of length >= 2\n"
    )

    expect_equal(
        .checkCoverage(
            msg = c(), coverage = test_bw,
            linecol = c(a = "blue"), type = "l",
            annotation = GRangesList(), annotcol = NULL,
            ylim = 0
        ),
        "ylim should be passed as a numeric vector of length >= 2"
    )
})

test_that(".checkGenes performs correct checks", {

    expect_equal(
        .checkGenes(
            c(), GRangesList(genes), genecol = NULL, collapseTranscripts = FALSE
        ),
        "All elements of the 'genes' GRangesList must be named \n"
    )

    expect_equal(
        .checkGenes(
            c(), GRangesList(up = genes), genecol = c(a = "blue", b = "red"),
            collapseTranscripts = FALSE
        ),
        "All elements of the 'genes' GRangesList should have a named colour in 'genecol'\n"
    )

})

test_that("Malformed cytobands are caught", {
    expect_true(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            cytobands = cyto_df
        )
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            cytobands = cyto_df[1]
        ),
        "Columns in cytobands must be exactly"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            cytobands = cyto_df[c(),]
        ),
        "seqnames not found in cytobands"
    )
    expect_message(
        .checkHFGCArgs(
            gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
            cytobands = c()
        ),
        "cytobands must be supplied as a data.frame"
    )
})

test_that("Ideogram forms correctly", {
    ideo <- .makeIdeoTrack(gr1, cyto_df, 12)
    expect_true(is(ideo, "IdeogramTrack"))
    expect_error(
        .makeIdeoTrack(gr1, "", .fontsize = 12),
        "'bands' must be an object of class 'data.frame'"
    )
})

test_that("HiC Track forms correctly", {
    hic_track <- .makeHiCTrack(
        hic, gr1, .tracksize = 2, .fontsize = 10, .cex = 1, .rot = 90,
        .col = hiccol, .name = "HiC", .col.title = "black", .bg.title = "white"
    )
    expect_true(is(hic_track, "InteractionTrack"))
    expect_equal("HiC", hic_track@name)
    expect_equal(hic_track@dp@pars$col.anchors.line, hiccol$anchors)
    expect_equal(hic_track@dp@pars$col.anchors.fill, hiccol$anchors)
    expect_equal(hic_track@dp@pars$col.interactions, hiccol$interactions)
    expect_equal(hic_track@dp@pars$fontsize, 10)
    expect_equal(hic_track@dp@pars$fontcolor.title, "black")
    expect_equal(hic_track@dp@pars$background.title, "white")
    expect_equal(hic_track@dp@pars$size, 2)
    expect_equal(hic_track@dp@pars$rotation.title, 90)
    expect_null(.makeHiCTrack(hic, GRanges()))
    expect_null(.makeHiCTrack())
})

test_that("Feature Track forms correctly", {
    feat_track <- .makeFeatureTrack(
        feat, gr1, .fontsize = 12, list(a = "blue"), .cex = 1, .tracksize = 1,
        .rot = 0, .name = "Features", .stacking = "full", .col.title = "black",
        .bg.title = "white"
    )
    expect_true(is(feat_track, "AnnotationTrack"))
    expect_equal(feat_track@dp@pars$fill, c(a = "blue"))
    expect_equal(feat_track@dp@pars$background.title, "white")
    expect_equal(feat_track@dp@pars$fontcolor.title, "black")
    expect_equal("Features", feat_track@name)
    expect_null(
        .makeFeatureTrack(
            feat, gr2, .fontsize = 12, list(a = "blue"), .cex = 1, .tracksize = 1,
            .rot = 0, .name = "Features", .stacking = "full"
        )
    )
    expect_null(.makeFeatureTrack())
    ## The automatic colour assignment
    feat_track <- .makeFeatureTrack(
        feat, gr1, .fontsize = 12, .cex = 1, .tracksize = 1,
        .rot = 0, .name = "Features", .stacking = "full", .col.title = "black",
        .bg.title = "white"
    )
    expect_equal(
        feat_track@dp@pars$fill, setNames(brewer.pal(8, "Set2")[[1]], "a")
    )
})

test_that("Genes Track forms correctly", {
    expect_equal(.makeGeneTracks(), list(NULL))
    gene_track <- .makeGeneTracks(
        genes, gr1, "meta", .tracksize = 1, .cex = 1, .rot = 0, .fontsize = 12,
        .col.title = "black", .bg.title = "white"
    )
    expect_true(is(gene_track, "GeneRegionTrack"))
    expect_equal(length(gene_track), 1)
    expect_equal(gene_track@dp@pars$fill, "#FFD58A")
    expect_equal(gene_track@dp@pars$fontsize, 12)
    expect_equal(gene_track@dp@pars$fontcolor.title, "black")
    expect_equal(gene_track@dp@pars$background.title, "white")
    expect_equal(gene_track@dp@pars$cex, 1)
    expect_equal(gene_track@dp@pars$collapseTranscripts, "meta")
    expect_equal(gene_track@dp@pars$rotation, 0)
    expect_equal(
        .makeGeneTracks(genes, gr2, "meta", .tracksize = 1, .cex = 1, .rot = 0),
        list(NULL)
    )

    gene_track <- .makeGeneTracks(
        GRangesList(a = genes),
        gr1, "meta", .tracksize = 1, .cex = 1, .rot = 0, .fontsize = 12,
        .col.title = "black", .bg.title = "white"
    )
    expect_true(is(gene_track, "list"))
    expect_true(is(gene_track[[1]], "GeneRegionTrack"))
    expect_equal(names(gene_track[[1]]), str_pad("A", 5))
    expect_equal(gene_track[[1]]@dp@pars$fill, "#E41A1C")

})

test_that("Coverage Track forms correctly", {
    cov_track <- .makeCoverageTracks(
        .coverage = test_bw, .gr = GRanges("chr2:1-10"),
        .fontsize = 12, .type = "l", .gradient = "blue", .tracksize = 1, .cex = 1,
        .rot = 0, .linecol = c(a = "red"), .ylim = c(0, 1),
        .col.title = "black", .bg.title = "white"
    )
    expect_true(is(cov_track, "list"))
    expect_true(is(cov_track[[1]], "DataTrack"))
    expect_equal(length(cov_track), 1)
    expect_equal(cov_track[[1]]@dp@pars$gradient, "blue")
    expect_equal(cov_track[[1]]@dp@pars$size, 1)
    expect_equal(cov_track[[1]]@dp@pars$rotation, 0)
    expect_equal(cov_track[[1]]@dp@pars$col, "red")
    expect_equal(cov_track[[1]]@dp@pars$fontcolor.title, "black")
    expect_equal(cov_track[[1]]@dp@pars$background.title, "white")
    expect_equal(levels(cov_track[[1]]@dp@pars$groups), "a")
    expect_error(
        suppressWarnings(
            .makeCoverageTracks(
                .coverage = test_bw, .gr = gr1, .linecol = c(), .ylim = c(),
                .col.title = "black", .bg.title = "white"
            )
        ),
        "length\\(.gr\\) > 0 is not TRUE"
    )
    expect_equal(.makeCoverageTracks(), list(NULL))
    expect_equal(.makeCoverageTracks(.coverage = NULL), list(NULL))
    cov_heat <- .makeCoverageTracks(
        .coverage = test_bw, .gr = GRanges("chr2:1-10"), .type = "heatmap",
        .fontsize = 12, .gradient = "blue", .tracksize = 1, .cex = 1, .rot = 0,
        .linecol = c(), .ylim = c(), .col.title = "black", .bg.title = "white"
    )
    expect_null(cov_heat[[1]]@dp@pars$col)

    ann_tr <- .addAnnotations(
        .annotation = list(a = GRangesList(Up = "chr2:5")),
        .gr = GRanges("chr2:1-10"), .cov_tracks = cov_track,
        .coverage = list(a = test_bw), .col.title = "black",
        .bg.title = "white", .size = 1
    )
    expect_equal(length(ann_tr), 2L)
    expect_equal(
        vapply(ann_tr, function(x) is(x)[[1]], character(1)),
        c("AnnotationTrack", "DataTrack")
    )

})

test_that("ylim is parsed correctly", {
    expect_equal(.assignLimits(test_bw, c()), list(a = c()))
    expect_equal(.assignLimits(test_bw, c(0, 1)), list(a = c(0, 1)))
    expect_equal(
        .assignLimits(list(test = test_bw), c(0, 1)),
        list(test = c(0, 1))
    )
    expect_equal(
        .assignLimits(
            list(test1 = test_bw, test2 = test_bw),
            list(test1 = c(0, 1), test2 = c(1,5))
        ),
        list(
            test1 = c(0, 1), test2 = c(1, 5)
        )
    )
    expect_equal(.assignLimits(test_bw, list(c(0, 1))), list(a = c(0, 1)))
})

test_that("linecol is parsed correctly",{
    expect_equal(.assignColours(test_bw, c(), "heatmap"), list(a = NULL))
    expect_equal(.assignColours(test_bw, c(), "l"), list(a = "#0080ff"))
    expect_equal(
        .assignColours(list(test = c(test_bw, test_bw)), c(), "l"),
        list(test = c("#0080ff", "#ff00ff"))
    )
    expect_equal(.assignColours(test_bw, "red", "l"), list(a = "red"))
    expect_equal(
        .assignColours(list(a = test_bw, b = test_bw), c("red", "blue"), "l"),
        list(a = "red", b = "blue")
    )
    expect_equal(
        .assignColours(
            setNames(c(test_bw, test_bw), c("a", "b")), c("red", "blue"), "l"
        ),
        list(a = "red", b = "blue")
    )
})

test_that(".addAnnotations behaves as expected", {
    expect_null(.addAnnotations(.gr = gr1, .cov_tracks = NULL))
})

test_that("Poor choices for max are caught", {
    expect_warning(
        plotHFGC(gr1, cytobands = cyto_df, max = 5),
        "Provided range is wider.+"
    )
})
