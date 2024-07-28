library(rtracklayer)
library(Rsamtools)
bw <- BigWigFile(system.file("extdata/bigwig/ex1.bw", package = "extraChIPs"))
peaks <- importPeaks(
    system.file("extdata/peaks.bed.gz", package = "extraChIPs"), type = "bed"
)[[1]]
bf <- BamFile(system.file("extdata/bam/ex1.bam", package = "extraChIPs"))

test_that("All bw methods work as expected", {
    gr <- centrePeaks(peaks[1], bw)
    expect_true(is(gr, "GRanges"))
    expect_true(length(gr) == 1)

    gr <- centrePeaks(peaks[1], bw, f = "mean")
    expect_true(is(gr, "GRanges"))
    expect_true(length(gr) == 1)

    gr <- centrePeaks(peaks[1], bw, f = "med")
    expect_true(is(gr, "GRanges"))
    expect_true(length(gr) == 1)
})

test_that("All bam methods work as expected", {
    gr <- centrePeaks(peaks[1], path(bf))
    expect_true(is(gr, "GRanges"))
    expect_true(length(gr) == 1)
})

test_that("Character methods fail on mixed files", {
    f <- c(path(bw), path(bf))
    expect_error(centrePeaks(gr, f))
})
