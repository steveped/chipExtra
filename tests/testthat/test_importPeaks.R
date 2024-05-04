test_that("narrowPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )

  grl <- importPeaks(fl, type = "narrow")
  gr <- unlist(grl)

  expect_equal(length(grl), 1L)
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 5L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))

  expect_true(!"centre" %in% .mcolnames(gr))
  gr <- unlist(importPeaks(fl, type = "narrow", nameRanges = FALSE, centre = TRUE))
  expect_true(all(c("centre", "name") %in% .mcolnames(gr)))


})

test_that("broadPeak files parse correctly", {
  fl <- system.file(
    "extdata/testFiles", "test.broadPeak", package = "extraChIPs"
  )

  grl <- importPeaks(fl, type = "broad")
  gr <- unlist(grl)

  expect_equal(length(grl), 1L)
  expect_equal(length(gr), 2L)
  expect_equal(width(gr), c(171, 343)) # Checks for 0-based ranges
  expect_equal(ncol(mcols(gr)), 4L)
  expect_true(all(vapply(mcols(gr), is.numeric, logical(1))))
})

test_that("bed files parse correctly", {

  fl <- system.file("extdata", "peaks.bed.gz", package = "extraChIPs")
  grl <- importPeaks(fl, type = "bed")
  gr <- unlist(grl)

  expect_equal(length(grl), 1L)
  expect_equal(length(gr), 5L)
  expect_equal(.mcolnames(gr), "score")

  grl <- importPeaks(fl, type = "bed", nameRanges = FALSE)
  gr <- unlist(grl)
  expect_equal(.mcolnames(gr), c("name", "score"))


})

test_that("seqinfo objects behave correctly for narrowPeak files", {
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )
  ## Succeed
  expect_type(
    importPeaks(
      fl, type = "narrow",
      seqinfo = Seqinfo(seqnames = "chr1", seqlengths = 8100000)
    ),
    "S4"
  )
  ## Empty GRanges
  expect_message(
    importPeaks(
      fl, type = "narrow", seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "coarse"
    ),
    "No ranges match the supplied seqinfo object"
  )
  ## Error
  expect_error(
    importPeaks(
      fl, type = "narrow",
      seqinfo = Seqinfo(seqnames = "chr2"),
      pruning.mode = "error"
    )
  )
  ## Error
  expect_error(importPeaks(fl, "narrow", seqinfo = NULL))
  ## Error
  expect_error(
    importPeaks(fl, "narrow", glueNames = "{1:5}"), "length.+is not TRUE"
  )

})

test_that("blacklists behave correctly",{
  fl <- system.file(
    "extdata/testFiles", "test.narrowPeak", package = "extraChIPs"
  )
  expect_error(
    importPeaks(gr, "narrow", blacklist = NULL)
  )
  expect_equal(
    length(importPeaks(fl, blacklist = GRanges("chr1:1299700"))),
    1L
  )
})


test_that("Empty files parse empty GRanges", {
  ## This now catches empty files that are compressed too!
  fl <- system.file("extdata/empty.bed.gz", package = "extraChIPs")
  np1 <- importPeaks(fl, type = "narrow", setNames = FALSE)
  bp1 <- importPeaks(fl, type = "broad", setNames = FALSE)
  expect_true(is(c(np1, bp1), "GRangesList"))
  expect_equal(vapply(c(np1, bp1), length, integer(1)), c(0, 0))

  # And with a seqinfo
  sq <- Seqinfo(seqnames = "chr1")
  np2 <- importPeaks(fl, type = "narrow", seqinfo = sq, setNames = FALSE)
  bp2 <- importPeaks(fl, type = "broad", seqinfo = sq, setNames = FALSE)
  expect_true(is(c(np2, bp2), "GRangesList"))
  expect_equal(vapply(c(np1, bp1), length, integer(1)), c(0, 0))
})
