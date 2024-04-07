# Changes in version 0.99.0 (2022-03-28)

- Submitted to Bioconductor

# Changes in version 1.1.1 (2022-06-01)

- Bugfix so `as_tibble()` respects original column names

# Changes in version 1.1.2 (2022-07-22)

- Added `collapseTranscripts = "auto"` as the default for `plotHFGC()`
- `getProfileData()` now returns log2 transformed data by default 
- `getProfileData()` now uses `bplapply()` internally
- Bugfixes for `plotPie()`, `distinctMC()` `colToRanges()` and `stitchRanges()`

# Changes in version 1.1.3 (2022-08-01)

- Added `plotOverlaps()` for generation of Venn Diagrams and ComplexUpset plots

# Changes in version 1.1.5 (2022-10-11)

- Added `makeConsensus()` and updated vignette

# Changes in version 1.3.4 (22-12-20)

- Added `mergeBySig()`

# Changes in version 1.3.5 (23-01-31)

- Added `plotSplitDonut()`
- Bugfix in `plotAssayDensities()` and `plotAssayRle()` along with enabling plotting by group

# Changes in 1.3.6 (23-03-12)

- Expanded arguments for `plotSpliDonut()`
- Added `mergeByHMP()` for merging overlapping windows using the harmonic mean p

# Changes in 1.3.7 (23-03-20)

- Added `plotAssayHeatmap()`
- Added `fitAssayDiff()` and added coercion of `TopTags` objects

# Changes in 1.3.8 (23-03-20)

- Added using coverage to makeConsensus()

# Changes in 1.3.9 (23-04-10)

- Changed labelling strategy for `plotPie()`

# Changes in 1.4.2 (23-06-20)

- Added Fixed-width vignette & edited sliding window
- Added `se` and `peaks` as example data for man pages and vignettes
- Included `defineRegions()`
- Enabled `plotHFCG()` without Ideogram tracks
- Enabled use of offsets for normalisation in `fitAssayDiff()`

# Changes in 1.5.5 (23-06-28)

- Added `plotPairwise()` and `addDiffStatus()`

# Changes in 1.5.6 (23-07-01)

- Added `mapGrlCols()`

# Changes in 1.5.7 (23-07-08)

- Matched DiffBind and csaw settings for `fitAssayDiff()`
- Added `min_win` to all merging functions
- Added `n_max` to `getProfileData()`

# Changes in 1.5.8 (23-07-14)

- Added control of side-axis label position for `plotProfileHeatmap()`
- Added option to return merged key-value ranges for `mergeByHMP()`

# Changes in 1.5.10 (23-08-11)

- Added bed format to `importPeaks()`

# Changes in 1.5.11 (23-08-26)

- Added `defineSeqinfo()`

# Changes in 1.5.12 (23-08-30)

- Added `plotGrlCol()`

# Changes in 1.5.13 (23-09-24)

- Added handling of unquoted column names to most plotting functions
- Added passing of specific columns to dualFilter
- Added `drop` to addDiffStatus

# Changes in 1.5.14 (23-09-24)

- Added p_mu0 to output of `fitAssayDiff()`
- Added `respectLevels` and filtering to `plotProfileHeatmap()`

# Changes in 1.7.6

- Added the DESeq2 Wald statistic to options for `fitAssayDiff()`

# Changes in 1.7.7

- Added `merge_within` to `makeConsensus()` for better handling when `method = "coverage"`
