test_that("centwave_param() passes the settings on", {
  param <- centwave_param(
    ppm = 15,
    peakwidth = c(5, 40),
    snthresh = 8,
    prefilter = c(4, 2000),
    noise = 100,
    integrate = 2L,
    fitgauss = TRUE
  )

  expect_s4_class(param, "CentWaveParam")
  expect_equal(param@ppm, 15)
  expect_equal(param@peakwidth, c(5, 40))
  expect_equal(param@snthresh, 8)
  expect_equal(param@prefilter, c(4, 2000))
  expect_equal(param@noise, 100)
  expect_equal(param@integrate, 2L)
  expect_true(param@fitgauss)
})

test_that("merge_peaks_param() passes the settings on", {
  param <- merge_peaks_param(expandRt = 3, expandMz = 0.01, ppm = 5, minProp = 0.5)

  expect_s4_class(param, "MergeNeighboringPeaksParam")
  expect_equal(param@expandRt, 3)
  expect_equal(param@expandMz, 0.01)
  expect_equal(param@ppm, 5)
  expect_equal(param@minProp, 0.5)
})

test_that("parallel_param() uses serial processing for a single worker", {
  expect_s4_class(parallel_param(1), "SerialParam")
  expect_s4_class(parallel_param(0), "SerialParam")
  expect_s4_class(parallel_param(2), "BiocParallelParam")
  # R CMD check limits the number of cores, so do not ask for more than two.
  expect_equal(BiocParallel::bpnworkers(parallel_param(2)), 2L)
})
