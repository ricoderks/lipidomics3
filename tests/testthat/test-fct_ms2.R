test_that("ms2_spectra_table() returns an empty table for no spectra", {
  empty <- ms2_spectra_table(sps = Spectra::Spectra(), x = NULL)

  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0)
  expect_true(
    all(
      c("spectrum", "peak_id", "sample_name", "precursor_mz", "rt") %in%
        colnames(empty)
    )
  )
})

test_that("plot_ms2_spectrum() labels the most intense fragments", {
  peaks <- data.frame(
    mz = c(100.1, 200.2, 300.3),
    intensity = c(10, 1000, 500)
  )

  plot <- plot_ms2_spectrum(peaks, title = "test", n_labels = 2)

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$layers[[2]]$data), 2)
  expect_equal(plot$layers[[2]]$data$mz, c(200.2, 300.3))
})

test_that("plot_ms2_spectrum() only labels the fragments above the cut off", {
  peaks <- data.frame(
    mz = c(100.1, 200.2, 300.3),
    intensity = c(10, 1000, 500)
  )

  plot <- plot_ms2_spectrum(peaks, n_labels = 10, min_rel_intensity = 0.05)

  expect_equal(nrow(plot$layers[[2]]$data), 2)
})

test_that("plot_ms2_spectrum() handles an empty spectrum", {
  plot <- plot_ms2_spectrum(data.frame(mz = numeric(0), intensity = numeric(0)))

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$layers[[2]]$data), 0)
})
