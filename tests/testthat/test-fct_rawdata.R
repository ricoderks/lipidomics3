test_that("build_sample_data() keeps the order of the raw data files", {
  meta_data <- data.frame(
    sample_name = c("s1", "s2"),
    file_name = c("a.mzML", "b.mzML"),
    sample_group = c("QC", "blank"),
    sample_type = c("qc", "blank"),
    stringsAsFactors = FALSE
  )
  map <- list(
    sample = "sample_name",
    file = "file_name",
    group = "sample_group",
    type = "sample_type"
  )
  file_info <- data.frame(
    name = c("b.mzML", "a.mzML"),
    path = c("/tmp/b.mzML", "/tmp/a.mzML"),
    meta_row = c(2L, 1L),
    stringsAsFactors = FALSE
  )

  sample_data <- build_sample_data(file_info, meta_data, map)

  expect_equal(sample_data$sample_name, c("s2", "s1"))
  expect_equal(sample_data$sample_group, c("blank", "QC"))
  expect_equal(sample_data$sample_type, c("blank", "qc"))
  expect_equal(sample_data$file_name, c("b.mzML", "a.mzML"))
})

test_that("build_sample_data() falls back to the file name", {
  meta_data <- data.frame(
    sample_name = "s1",
    file_name = "a.mzML",
    stringsAsFactors = FALSE
  )
  map <- list(
    sample = "sample_name",
    file = "file_name",
    group = NA_character_,
    type = NA_character_
  )
  file_info <- data.frame(
    name = c("a.mzML", "unknown.mzML"),
    path = c("/tmp/a.mzML", "/tmp/unknown.mzML"),
    meta_row = c(1L, NA_integer_),
    stringsAsFactors = FALSE
  )

  sample_data <- build_sample_data(file_info, meta_data, map)

  expect_equal(sample_data$sample_name, c("s1", "unknown"))
  expect_equal(sample_data$sample_group, c("all", "all"))
  expect_equal(sample_data$sample_type, c("sample", "sample"))
})

test_that("build_sample_data() names a missing sample type", {
  meta_data <- data.frame(
    sample_name = c("s1", "s2"),
    file_name = c("a.mzML", "b.mzML"),
    type = c("qc", NA_character_),
    stringsAsFactors = FALSE
  )
  map <- list(
    sample = "sample_name",
    file = "file_name",
    group = NA_character_,
    type = "type"
  )
  file_info <- data.frame(
    name = c("a.mzML", "b.mzML"),
    path = c("/tmp/a.mzML", "/tmp/b.mzML"),
    meta_row = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  sample_data <- build_sample_data(file_info, meta_data, map)

  expect_equal(sample_data$sample_type, c("qc", "unknown"))
})

test_that("mz_window() adds both tolerances", {
  window <- mz_window(mz = 700, tolerance = 0.01, ppm = 0)

  expect_true(is.matrix(window))
  expect_equal(dim(window), c(1L, 2L))
  expect_equal(as.numeric(window), c(699.99, 700.01))

  expect_equal(
    as.numeric(mz_window(mz = 1000, tolerance = 0, ppm = 10)),
    c(999.99, 1000.01)
  )
  expect_equal(
    as.numeric(mz_window(mz = 1000, tolerance = 0.01, ppm = 10)),
    c(999.98, 1000.02)
  )
})

test_that("mz_window() keeps a missing m/z missing", {
  expect_true(all(is.na(mz_window(mz = NA_real_))))
})

test_that("chromatograms_data() puts the samples below each other", {
  skip_if_not_installed("MSnbase")

  chroms <- MSnbase::MChromatograms(
    data = list(
      MSnbase::Chromatogram(rtime = c(60, 120), intensity = c(3, 9)),
      MSnbase::Chromatogram(rtime = c(60, 120), intensity = c(4, NA))
    ),
    nrow = 1
  )
  sample_data <- data.frame(
    sample_name = c("s1", "s2"),
    sample_group = c("QC", "blank"),
    stringsAsFactors = FALSE
  )

  chrom_data <- chromatograms_data(chroms = chroms, sample_data = sample_data)

  expect_equal(nrow(chrom_data), 4)
  expect_equal(chrom_data$rtime, c(1, 2, 1, 2))
  # A spectrum without a peak in the m/z window has no signal, not a missing
  # value.
  expect_equal(chrom_data$intensity, c(3, 9, 4, 0))
  expect_equal(chrom_data$sample_name, c("s1", "s1", "s2", "s2"))
  expect_equal(chrom_data$sample_group, c("QC", "QC", "blank", "blank"))
})

test_that("plot_chromatograms() draws one line per sample", {
  chrom_data <- data.frame(
    rtime = c(1, 2, 3, 1, 2, 3),
    intensity = c(10, 50, 20, 5, 8, 3),
    sample_name = rep(c("s1", "s2"), each = 3),
    sample_group = rep(c("QC", "blank"), each = 3),
    stringsAsFactors = FALSE
  )

  plot <- plot_chromatograms(chrom_data = chrom_data, y_label = "TIC")
  built <- plotly::plotly_build(plot)$x$data

  expect_s3_class(plot, "plotly")
  expect_equal(length(built), 2)
  expect_equal(as.numeric(built[[1]]$x), c(1, 2, 3))
  expect_equal(as.numeric(built[[1]]$y), c(10, 50, 20))
  expect_equal(built[[1]]$name, "s1")
  expect_equal(as.numeric(built[[2]]$y), c(5, 8, 3))
  expect_equal(built[[2]]$name, "s2")
})

test_that("plot_chromatograms() handles an empty chromatogram", {
  chrom_data <- data.frame(
    rtime = numeric(0),
    intensity = numeric(0),
    sample_name = character(0),
    sample_group = character(0),
    stringsAsFactors = FALSE
  )

  plot <- plot_chromatograms(chrom_data = chrom_data)
  built <- suppressWarnings(plotly::plotly_build(plot)$x$data)

  expect_s3_class(plot, "plotly")
  expect_true(all(vapply(built, function(z) length(z[["x"]]) == 0, logical(1))))
})
