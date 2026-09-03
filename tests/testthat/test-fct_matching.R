peaks <- function(mz, intensity) {
  matrix(
    c(mz, intensity),
    ncol = 2L,
    dimnames = list(NULL, c("mz", "intensity"))
  )
}

test_that("align_peaks() pairs the peaks that are within the tolerance", {
  x <- peaks(c(100, 200, 300), c(10, 20, 30))
  y <- peaks(c(100.004, 300), c(1, 3))

  aligned <- align_peaks(x, y, tolerance = 0.01, ppm = 0)

  expect_equal(nrow(aligned$x), 3)
  expect_equal(aligned$matched, c(TRUE, FALSE, TRUE))
  expect_equal(aligned$y[, "intensity"], c(1, NA, 3))
})

test_that("align_peaks() separates peaks outside the tolerance", {
  x <- peaks(100, 10)
  y <- peaks(100.05, 10)

  aligned <- align_peaks(x, y, tolerance = 0.01, ppm = 0)

  expect_equal(nrow(aligned$x), 2)
  expect_false(any(aligned$matched))
})

test_that("align_peaks() copes with an empty spectrum", {
  aligned <- align_peaks(peaks(numeric(0), numeric(0)), peaks(100, 10))

  expect_equal(nrow(aligned$x), 0)
  expect_length(aligned$matched, 0)
})

test_that("identical spectra score one on all three scores", {
  x <- peaks(c(100, 200, 300), c(10, 100, 55))
  aligned <- align_peaks(x, x, tolerance = 0.01, ppm = 0)

  expect_equal(dot_product(aligned), 1)
  expect_equal(weighted_dot_product(aligned), 1)
  expect_equal(reverse_dot_product(aligned), 1)
})

test_that("a peak only the query has lowers the dot product but not the reverse", {
  query <- peaks(c(100, 200, 300), c(100, 50, 10))
  reference <- peaks(c(100, 200), c(100, 50))
  aligned <- align_peaks(query, reference, tolerance = 0.01, ppm = 0)

  # 12500^2 / (12600 * 12500)
  expect_equal(dot_product(aligned), 12500^2 / (12600 * 12500))
  expect_equal(reverse_dot_product(aligned), 1)
  expect_equal(sum(aligned$matched), 2)
})

test_that("a peak only the reference has lowers both the dot and the reverse", {
  query <- peaks(c(100, 200), c(100, 50))
  reference <- peaks(c(100, 200, 300), c(100, 50, 10))
  aligned <- align_peaks(query, reference, tolerance = 0.01, ppm = 0)

  expected <- 12500^2 / (12500 * 12600)

  expect_equal(dot_product(aligned), expected)
  expect_equal(reverse_dot_product(aligned), expected)
})

test_that("the weighted dot product without weights is the dot product", {
  query <- peaks(c(100, 200, 300), c(100, 50, 10))
  reference <- peaks(c(100, 200), c(100, 50))
  aligned <- align_peaks(query, reference, tolerance = 0.01, ppm = 0)

  expect_equal(
    weighted_dot_product(aligned, m = 0, n = 1),
    dot_product(aligned)
  )
})

test_that("the weighted dot product weights by m/z and intensity", {
  query <- peaks(c(100, 200), c(100, 50))
  reference <- peaks(c(100, 200), c(50, 100))
  aligned <- align_peaks(query, reference, tolerance = 0.01, ppm = 0)

  wx <- c(100, 200)^3 * c(100, 50)^0.6
  wy <- c(100, 200)^3 * c(50, 100)^0.6

  expect_equal(
    weighted_dot_product(aligned, m = 3, n = 0.6),
    sum(wx * wy)^2 / (sum(wx^2) * sum(wy^2))
  )
})

test_that("spectrum_scores() returns all three scores and the match count", {
  query <- peaks(c(100, 200, 300), c(100, 50, 10))
  reference <- peaks(c(100, 200), c(100, 50))

  scores <- spectrum_scores(query, reference, tolerance = 0.01, ppm = 0)

  expect_named(scores, c("dot", "weighted_dot", "reverse_dot", "n_matched"))
  expect_equal(unname(scores["n_matched"]), 2)
  expect_equal(unname(scores["reverse_dot"]), 1)
})

test_that("the scores are NA when a spectrum has no peaks", {
  aligned <- align_peaks(peaks(numeric(0), numeric(0)), peaks(100, 10))

  expect_true(is.na(dot_product(aligned)))
  expect_true(is.na(reverse_dot_product(aligned)))
})

test_that("match_spectrum() keeps the best hits above the match count", {
  encode <- function(p) {
    writeBin(as.vector(t(p)), raw(), size = 4L, endian = "little")
  }

  query <- peaks(c(100, 200, 300), c(100, 50, 10))
  candidates <- data.frame(
    id = 1:3,
    name = c("perfect", "partial", "one peak"),
    lipid_class = "PC",
    precursor_mz = 500,
    adduct = "[M+H]+",
    retention_time = 5,
    n_peaks = c(3L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  candidates$peaks <- list(
    encode(peaks(c(100, 200, 300), c(100, 50, 10))),
    encode(peaks(c(100, 200), c(100, 50))),
    encode(peaks(100, 100))
  )

  hits <- match_spectrum(
    query = query,
    candidates = candidates,
    tolerance = 0.01,
    ppm = 0,
    min_matched = 2,
    top_n = 5
  )

  # The single peak reference matches only one fragment and is dropped.
  expect_equal(nrow(hits), 2)
  expect_equal(hits$name[1], "perfect")
  expect_equal(hits$dot[1], 1)
  expect_true(all(hits$n_matched >= 2))
})

test_that("match_spectrum() returns NULL without candidates", {
  expect_null(
    match_spectrum(peaks(100, 10), candidates = data.frame())
  )
})

test_that("empty_match_table() has the columns of a real hit table", {
  empty <- empty_match_table()

  expect_equal(nrow(empty), 0)
  expect_true(
    all(
      c("peak_id", "name", "dot", "weighted_dot", "reverse_dot", "n_matched") %in%
        colnames(empty)
    )
  )
})

test_that("mirror_spectrum_data() puts the reference below the axis", {
  query <- peaks(c(100, 200), c(100, 50))
  reference <- peaks(c(100, 200), c(50, 100))

  plot_data <- mirror_spectrum_data(query, reference, tolerance = 0.01, ppm = 0)

  expect_equal(nrow(plot_data), 4)
  expect_true(all(plot_data$intensity[plot_data$spectrum == "query"] > 0))
  expect_true(all(plot_data$intensity[plot_data$spectrum == "reference"] < 0))
  expect_true(all(plot_data$matched))
})

test_that("mirror_spectrum_data() scales both spectra to their base peak", {
  query <- peaks(c(100, 200), c(20, 10))
  reference <- peaks(c(100, 200), c(4000, 1000))

  plot_data <- mirror_spectrum_data(query, reference, tolerance = 0.01, ppm = 0)

  expect_equal(plot_data$intensity[plot_data$spectrum == "query"], c(100, 50))
  expect_equal(plot_data$intensity[plot_data$spectrum == "reference"], c(-100, -25))
})

test_that("mirror_spectrum_data() marks the peaks that only one spectrum has", {
  query <- peaks(c(100, 300), c(100, 50))
  reference <- peaks(100, 100)

  plot_data <- mirror_spectrum_data(query, reference, tolerance = 0.01, ppm = 0)

  # two query peaks and one reference peak, the peak at 300 is unmatched
  expect_equal(nrow(plot_data), 3)
  expect_equal(sum(!plot_data$matched), 1)
  expect_equal(plot_data$mz[!plot_data$matched], 300)
})

test_that("mirror_spectrum_data() removes the precursor from both spectra", {
  query <- peaks(c(100, 800), c(50, 100))
  reference <- peaks(c(100, 800), c(50, 100))

  plot_data <- mirror_spectrum_data(
    query, reference,
    tolerance = 0.01, ppm = 0, precursor_mz = 800
  )

  expect_equal(unique(plot_data$mz), 100)
})

test_that("mirror_spectrum_data() copes with an empty spectrum", {
  plot_data <- mirror_spectrum_data(
    peaks(numeric(0), numeric(0)),
    peaks(100, 10)
  )

  expect_s3_class(plot_data, "data.frame")
  expect_equal(nrow(plot_data), 0)
})

test_that("plot_mirror_spectrum() returns an interactive plot", {
  query <- peaks(c(100, 200), c(100, 50))
  reference <- peaks(c(100, 300), c(50, 100))

  plot <- plot_mirror_spectrum(
    query, reference,
    tolerance = 0.01, ppm = 0,
    title = "PC 16:0_18:1", subtitle = "dot 0.500"
  )

  expect_s3_class(plot, "plotly")

  built <- plotly::plotly_build(plot)
  # one segment trace and one hover marker trace for matched and unmatched
  expect_equal(length(built$x$data), 4)
  expect_true(any(grepl("PC 16:0_18:1", built$x$layout$title$text, fixed = TRUE)))
})

test_that("plot_mirror_spectrum() survives an empty spectrum", {
  plot <- plot_mirror_spectrum(
    peaks(numeric(0), numeric(0)),
    peaks(numeric(0), numeric(0))
  )

  expect_s3_class(plot, "plotly")
})

test_that("drop_precursor_peaks() removes the precursor and the window below it", {
  x <- peaks(c(100, 200, 798.9, 800), c(10, 20, 30, 40))

  kept <- drop_precursor_peaks(x, precursor_mz = 800, window = 1.5)

  expect_equal(kept[, "mz"], c(100, 200))
  expect_equal(nrow(drop_precursor_peaks(x, precursor_mz = NA)), 4)
  expect_equal(nrow(drop_precursor_peaks(x, precursor_mz = NULL)), 4)
})

test_that("removing the precursor stops it from carrying the score", {
  # The reference disagrees with the query on every real fragment, they only
  # share the precursor at m/z 800.
  query <- peaks(c(100, 200, 800), c(999, 500, 300))
  reference <- peaks(c(150, 250, 800), c(999, 500, 300))

  with_precursor <- spectrum_scores(
    query, reference, tolerance = 0.01, ppm = 0, precursor_mz = NA_real_
  )
  without_precursor <- spectrum_scores(
    query, reference, tolerance = 0.01, ppm = 0, precursor_mz = 800
  )

  expect_equal(unname(with_precursor["n_matched"]), 1)
  expect_gt(with_precursor["weighted_dot"], 0.99)

  expect_equal(unname(without_precursor["n_matched"]), 0)
  expect_equal(unname(without_precursor["weighted_dot"]), 0)
})
