example_matches <- function() {
  data.frame(
    spectrum = 1:4,
    peak_id = c("CP0002", "CP0002", "CP0001", "CP0001"),
    sample_name = "sample_a",
    precursor_mz = c(786.6343, 786.6343, 720.5530, 720.5530),
    rt = c(8.203, 8.310, 7.402, 7.512),
    peak_mz = c(786.6343, 786.6343, 720.5530, 720.5530),
    peak_rt = c(8.21, 8.21, 7.45, 7.45),
    library_id = 1:4,
    name = c("PE 18:0_21:3", "PE O-19:1_21:1", "PE 16:0_18:0", "PE 16:0_18:0"),
    lipid_class = c("PE", "EtherPE", "PE", "PE"),
    adduct = "[M+H]+",
    library_mz = c(786.6343, 786.6343, 720.5530, 720.5530),
    library_rt = c(8.1, 8.4, 7.3, 7.3),
    dot = c(0.60, 0.80, 0.40, 0.70),
    weighted_dot = c(0.70, 0.90, 0.50, 0.80),
    reverse_dot = c(0.75, 0.95, 0.55, 0.85),
    n_matched = 2L,
    n_library_peaks = 10L,
    stringsAsFactors = FALSE
  )
}


test_that("peak_match_summary() collapses the hits to one row per peak", {
  summary <- peak_match_summary(example_matches())

  expect_equal(nrow(summary), 2)
  expect_equal(summary$peak_id, c("CP0001", "CP0002"))
  expect_equal(summary$peak_rt, c(7.45, 8.21))
  expect_equal(summary$n_hits, c(2L, 2L))
  expect_equal(summary$n_spectra, c(2L, 2L))
})


test_that("peak_match_summary() counts how ambiguous a peak is", {
  summary <- peak_match_summary(example_matches())

  # Both hits of CP0001 are the same lipid, CP0002 has two different lipids
  # of two different classes.
  expect_equal(summary$n_lipids, c(1L, 2L))
  expect_equal(summary$n_classes, c(1L, 2L))
})


test_that("peak_match_summary() reports the best hit of the chosen score", {
  by_weighted <- peak_match_summary(example_matches())
  by_dot <- peak_match_summary(example_matches(), score_column = "dot")

  expect_equal(by_weighted$best_score, c(0.80, 0.90))
  expect_equal(by_weighted$best_name, c("PE 16:0_18:0", "PE O-19:1_21:1"))
  expect_equal(by_dot$best_score, c(0.70, 0.80))
})


test_that("peak_match_summary() returns an empty summary without hits", {
  expect_equal(nrow(peak_match_summary(NULL)), 0)
  expect_equal(nrow(peak_match_summary(empty_match_table())), 0)
  expect_equal(
    colnames(peak_match_summary(NULL)),
    colnames(peak_match_summary(example_matches()))
  )
})


test_that("peak_coordinate() falls back on the spectrum", {
  matches <- example_matches()
  matches$peak_rt[1] <- NA_real_

  expect_equal(peak_coordinate(matches, "peak_rt", "rt")[1], 8.203)
  expect_equal(peak_coordinate(matches, "peak_rt", "rt")[2], 8.21)

  matches$peak_rt <- NULL

  expect_equal(peak_coordinate(matches, "peak_rt", "rt"), matches$rt)
})


test_that("peak_group_label() names the peak and counts its hits", {
  labels <- peak_group_label(example_matches())

  expect_equal(length(labels), 4)
  expect_equal(length(unique(labels)), 2)
  expect_equal(
    labels[1],
    "CP0002 | m/z 786.6343 | 8.21 min | sample_a | 2 hit(s)"
  )
  expect_equal(peak_group_label(NULL), character(0))
})


test_that("group_matches_by_peak() sorts the peaks by retention time", {
  grouped <- group_matches_by_peak(example_matches())

  expect_equal(colnames(grouped)[1:2], c("peak_order", "peak"))
  expect_equal(grouped$peak_id, c("CP0001", "CP0001", "CP0002", "CP0002"))
  expect_equal(grouped$peak_order, c(1L, 1L, 2L, 2L))
})


test_that("group_matches_by_peak() puts the best hit of a peak first", {
  grouped <- group_matches_by_peak(example_matches())

  expect_equal(grouped$weighted_dot, c(0.80, 0.50, 0.90, 0.70))

  by_dot <- group_matches_by_peak(example_matches(), score_column = "dot")

  expect_equal(by_dot$dot, c(0.70, 0.40, 0.80, 0.60))
})


test_that("group_matches_by_peak() keeps its columns without hits", {
  grouped <- group_matches_by_peak(empty_match_table())

  expect_equal(nrow(grouped), 0)
  expect_equal(
    colnames(grouped),
    colnames(group_matches_by_peak(example_matches()))
  )
})


test_that("plot_peak_map() draws the peaks and the selection", {
  peaks <- peak_match_summary(example_matches())
  plot <- plot_peak_map(peaks = peaks, selected = "CP0002")

  expect_s3_class(plot, "plotly")

  built <- plotly::plotly_build(plot)$x$data

  expect_equal(length(built), 2)
  expect_equal(as.numeric(built[[1]]$x), peaks$peak_rt)
  expect_equal(as.numeric(built[[1]]$y), peaks$peak_mz)
  expect_equal(as.numeric(built[[2]]$x), 8.21)
  expect_equal(as.numeric(built[[2]]$y), 786.6343)
  expect_true(built[[2]]$visible)
})


test_that("plot_peak_map() colours the peaks by the chosen column", {
  peaks <- peak_match_summary(example_matches())

  built <- plotly::plotly_build(
    plot_peak_map(peaks = peaks, colour_by = "n_lipids")
  )$x$data

  expect_equal(as.integer(built[[1]]$marker$color), c(1L, 2L))

  # An unknown column falls back on the best score, and nothing is marked.
  built <- plotly::plotly_build(
    plot_peak_map(peaks = peaks, colour_by = "does_not_exist")
  )$x$data

  expect_equal(as.numeric(built[[1]]$marker$color), peaks$best_score)
  expect_false(built[[2]]$visible)
})


test_that("map_colour_label() labels the colour bar", {
  expect_equal(map_colour_label("n_lipids"), "lipids")
  expect_equal(map_colour_label("best_score"), "best score")
  expect_equal(map_colour_label("nonsense"), "best score")
})
