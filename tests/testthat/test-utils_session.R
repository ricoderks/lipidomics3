test_that("has_unsaved_work() is FALSE for an empty session", {
  r <- list(
    metadata = NULL,
    metadata_map = NULL,
    raw_files = NULL,
    ms_data = NULL,
    xcms_data = NULL,
    ms2_spectra = NULL,
    matches = NULL
  )

  expect_false(has_unsaved_work(r))
})

test_that("has_unsaved_work() finds the result of every step", {
  steps <- c(
    "metadata", "raw_files", "ms_data", "xcms_data", "ms2_spectra", "matches"
  )

  for (step in steps) {
    r <- list()
    r[[step]] <- "a result"

    expect_true(has_unsaved_work(r), label = step)
  }
})

test_that("has_unsaved_work() ignores the meta data column mapping", {
  # The mapping is derived from the meta data and is nothing on its own.
  expect_false(has_unsaved_work(list(metadata_map = list(sample = "name"))))
})
