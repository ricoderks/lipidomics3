test_that("decode_library_peaks() decodes the 4 byte float pairs", {
  original <- matrix(
    c(59.01385, 500, 113.09719, 800),
    ncol = 2L,
    byrow = TRUE
  )
  blob <- writeBin(as.vector(t(original)), raw(), size = 4L, endian = "little")

  decoded <- decode_library_peaks(blob, n_peaks = 2L)

  expect_equal(colnames(decoded), c("mz", "intensity"))
  expect_equal(nrow(decoded), 2)
  expect_equal(decoded[, "mz"], original[, 1], tolerance = 1e-5)
  expect_equal(decoded[, "intensity"], original[, 2], tolerance = 1e-5)
})

test_that("decode_library_peaks() returns an empty matrix for a broken blob", {
  expect_equal(nrow(decode_library_peaks(raw(0), n_peaks = 0L)), 0)
  expect_equal(nrow(decode_library_peaks(raw(4), n_peaks = 3L)), 0)
  expect_equal(nrow(decode_library_peaks(raw(8), n_peaks = NA)), 0)
})

test_that("lipid_db_connect() reports a missing database", {
  expect_error(
    lipid_db_connect(file.path(tempdir(), "does-not-exist.sqlite")),
    "was not found"
  )
})

test_that("polarity_to_ion_mode() translates the mzML polarity", {
  expect_equal(
    polarity_to_ion_mode(c(1L, 0L, NA_integer_, 9L)),
    c("Positive", "Negative", NA, NA)
  )
  expect_length(polarity_to_ion_mode(integer(0)), 0)
})
