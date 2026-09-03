test_that("build_sample_data() keeps the order of the raw data files", {
  meta_data <- data.frame(
    sample_name = c("s1", "s2"),
    file_name = c("a.mzML", "b.mzML"),
    sample_group = c("QC", "blank"),
    stringsAsFactors = FALSE
  )
  map <- list(sample = "sample_name", file = "file_name", group = "sample_group")
  file_info <- data.frame(
    name = c("b.mzML", "a.mzML"),
    path = c("/tmp/b.mzML", "/tmp/a.mzML"),
    meta_row = c(2L, 1L),
    stringsAsFactors = FALSE
  )

  sample_data <- build_sample_data(file_info, meta_data, map)

  expect_equal(sample_data$sample_name, c("s2", "s1"))
  expect_equal(sample_data$sample_group, c("blank", "QC"))
  expect_equal(sample_data$file_name, c("b.mzML", "a.mzML"))
})

test_that("build_sample_data() falls back to the file name", {
  meta_data <- data.frame(
    sample_name = "s1",
    file_name = "a.mzML",
    stringsAsFactors = FALSE
  )
  map <- list(sample = "sample_name", file = "file_name", group = NA_character_)
  file_info <- data.frame(
    name = c("a.mzML", "unknown.mzML"),
    path = c("/tmp/a.mzML", "/tmp/unknown.mzML"),
    meta_row = c(1L, NA_integer_),
    stringsAsFactors = FALSE
  )

  sample_data <- build_sample_data(file_info, meta_data, map)

  expect_equal(sample_data$sample_name, c("s1", "unknown"))
  expect_equal(sample_data$sample_group, c("all", "all"))
})
