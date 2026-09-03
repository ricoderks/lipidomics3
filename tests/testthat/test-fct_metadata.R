test_that("metadata_file_type() recognises the supported file types", {
  expect_equal(metadata_file_type("a/b/meta.xlsx"), "excel")
  expect_equal(metadata_file_type("meta.XLS"), "excel")
  expect_equal(metadata_file_type("meta.csv"), "csv")
  expect_equal(metadata_file_type("meta.tsv"), "tsv")
  expect_equal(metadata_file_type("meta.txt"), "tsv")
  expect_equal(metadata_file_type("meta.mzML"), "unknown")
})

test_that("metadata_sheets() returns nothing for text files", {
  expect_equal(metadata_sheets("meta.csv"), character(0))
})

test_that("guess_delimiter() picks the most frequent candidate", {
  comma <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("a,b,c", "1,2,3"), comma)
  expect_equal(guess_delimiter(comma), ",")

  semi <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("a;b;c", "1;2;3"), semi)
  expect_equal(guess_delimiter(semi), ";")

  tab <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("a\tb\tc", "1\t2\t3"), tab)
  expect_equal(guess_delimiter(tab), "\t")
})

test_that("read_metadata() reads delimited files and refuses other types", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("sample name;file name", "s1;a.mzML"), path)

  meta_data <- read_metadata(path)
  expect_s3_class(meta_data, "data.frame")
  expect_equal(nrow(meta_data), 1)

  expect_error(read_metadata("meta.mzML"), "Unsupported meta data file type")
})

test_that("guess_metadata_column() follows the order of the patterns", {
  columns <- c("sample", "sample_name", "file_name")

  expect_equal(
    guess_metadata_column(columns, c("^sample.?name$", "sample")),
    "sample_name"
  )
  expect_equal(guess_metadata_column(columns, c("nothing")), "")
})

test_that("check_metadata() finds the problems of the meta data", {
  meta_data <- data.frame(
    sample_name = c("s1", "s2"),
    file_name = c("a.mzML", "b.mzML"),
    stringsAsFactors = FALSE
  )

  expect_equal(check_metadata(meta_data, "sample_name", "file_name"), character(0))
  expect_equal(check_metadata(NULL, "sample_name", "file_name"), "No meta data has been read yet.")
  expect_match(
    check_metadata(meta_data, "does_not_exist", "file_name"),
    "is not available",
    all = FALSE
  )

  duplicated_data <- data.frame(
    sample_name = c("s1", "s1"),
    file_name = c("a.mzML", "a.mzML"),
    stringsAsFactors = FALSE
  )
  expect_length(check_metadata(duplicated_data, "sample_name", "file_name"), 2)

  empty_data <- data.frame(
    sample_name = c("s1", "s2"),
    file_name = c("a.mzML", ""),
    stringsAsFactors = FALSE
  )
  expect_match(
    check_metadata(empty_data, "sample_name", "file_name"),
    "empty file names",
    all = FALSE
  )
})

test_that("link_metadata_files() ignores the path and the extension", {
  meta_data <- data.frame(
    file_name = c("QC_pool_006.mzML", "Blank_002"),
    stringsAsFactors = FALSE
  )

  link <- link_metadata_files(
    meta_data = meta_data,
    file_column = "file_name",
    file_names = c("Blank_002.mzML", "QC_pool_006.mzML", "other.mzML")
  )

  expect_equal(link$meta_row, c(2L, 1L, NA_integer_))
})

test_that("normalise_file_key() only keeps letters and digits", {
  expect_equal(normalise_file_key("QC pool_006"), "qc_pool_006")
  expect_equal(normalise_file_key("QC-pool 006"), "qc_pool_006")
  expect_equal(normalise_file_key(" sample (1) "), "sample_1")
})

test_that("link_metadata_files() falls back to a normalised comparison", {
  meta_data <- data.frame(
    file_name = c("QC pool_006", "QC pool_012"),
    stringsAsFactors = FALSE
  )

  link <- link_metadata_files(
    meta_data = meta_data,
    file_column = "file_name",
    file_names = c("QC_pool_012.mzML", "QC_pool_006.mzML", "QC_pool_099.mzML")
  )

  expect_equal(link$meta_row, c(2L, 1L, NA_integer_))
})

test_that("link_metadata_files() prefers an exact match", {
  meta_data <- data.frame(
    file_name = c("a b", "a_b"),
    stringsAsFactors = FALSE
  )

  link <- link_metadata_files(
    meta_data = meta_data,
    file_column = "file_name",
    file_names = "a_b.mzML"
  )

  expect_equal(link$meta_row, 2L)
})
