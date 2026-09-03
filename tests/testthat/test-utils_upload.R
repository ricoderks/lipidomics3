test_that("stage_upload() restores the original file names", {
  upload_dir <- withr::local_tempdir()
  staging_dir <- withr::local_tempdir()

  datapath <- file.path(upload_dir, c("0.mzML", "1.mzML"))
  writeLines("a", datapath[1])
  writeLines("b", datapath[2])

  staged <- stage_upload(
    upload = data.frame(
      name = c("QC pool_006.mzML", "Blank_002.mzML"),
      datapath = datapath,
      stringsAsFactors = FALSE
    ),
    dir = staging_dir
  )

  expect_equal(staged$name, c("QC pool_006.mzML", "Blank_002.mzML"))
  expect_true(all(file.exists(staged$path)))
  expect_equal(readLines(staged$path[1]), "a")
})

test_that("stage_upload() strips a directory from the uploaded name", {
  upload_dir <- withr::local_tempdir()
  staging_dir <- withr::local_tempdir()

  datapath <- file.path(upload_dir, "0.mzML")
  writeLines("a", datapath)

  staged <- stage_upload(
    upload = data.frame(
      name = "../escaped.mzML",
      datapath = datapath,
      stringsAsFactors = FALSE
    ),
    dir = staging_dir
  )

  expect_equal(staged$name, "escaped.mzML")
  expect_equal(dirname(staged$path), staging_dir)
})

test_that("stage_upload() creates the staging directory", {
  parent <- withr::local_tempdir()
  staging_dir <- file.path(parent, "does", "not", "exist")

  datapath <- file.path(parent, "0.mzML")
  writeLines("a", datapath)

  staged <- stage_upload(
    upload = data.frame(name = "a.mzML", datapath = datapath, stringsAsFactors = FALSE),
    dir = staging_dir
  )

  expect_true(dir.exists(staging_dir))
  expect_true(file.exists(staged$path))
})

test_that("format_file_size() formats bytes", {
  expect_equal(format_file_size(1024), "1 Kb")
  expect_length(format_file_size(c(1024, 1048576)), 2)
})
