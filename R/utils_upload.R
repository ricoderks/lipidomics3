#' Restore the names of uploaded files
#'
#' Moves files uploaded with [shiny::fileInput()] to a directory, using their
#' original file names. Shiny stores uploads under generated names, which
#' breaks the detection of the file type and the link between the raw data
#' files and the meta data.
#'
#' @param upload A `data.frame` as returned by [shiny::fileInput()], with at
#'   least the columns `name` and `datapath`.
#' @param dir Character(1), directory to move the files to. It is created when
#'   it does not exist yet.
#'
#' @returns A `data.frame` with the columns `name` and `path`, holding the
#'   original file name and the path of the staged file.
#'
#' @noRd
stage_upload <- function(upload, dir) {
  if (!dir.exists(dir)) {
    dir.create(path = dir, recursive = TRUE)
  }

  names <- basename(as.character(upload$name))
  paths <- file.path(dir, names)

  moved <- suppressWarnings(
    file.rename(from = upload$datapath, to = paths)
  )

  # file.rename() fails when the temporary upload directory is on another file
  # system than the staging directory.
  if (any(!moved)) {
    file.copy(
      from = upload$datapath[!moved],
      to = paths[!moved],
      overwrite = TRUE
    )
  }

  data.frame(
    name = names,
    path = paths,
    stringsAsFactors = FALSE
  )
}


#' Directory to stage the uploads of a session in
#'
#' Creates the path of the directory in which the files uploaded during a
#' Shiny session are stored.
#'
#' @param session The Shiny `session` object.
#' @param sub_dir Character(1), name of the sub directory to use.
#'
#' @returns Character(1), the path of the staging directory.
#'
#' @noRd
session_upload_dir <- function(session, sub_dir) {
  file.path(tempdir(), "lipidomics3", session$token, sub_dir)
}


#' Format a file size
#'
#' Formats a file size in bytes as a human readable string.
#'
#' @param size Numeric vector with the file sizes in bytes.
#'
#' @returns A character vector with the formatted file sizes.
#'
#' @noRd
format_file_size <- function(size) {
  vapply(
    X = size,
    FUN = function(x) format(structure(x, class = "object_size"), units = "auto"),
    FUN.VALUE = character(1)
  )
}
