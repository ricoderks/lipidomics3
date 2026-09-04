#' Is there work that would be lost
#'
#' Tells whether the session holds a result that the user would lose when the
#' page is left, since everything the application creates lives in the memory
#' of the R session of the browser tab and is never written to disk.
#'
#' @param r A `reactiveValues` object with the state that is shared between the
#'   modules of the application.
#'
#' @returns `TRUE` when at least one step of the workflow has a result.
#'
#' @noRd
has_unsaved_work <- function(r) {
  results <- c(
    "metadata",
    "raw_files",
    "ms_data",
    "xcms_data",
    "ms2_spectra",
    "matches"
  )

  for (result in results) {
    if (!is.null(r[[result]])) {
      return(TRUE)
    }
  }

  FALSE
}
