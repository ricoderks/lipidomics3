#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#'
#' @returns Nothing, the function is called for its side effects.
#'
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # The state that is shared between the modules. Every step of the workflow
  # writes its result here and invalidates the results of the steps after it.
  r <- shiny::reactiveValues(
    metadata = NULL,
    metadata_map = NULL,
    raw_files = NULL,
    ms_data = NULL,
    xcms_data = NULL,
    ms2_spectra = NULL,
    matches = NULL
  )

  mod_metadata_server(id = "metadata_1", r = r)
  mod_rawdata_server(id = "rawdata_1", r = r)
  mod_peakpicking_server(id = "peakpicking_1", r = r)
  mod_ms2_server(id = "ms2_1", r = r)
  mod_identification_server(id = "identification_1", r = r)

  # The uploaded files can be large, so they are removed when the user leaves.
  upload_dir <- file.path(tempdir(), "lipidomics3", session$token)

  session$onSessionEnded(function() {
    unlink(x = upload_dir, recursive = TRUE, force = TRUE)
  })
}
