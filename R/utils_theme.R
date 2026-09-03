#' Application theme
#'
#' Creates the Bootstrap 5 theme that is used throughout the application.
#'
#' @returns A `bs_theme` object as created by [bslib::bs_theme()].
#'
#' @importFrom bslib bs_theme
#' @noRd
app_theme <- function() {
  bslib::bs_theme(
    version = 5,
    preset = "shiny",
    primary = "#1b5e7e",
    success = "#3f8f5b",
    warning = "#d9922b",
    danger = "#b4453c"
  )
}


#' Show an error notification
#'
#' Shows the message of a condition as a Shiny error notification. Used to
#' report failures of the long running mass spectrometry steps without
#' crashing the session.
#'
#' @param title Character(1), short title shown in bold.
#' @param cnd A condition object, i.e. the error that was caught.
#'
#' @returns The id of the notification, invisibly.
#'
#' @importFrom shiny showNotification
#' @importFrom htmltools tags tagList
#' @noRd
notify_error <- function(title, cnd) {
  shiny::showNotification(
    ui = htmltools::tagList(
      htmltools::tags$b(title),
      htmltools::tags$br(),
      conditionMessage(cnd)
    ),
    type = "error",
    duration = NULL
  )
}
