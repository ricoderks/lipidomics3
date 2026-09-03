#' Run the Shiny Application
#'
#' Starts the lipidomics3 application. The maximum size of the files that can
#' be uploaded is taken from the `max_upload_size_mb` entry of the golem
#' configuration file, because mzML files are much larger than the Shiny
#' default of 5 MB.
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @returns An object that represents the application, see [shiny::shinyApp()].
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
	onStart = NULL,
	options = list(),
	enableBookmarking = NULL,
	uiPattern = "/",
	...
) {
  max_upload_size_mb <- get_golem_config("max_upload_size_mb")

  if (!is.null(max_upload_size_mb)) {
    base::options(shiny.maxRequestSize = max_upload_size_mb * 1024^2)
  }

	with_golem_options(
		app = shinyApp(
			ui = app_ui,
			server = app_server,
			onStart = onStart,
			options = options,
			enableBookmarking = enableBookmarking,
			uiPattern = uiPattern
		),
		golem_opts = list(...)
	)
}
