#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#'
#' @returns A `shiny.tag.list` with the user interface of the application.
#'
#' @import shiny
#' @importFrom bslib page_navbar nav_panel
#' @importFrom bsicons bs_icon
#' @importFrom utils packageVersion
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    bslib::page_navbar(
      id = "main_nav",
      title = sprintf("lipidomics3 %s", utils::packageVersion("lipidomics3")),
      theme = app_theme(),
      fillable = FALSE,
      bslib::nav_panel(
        title = "Meta data",
        icon = bsicons::bs_icon("table"),
        mod_metadata_ui("metadata_1")
      ),
      bslib::nav_panel(
        title = "Raw data",
        icon = bsicons::bs_icon("file-earmark-binary"),
        mod_rawdata_ui("rawdata_1")
      ),
      bslib::nav_panel(
        title = "Peak picking",
        icon = bsicons::bs_icon("graph-up"),
        mod_peakpicking_ui("peakpicking_1")
      ),
      bslib::nav_panel(
        title = "MS/MS spectra",
        icon = bsicons::bs_icon("soundwave"),
        mod_ms2_ui("ms2_1")
      ),
      bslib::nav_panel(
        title = "Identification",
        icon = bsicons::bs_icon("search"),
        mod_identification_ui("identification_1")
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @returns A `shiny.tag` with the external resources of the application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
	add_resource_path(
		"www",
		app_sys("app/www")
	)

	tags$head(
		favicon(),
		bundle_resources(
			path = app_sys("app/www"),
			app_title = "lipidomics3"
		)
		# Add here other external resources
		# for example, you can add shinyalert::useShinyalert()
	)
}
