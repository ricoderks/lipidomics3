#' Meta data UI Function
#'
#' The user interface of the module that reads the sample meta data and lets
#' the user point out which columns hold the sample name, the raw data file
#' name and the sample group.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#'
#' @returns A `shiny.tag.list` with the user interface of the module.
#'
#' @importFrom shiny NS fileInput selectInput uiOutput
#' @importFrom bslib layout_sidebar sidebar card card_header layout_column_wrap
#'   value_box
#' @importFrom bsicons bs_icon
#' @importFrom DT DTOutput
#' @noRd
mod_metadata_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      title = "Meta data file",
      shiny::fileInput(
        inputId = ns("file"),
        label = "Select a meta data file",
        accept = c(".xlsx", ".xlsm", ".xls", ".csv", ".tsv", ".txt"),
        width = "100%"
      ),
      shiny::selectInput(
        inputId = ns("sheet"),
        label = "Worksheet",
        choices = character(0),
        width = "100%"
      ),
      shiny::selectInput(
        inputId = ns("sample_column"),
        label = "Sample name column",
        choices = character(0),
        width = "100%"
      ),
      shiny::selectInput(
        inputId = ns("file_column"),
        label = "Raw data file name column",
        choices = character(0),
        width = "100%"
      ),
      shiny::selectInput(
        inputId = ns("group_column"),
        label = "Sample group column",
        choices = c("None" = ""),
        width = "100%"
      )
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      fill = FALSE,
      bslib::value_box(
        title = "Samples",
        value = shiny::textOutput(outputId = ns("n_samples")),
        showcase = bsicons::bs_icon("droplet-half"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Sample groups",
        value = shiny::textOutput(outputId = ns("n_groups")),
        showcase = bsicons::bs_icon("diagram-3"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Status",
        value = shiny::textOutput(outputId = ns("status")),
        showcase = bsicons::bs_icon("clipboard-check"),
        theme = "primary"
      )
    ),
    shiny::uiOutput(outputId = ns("problems")),
    bslib::card(
      full_screen = TRUE,
      height = 520,
      bslib::card_header("Sample meta data"),
      DT::DTOutput(outputId = ns("table"))
    )
  )
}


#' Meta data Server Functions
#'
#' The server side of the meta data module. Reads the meta data file, guesses
#' the meaning of the columns and stores the validated meta data in the shared
#' reactive values object.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#' @param r A `reactiveValues` object with the state that is shared between the
#'   modules of the application.
#'
#' @returns Nothing, the module is called for its side effects.
#'
#' @importFrom shiny moduleServer reactiveValues observeEvent observe req
#'   updateSelectInput renderText renderUI
#' @importFrom htmltools tags
#' @importFrom DT renderDT datatable
#' @noRd
mod_metadata_server <- function(id, r) {
  shiny::moduleServer(id, function(input, output, session) {
    local_r <- shiny::reactiveValues(
      meta_data = NULL,
      problems = character(0)
    )

    # A newly uploaded file is staged under its original name, so that both the
    # file type detection and the sheet names work as expected.
    staged_file <- shiny::reactive({
      shiny::req(input$file)

      stage_upload(
        upload = input$file,
        dir = session_upload_dir(session = session, sub_dir = "metadata")
      )
    })

    shiny::observeEvent(staged_file(), {
      sheets <- tryCatch(
        expr = metadata_sheets(path = staged_file()$path[1]),
        error = function(e) {
          notify_error(title = "Could not read the meta data file", cnd = e)
          character(0)
        }
      )

      shiny::updateSelectInput(
        session = session,
        inputId = "sheet",
        choices = if (length(sheets) > 0) sheets else c("Not applicable" = ""),
        selected = if (length(sheets) > 0) sheets[1] else ""
      )
    })

    shiny::observeEvent(list(staged_file(), input$sheet), {
      shiny::req(staged_file())

      meta_data <- tryCatch(
        expr = read_metadata(
          path = staged_file()$path[1],
          sheet = if (isTRUE(nzchar(input$sheet))) input$sheet else NULL
        ),
        error = function(e) {
          notify_error(title = "Could not read the meta data", cnd = e)
          NULL
        }
      )

      local_r$meta_data <- meta_data
      shiny::req(meta_data)

      columns <- colnames(meta_data)

      shiny::updateSelectInput(
        session = session,
        inputId = "sample_column",
        choices = columns,
        selected = guess_metadata_column(
          column_names = columns,
          patterns = c("^sample.?name$", "^sample.?id$", "^name$", "sample")
        )
      )
      shiny::updateSelectInput(
        session = session,
        inputId = "file_column",
        choices = columns,
        selected = guess_metadata_column(
          column_names = columns,
          patterns = c("^file.?name$", "mzml", "^file$", "raw", "file")
        )
      )
      shiny::updateSelectInput(
        session = session,
        inputId = "group_column",
        choices = c("None" = "", columns),
        selected = guess_metadata_column(
          column_names = columns,
          patterns = c("^sample.?group$", "^group$", "^class$", "^type$", "group")
        )
      )
    })

    # Only meta data that passes the checks is shared with the other modules.
    shiny::observe({
      meta_data <- local_r$meta_data

      if (is.null(meta_data)) {
        local_r$problems <- character(0)
        r$metadata <- NULL
        r$metadata_map <- NULL
        return()
      }

      problems <- check_metadata(
        meta_data = meta_data,
        sample_column = input$sample_column,
        file_column = input$file_column
      )

      local_r$problems <- problems

      if (length(problems) > 0) {
        r$metadata <- NULL
        r$metadata_map <- NULL
        return()
      }

      r$metadata <- meta_data
      r$metadata_map <- list(
        sample = input$sample_column,
        file = input$file_column,
        group = if (isTRUE(nzchar(input$group_column))) {
          input$group_column
        } else {
          NA_character_
        }
      )
    })

    output$n_samples <- shiny::renderText({
      if (is.null(local_r$meta_data)) "-" else nrow(local_r$meta_data)
    })

    output$n_groups <- shiny::renderText({
      map <- r$metadata_map

      if (is.null(map) || is.na(map$group)) {
        return("-")
      }

      length(unique(r$metadata[[map$group]]))
    })

    output$status <- shiny::renderText({
      if (is.null(local_r$meta_data)) {
        "No file"
      } else if (length(local_r$problems) > 0) {
        "Invalid"
      } else {
        "Ready"
      }
    })

    output$problems <- shiny::renderUI({
      shiny::req(length(local_r$problems) > 0)

      htmltools::tags$div(
        class = "alert alert-warning",
        htmltools::tags$b("The meta data can not be used yet:"),
        htmltools::tags$ul(
          lapply(X = local_r$problems, FUN = htmltools::tags$li)
        )
      )
    })

    output$table <- DT::renderDT({
      shiny::req(local_r$meta_data)

      DT::datatable(
        data = local_r$meta_data,
        rownames = FALSE,
        filter = "top",
        fillContainer = TRUE,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "tip"
        )
      )
    })
  })
}
