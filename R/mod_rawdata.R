#' Raw data UI Function
#'
#' The user interface of the module that reads the mzML files and links them to
#' the sample meta data.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#'
#' @returns A `shiny.tag.list` with the user interface of the module.
#'
#' @importFrom shiny NS fileInput actionButton radioButtons helpText uiOutput
#'   textOutput plotOutput
#' @importFrom bslib layout_sidebar sidebar card card_header card_body
#'   layout_column_wrap value_box
#' @importFrom bsicons bs_icon
#' @importFrom htmltools tags
#' @importFrom DT DTOutput
#' @noRd
mod_rawdata_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      title = "Raw data files",
      shiny::fileInput(
        inputId = ns("files"),
        label = "Select mzML files",
        multiple = TRUE,
        accept = c(".mzML", ".mzml"),
        width = "100%"
      ),
      shiny::helpText(
        "The files are uploaded to the server before they can be read.",
        "For large files this takes a while."
      ),
      shiny::actionButton(
        inputId = ns("read"),
        label = "Read raw data",
        icon = shiny::icon("play"),
        class = "btn-primary",
        width = "100%"
      ),
      htmltools::tags$hr(),
      shiny::radioButtons(
        inputId = ns("chrom_type"),
        label = "Chromatogram",
        choices = c(
          "Total ion chromatogram" = "tic",
          "Base peak chromatogram" = "bpc"
        ),
        selected = "tic"
      ),
      shiny::helpText(
        "The total ion chromatogram is read from the file headers.",
        "The base peak chromatogram reads all spectra and is slow."
      )
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      fill = FALSE,
      bslib::value_box(
        title = "mzML files",
        value = shiny::textOutput(outputId = ns("n_files")),
        showcase = bsicons::bs_icon("file-earmark-binary"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Linked to meta data",
        value = shiny::textOutput(outputId = ns("n_matched")),
        showcase = bsicons::bs_icon("link-45deg"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Spectra",
        value = shiny::textOutput(outputId = ns("n_spectra")),
        showcase = bsicons::bs_icon("bar-chart-steps"),
        theme = "primary"
      )
    ),
    shiny::uiOutput(outputId = ns("messages")),
    bslib::card(
      full_screen = TRUE,
      height = 340,
      bslib::card_header("Raw data files"),
      DT::DTOutput(outputId = ns("table"))
    ),
    bslib::card(
      full_screen = TRUE,
      height = 440,
      bslib::card_header("Chromatograms"),
      shiny::plotOutput(outputId = ns("chromatogram"), height = "100%")
    )
  )
}


#' Raw data Server Functions
#'
#' The server side of the raw data module. Stages the uploaded mzML files,
#' links them to the sample meta data and reads them into an `MsExperiment`
#' object.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#' @param r A `reactiveValues` object with the state that is shared between the
#'   modules of the application.
#'
#' @returns Nothing, the module is called for its side effects.
#'
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req
#'   renderText renderUI renderPlot showNotification withProgress incProgress
#' @importFrom htmltools tags
#' @importFrom DT renderDT datatable
#' @importFrom MsExperiment spectra
#' @noRd
mod_rawdata_server <- function(id, r) {
  shiny::moduleServer(id, function(input, output, session) {
    local_r <- shiny::reactiveValues(summary = NULL)

    staged_files <- shiny::reactive({
      shiny::req(input$files)

      stage_upload(
        upload = input$files,
        dir = session_upload_dir(session = session, sub_dir = "mzml")
      )
    })

    # Overview of the uploaded files and the meta data row they belong to.
    file_info <- shiny::reactive({
      files <- staged_files()

      if (is.null(r$metadata)) {
        files$meta_row <- NA_integer_
        return(files)
      }

      link <- link_metadata_files(
        meta_data = r$metadata,
        file_column = r$metadata_map$file,
        file_names = files$name
      )
      files$meta_row <- link$meta_row

      files
    })

    shiny::observeEvent(input$read, {
      shiny::req(staged_files())

      if (is.null(r$metadata)) {
        shiny::showNotification(
          ui = "Read a valid meta data file first.",
          type = "warning"
        )
        return()
      }

      files <- file_info()

      if (all(is.na(files$meta_row))) {
        shiny::showNotification(
          ui = paste(
            "None of the mzML files could be linked to the meta data.",
            "Check the raw data file name column of the meta data."
          ),
          type = "error",
          duration = NULL
        )
        return()
      }

      shiny::withProgress(message = "Reading the raw data", value = 0.1, {
        ms_data <- tryCatch(
          expr = {
            sample_data <- build_sample_data(
              file_info = files,
              meta_data = r$metadata,
              map = r$metadata_map
            )

            shiny::incProgress(amount = 0.4)

            create_ms_experiment(
              paths = files$path,
              sample_data = sample_data
            )
          },
          error = function(e) {
            notify_error(title = "Could not read the raw data", cnd = e)
            NULL
          }
        )

        shiny::req(ms_data)
        shiny::incProgress(amount = 0.4, message = "Summarising the raw data")

        r$ms_data <- ms_data
        r$raw_files <- files

        # The results further down the workflow are no longer valid.
        r$xcms_data <- NULL
        r$ms2_spectra <- NULL

        local_r$summary <- ms_experiment_summary(ms_data)
      })

      shiny::showNotification(
        ui = sprintf("Read %d raw data file(s).", nrow(files)),
        type = "message"
      )
    })

    output$n_files <- shiny::renderText({
      if (is.null(input$files)) "-" else nrow(staged_files())
    })

    output$n_matched <- shiny::renderText({
      if (is.null(input$files)) {
        return("-")
      }

      files <- file_info()

      sprintf("%d / %d", sum(!is.na(files$meta_row)), nrow(files))
    })

    output$n_spectra <- shiny::renderText({
      if (is.null(r$ms_data)) {
        return("-")
      }

      format(length(MsExperiment::spectra(r$ms_data)), big.mark = ",")
    })

    output$messages <- shiny::renderUI({
      if (is.null(r$metadata)) {
        return(
          htmltools::tags$div(
            class = "alert alert-info",
            "Read a valid meta data file on the ",
            htmltools::tags$b("Meta data"),
            " page first."
          )
        )
      }

      shiny::req(input$files)
      files <- file_info()
      unmatched <- files$name[is.na(files$meta_row)]

      shiny::req(length(unmatched) > 0)

      htmltools::tags$div(
        class = "alert alert-warning",
        htmltools::tags$b("Not in the meta data:"),
        htmltools::tags$ul(lapply(X = unmatched, FUN = htmltools::tags$li)),
        "These files will be given the name of the file as sample name."
      )
    })

    output$table <- DT::renderDT({
      if (!is.null(local_r$summary)) {
        return(
          DT::datatable(
            data = local_r$summary,
            rownames = FALSE,
            fillContainer = TRUE,
            options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
          )
        )
      }

      shiny::req(input$files)
      files <- file_info()

      DT::datatable(
        data = data.frame(
          file_name = files$name,
          in_meta_data = ifelse(is.na(files$meta_row), "no", "yes"),
          stringsAsFactors = FALSE
        ),
        rownames = FALSE,
        fillContainer = TRUE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )
    })

    output$chromatogram <- shiny::renderPlot({
      shiny::req(r$ms_data)

      if (identical(input$chrom_type, "bpc")) {
        chrom_data <- shiny::withProgress(
          message = "Reading all spectra for the base peak chromatogram",
          value = 0.5,
          expr = bpc_data(x = r$ms_data, ms_level = 1L)
        )

        return(
          plot_chromatograms(
            chrom_data = chrom_data,
            y_label = "Base peak intensity"
          )
        )
      }

      plot_chromatograms(
        chrom_data = tic_data(x = r$ms_data, ms_level = 1L),
        y_label = "Total ion current"
      )
    })
  })
}
