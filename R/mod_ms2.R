#' MS/MS spectra UI Function
#'
#' The user interface of the module that extracts the MS/MS spectra of the
#' chromatographic peaks.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#'
#' @returns A `shiny.tag.list` with the user interface of the module.
#'
#' @importFrom shiny NS selectInput numericInput actionButton icon uiOutput
#'   textOutput plotOutput helpText
#' @importFrom bslib layout_sidebar sidebar card card_header
#'   layout_column_wrap value_box
#' @importFrom bsicons bs_icon
#' @importFrom DT DTOutput
#' @noRd
mod_ms2_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      title = "MS/MS spectra",
      shiny::selectInput(
        inputId = ns("method"),
        label = "Spectra per chromatographic peak",
        choices = c(
          "All" = "all",
          "Closest in retention time" = "closest_rt",
          "Closest in m/z" = "closest_mz",
          "Largest total ion current" = "largest_tic",
          "Largest base peak" = "largest_bpi"
        ),
        selected = "all",
        width = "100%"
      ),
      shiny::numericInput(
        inputId = ns("expand_rt"),
        label = "Expand retention time range [s]",
        value = 0,
        min = 0,
        step = 1
      ),
      shiny::numericInput(
        inputId = ns("expand_mz"),
        label = "Expand m/z range",
        value = 0,
        min = 0,
        step = 0.001
      ),
      shiny::numericInput(
        inputId = ns("ppm"),
        label = "Expand m/z range [ppm]",
        value = 0,
        min = 0,
        step = 1
      ),
      shiny::numericInput(
        inputId = ns("workers"),
        label = "Number of workers",
        value = 2,
        min = 1,
        step = 1
      ),
      shiny::actionButton(
        inputId = ns("extract"),
        label = "Extract MS/MS spectra",
        icon = shiny::icon("play"),
        class = "btn-primary",
        width = "100%"
      ),
      shiny::helpText(
        "Select a row in the table to show the corresponding spectrum."
      )
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      fill = FALSE,
      bslib::value_box(
        title = "MS/MS spectra",
        value = shiny::textOutput(outputId = ns("n_spectra")),
        showcase = bsicons::bs_icon("soundwave"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Peaks with MS/MS",
        value = shiny::textOutput(outputId = ns("n_peaks_with_ms2")),
        showcase = bsicons::bs_icon("check2-circle"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Coverage",
        value = shiny::textOutput(outputId = ns("coverage")),
        showcase = bsicons::bs_icon("pie-chart"),
        theme = "primary"
      )
    ),
    shiny::uiOutput(outputId = ns("messages")),
    bslib::card(
      full_screen = TRUE,
      height = 460,
      bslib::card_header("Extracted MS/MS spectra"),
      DT::DTOutput(outputId = ns("table"))
    ),
    bslib::card(
      full_screen = TRUE,
      height = 440,
      bslib::card_header("Selected spectrum"),
      shiny::plotOutput(outputId = ns("spectrum"), height = "100%")
    )
  )
}


#' MS/MS spectra Server Functions
#'
#' The server side of the MS/MS module. Extracts the MS/MS spectra that belong
#' to the chromatographic peaks and shows the spectrum of the selected row.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#' @param r A `reactiveValues` object with the state that is shared between the
#'   modules of the application.
#'
#' @returns Nothing, the module is called for its side effects.
#'
#' @importFrom shiny moduleServer reactiveValues observeEvent req renderText
#'   renderUI renderPlot showNotification withProgress incProgress
#' @importFrom htmltools tags
#' @importFrom DT renderDT datatable
#' @noRd
mod_ms2_server <- function(id, r) {
  shiny::moduleServer(id, function(input, output, session) {
    local_r <- shiny::reactiveValues(
      table = NULL,
      coverage = NULL
    )

    # New chromatographic peaks invalidate the extracted spectra.
    shiny::observeEvent(r$xcms_data, {
      local_r$table <- NULL
      local_r$coverage <- NULL
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$extract, {
      if (is.null(r$xcms_data)) {
        shiny::showNotification(
          ui = "Pick the chromatographic peaks first.",
          type = "warning"
        )
        return()
      }

      shiny::withProgress(message = "Extracting the MS/MS spectra", value = 0.1, {
        spectra <- tryCatch(
          expr = extract_ms2_spectra(
            x = r$xcms_data,
            method = input$method,
            ms_level = 2L,
            expand_rt = input$expand_rt,
            expand_mz = input$expand_mz,
            ppm = input$ppm,
            bpparam = parallel_param(workers = input$workers)
          ),
          error = function(e) {
            notify_error(title = "Could not extract the MS/MS spectra", cnd = e)
            NULL
          }
        )

        shiny::req(spectra)
        shiny::incProgress(amount = 0.8, message = "Collecting the spectra")

        r$ms2_spectra <- spectra

        # The hits belong to the previous set of MS/MS spectra.
        r$matches <- NULL

        local_r$table <- ms2_spectra_table(sps = spectra, x = r$xcms_data)
        local_r$coverage <- ms2_coverage(sps = spectra, x = r$xcms_data)
      })

      shiny::showNotification(
        ui = sprintf(
          "Extracted %d MS/MS spectra for %d chromatographic peaks.",
          local_r$coverage$n_spectra,
          local_r$coverage$n_peaks_with_ms2
        ),
        type = "message"
      )
    })

    output$n_spectra <- shiny::renderText({
      if (is.null(local_r$coverage)) {
        return("-")
      }

      format(local_r$coverage$n_spectra, big.mark = ",")
    })

    output$n_peaks_with_ms2 <- shiny::renderText({
      if (is.null(local_r$coverage)) {
        return("-")
      }

      sprintf(
        "%s / %s",
        format(local_r$coverage$n_peaks_with_ms2, big.mark = ","),
        format(local_r$coverage$n_peaks, big.mark = ",")
      )
    })

    output$coverage <- shiny::renderText({
      if (is.null(local_r$coverage)) {
        return("-")
      }

      sprintf("%.1f%%", 100 * local_r$coverage$coverage)
    })

    output$messages <- shiny::renderUI({
      if (is.null(r$xcms_data)) {
        return(
          htmltools::tags$div(
            class = "alert alert-info",
            "Pick the chromatographic peaks on the ",
            htmltools::tags$b("Peak picking"),
            " page first."
          )
        )
      }

      shiny::req(local_r$coverage)
      shiny::req(local_r$coverage$n_spectra == 0)

      htmltools::tags$div(
        class = "alert alert-warning",
        "No MS/MS spectra were found for the chromatographic peaks. ",
        "Check that the mzML files contain MS2 spectra and try to expand ",
        "the retention time or m/z range."
      )
    })

    output$table <- DT::renderDT({
      shiny::req(local_r$table)

      DT::datatable(
        data = local_r$table,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        fillContainer = TRUE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "tip")
      )
    })

    output$spectrum <- shiny::renderPlot({
      shiny::req(r$ms2_spectra, local_r$table, input$table_rows_selected)

      selected <- local_r$table[input$table_rows_selected, , drop = FALSE]

      plot_ms2_spectrum(
        peaks = ms2_peaks_data(
          sps = r$ms2_spectra,
          index = selected$spectrum
        ),
        title = sprintf(
          "%s | precursor m/z %.4f | rt %.2f min | %s",
          selected$peak_id,
          selected$precursor_mz,
          selected$rt,
          selected$sample_name
        )
      )
    })
  })
}
