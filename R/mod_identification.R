#' Identification UI Function
#'
#' The user interface of the module that matches the extracted MS/MS spectra
#' against the reference spectra of the lipid database.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#'
#' @returns A `shiny.tag.list` with the user interface of the module.
#'
#' @importFrom shiny NS textInput numericInput selectInput checkboxInput
#'   actionButton icon uiOutput textOutput helpText
#' @importFrom plotly plotlyOutput
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel card
#'   card_header layout_column_wrap value_box
#' @importFrom bsicons bs_icon
#' @importFrom DT DTOutput
#' @noRd
mod_identification_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      title = "Identification",
      bslib::accordion(
        open = "database",
        bslib::accordion_panel(
          title = "Database",
          value = "database",
          icon = bsicons::bs_icon("database"),
          shiny::textInput(
            inputId = ns("db_path"),
            label = "Lipid database",
            value = get_golem_config("lipid_db"),
            width = "100%"
          ),
          shiny::numericInput(
            inputId = ns("precursor_ppm"),
            label = "Precursor tolerance [ppm]",
            value = 10,
            min = 0,
            step = 1
          )
        ),
        bslib::accordion_panel(
          title = "Fragment matching",
          value = "matching",
          icon = bsicons::bs_icon("bullseye"),
          shiny::numericInput(
            inputId = ns("tolerance"),
            label = "Fragment tolerance [Da]",
            value = 0.01,
            min = 0,
            step = 0.001
          ),
          shiny::numericInput(
            inputId = ns("ppm"),
            label = "Fragment tolerance [ppm]",
            value = 20,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("min_matched"),
            label = "Minimum matching fragments",
            value = 2,
            min = 1,
            step = 1
          ),
          shiny::checkboxInput(
            inputId = ns("remove_precursor"),
            label = "Ignore the precursor",
            value = TRUE
          ),
          shiny::numericInput(
            inputId = ns("precursor_window"),
            label = "Ignored window below the precursor [Da]",
            value = 1.5,
            min = 0,
            step = 0.5
          ),
          shiny::helpText(
            "Every reference spectrum contains its own precursor, and the",
            "candidates are selected on their precursor m/z, so that peak",
            "always matches and flatters the scores."
          )
        ),
        bslib::accordion_panel(
          title = "Scoring",
          value = "scoring",
          icon = bsicons::bs_icon("calculator"),
          shiny::numericInput(
            inputId = ns("weight_m"),
            label = "Weight: m/z exponent",
            value = 3,
            min = 0,
            step = 0.1
          ),
          shiny::numericInput(
            inputId = ns("weight_n"),
            label = "Weight: intensity exponent",
            value = 0.6,
            min = 0,
            step = 0.1
          ),
          shiny::selectInput(
            inputId = ns("rank_by"),
            label = "Sort the hits by",
            choices = c(
              "Weighted dot product" = "weighted_dot",
              "Dot product" = "dot",
              "Reverse dot product" = "reverse_dot"
            ),
            selected = "weighted_dot"
          ),
          shiny::numericInput(
            inputId = ns("top_n"),
            label = "Hits per spectrum",
            value = 5,
            min = 1,
            step = 1
          )
        )
      ),
      shiny::actionButton(
        inputId = ns("start"),
        label = "Match against the database",
        icon = shiny::icon("play"),
        class = "btn-primary",
        width = "100%"
      ),
      shiny::helpText(
        "Select a row in the table to compare the query spectrum with the",
        "reference spectrum."
      )
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      fill = FALSE,
      bslib::value_box(
        title = "Hits",
        value = shiny::textOutput(outputId = ns("n_hits")),
        showcase = bsicons::bs_icon("search"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Spectra with a hit",
        value = shiny::textOutput(outputId = ns("n_spectra")),
        showcase = bsicons::bs_icon("soundwave"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Peaks with a hit",
        value = shiny::textOutput(outputId = ns("n_peaks")),
        showcase = bsicons::bs_icon("graph-up"),
        theme = "primary"
      )
    ),
    shiny::uiOutput(outputId = ns("messages")),
    bslib::card(
      full_screen = TRUE,
      height = 460,
      bslib::card_header("Hits"),
      DT::DTOutput(outputId = ns("table"))
    ),
    bslib::card(
      full_screen = TRUE,
      height = 440,
      bslib::card_header("Query against reference spectrum"),
      plotly::plotlyOutput(outputId = ns("mirror"), height = "100%")
    )
  )
}


#' Identification Server Functions
#'
#' The server side of the identification module. Scores every extracted MS/MS
#' spectrum against the reference spectra of the lipid database that have a
#' matching precursor m/z.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#' @param r A `reactiveValues` object with the state that is shared between the
#'   modules of the application.
#'
#' @returns Nothing, the module is called for its side effects.
#'
#' @importFrom shiny moduleServer reactiveValues observeEvent req renderText
#'   renderUI showNotification withProgress setProgress
#' @importFrom plotly renderPlotly
#' @importFrom htmltools tags
#' @importFrom DT renderDT datatable
#' @importFrom DBI dbDisconnect
#' @noRd
mod_identification_server <- function(id, r) {
  shiny::moduleServer(id, function(input, output, session) {
    local_r <- shiny::reactiveValues(
      matches = NULL,
      db_info = NULL,
      elapsed = NULL
    )

    # New MS/MS spectra invalidate the hits that were found before.
    shiny::observeEvent(r$ms2_spectra, {
      local_r$matches <- NULL
      local_r$elapsed <- NULL
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$start, {
      if (is.null(r$ms2_spectra) || length(r$ms2_spectra) == 0) {
        shiny::showNotification(
          ui = "Extract the MS/MS spectra first.",
          type = "warning"
        )
        return()
      }

      started <- Sys.time()

      shiny::withProgress(message = "Matching against the lipid database", value = 0.02, {
        con <- tryCatch(
          expr = lipid_db_connect(path = input$db_path),
          error = function(e) {
            notify_error(title = "Could not open the lipid database", cnd = e)
            NULL
          }
        )

        shiny::req(con)
        on.exit(DBI::dbDisconnect(con), add = TRUE)

        local_r$db_info <- lipid_db_info(con)

        matches <- tryCatch(
          expr = match_ms2_spectra(
            sps = r$ms2_spectra,
            spectra_info = ms2_spectra_table(sps = r$ms2_spectra, x = r$xcms_data),
            con = con,
            precursor_ppm = input$precursor_ppm,
            tolerance = input$tolerance,
            ppm = input$ppm,
            m = input$weight_m,
            n = input$weight_n,
            min_matched = input$min_matched,
            top_n = input$top_n,
            rank_by = input$rank_by,
            remove_precursor = input$remove_precursor,
            precursor_window = input$precursor_window,
            progress = function(fraction) {
              shiny::setProgress(value = fraction)
            }
          ),
          error = function(e) {
            notify_error(title = "The database search failed", cnd = e)
            NULL
          }
        )

        shiny::req(matches)

        local_r$matches <- matches
        local_r$elapsed <- as.numeric(
          difftime(Sys.time(), started, units = "secs")
        )
        r$matches <- matches
      })

      shiny::showNotification(
        ui = sprintf(
          "Found %d hit(s) for %d MS/MS spectra in %.0f s.",
          nrow(local_r$matches),
          length(unique(local_r$matches$spectrum)),
          local_r$elapsed
        ),
        type = "message"
      )
    })

    output$n_hits <- shiny::renderText({
      if (is.null(local_r$matches)) "-" else format(nrow(local_r$matches), big.mark = ",")
    })

    output$n_spectra <- shiny::renderText({
      if (is.null(local_r$matches)) {
        return("-")
      }

      sprintf(
        "%s / %s",
        format(length(unique(local_r$matches$spectrum)), big.mark = ","),
        format(length(r$ms2_spectra), big.mark = ",")
      )
    })

    output$n_peaks <- shiny::renderText({
      if (is.null(local_r$matches)) "-" else format(
        length(unique(local_r$matches$peak_id)),
        big.mark = ","
      )
    })

    output$messages <- shiny::renderUI({
      if (is.null(r$ms2_spectra) || length(r$ms2_spectra) == 0) {
        return(
          htmltools::tags$div(
            class = "alert alert-info",
            "Extract the MS/MS spectra on the ",
            htmltools::tags$b("MS/MS spectra"),
            " page first."
          )
        )
      }

      if (!isTRUE(file.exists(input$db_path))) {
        return(
          htmltools::tags$div(
            class = "alert alert-warning",
            htmltools::tags$b("The lipid database was not found."),
            htmltools::tags$br(),
            sprintf(
              "Looked for '%s', relative to '%s'.",
              input$db_path,
              getwd()
            )
          )
        )
      }

      shiny::req(local_r$db_info)

      htmltools::tags$div(
        class = "alert alert-light border",
        sprintf(
          "Searched %s reference spectra of %s lipid classes.",
          format(local_r$db_info$n_spectra, big.mark = ","),
          format(local_r$db_info$n_classes, big.mark = ",")
        )
      )
    })

    output$table <- DT::renderDT({
      shiny::req(local_r$matches)

      DT::datatable(
        data = local_r$matches,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        fillContainer = TRUE,
        options = list(pageLength = 15, scrollX = TRUE, dom = "tip")
      )
    })

    output$mirror <- plotly::renderPlotly({
      shiny::req(
        r$ms2_spectra,
        local_r$matches,
        input$table_rows_selected
      )

      hit <- local_r$matches[input$table_rows_selected, , drop = FALSE]

      con <- lipid_db_connect(path = input$db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      reference <- lipid_db_spectrum(con = con, id = hit$library_id)
      shiny::req(reference)

      plot_mirror_spectrum(
        query = as.matrix(
          ms2_peaks_data(sps = r$ms2_spectra, index = hit$spectrum)
        ),
        reference = reference$peaks,
        tolerance = input$tolerance,
        ppm = input$ppm,
        precursor_mz = if (isTRUE(input$remove_precursor)) {
          hit$precursor_mz
        } else {
          NA_real_
        },
        precursor_window = input$precursor_window,
        title = paste(reference$name, reference$adduct),
        subtitle = sprintf(
          paste(
            "dot %.3f | weighted dot %.3f | reverse dot %.3f |",
            "%d matching fragments | %s at %.2f min"
          ),
          hit$dot,
          hit$weighted_dot,
          hit$reverse_dot,
          hit$n_matched,
          hit$peak_id,
          hit$rt
        )
      )
    })
  })
}
