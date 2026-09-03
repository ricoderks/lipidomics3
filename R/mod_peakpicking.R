#' Peak picking UI Function
#'
#' The user interface of the module that performs the chromatographic peak
#' picking with the centWave algorithm of `{xcms}`.
#'
#' @param id Character(1), internal parameter for `{shiny}`.
#'
#' @returns A `shiny.tag.list` with the user interface of the module.
#'
#' @importFrom shiny NS numericInput selectInput checkboxInput actionButton
#'   icon uiOutput textOutput plotOutput
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel card
#'   card_header layout_column_wrap value_box
#' @importFrom bsicons bs_icon
#' @importFrom DT DTOutput
#' @noRd
mod_peakpicking_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      title = "Peak picking",
      bslib::accordion(
        open = "centwave",
        bslib::accordion_panel(
          title = "centWave",
          value = "centwave",
          icon = bsicons::bs_icon("sliders"),
          shiny::numericInput(
            inputId = ns("ppm"),
            label = "ppm",
            value = 25,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("peakwidth_min"),
            label = "Minimum peak width [s]",
            value = 4,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("peakwidth_max"),
            label = "Maximum peak width [s]",
            value = 30,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("snthresh"),
            label = "Signal to noise threshold",
            value = 10,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("prefilter_k"),
            label = "Prefilter: number of scans",
            value = 3,
            min = 1,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("prefilter_i"),
            label = "Prefilter: intensity",
            value = 1000,
            min = 0,
            step = 100
          ),
          shiny::numericInput(
            inputId = ns("noise"),
            label = "Noise",
            value = 500,
            min = 0,
            step = 100
          ),
          shiny::numericInput(
            inputId = ns("mzdiff"),
            label = "Minimum m/z difference",
            value = -0.001,
            step = 0.001
          ),
          shiny::selectInput(
            inputId = ns("integrate"),
            label = "Integration method",
            choices = c(
              "Mexican hat filter" = "1",
              "Real peak boundaries" = "2"
            ),
            selected = "1"
          ),
          shiny::checkboxInput(
            inputId = ns("fitgauss"),
            label = "Fit a Gaussian to each peak",
            value = FALSE
          )
        ),
        bslib::accordion_panel(
          title = "Merge split peaks",
          value = "merge",
          icon = bsicons::bs_icon("union"),
          shiny::checkboxInput(
            inputId = ns("refine"),
            label = "Merge neighbouring peaks",
            value = TRUE
          ),
          shiny::numericInput(
            inputId = ns("expand_rt"),
            label = "Retention time window [s]",
            value = 2,
            min = 0,
            step = 0.5
          ),
          shiny::numericInput(
            inputId = ns("expand_mz"),
            label = "m/z window",
            value = 0,
            min = 0,
            step = 0.001
          ),
          shiny::numericInput(
            inputId = ns("merge_ppm"),
            label = "m/z window [ppm]",
            value = 10,
            min = 0,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("min_prop"),
            label = "Minimum relative intensity between peaks",
            value = 0.75,
            min = 0,
            max = 1,
            step = 0.05
          )
        ),
        bslib::accordion_panel(
          title = "Performance",
          value = "performance",
          icon = bsicons::bs_icon("cpu"),
          shiny::numericInput(
            inputId = ns("workers"),
            label = "Number of workers",
            value = 2,
            min = 1,
            step = 1
          ),
          shiny::numericInput(
            inputId = ns("chunk_size"),
            label = "Files in memory at the same time",
            value = 2,
            min = 1,
            step = 1
          )
        )
      ),
      shiny::actionButton(
        inputId = ns("start"),
        label = "Start peak picking",
        icon = shiny::icon("play"),
        class = "btn-primary",
        width = "100%"
      )
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      fill = FALSE,
      bslib::value_box(
        title = "Chromatographic peaks",
        value = shiny::textOutput(outputId = ns("n_peaks")),
        showcase = bsicons::bs_icon("graph-up"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Peaks per sample",
        value = shiny::textOutput(outputId = ns("mean_peaks")),
        showcase = bsicons::bs_icon("calculator"),
        theme = "primary"
      ),
      bslib::value_box(
        title = "Processing time",
        value = shiny::textOutput(outputId = ns("elapsed")),
        showcase = bsicons::bs_icon("stopwatch"),
        theme = "primary"
      )
    ),
    shiny::uiOutput(outputId = ns("messages")),
    bslib::layout_column_wrap(
      width = 1 / 2,
      bslib::card(
        full_screen = TRUE,
        height = 400,
        bslib::card_header("Peaks per sample"),
        DT::DTOutput(outputId = ns("per_sample"))
      ),
      bslib::card(
        full_screen = TRUE,
        height = 400,
        bslib::card_header("Peak map"),
        shiny::plotOutput(outputId = ns("peak_map"), height = "100%")
      )
    ),
    bslib::card(
      full_screen = TRUE,
      height = 520,
      bslib::card_header("Chromatographic peaks"),
      DT::DTOutput(outputId = ns("peaks"))
    )
  )
}


#' Peak picking Server Functions
#'
#' The server side of the peak picking module. Runs centWave on the raw data
#' and stores the resulting `XcmsExperiment` object in the shared reactive
#' values object.
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
mod_peakpicking_server <- function(id, r) {
  shiny::moduleServer(id, function(input, output, session) {
    local_r <- shiny::reactiveValues(
      peaks = NULL,
      per_sample = NULL,
      elapsed = NULL
    )

    # A new set of raw data invalidates the peaks that were found before.
    shiny::observeEvent(r$ms_data, {
      local_r$peaks <- NULL
      local_r$per_sample <- NULL
      local_r$elapsed <- NULL
    })

    shiny::observeEvent(input$start, {
      if (is.null(r$ms_data)) {
        shiny::showNotification(
          ui = "Read the raw data first.",
          type = "warning"
        )
        return()
      }

      started <- Sys.time()

      shiny::withProgress(message = "Picking chromatographic peaks", value = 0.1, {
        result <- tryCatch(
          expr = do_peak_picking(
            x = r$ms_data,
            param = centwave_param(
              ppm = input$ppm,
              peakwidth = c(input$peakwidth_min, input$peakwidth_max),
              snthresh = input$snthresh,
              prefilter = c(input$prefilter_k, input$prefilter_i),
              noise = input$noise,
              mzdiff = input$mzdiff,
              integrate = as.integer(input$integrate),
              fitgauss = input$fitgauss
            ),
            ms_level = 1L,
            chunk_size = input$chunk_size,
            refine_param = if (isTRUE(input$refine)) {
              merge_peaks_param(
                expandRt = input$expand_rt,
                expandMz = input$expand_mz,
                ppm = input$merge_ppm,
                minProp = input$min_prop
              )
            } else {
              NULL
            },
            bpparam = parallel_param(workers = input$workers)
          ),
          error = function(e) {
            notify_error(title = "The peak picking failed", cnd = e)
            NULL
          }
        )

        shiny::req(result)
        shiny::incProgress(amount = 0.8, message = "Collecting the peaks")

        r$xcms_data <- result

        # The MS/MS spectra belong to the previous set of peaks.
        r$ms2_spectra <- NULL

        local_r$peaks <- chrom_peaks_table(result)
        local_r$per_sample <- chrom_peaks_per_sample(result)
        local_r$elapsed <- as.numeric(
          difftime(Sys.time(), started, units = "secs")
        )
      })

      shiny::showNotification(
        ui = sprintf("Found %d chromatographic peaks.", nrow(local_r$peaks)),
        type = "message"
      )
    })

    output$n_peaks <- shiny::renderText({
      if (is.null(local_r$peaks)) "-" else format(nrow(local_r$peaks), big.mark = ",")
    })

    output$mean_peaks <- shiny::renderText({
      if (is.null(local_r$per_sample)) {
        return("-")
      }

      format(round(mean(local_r$per_sample$peaks)), big.mark = ",")
    })

    output$elapsed <- shiny::renderText({
      if (is.null(local_r$elapsed)) {
        return("-")
      }

      sprintf("%.0f s", local_r$elapsed)
    })

    output$messages <- shiny::renderUI({
      shiny::req(is.null(r$ms_data))

      htmltools::tags$div(
        class = "alert alert-info",
        "Read the mzML files on the ",
        htmltools::tags$b("Raw data"),
        " page first."
      )
    })

    output$per_sample <- DT::renderDT({
      shiny::req(local_r$per_sample)

      DT::datatable(
        data = local_r$per_sample,
        rownames = FALSE,
        fillContainer = TRUE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )
    })

    output$peak_map <- shiny::renderPlot({
      shiny::req(local_r$peaks)

      plot_chrom_peaks(peaks = local_r$peaks)
    })

    output$peaks <- DT::renderDT({
      shiny::req(local_r$peaks)

      DT::datatable(
        data = local_r$peaks,
        rownames = FALSE,
        filter = "top",
        fillContainer = TRUE,
        options = list(pageLength = 25, scrollX = TRUE, dom = "tip")
      )
    })
  })
}
