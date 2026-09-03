#' Summarise the hits per chromatographic peak
#'
#' Collapses the table of hits to one row per chromatographic peak. Several
#' MS/MS spectra can be recorded for one peak and every spectrum can have
#' several hits, so the table of hits usually holds many rows per peak. This
#' summary is what the peak map draws, and it is also the quickest way to see
#' how ambiguous the identification of a peak is.
#'
#' @param matches A `data.frame` with the hits, as created by
#'   [match_ms2_spectra()].
#' @param score_column Character(1), the score that is used to pick the best
#'   hit of a peak.
#'
#' @returns A `data.frame` with one row per chromatographic peak, sorted by
#'   retention time and m/z.
#'
#' @importFrom stats median
#' @noRd
peak_match_summary <- function(matches, score_column = "weighted_dot") {
  if (is.null(matches) || nrow(matches) == 0) {
    return(empty_peak_summary())
  }

  if (!score_column %in% colnames(matches)) {
    score_column <- "weighted_dot"
  }

  matches$peak_id <- as.character(matches$peak_id)
  rt <- peak_coordinate(matches, "peak_rt", "rt")
  mz <- peak_coordinate(matches, "peak_mz", "precursor_mz")
  score <- matches[[score_column]]

  rows <- split(x = seq_len(nrow(matches)), f = matches$peak_id)

  summary <- lapply(
    X = rows,
    FUN = function(i) {
      best <- i[which.max(replace(score[i], is.na(score[i]), -Inf))]

      data.frame(
        peak_id = matches$peak_id[i[1L]],
        sample_name = matches$sample_name[i[1L]],
        peak_rt = stats::median(rt[i], na.rm = TRUE),
        peak_mz = stats::median(mz[i], na.rm = TRUE),
        n_spectra = length(unique(matches$spectrum[i])),
        n_hits = length(i),
        n_lipids = length(unique(matches$name[i])),
        n_classes = length(unique(matches$lipid_class[i])),
        best_score = score[best],
        best_name = matches$name[best],
        stringsAsFactors = FALSE
      )
    }
  )

  summary <- do.call(what = rbind, args = summary)
  rownames(summary) <- NULL

  summary[order(summary$peak_rt, summary$peak_mz), , drop = FALSE]
}


#' The coordinate of a chromatographic peak
#'
#' Reads the retention time or the m/z of the chromatographic peak from the
#' table of hits and falls back on the value of the MS/MS spectrum itself when
#' the raw files do not provide the peak coordinates.
#'
#' @param matches A `data.frame` with the hits.
#' @param column Character(1), the preferred column.
#' @param fallback Character(1), the column that is used where `column` is
#'   missing or `NA`.
#'
#' @returns A `numeric` vector with one value per row of `matches`.
#'
#' @noRd
peak_coordinate <- function(matches, column, fallback) {
  value <- if (column %in% colnames(matches)) {
    as.numeric(matches[[column]])
  } else {
    rep(NA_real_, nrow(matches))
  }

  if (!fallback %in% colnames(matches)) {
    return(value)
  }

  missing <- is.na(value)
  value[missing] <- as.numeric(matches[[fallback]])[missing]

  value
}


#' An empty summary of the hits per peak
#'
#' Creates the summary per chromatographic peak with zero rows, so that the
#' peak map always gets the same columns.
#'
#' @returns A `data.frame` with zero rows.
#'
#' @noRd
empty_peak_summary <- function() {
  data.frame(
    peak_id = character(0),
    sample_name = character(0),
    peak_rt = numeric(0),
    peak_mz = numeric(0),
    n_spectra = integer(0),
    n_hits = integer(0),
    n_lipids = integer(0),
    n_classes = integer(0),
    best_score = numeric(0),
    best_name = character(0),
    stringsAsFactors = FALSE
  )
}


#' Label the chromatographic peak of every hit
#'
#' Creates the label that the table of hits is grouped by, so that all hits of
#' one chromatographic peak are shown together under one heading.
#'
#' @param matches A `data.frame` with the hits, as created by
#'   [match_ms2_spectra()].
#'
#' @returns A `character` vector with one label per row of `matches`.
#'
#' @noRd
peak_group_label <- function(matches) {
  if (is.null(matches) || nrow(matches) == 0) {
    return(character(0))
  }

  rt <- peak_coordinate(matches, "peak_rt", "rt")
  mz <- peak_coordinate(matches, "peak_mz", "precursor_mz")
  n_hits <- as.integer(table(matches$peak_id)[as.character(matches$peak_id)])

  sprintf(
    "%s | m/z %.4f | %.2f min | %s | %d hit(s)",
    matches$peak_id,
    mz,
    rt,
    matches$sample_name,
    n_hits
  )
}


#' Group the hits by chromatographic peak
#'
#' Adds the label of the chromatographic peak to the table of hits and sorts
#' the table so that the hits of one peak follow each other, best hit first.
#'
#' The column `peak_order` is the rank of the chromatographic peak by retention
#' time and m/z. The table of the user interface sorts on it and groups on
#' `peak`, so that the groups stay in the order of the peak map even when the
#' user sorts the hits within a group on another column.
#'
#' @param matches A `data.frame` with the hits, as created by
#'   [match_ms2_spectra()].
#' @param score_column Character(1), the score that the hits of a peak are
#'   sorted by.
#'
#' @returns The `data.frame` with the columns `peak_order` and `peak` added as
#'   its first two columns.
#'
#' @noRd
group_matches_by_peak <- function(matches, score_column = "weighted_dot") {
  if (is.null(matches) || nrow(matches) == 0) {
    return(
      cbind(
        peak_order = integer(0),
        peak = character(0),
        empty_match_table()
      )
    )
  }

  if (!score_column %in% colnames(matches)) {
    score_column <- "weighted_dot"
  }

  label <- peak_group_label(matches)
  rt <- peak_coordinate(matches, "peak_rt", "rt")
  mz <- peak_coordinate(matches, "peak_mz", "precursor_mz")
  peak_order <- match(label, unique(label[order(rt, mz, label)]))

  grouped <- cbind(
    peak_order = peak_order,
    peak = label,
    matches,
    stringsAsFactors = FALSE
  )

  grouped[
    order(grouped$peak_order, -grouped[[score_column]]), ,
    drop = FALSE
  ]
}


#' Draw the map of the chromatographic peaks with a hit
#'
#' Draws every chromatographic peak that has at least one hit as a dot, with
#' the retention time on the x axis and the m/z on the y axis. Clicking a dot
#' selects the peak, which limits the table of hits to that peak.
#'
#' @param peaks A `data.frame` as created by [peak_match_summary()].
#' @param colour_by Character(1), the column of `peaks` that the dots are
#'   coloured by.
#' @param colour_label Character(1), the title of the colour bar.
#' @param selected Character, the identifiers of the chromatographic peaks that
#'   are marked as selected, or `NULL` when nothing is selected.
#' @param source Character(1), the identifier that the click events of the plot
#'   are registered under.
#'
#' @returns A `plotly` object.
#'
#' @importFrom plotly plot_ly add_markers layout config event_register
#' @noRd
plot_peak_map <- function(peaks,
                          colour_by = "best_score",
                          colour_label = NULL,
                          selected = NULL,
                          source = "peak_map") {
  peaks <- as.data.frame(peaks)

  if (!colour_by %in% colnames(peaks)) {
    colour_by <- "best_score"
  }

  if (is.null(colour_label)) {
    colour_label <- colour_by
  }

  hover <- sprintf(
    paste(
      "%s<br>m/z %.4f | %.2f min<br>%s<br>",
      "%d spectrum(s) | %d hit(s)<br>%d lipid(s) in %d class(es)<br>",
      "best: %s (%.3f)"
    ),
    peaks$peak_id,
    peaks$peak_mz,
    peaks$peak_rt,
    peaks$sample_name,
    peaks$n_spectra,
    peaks$n_hits,
    peaks$n_lipids,
    peaks$n_classes,
    peaks$best_name,
    peaks$best_score
  )

  plot <- plotly::plot_ly(source = source)

  plot <- plotly::add_markers(
    p = plot,
    x = peaks$peak_rt,
    y = peaks$peak_mz,
    customdata = peaks$peak_id,
    text = hover,
    hoverinfo = "text",
    name = "peaks",
    showlegend = FALSE,
    marker = list(
      size = 9,
      color = peaks[[colour_by]],
      colorscale = "Viridis",
      showscale = TRUE,
      colorbar = list(title = list(text = colour_label), thickness = 12),
      line = list(width = 0.5, color = "rgba(0, 0, 0, 0.35)")
    )
  )

  # The trace of the selected peak is updated through a proxy when the user
  # clicks a dot, so that the selection can be shown without redrawing the map
  # and losing the zoom. A trace without peaks is not allowed, so the trace is
  # hidden instead when nothing is selected.
  marked <- peaks[peaks$peak_id %in% selected, , drop = FALSE]
  has_selection <- nrow(marked) > 0

  plot <- plotly::add_markers(
    p = plot,
    x = if (has_selection) marked$peak_rt else first_or_zero(peaks$peak_rt),
    y = if (has_selection) marked$peak_mz else first_or_zero(peaks$peak_mz),
    visible = has_selection,
    name = "selected",
    hoverinfo = "skip",
    showlegend = FALSE,
    marker = list(
      size = 18,
      color = "rgba(0, 0, 0, 0)",
      line = list(width = 2.5, color = "#d62728")
    )
  )

  plot <- plotly::layout(
    p = plot,
    xaxis = list(title = "Retention time [min]", zeroline = FALSE),
    yaxis = list(title = "m/z", zeroline = FALSE),
    hovermode = "closest",
    margin = list(t = 20)
  )

  plot <- plotly::config(
    p = plot,
    displaylogo = FALSE,
    modeBarButtonsToRemove = list("select2d", "lasso2d", "autoScale2d")
  )

  plot <- plotly::event_register(p = plot, event = "plotly_click")

  plotly::event_register(p = plot, event = "plotly_doubleclick")
}


#' The title of the colour bar of the peak map
#'
#' Translates the name of the column that the peaks are coloured by into the
#' label that is shown next to the colour bar.
#'
#' @param colour_by Character(1), the name of the column.
#'
#' @returns A `character(1)` with the label.
#'
#' @noRd
map_colour_label <- function(colour_by) {
  labels <- c(
    best_score = "best score",
    n_hits = "hits",
    n_lipids = "lipids",
    n_classes = "classes",
    n_spectra = "spectra"
  )

  if (!isTRUE(colour_by %in% names(labels))) {
    return("best score")
  }

  unname(labels[colour_by])
}


#' The first value of a vector
#'
#' Returns the first value of a vector, or zero when the vector is empty. It is
#' used to give the hidden trace of the peak map a coordinate, since `{plotly}`
#' refuses a trace without any peak.
#'
#' @param x A `numeric` vector.
#'
#' @returns A `numeric(1)`.
#'
#' @noRd
first_or_zero <- function(x) {
  if (length(x) == 0) 0 else x[1L]
}
