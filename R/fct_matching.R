#' Remove the precursor from a peak list
#'
#' Removes the peaks at and above the precursor m/z. Every reference spectrum
#' of the lipid database contains its own precursor as a peak, and the
#' candidates are selected on their precursor m/z, so that peak matches by
#' construction. Leaving it in inflates all three scores, and the weighted dot
#' product most of all, because the precursor is the fragment with the highest
#' m/z and therefore gets the largest weight.
#'
#' @param x A `matrix` with the columns `mz` and `intensity`.
#' @param precursor_mz Numeric(1), the precursor m/z of the spectrum.
#' @param window Numeric(1), peaks within this distance below the precursor m/z
#'   are removed as well, since a loss of less than about one dalton is not a
#'   real fragment.
#'
#' @returns A `matrix` with the columns `mz` and `intensity`.
#'
#' @noRd
drop_precursor_peaks <- function(x, precursor_mz, window = 1.5) {
  if (nrow(x) == 0 || is.null(precursor_mz) || is.na(precursor_mz)) {
    return(x)
  }

  x[x[, 1L] <= precursor_mz - window, , drop = FALSE]
}


#' Align the peaks of two spectra
#'
#' Matches the peaks of a query spectrum to the peaks of a reference spectrum.
#' Both peak lists are returned with one row per aligned position, where a row
#' of `NA` means that the spectrum has no peak at that m/z. This alignment is
#' the input of all three similarity scores, so that they are all calculated on
#' exactly the same set of peaks.
#'
#' @param x A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum.
#' @param y A `matrix` with the columns `mz` and `intensity` of the reference
#'   spectrum.
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks, in
#'   ppm. It is added to `tolerance`.
#'
#' @returns A named `list` with the aligned peak lists `x` and `y` and the
#'   logical vector `matched`, which is `TRUE` for the positions where both
#'   spectra have a peak.
#'
#' @importFrom MsCoreUtils join
#' @noRd
align_peaks <- function(x, y, tolerance = 0.01, ppm = 20) {
  empty <- matrix(
    numeric(0),
    ncol = 2L,
    dimnames = list(NULL, c("mz", "intensity"))
  )

  if (nrow(x) == 0 || nrow(y) == 0) {
    return(list(x = empty, y = empty, matched = logical(0)))
  }

  x <- x[order(x[, 1L]), , drop = FALSE]
  y <- y[order(y[, 1L]), , drop = FALSE]

  index <- MsCoreUtils::join(
    x = x[, 1L],
    y = y[, 1L],
    tolerance = tolerance,
    ppm = ppm,
    type = "outer"
  )

  list(
    x = x[index$x, , drop = FALSE],
    y = y[index$y, , drop = FALSE],
    matched = !is.na(index$x) & !is.na(index$y)
  )
}


#' Normalised dot product of two aligned peak lists
#'
#' Calculates the normalised dot product, i.e. the squared cosine of the angle
#' between the two intensity vectors. Every peak is weighted by
#' `mz^m * intensity^n`, so that `m = 0` and `n = 1` give the plain dot product
#' of the intensities. Positions where one of the spectra has no peak count
#' towards the length of the vector of the other spectrum, and therefore lower
#' the score.
#'
#' @param x A `matrix` with the columns `mz` and `intensity`, aligned to `y`.
#' @param y A `matrix` with the columns `mz` and `intensity`, aligned to `x`.
#' @param m Numeric(1), the exponent of the m/z in the weight.
#' @param n Numeric(1), the exponent of the intensity in the weight.
#'
#' @returns Numeric(1) between 0 and 1, or `NA` when one of the spectra has no
#'   intensity at all.
#'
#' @noRd
normalised_dot <- function(x, y, m = 0, n = 1) {
  if (nrow(x) == 0) {
    return(NA_real_)
  }

  weight_x <- x[, 1L]^m * x[, 2L]^n
  weight_y <- y[, 1L]^m * y[, 2L]^n

  denominator <- sum(weight_x^2, na.rm = TRUE) * sum(weight_y^2, na.rm = TRUE)

  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  sum(weight_x * weight_y, na.rm = TRUE)^2 / denominator
}


#' Dot product of a query and a reference spectrum
#'
#' The normalised dot product of the intensities, without any weighting. All
#' peaks of both spectra are taken into account, so peaks that only one of the
#' two spectra has lower the score.
#'
#' @param aligned A `list` as created by [align_peaks()].
#'
#' @returns Numeric(1) between 0 and 1.
#'
#' @noRd
dot_product <- function(aligned) {
  normalised_dot(x = aligned$x, y = aligned$y, m = 0, n = 1)
}


#' Weighted dot product of a query and a reference spectrum
#'
#' The normalised dot product of the intensities, where every peak is weighted
#' by `mz^m * intensity^n`. Weighting the intensity by a power smaller than one
#' reduces the influence of the most intense peaks, and weighting by the m/z
#' gives the fragments of a higher mass, which are more specific for a
#' structure, more influence. The default `m = 3` and `n = 0.6` are the
#' exponents of Stein and Scott.
#'
#' @param aligned A `list` as created by [align_peaks()].
#' @param m Numeric(1), the exponent of the m/z in the weight.
#' @param n Numeric(1), the exponent of the intensity in the weight.
#'
#' @returns Numeric(1) between 0 and 1.
#'
#' @noRd
weighted_dot_product <- function(aligned, m = 3, n = 0.6) {
  normalised_dot(x = aligned$x, y = aligned$y, m = m, n = n)
}


#' Reverse dot product of a query and a reference spectrum
#'
#' The normalised dot product calculated over the peaks of the reference
#' spectrum only. Peaks that the query spectrum has and the reference spectrum
#' does not are ignored, peaks that the reference spectrum has and the query
#' spectrum does not still lower the score. This makes the score insensitive to
#' co-eluting compounds that add peaks to the query spectrum.
#'
#' @param aligned A `list` as created by [align_peaks()].
#' @param m Numeric(1), the exponent of the m/z in the weight.
#' @param n Numeric(1), the exponent of the intensity in the weight.
#'
#' @returns Numeric(1) between 0 and 1.
#'
#' @noRd
reverse_dot_product <- function(aligned, m = 0, n = 1) {
  if (nrow(aligned$y) == 0) {
    return(NA_real_)
  }

  keep <- !is.na(aligned$y[, 1L])

  normalised_dot(
    x = aligned$x[keep, , drop = FALSE],
    y = aligned$y[keep, , drop = FALSE],
    m = m,
    n = n
  )
}


#' All similarity scores of a query and a reference spectrum
#'
#' Aligns the two peak lists once and calculates the dot product, the weighted
#' dot product and the reverse dot product on that alignment.
#'
#' @param query A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum.
#' @param reference A `matrix` with the columns `mz` and `intensity` of the
#'   reference spectrum.
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param m Numeric(1), the exponent of the m/z in the weight of the weighted
#'   dot product.
#' @param n Numeric(1), the exponent of the intensity in the weight of the
#'   weighted dot product.
#' @param precursor_mz Numeric(1), the precursor m/z, or `NA` when the
#'   precursor should not be removed. See [drop_precursor_peaks()].
#' @param precursor_window Numeric(1), the window below the precursor m/z that
#'   is removed together with the precursor.
#'
#' @returns A named numeric vector with the elements `dot`, `weighted_dot`,
#'   `reverse_dot` and `n_matched`.
#'
#' @noRd
spectrum_scores <- function(query,
                            reference,
                            tolerance = 0.01,
                            ppm = 20,
                            m = 3,
                            n = 0.6,
                            precursor_mz = NA_real_,
                            precursor_window = 1.5) {
  query <- drop_precursor_peaks(query, precursor_mz, precursor_window)
  reference <- drop_precursor_peaks(reference, precursor_mz, precursor_window)

  aligned <- align_peaks(x = query, y = reference, tolerance = tolerance, ppm = ppm)

  c(
    dot = dot_product(aligned),
    weighted_dot = weighted_dot_product(aligned, m = m, n = n),
    reverse_dot = reverse_dot_product(aligned),
    n_matched = sum(aligned$matched)
  )
}


#' Match one query spectrum against the candidate reference spectra
#'
#' Scores a query spectrum against every candidate reference spectrum and
#' returns the best hits.
#'
#' @param query A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum.
#' @param candidates A `data.frame` as returned by [lipid_db_candidates()].
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param m Numeric(1), the exponent of the m/z in the weight of the weighted
#'   dot product.
#' @param n Numeric(1), the exponent of the intensity in the weight of the
#'   weighted dot product.
#' @param min_matched Integer(1), the minimum number of peaks a query and a
#'   reference spectrum must have in common before the hit is reported. The
#'   reference spectra are small, so a single matching fragment already gives a
#'   high score.
#' @param top_n Integer(1), the number of best hits to report.
#' @param rank_by Character(1), the score to sort the hits by, one of `"dot"`,
#'   `"weighted_dot"` or `"reverse_dot"`.
#' @param precursor_mz Numeric(1), the precursor m/z, or `NA` when the
#'   precursor should not be removed before scoring.
#' @param precursor_window Numeric(1), the window below the precursor m/z that
#'   is removed together with the precursor.
#'
#' @returns A `data.frame` with at most `top_n` rows, ordered from the best to
#'   the worst hit.
#'
#' @noRd
match_spectrum <- function(query,
                           candidates,
                           tolerance = 0.01,
                           ppm = 20,
                           m = 3,
                           n = 0.6,
                           min_matched = 2,
                           top_n = 5,
                           rank_by = "weighted_dot",
                           precursor_mz = NA_real_,
                           precursor_window = 1.5) {
  if (nrow(candidates) == 0 || nrow(query) == 0) {
    return(NULL)
  }

  scores <- vapply(
    X = seq_len(nrow(candidates)),
    FUN = function(i) {
      spectrum_scores(
        query = query,
        reference = decode_library_peaks(
          blob = candidates$peaks[[i]],
          n_peaks = candidates$n_peaks[i]
        ),
        tolerance = tolerance,
        ppm = ppm,
        m = m,
        n = n,
        precursor_mz = precursor_mz,
        precursor_window = precursor_window
      )
    },
    FUN.VALUE = numeric(4)
  )

  keep <- which(scores["n_matched", ] >= min_matched & !is.na(scores[rank_by, ]))

  if (length(keep) == 0) {
    return(NULL)
  }

  keep <- keep[order(scores[rank_by, keep], decreasing = TRUE)]
  keep <- utils::head(keep, top_n)

  data.frame(
    library_id = candidates$id[keep],
    name = candidates$name[keep],
    lipid_class = candidates$lipid_class[keep],
    adduct = candidates$adduct[keep],
    library_mz = candidates$precursor_mz[keep],
    library_rt = candidates$retention_time[keep],
    dot = scores["dot", keep],
    weighted_dot = scores["weighted_dot", keep],
    reverse_dot = scores["reverse_dot", keep],
    n_matched = as.integer(scores["n_matched", keep]),
    n_library_peaks = as.integer(candidates$n_peaks[keep]),
    stringsAsFactors = FALSE
  )
}


#' Read one value from the table of spectra
#'
#' Reads a single value from the table with the information on the extracted
#' MS/MS spectra, and returns `NA` when the column is not there. The columns of
#' that table depend on the spectra variables that the raw files provide.
#'
#' @param spectra_info A `data.frame` as created by [ms2_spectra_table()].
#' @param column Character(1), the name of the column.
#' @param i Integer(1), the row to read.
#'
#' @returns The value of the column, or `NA` when the column is missing.
#'
#' @noRd
spectra_info_column <- function(spectra_info, column, i) {
  if (!column %in% colnames(spectra_info)) {
    return(NA_real_)
  }

  spectra_info[[column]][i]
}


#' Match the extracted MS/MS spectra against the lipid database
#'
#' Looks up, for every extracted MS/MS spectrum, the reference spectra with a
#' matching precursor m/z and scores them. The reference spectra are read from
#' the database one query spectrum at a time, so that the memory use does not
#' depend on the size of the database.
#'
#' @param sps A `Spectra` object with the extracted MS/MS spectra.
#' @param spectra_info A `data.frame` as created by [ms2_spectra_table()], with
#'   one row per spectrum of `sps`.
#' @param con A `SQLiteConnection` object, see [lipid_db_connect()].
#' @param precursor_ppm Numeric(1), the precursor m/z tolerance in ppm used to
#'   select the candidate reference spectra.
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param m Numeric(1), the exponent of the m/z in the weight of the weighted
#'   dot product.
#' @param n Numeric(1), the exponent of the intensity in the weight of the
#'   weighted dot product.
#' @param min_matched Integer(1), the minimum number of peaks a query and a
#'   reference spectrum must have in common.
#' @param top_n Integer(1), the number of best hits per query spectrum.
#' @param rank_by Character(1), the score to sort the hits by.
#' @param remove_precursor Logical(1), whether the precursor should be removed
#'   from both spectra before they are scored. See [drop_precursor_peaks()].
#' @param precursor_window Numeric(1), the window below the precursor m/z that
#'   is removed together with the precursor.
#' @param progress A `function` that is called with the fraction of the work
#'   that is done, or `NULL` when no progress should be reported.
#'
#' @returns A `data.frame` with one row per hit.
#'
#' @importFrom Spectra peaksData polarity precursorMz
#' @noRd
match_ms2_spectra <- function(sps,
                              spectra_info,
                              con,
                              precursor_ppm = 10,
                              tolerance = 0.01,
                              ppm = 20,
                              m = 3,
                              n = 0.6,
                              min_matched = 2,
                              top_n = 5,
                              rank_by = "weighted_dot",
                              remove_precursor = TRUE,
                              precursor_window = 1.5,
                              progress = NULL) {
  if (length(sps) == 0) {
    return(empty_match_table())
  }

  peaks <- Spectra::peaksData(sps)
  precursor <- Spectra::precursorMz(sps)
  ion_mode <- polarity_to_ion_mode(Spectra::polarity(sps))

  hits <- vector(mode = "list", length = length(sps))
  report_every <- max(1L, floor(length(sps) / 20))

  for (i in seq_along(sps)) {
    if (!is.null(progress) && i %% report_every == 0) {
      progress(i / length(sps))
    }

    if (is.na(precursor[i]) || is.na(ion_mode[i])) {
      next
    }

    hit <- match_spectrum(
      query = as.matrix(peaks[[i]]),
      candidates = lipid_db_candidates(
        con = con,
        precursor_mz = precursor[i],
        ion_mode = ion_mode[i],
        ppm = precursor_ppm
      ),
      tolerance = tolerance,
      ppm = ppm,
      m = m,
      n = n,
      min_matched = min_matched,
      top_n = top_n,
      rank_by = rank_by,
      precursor_mz = if (isTRUE(remove_precursor)) precursor[i] else NA_real_,
      precursor_window = precursor_window
    )

    if (is.null(hit)) {
      next
    }

    hits[[i]] <- cbind(
      data.frame(
        spectrum = spectra_info$spectrum[i],
        peak_id = spectra_info$peak_id[i],
        sample_name = spectra_info$sample_name[i],
        precursor_mz = spectra_info$precursor_mz[i],
        rt = spectra_info$rt[i],
        peak_mz = spectra_info_column(spectra_info, "peak_mz", i),
        peak_rt = spectra_info_column(spectra_info, "peak_rt", i),
        stringsAsFactors = FALSE
      ),
      hit
    )
  }

  hits <- do.call(what = rbind, args = hits)

  if (is.null(hits)) {
    return(empty_match_table())
  }

  hits$library_mz <- round(hits$library_mz, 4)
  hits$library_rt <- round(hits$library_rt, 2)
  hits$dot <- round(hits$dot, 4)
  hits$weighted_dot <- round(hits$weighted_dot, 4)
  hits$reverse_dot <- round(hits$reverse_dot, 4)

  hits
}


#' An empty table of hits
#'
#' Creates the table of hits with zero rows, so that the user interface always
#' gets the same columns.
#'
#' @returns A `data.frame` with zero rows.
#'
#' @noRd
empty_match_table <- function() {
  data.frame(
    spectrum = integer(0),
    peak_id = character(0),
    sample_name = character(0),
    precursor_mz = numeric(0),
    rt = numeric(0),
    peak_mz = numeric(0),
    peak_rt = numeric(0),
    library_id = integer(0),
    name = character(0),
    lipid_class = character(0),
    adduct = character(0),
    library_mz = numeric(0),
    library_rt = numeric(0),
    dot = numeric(0),
    weighted_dot = numeric(0),
    reverse_dot = numeric(0),
    n_matched = integer(0),
    n_library_peaks = integer(0),
    stringsAsFactors = FALSE
  )
}


#' Data for the mirror plot of two spectra
#'
#' Prepares the peaks of a query and a reference spectrum for a mirror plot.
#' Both spectra are scaled to their most intense peak, and the intensities of
#' the reference spectrum are made negative so that it is drawn downwards. The
#' same peaks are removed as the scores removed, so that the plot shows what
#' was scored.
#'
#' @param query A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum.
#' @param reference A `matrix` with the columns `mz` and `intensity` of the
#'   reference spectrum.
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param precursor_mz Numeric(1), the precursor m/z, or `NA` when the
#'   precursor should be kept.
#' @param precursor_window Numeric(1), the window below the precursor m/z that
#'   is removed together with the precursor.
#'
#' @returns A `data.frame` with the columns `mz`, `intensity`, `spectrum` and
#'   `matched`, with one row per peak of both spectra.
#'
#' @noRd
mirror_spectrum_data <- function(query,
                                 reference,
                                 tolerance = 0.01,
                                 ppm = 20,
                                 precursor_mz = NA_real_,
                                 precursor_window = 1.5) {
  query <- drop_precursor_peaks(query, precursor_mz, precursor_window)
  reference <- drop_precursor_peaks(reference, precursor_mz, precursor_window)

  aligned <- align_peaks(x = query, y = reference, tolerance = tolerance, ppm = ppm)

  relative <- function(x) {
    if (nrow(x) == 0 || all(is.na(x[, 2L]))) {
      return(x)
    }

    x[, 2L] <- 100 * x[, 2L] / max(x[, 2L], na.rm = TRUE)
    x
  }

  query_peaks <- relative(aligned$x)
  reference_peaks <- relative(aligned$y)

  plot_data <- rbind(
    data.frame(
      mz = query_peaks[, 1L],
      intensity = query_peaks[, 2L],
      spectrum = rep("query", nrow(query_peaks)),
      matched = aligned$matched,
      stringsAsFactors = FALSE
    ),
    data.frame(
      mz = reference_peaks[, 1L],
      intensity = -reference_peaks[, 2L],
      spectrum = rep("reference", nrow(reference_peaks)),
      matched = aligned$matched,
      stringsAsFactors = FALSE
    )
  )

  plot_data[!is.na(plot_data$mz) & !is.na(plot_data$intensity), , drop = FALSE]
}


#' Plot a query spectrum against a reference spectrum
#'
#' Plots the query spectrum upwards and the reference spectrum downwards, both
#' scaled to their most intense peak, so that the two can be compared. The
#' peaks that the two spectra have in common are coloured. The plot is
#' interactive, hovering over a peak shows its m/z and relative intensity.
#'
#' @param query A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum.
#' @param reference A `matrix` with the columns `mz` and `intensity` of the
#'   reference spectrum.
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param precursor_mz Numeric(1), the precursor m/z, or `NA` when the
#'   precursor should be kept. The plot removes the same peaks as the scores
#'   did, so that the plot shows what was scored.
#' @param precursor_window Numeric(1), the window below the precursor m/z that
#'   is removed together with the precursor.
#' @param title Character(1), the title of the plot.
#' @param subtitle Character(1), the subtitle of the plot. Keeping the scores
#'   out of the title stops a long lipid name from being clipped.
#'
#' @returns A `plotly` object.
#'
#' @importFrom plotly plot_ly add_segments add_markers layout config
#' @importFrom htmltools htmlEscape
#' @noRd
plot_mirror_spectrum <- function(query,
                                 reference,
                                 tolerance = 0.01,
                                 ppm = 20,
                                 precursor_mz = NA_real_,
                                 precursor_window = 1.5,
                                 title = "",
                                 subtitle = "") {
  plot_data <- mirror_spectrum_data(
    query = query,
    reference = reference,
    tolerance = tolerance,
    ppm = ppm,
    precursor_mz = precursor_mz,
    precursor_window = precursor_window
  )

  plot_data$label <- ifelse(plot_data$matched, "matched", "not matched")
  plot_data$base <- rep(0, nrow(plot_data))
  plot_data$hover <- sprintf(
    "%s<br>m/z %.4f<br>intensity %.2f%% of the base peak<br>%s",
    plot_data$spectrum,
    plot_data$mz,
    abs(plot_data$intensity),
    plot_data$label
  )

  colours <- c("matched" = "#1b5e7e", "not matched" = "#b0b0b0")

  plot <- plotly::plot_ly()

  for (label in names(colours)) {
    peaks <- plot_data[plot_data$label == label, , drop = FALSE]

    if (nrow(peaks) == 0) {
      next
    }

    plot <- plotly::add_segments(
      p = plot,
      data = peaks,
      x = ~mz,
      xend = ~mz,
      y = ~base,
      yend = ~intensity,
      line = list(color = unname(colours[label]), width = 1.3),
      name = label,
      legendgroup = label,
      hoverinfo = "none"
    )

    # Segments are only hoverable at their ends, so a transparent marker on
    # top of every peak carries the tooltip.
    plot <- plotly::add_markers(
      p = plot,
      data = peaks,
      x = ~mz,
      y = ~intensity,
      marker = list(color = unname(colours[label]), size = 8, opacity = 0),
      name = label,
      legendgroup = label,
      showlegend = FALSE,
      hoverinfo = "text",
      text = ~hover
    )
  }

  heading <- htmltools::htmlEscape(title)

  if (nzchar(subtitle)) {
    heading <- paste0(
      heading,
      "<br><sub>",
      htmltools::htmlEscape(subtitle),
      "</sub>"
    )
  }

  plot <- plotly::layout(
    p = plot,
    title = list(text = heading, x = 0, xanchor = "left", font = list(size = 14)),
    xaxis = list(title = "m/z", zeroline = FALSE),
    yaxis = list(
      title = "Query (up) and reference (down) intensity [%]",
      zeroline = TRUE,
      zerolinecolor = "#666666",
      zerolinewidth = 1
    ),
    hovermode = "closest",
    legend = list(orientation = "h", x = 0, y = -0.18),
    margin = list(t = 60)
  )

  plotly::config(
    p = plot,
    displaylogo = FALSE,
    modeBarButtonsToRemove = list("select2d", "lasso2d", "autoScale2d")
  )
}
