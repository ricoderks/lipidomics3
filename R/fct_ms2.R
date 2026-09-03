#' Extract the MS/MS spectra of the chromatographic peaks
#'
#' Extracts the MS/MS spectra whose precursor m/z and retention time fall
#' within the m/z and retention time range of a chromatographic peak. The
#' spectra are linked to the peaks through the `chrom_peak_id` spectra
#' variable.
#'
#' @param x An `XcmsExperiment` object with chromatographic peaks.
#' @param method Character(1), which spectra to return per chromatographic
#'   peak. One of `"all"`, `"closest_rt"`, `"closest_mz"`, `"largest_tic"` or
#'   `"largest_bpi"`.
#' @param ms_level Integer(1), the MS level of the spectra to extract.
#' @param expand_rt Numeric(1), the retention time range of the peaks is
#'   expanded by this value, in seconds, on both sides.
#' @param expand_mz Numeric(1), the m/z range of the peaks is expanded by this
#'   value on both sides.
#' @param ppm Numeric(1), the m/z range of the peaks is expanded by this many
#'   ppm on both sides.
#' @param bpparam A `BiocParallelParam` object describing how to parallelise.
#'
#' @returns A `Spectra` object with the MS/MS spectra.
#'
#' @importFrom xcms chromPeakSpectra
#' @importFrom BiocParallel bpparam
#' @noRd
extract_ms2_spectra <- function(x,
                                method = "all",
                                ms_level = 2L,
                                expand_rt = 0,
                                expand_mz = 0,
                                ppm = 0,
                                bpparam = BiocParallel::bpparam()) {
  xcms::chromPeakSpectra(
    object = x,
    method = method,
    msLevel = as.integer(ms_level),
    expandRt = expand_rt,
    expandMz = expand_mz,
    ppm = ppm,
    return.type = "Spectra",
    BPPARAM = bpparam
  )
}


#' Table with the extracted MS/MS spectra
#'
#' Creates a table with one row per extracted MS/MS spectrum. All information
#' comes from the file headers, the peak lists themselves are not read.
#'
#' @param sps A `Spectra` object as created by [extract_ms2_spectra()].
#' @param x An `XcmsExperiment` object, used to look up the sample names.
#'
#' @returns A `data.frame` with one row per MS/MS spectrum.
#'
#' @importFrom Spectra spectraData spectraVariables dataOrigin
#' @importFrom MsExperiment sampleData
#' @noRd
ms2_spectra_table <- function(sps, x) {
  if (length(sps) == 0) {
    return(
      data.frame(
        spectrum = integer(0),
        peak_id = character(0),
        sample_name = character(0),
        precursor_mz = numeric(0),
        rt = numeric(0),
        collision_energy = numeric(0),
        peak_mz = numeric(0),
        peak_rt = numeric(0),
        n_peaks = integer(0),
        tic = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }

  wanted <- c(
    "chrom_peak_id", "chrom_peak_mz", "chrom_peak_rt", "rtime",
    "precursorMz", "collisionEnergy", "peaksCount", "totIonCurrent"
  )
  available <- intersect(wanted, Spectra::spectraVariables(sps))
  spectra_data <- as.data.frame(Spectra::spectraData(sps, columns = available))

  sample_data <- as.data.frame(MsExperiment::sampleData(x))
  index <- match(basename(Spectra::dataOrigin(sps)), sample_data$file_name)

  column <- function(name) {
    if (name %in% colnames(spectra_data)) spectra_data[[name]] else NA
  }

  data.frame(
    spectrum = seq_along(sps),
    peak_id = column("chrom_peak_id"),
    sample_name = sample_data$sample_name[index],
    precursor_mz = round(column("precursorMz"), 4),
    rt = round(column("rtime") / 60, 3),
    collision_energy = column("collisionEnergy"),
    peak_mz = round(column("chrom_peak_mz"), 4),
    peak_rt = round(column("chrom_peak_rt") / 60, 3),
    n_peaks = as.integer(column("peaksCount")),
    tic = signif(column("totIonCurrent"), 4),
    stringsAsFactors = FALSE
  )
}


#' Summarise the MS/MS coverage of the chromatographic peaks
#'
#' Counts for how many of the chromatographic peaks at least one MS/MS spectrum
#' was found.
#'
#' @param sps A `Spectra` object as created by [extract_ms2_spectra()].
#' @param x An `XcmsExperiment` object with the chromatographic peaks.
#'
#' @returns A named `list` with the number of spectra (`n_spectra`), the number
#'   of chromatographic peaks with at least one MS/MS spectrum
#'   (`n_peaks_with_ms2`), the total number of chromatographic peaks
#'   (`n_peaks`) and the coverage as a fraction (`coverage`).
#'
#' @importFrom Spectra spectraData
#' @importFrom xcms chromPeaks
#' @noRd
ms2_coverage <- function(sps, x) {
  n_peaks <- nrow(xcms::chromPeaks(x))

  n_with_ms2 <- if (length(sps) == 0) {
    0L
  } else {
    length(unique(as.data.frame(Spectra::spectraData(
      sps,
      columns = "chrom_peak_id"
    ))$chrom_peak_id))
  }

  list(
    n_spectra = length(sps),
    n_peaks_with_ms2 = n_with_ms2,
    n_peaks = n_peaks,
    coverage = if (n_peaks > 0) n_with_ms2 / n_peaks else 0
  )
}


#' Peak list of a single MS/MS spectrum
#'
#' Reads the m/z and intensity values of a single MS/MS spectrum from disk.
#'
#' @param sps A `Spectra` object as created by [extract_ms2_spectra()].
#' @param index Integer(1), the index of the spectrum to read.
#'
#' @returns A `data.frame` with the columns `mz` and `intensity`.
#'
#' @importFrom Spectra peaksData
#' @noRd
ms2_peaks_data <- function(sps, index) {
  peaks <- as.data.frame(Spectra::peaksData(sps[index])[[1]])

  if (nrow(peaks) == 0) {
    return(data.frame(mz = numeric(0), intensity = numeric(0)))
  }

  peaks
}


#' Plot a single MS/MS spectrum
#'
#' Plots the peak list of an MS/MS spectrum, with the m/z values of the most
#' intense fragments labelled.
#'
#' @param peaks A `data.frame` as created by [ms2_peaks_data()].
#' @param title Character(1), the title of the plot.
#' @param n_labels Integer(1), the maximum number of fragments to label.
#' @param min_rel_intensity Numeric(1), only fragments with at least this
#'   intensity, relative to the most intense fragment, are labelled. This keeps
#'   the labels of the small fragments from overlapping.
#'
#' @returns A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_text labs theme_bw
#'   expansion scale_y_continuous
#' @importFrom rlang .data
#' @noRd
plot_ms2_spectrum <- function(peaks,
                              title = "",
                              n_labels = 10,
                              min_rel_intensity = 0.05) {
  labelled <- peaks[0, , drop = FALSE]

  if (nrow(peaks) > 0) {
    labelled <- peaks[
      peaks$intensity >= min_rel_intensity * max(peaks$intensity), ,
      drop = FALSE
    ]
    labelled <- labelled[
      order(labelled$intensity, decreasing = TRUE)[
        seq_len(min(n_labels, nrow(labelled)))
      ], ,
      drop = FALSE
    ]
  }

  ggplot2::ggplot(
    data = peaks,
    mapping = ggplot2::aes(x = .data$mz, y = .data$intensity)
  ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(xend = .data$mz, yend = 0),
      linewidth = 0.4,
      colour = "#1b5e7e"
    ) +
    ggplot2::geom_text(
      data = labelled,
      mapping = ggplot2::aes(label = round(.data$mz, 4)),
      vjust = -0.4,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(x = "m/z", y = "Intensity", title = title) +
    ggplot2::theme_bw()
}
