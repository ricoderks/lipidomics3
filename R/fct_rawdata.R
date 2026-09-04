#' Build the sample data of the raw data files
#'
#' Combines the uploaded mzML files with the sample meta data into the sample
#' data `data.frame` that describes the samples of an `MsExperiment` object.
#' The rows are in the same order as `file_info`, i.e. in the same order as the
#' raw data files.
#'
#' @param file_info A `data.frame` with the columns `name`, `path` and
#'   `meta_row`, as created by [link_metadata_files()].
#' @param meta_data A `data.frame` with the sample meta data.
#' @param map A named `list` with the elements `sample`, `file`, `group` and
#'   `type`, holding the names of the meta data columns.
#'
#' @returns A `data.frame` with one row per raw data file, with at least the
#'   columns `sample_name`, `sample_group`, `sample_type` and `file_name`,
#'   followed by all columns of the meta data.
#'
#' @importFrom tools file_path_sans_ext
#' @noRd
build_sample_data <- function(file_info, meta_data, map) {
  matched <- meta_data[file_info$meta_row, , drop = FALSE]
  rownames(matched) <- NULL

  sample_name <- as.character(matched[[map$sample]])
  missing_name <- is.na(sample_name)
  sample_name[missing_name] <-
    tools::file_path_sans_ext(file_info$name[missing_name])

  sample_group <- if (is.na(map$group)) {
    rep("all", nrow(file_info))
  } else {
    as.character(matched[[map$group]])
  }
  sample_group[is.na(sample_group)] <- "unknown"

  # Without a sample type column every file is an ordinary sample, so that the
  # steps further down the workflow can always count on the column.
  sample_type <- if (is.null(map$type) || is.na(map$type)) {
    rep("sample", nrow(file_info))
  } else {
    as.character(matched[[map$type]])
  }
  sample_type[is.na(sample_type)] <- "unknown"

  cbind(
    data.frame(
      sample_name = sample_name,
      sample_group = sample_group,
      sample_type = sample_type,
      file_name = file_info$name,
      stringsAsFactors = FALSE
    ),
    matched
  )
}


#' Read the raw data files
#'
#' Reads the mzML files into an `MsExperiment` object. The spectra are not read
#' into memory, only the file headers are, so that also large data sets can be
#' handled.
#'
#' @param paths Character vector with the paths of the mzML files.
#' @param sample_data A `data.frame` with one row per file, as created by
#'   [build_sample_data()].
#'
#' @returns An `MsExperiment` object.
#'
#' @importFrom MsExperiment readMsExperiment
#' @noRd
create_ms_experiment <- function(paths, sample_data) {
  MsExperiment::readMsExperiment(
    spectraFiles = paths,
    sampleData = sample_data
  )
}


#' Summarise the raw data
#'
#' Creates a per sample overview of the raw data, with the number of spectra
#' per MS level and the retention time range. All information is taken from the
#' file headers, no spectra are read.
#'
#' @param x An `MsExperiment` or `XcmsExperiment` object.
#'
#' @returns A `data.frame` with one row per sample.
#'
#' @importFrom MsExperiment spectra sampleData
#' @importFrom Spectra dataOrigin rtime msLevel
#' @noRd
ms_experiment_summary <- function(x) {
  sps <- MsExperiment::spectra(x)
  sample_data <- as.data.frame(MsExperiment::sampleData(x))

  origin <- factor(
    x = basename(Spectra::dataOrigin(sps)),
    levels = sample_data$file_name
  )
  level <- Spectra::msLevel(sps)
  retention_time <- Spectra::rtime(sps)

  data.frame(
    sample_name = sample_data$sample_name,
    sample_group = sample_data$sample_group,
    sample_type = sample_data$sample_type,
    file_name = sample_data$file_name,
    ms1_spectra = as.integer(tapply(level, origin, function(y) sum(y == 1L))),
    ms2_spectra = as.integer(tapply(level, origin, function(y) sum(y == 2L))),
    rt_min = round(as.numeric(tapply(retention_time, origin, min)) / 60, 2),
    rt_max = round(as.numeric(tapply(retention_time, origin, max)) / 60, 2),
    stringsAsFactors = FALSE
  )
}


#' Total ion chromatograms of the raw data
#'
#' Creates the total ion chromatograms of all samples. The intensities are read
#' from the file headers, which makes this fast even for large files.
#'
#' @param x An `MsExperiment` or `XcmsExperiment` object.
#' @param ms_level Integer(1), the MS level to show.
#'
#' @returns A `data.frame` with the columns `rtime` (in minutes), `intensity`,
#'   `sample_name` and `sample_group`.
#'
#' @importFrom MsExperiment spectra sampleData
#' @importFrom Spectra filterMsLevel dataOrigin rtime tic
#' @noRd
tic_data <- function(x, ms_level = 1L) {
  sps <- Spectra::filterMsLevel(
    MsExperiment::spectra(x),
    msLevel = as.integer(ms_level)
  )
  sample_data <- as.data.frame(MsExperiment::sampleData(x))

  index <- match(basename(Spectra::dataOrigin(sps)), sample_data$file_name)

  data.frame(
    rtime = Spectra::rtime(sps) / 60,
    intensity = Spectra::tic(sps, initial = TRUE),
    sample_name = sample_data$sample_name[index],
    sample_group = sample_data$sample_group[index],
    stringsAsFactors = FALSE
  )
}


#' Base peak chromatograms of the raw data
#'
#' Creates the base peak chromatograms of all samples. In contrast to
#' [tic_data()] this reads all spectra of the requested MS level from disk,
#' which is slow for large data sets.
#'
#' @param x An `MsExperiment` or `XcmsExperiment` object.
#' @param ms_level Integer(1), the MS level to show.
#' @param bpparam A `BiocParallelParam` object describing how to parallelise.
#'
#' @returns A `data.frame` with the columns `rtime` (in minutes), `intensity`,
#'   `sample_name` and `sample_group`.
#'
#' @importFrom MsExperiment sampleData
#' @importFrom xcms chromatogram
#' @importFrom BiocParallel bpparam
#' @noRd
bpc_data <- function(x, ms_level = 1L, bpparam = BiocParallel::bpparam()) {
  chromatograms_data(
    chroms = xcms::chromatogram(
      x,
      aggregationFun = "max",
      msLevel = as.integer(ms_level),
      BPPARAM = bpparam
    ),
    sample_data = as.data.frame(MsExperiment::sampleData(x))
  )
}


#' Extracted ion chromatograms of the raw data
#'
#' Creates the extracted ion chromatogram of a single m/z of all samples, by
#' summing the intensities of all peaks that fall within the m/z window of
#' every spectrum. Like [bpc_data()] this reads all spectra of the requested MS
#' level from disk, which is slow for large data sets.
#'
#' @param x An `MsExperiment` or `XcmsExperiment` object.
#' @param mz Numeric(1), the m/z to extract.
#' @param tolerance Numeric(1), the absolute half width of the m/z window.
#' @param ppm Numeric(1), the relative half width of the m/z window.
#' @param ms_level Integer(1), the MS level to show.
#' @param bpparam A `BiocParallelParam` object describing how to parallelise.
#'
#' @returns A `data.frame` with the columns `rtime` (in minutes), `intensity`,
#'   `sample_name` and `sample_group`.
#'
#' @importFrom MsExperiment sampleData
#' @importFrom xcms chromatogram
#' @importFrom BiocParallel bpparam
#' @noRd
eic_data <- function(x,
                     mz,
                     tolerance = 0.01,
                     ppm = 0,
                     ms_level = 1L,
                     bpparam = BiocParallel::bpparam()) {
  chromatograms_data(
    chroms = xcms::chromatogram(
      x,
      mz = mz_window(mz = mz, tolerance = tolerance, ppm = ppm),
      aggregationFun = "sum",
      msLevel = as.integer(ms_level),
      BPPARAM = bpparam
    ),
    sample_data = as.data.frame(MsExperiment::sampleData(x))
  )
}


#' The m/z window around an m/z value
#'
#' Creates the m/z range that is extracted for an extracted ion chromatogram.
#' The half width of the window is the sum of the absolute and the relative
#' tolerance, so that both can be used together or on their own.
#'
#' @param mz Numeric(1), the m/z to extract.
#' @param tolerance Numeric(1), the absolute half width of the m/z window.
#' @param ppm Numeric(1), the relative half width of the m/z window.
#'
#' @returns A `matrix` with one row and the two columns of the m/z range, which
#'   is the shape that [xcms::chromatogram()] expects.
#'
#' @noRd
mz_window <- function(mz, tolerance = 0.01, ppm = 0) {
  half_width <- tolerance + ppm * mz / 1e6

  matrix(c(mz - half_width, mz + half_width), nrow = 1, ncol = 2)
}


#' Convert chromatograms into a data frame
#'
#' Puts the chromatograms of all samples below each other in a single
#' `data.frame`, together with the name and the group of the sample they belong
#' to.
#'
#' @param chroms An `MChromatograms` object with one row and one column per
#'   sample, as returned by [xcms::chromatogram()].
#' @param sample_data A `data.frame` with one row per sample, in the same order
#'   as the columns of `chroms`.
#'
#' @returns A `data.frame` with the columns `rtime` (in minutes), `intensity`,
#'   `sample_name` and `sample_group`.
#'
#' @importFrom Spectra rtime intensity
#' @noRd
chromatograms_data <- function(chroms, sample_data) {
  do.call(
    what = rbind,
    args = lapply(
      X = seq_len(ncol(chroms)),
      FUN = function(i) {
        chrom <- chroms[1, i]
        intensity <- Spectra::intensity(chrom)

        # A spectrum without a peak in the m/z window has no intensity, which
        # is no signal and not a missing value, so the line is drawn through
        # zero instead of being interrupted.
        intensity[is.na(intensity)] <- 0

        data.frame(
          rtime = Spectra::rtime(chrom) / 60,
          intensity = intensity,
          sample_name = sample_data$sample_name[i],
          sample_group = sample_data$sample_group[i],
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


#' Plot chromatograms
#'
#' Plots the chromatograms of all samples in a single panel, coloured by
#' sample. Clicking a sample in the legend hides its chromatogram, which makes
#' it easier to compare a few samples of a large data set.
#'
#' @param chrom_data A `data.frame` as created by [tic_data()], [bpc_data()] or
#'   [eic_data()].
#' @param y_label Character(1), the label of the y axis.
#' @param title Character(1), the title of the plot.
#'
#' @returns A `plotly` object.
#'
#' @importFrom plotly plot_ly add_lines layout config
#' @importFrom htmltools htmlEscape
#' @noRd
plot_chromatograms <- function(chrom_data, y_label = "Intensity", title = "") {
  plot <- plotly::plot_ly()

  # A trace without any point is not allowed, so an empty data set is shown as
  # an empty pair of axes.
  if (nrow(chrom_data) > 0) {
    chrom_data$hover <- sprintf(
      "%s (%s)<br>%.3f min<br>intensity %s",
      chrom_data$sample_name,
      chrom_data$sample_group,
      chrom_data$rtime,
      # `formatC()` formats every intensity on its own, so a chromatogram with
      # a large peak does not push the small ones into scientific notation.
      trimws(
        formatC(
          x = chrom_data$intensity,
          format = "fg",
          digits = 4,
          big.mark = ",",
          drop0trailing = TRUE
        )
      )
    )

    plot <- plotly::add_lines(
      p = plot,
      data = chrom_data,
      x = ~rtime,
      y = ~intensity,
      color = ~sample_name,
      line = list(width = 1.2),
      hoverinfo = "text",
      text = ~hover
    )
  }

  plot <- plotly::layout(
    p = plot,
    title = list(
      text = htmltools::htmlEscape(title),
      x = 0,
      xanchor = "left",
      font = list(size = 14)
    ),
    xaxis = list(title = "Retention time [min]", zeroline = FALSE),
    yaxis = list(title = y_label, zeroline = FALSE),
    hovermode = "closest",
    legend = list(orientation = "h", x = 0, y = -0.18, title = list(text = "")),
    margin = list(t = 40)
  )

  plotly::config(
    p = plot,
    displaylogo = FALSE,
    modeBarButtonsToRemove = list("select2d", "lasso2d", "autoScale2d")
  )
}
