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
#' @param map A named `list` with the elements `sample`, `file` and `group`,
#'   holding the names of the meta data columns.
#'
#' @returns A `data.frame` with one row per raw data file, with at least the
#'   columns `sample_name`, `sample_group` and `file_name`, followed by all
#'   columns of the meta data.
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

  cbind(
    data.frame(
      sample_name = sample_name,
      sample_group = sample_group,
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
#' @importFrom Spectra rtime intensity
#' @importFrom BiocParallel bpparam
#' @noRd
bpc_data <- function(x, ms_level = 1L, bpparam = BiocParallel::bpparam()) {
  chroms <- xcms::chromatogram(
    x,
    aggregationFun = "max",
    msLevel = as.integer(ms_level),
    BPPARAM = bpparam
  )
  sample_data <- as.data.frame(MsExperiment::sampleData(x))

  do.call(
    what = rbind,
    args = lapply(
      X = seq_len(ncol(chroms)),
      FUN = function(i) {
        chrom <- chroms[1, i]

        data.frame(
          rtime = Spectra::rtime(chrom) / 60,
          intensity = Spectra::intensity(chrom),
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
#' sample.
#'
#' @param chrom_data A `data.frame` as created by [tic_data()] or [bpc_data()].
#' @param y_label Character(1), the label of the y axis.
#'
#' @returns A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw theme element_text
#' @importFrom rlang .data
#' @noRd
plot_chromatograms <- function(chrom_data, y_label = "Intensity") {
  ggplot2::ggplot(
    data = chrom_data,
    mapping = ggplot2::aes(
      x = .data$rtime,
      y = .data$intensity,
      colour = .data$sample_name
    )
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(
      x = "Retention time [min]",
      y = y_label,
      colour = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
