#' Create the parameters for the centWave peak picking
#'
#' Collects the settings of the user interface into a `CentWaveParam` object.
#'
#' @param ppm Numeric(1), maximum expected deviation of the m/z values of the
#'   centroids of a chromatographic peak, in ppm.
#' @param peakwidth Numeric(2), expected minimum and maximum peak width in
#'   seconds.
#' @param snthresh Numeric(1), signal to noise ratio cut off.
#' @param prefilter Numeric(2), the number of consecutive scans and the
#'   intensity a mass trace must have to be considered.
#' @param noise Numeric(1), centroids with an intensity below this value are
#'   ignored.
#' @param mzdiff Numeric(1), minimum difference in m/z for peaks with
#'   overlapping retention times.
#' @param integrate Integer(1), peak integration method, `1` uses the mexican
#'   hat filter and `2` the real peak boundaries.
#' @param fitgauss Logical(1), whether a Gaussian should be fitted to each
#'   peak.
#'
#' @returns A `CentWaveParam` object.
#'
#' @importFrom xcms CentWaveParam
#' @noRd
centwave_param <- function(ppm = 25,
                           peakwidth = c(4, 30),
                           snthresh = 10,
                           prefilter = c(3, 1000),
                           noise = 500,
                           mzdiff = -0.001,
                           integrate = 1L,
                           fitgauss = FALSE) {
  xcms::CentWaveParam(
    ppm = ppm,
    peakwidth = peakwidth,
    snthresh = snthresh,
    prefilter = prefilter,
    noise = noise,
    mzdiff = mzdiff,
    integrate = as.integer(integrate),
    fitgauss = fitgauss
  )
}


#' Create the parameters for merging neighbouring peaks
#'
#' Collects the settings of the user interface into a
#' `MergeNeighboringPeaksParam` object. Merging neighbouring peaks repairs
#' chromatographic peaks that centWave split into several parts.
#'
#' @param expandRt Numeric(1), peaks that are closer than twice this value in
#'   retention time, in seconds, are considered for merging.
#' @param expandMz Numeric(1), value by which the m/z range of a peak is
#'   expanded before it is compared to its neighbours.
#' @param ppm Numeric(1), the m/z range of a peak is expanded by this many ppm
#'   before it is compared to its neighbours.
#' @param minProp Numeric(1), minimum intensity between two peaks, relative to
#'   the smallest of the two apex intensities, for them to be merged.
#'
#' @returns A `MergeNeighboringPeaksParam` object.
#'
#' @importFrom xcms MergeNeighboringPeaksParam
#' @noRd
merge_peaks_param <- function(expandRt = 2,
                              expandMz = 0,
                              ppm = 10,
                              minProp = 0.75) {
  xcms::MergeNeighboringPeaksParam(
    expandRt = expandRt,
    expandMz = expandMz,
    ppm = ppm,
    minProp = minProp
  )
}


#' Parallel processing settings
#'
#' Creates the object that describes how the mass spectrometry data should be
#' processed in parallel. A single worker results in serial processing, which
#' makes debugging a lot easier.
#'
#' @param workers Integer(1), the number of workers to use.
#'
#' @returns A `BiocParallelParam` object.
#'
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#' @noRd
parallel_param <- function(workers = 1L) {
  workers <- max(1L, as.integer(workers))

  if (workers == 1L) {
    return(BiocParallel::SerialParam())
  }

  if (.Platform$OS.type == "windows") {
    return(BiocParallel::SnowParam(workers = workers))
  }

  BiocParallel::MulticoreParam(workers = workers)
}


#' Perform the chromatographic peak picking
#'
#' Runs the centWave peak detection on the raw data and, when requested, merges
#' neighbouring peaks afterwards.
#'
#' @param x An `MsExperiment` object with the raw data.
#' @param param A `CentWaveParam` object with the peak picking settings.
#' @param ms_level Integer(1), the MS level to pick peaks in.
#' @param chunk_size Integer(1), the number of files that are loaded into
#'   memory at the same time. Lower this for large files.
#' @param refine_param A `MergeNeighboringPeaksParam` object, or `NULL` when
#'   neighbouring peaks should not be merged.
#' @param bpparam A `BiocParallelParam` object describing how to parallelise.
#'
#' @returns An `XcmsExperiment` object with the chromatographic peaks.
#'
#' @importFrom xcms findChromPeaks refineChromPeaks
#' @importFrom BiocParallel bpparam
#' @noRd
do_peak_picking <- function(x,
                            param,
                            ms_level = 1L,
                            chunk_size = 2L,
                            refine_param = NULL,
                            bpparam = BiocParallel::bpparam()) {
  result <- xcms::findChromPeaks(
    object = x,
    param = param,
    msLevel = as.integer(ms_level),
    chunkSize = as.integer(chunk_size),
    BPPARAM = bpparam
  )

  if (!is.null(refine_param)) {
    result <- xcms::refineChromPeaks(
      object = result,
      param = refine_param,
      msLevel = as.integer(ms_level),
      chunkSize = as.integer(chunk_size),
      BPPARAM = bpparam
    )
  }

  result
}


#' Table with the chromatographic peaks
#'
#' Creates a table with the chromatographic peaks that were found, with the
#' sample name added and the retention times converted to minutes.
#'
#' @param x An `XcmsExperiment` object.
#'
#' @returns A `data.frame` with one row per chromatographic peak.
#'
#' @importFrom xcms chromPeaks
#' @importFrom MsExperiment sampleData
#' @noRd
chrom_peaks_table <- function(x) {
  peaks <- as.data.frame(xcms::chromPeaks(x))
  sample_data <- as.data.frame(MsExperiment::sampleData(x))

  data.frame(
    peak_id = rownames(peaks),
    sample_name = sample_data$sample_name[peaks$sample],
    mz = round(peaks$mz, 4),
    mzmin = round(peaks$mzmin, 4),
    mzmax = round(peaks$mzmax, 4),
    rt = round(peaks$rt / 60, 3),
    rtmin = round(peaks$rtmin / 60, 3),
    rtmax = round(peaks$rtmax / 60, 3),
    into = signif(peaks$into, 4),
    maxo = signif(peaks$maxo, 4),
    sn = round(peaks$sn, 1),
    stringsAsFactors = FALSE
  )
}


#' Number of chromatographic peaks per sample
#'
#' Counts how many chromatographic peaks were found in each sample.
#'
#' @param x An `XcmsExperiment` object.
#'
#' @returns A `data.frame` with the columns `sample_name`, `sample_group` and
#'   `peaks`.
#'
#' @importFrom xcms chromPeaks
#' @importFrom MsExperiment sampleData
#' @noRd
chrom_peaks_per_sample <- function(x) {
  peaks <- as.data.frame(xcms::chromPeaks(x))
  sample_data <- as.data.frame(MsExperiment::sampleData(x))

  data.frame(
    sample_name = sample_data$sample_name,
    sample_group = sample_data$sample_group,
    peaks = as.integer(table(
      factor(peaks$sample, levels = seq_len(nrow(sample_data)))
    )),
    stringsAsFactors = FALSE
  )
}


#' Plot the chromatographic peaks
#'
#' Plots the chromatographic peaks of all samples in a retention time versus
#' m/z map, with the intensity shown as the size of the points.
#'
#' @param peaks A `data.frame` as created by [chrom_peaks_table()].
#'
#' @returns A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous labs
#'   theme_bw theme
#' @importFrom rlang .data
#' @noRd
plot_chrom_peaks <- function(peaks) {
  ggplot2::ggplot(
    data = peaks,
    mapping = ggplot2::aes(
      x = .data$rt,
      y = .data$mz,
      size = .data$maxo,
      colour = .data$sample_name
    )
  ) +
    ggplot2::geom_point(alpha = 0.3, stroke = 0) +
    ggplot2::scale_size_continuous(transform = "log10", guide = "none") +
    ggplot2::labs(
      x = "Retention time [min]",
      y = "m/z",
      colour = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
