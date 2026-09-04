#' Remove the peaks without a signal
#'
#' Removes the entries of a peak list that carry no signal. The mzML files that
#' `msconvert` writes for the SCIEX data keep a zero intensity point on either
#' side of every centroid, so two thirds of the peaks of a query spectrum are
#' empty. Those points sit within the matching tolerance of the centroid they
#' belong to, and since they are closer to it than half a tolerance they are
#' matched to a reference peak in its place. The match then contributes a zero
#' to every score, which is why all scores of a spectrum can be zero while it
#' has matching fragments.
#'
#' @param x A `matrix` with the columns `mz` and `intensity`.
#'
#' @returns A `matrix` with the columns `mz` and `intensity`.
#'
#' @noRd
drop_empty_peaks <- function(x) {
  if (nrow(x) == 0) {
    return(x)
  }

  keep <- !is.na(x[, 1L]) & !is.na(x[, 2L]) & x[, 2L] > 0

  x[keep, , drop = FALSE]
}


#' Prepare a peak list for the alignment
#'
#' Removes the peaks without a signal and sorts the remaining peaks by their
#' m/z, which is what [align_peaks()] and the scores expect. A query spectrum
#' is compared with hundreds of reference spectra, so it is prepared once by
#' [match_spectrum()] instead of again for every candidate.
#'
#' @param x A `matrix` with the columns `mz` and `intensity`.
#'
#' @returns A `matrix` with the columns `mz` and `intensity`, sorted by m/z.
#'
#' @noRd
prepare_peaks <- function(x) {
  x <- drop_empty_peaks(x)

  # Both the raw data files and the database store their peaks by increasing
  # m/z, so the sort is nearly always unnecessary.
  if (nrow(x) < 2L || !is.unsorted(x[, 1L])) {
    return(x)
  }

  x[order(x[, 1L]), , drop = FALSE]
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
#'   spectrum, prepared by [prepare_peaks()].
#' @param y A `matrix` with the columns `mz` and `intensity` of the reference
#'   spectrum, prepared by [prepare_peaks()].
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

  weight_x <- peak_weights(x = x, m = m, n = n)
  weight_y <- peak_weights(x = y, m = m, n = n)

  denominator <- sum(weight_x^2, na.rm = TRUE) * sum(weight_y^2, na.rm = TRUE)

  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  sum(weight_x * weight_y, na.rm = TRUE)^2 / denominator
}


#' The weights of the peaks of a spectrum
#'
#' Calculates `mz^m * intensity^n`, the weight of every peak in the normalised
#' dot product. The plain and the reverse dot product use `m = 0` and `n = 1`,
#' where both powers are the identity, and skipping them saves two thirds of
#' the work of weighting a spectrum.
#'
#' @param x A `matrix` with the columns `mz` and `intensity`.
#' @param m Numeric(1), the exponent of the m/z in the weight.
#' @param n Numeric(1), the exponent of the intensity in the weight.
#'
#' @returns A numeric vector with one weight per peak.
#'
#' @noRd
peak_weights <- function(x, m, n) {
  weight <- if (n == 1) x[, 2L] else x[, 2L]^n

  if (m == 0) {
    return(weight)
  }

  x[, 1L]^m * weight
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
#' @param min_matched Integer(1), the number of matching peaks below which the
#'   three scores are not worth calculating. See [prepared_spectrum_scores()].
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
                            min_matched = 0) {
  prepared_spectrum_scores(
    query = prepare_peaks(query),
    reference = prepare_peaks(reference),
    tolerance = tolerance,
    ppm = ppm,
    m = m,
    n = n,
    min_matched = min_matched
  )
}


#' All similarity scores of two prepared spectra
#'
#' The part of [spectrum_scores()] that does the work, for peak lists that were
#' already prepared by [prepare_peaks()].
#'
#' A candidate that shares too few peaks with the query is dropped by
#' [match_spectrum()] whatever its scores are, and four out of five candidates
#' are, so the scores of such a candidate are not calculated at all and
#' returned as `NA`.
#'
#' @param query A `matrix` with the columns `mz` and `intensity` of the query
#'   spectrum, prepared by [prepare_peaks()].
#' @param reference A `matrix` with the columns `mz` and `intensity` of the
#'   reference spectrum, prepared by [prepare_peaks()].
#' @param tolerance Numeric(1), absolute m/z tolerance for matching two peaks.
#' @param ppm Numeric(1), relative m/z tolerance for matching two peaks.
#' @param m Numeric(1), the exponent of the m/z in the weight of the weighted
#'   dot product.
#' @param n Numeric(1), the exponent of the intensity in the weight of the
#'   weighted dot product.
#' @param min_matched Integer(1), the number of matching peaks below which the
#'   three scores are returned as `NA` instead of being calculated.
#'
#' @returns A named numeric vector with the elements `dot`, `weighted_dot`,
#'   `reverse_dot` and `n_matched`.
#'
#' @noRd
prepared_spectrum_scores <- function(query,
                                     reference,
                                     tolerance = 0.01,
                                     ppm = 20,
                                     m = 3,
                                     n = 0.6,
                                     min_matched = 0) {
  aligned <- align_peaks(x = query, y = reference, tolerance = tolerance, ppm = ppm)
  n_matched <- sum(aligned$matched)

  if (n_matched < min_matched) {
    return(
      c(
        dot = NA_real_,
        weighted_dot = NA_real_,
        reverse_dot = NA_real_,
        n_matched = n_matched
      )
    )
  }

  c(
    dot = dot_product(aligned),
    weighted_dot = weighted_dot_product(aligned, m = m, n = n),
    reverse_dot = reverse_dot_product(aligned),
    n_matched = n_matched
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
                           rank_by = "weighted_dot") {
  if (nrow(candidates) == 0 || nrow(query) == 0) {
    return(NULL)
  }

  # The query is compared with every candidate, so it is prepared once.
  query <- prepare_peaks(query)

  scores <- vapply(
    X = seq_len(nrow(candidates)),
    FUN = function(i) {
      prepared_spectrum_scores(
        query = query,
        reference = prepare_peaks(
          decode_library_peaks(
            blob = candidates$peaks[[i]],
            n_peaks = candidates$n_peaks[i]
          )
        ),
        tolerance = tolerance,
        ppm = ppm,
        m = m,
        n = n,
        min_matched = min_matched
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


#' Split the query spectra into chunks
#'
#' Deals the query spectra out over the chunks one by one, instead of cutting
#' them into blocks. A chunk is the unit of work of a parallel worker and the
#' step in which the progress is reported.
#'
#' The number of candidates of a spectrum runs from a handful to a few thousand
#' and neighbouring spectra elute at the same time, so blocks of neighbours
#' differ by a factor of four in the amount of work they hold, while the chunks
#' of a round robin differ by less than two.
#'
#' @param n Integer(1), the number of query spectra.
#' @param workers Integer(1), the number of parallel workers.
#'
#' @returns A `list` of integer vectors with the indices of the spectra.
#'
#' @noRd
match_chunks <- function(n, workers = 1) {
  # Twenty chunks keep the progress bar moving, four chunks per worker keep the
  # workers busy without paying for a database connection too often, and a
  # small number of spectra gives a chunk per spectrum.
  n_chunks <- min(as.integer(n), max(20L, 4L * as.integer(workers)))

  unname(split(seq_len(n), rep(seq_len(n_chunks), length.out = n)))
}


#' The parallel back end of the database search
#'
#' Describes how the database search should be parallelised, by taking the back
#' end that `{BiocParallel}` uses on this machine and giving it the requested
#' number of workers.
#'
#' @param workers Integer(1), the number of parallel workers.
#'
#' @returns A `BiocParallelParam` object.
#'
#' @importFrom BiocParallel bpparam bpworkers<-
#' @noRd
matching_bpparam <- function(workers) {
  param <- BiocParallel::bpparam()
  BiocParallel::bpworkers(param) <- as.integer(workers)

  param
}


#' A chunk of query spectra
#'
#' Collects everything that is needed to match a chunk of query spectra into a
#' single object. A parallel worker is sent the chunk it has to do and nothing
#' else, so the peak lists of the other chunks are not copied to it.
#'
#' @param indices Integer vector with the query spectra of this chunk.
#' @param peaks A `list` with the peak matrix of every query spectrum.
#' @param precursor Numeric vector with the precursor m/z of every query
#'   spectrum.
#' @param ion_mode Character vector with the ion mode of every query spectrum.
#' @param spectra_info A `data.frame` with one row per query spectrum.
#'
#' @returns A named `list` with the elements `peaks`, `precursor`, `ion_mode`
#'   and `spectra_info` of the spectra of this chunk.
#'
#' @noRd
spectra_chunk <- function(indices, peaks, precursor, ion_mode, spectra_info) {
  list(
    peaks = peaks[indices],
    precursor = precursor[indices],
    ion_mode = ion_mode[indices],
    spectra_info = spectra_info[indices, , drop = FALSE]
  )
}


#' Match a chunk of query spectra against the lipid database
#'
#' Looks up and scores the candidate reference spectra of every query spectrum
#' of one chunk. This is the work that a parallel worker does.
#'
#' @param chunk A `list` as created by [spectra_chunk()].
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
#'
#' @returns A `list` with a `data.frame` of hits per query spectrum, or `NULL`
#'   for the spectra without a hit.
#'
#' @noRd
match_spectra_chunk <- function(chunk,
                                con,
                                precursor_ppm = 10,
                                tolerance = 0.01,
                                ppm = 20,
                                m = 3,
                                n = 0.6,
                                min_matched = 2,
                                top_n = 5,
                                rank_by = "weighted_dot") {
  precursor <- chunk$precursor
  ion_mode <- chunk$ion_mode
  spectra_info <- chunk$spectra_info

  lapply(
    X = seq_along(precursor),
    FUN = function(i) {
      if (is.na(precursor[i]) || is.na(ion_mode[i])) {
        return(NULL)
      }

      hit <- match_spectrum(
        query = as.matrix(chunk$peaks[[i]]),
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
        rank_by = rank_by
      )

      if (is.null(hit)) {
        return(NULL)
      }

      cbind(
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
  )
}


#' Match the query spectra in one process
#'
#' Works through the chunks of query spectra one after the other, over a single
#' connection to the lipid database.
#'
#' @param chunks A `list` of integer vectors, as created by [match_chunks()].
#' @param arguments A named `list` with the other arguments of
#'   [match_spectra_chunk()].
#' @param db_path Character(1), path to the lipid database.
#' @param report A `function` that is called with the number of query spectra
#'   that are done.
#'
#' @returns A `list` with a `data.frame` of hits per query spectrum.
#'
#' @importFrom DBI dbDisconnect
#' @noRd
serial_matches <- function(chunks, arguments, db_path, report) {
  con <- lipid_db_connect(path = db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  done <- 0L
  hits <- vector(mode = "list", length = length(chunks))

  for (k in seq_along(chunks)) {
    hits[[k]] <- do.call(
      what = match_spectra_chunk,
      args = c(list(chunk = chunks[[k]], con = con), arguments)
    )

    done <- done + length(chunks[[k]]$precursor)
    report(done)
  }

  do.call(what = c, args = hits)
}


#' Match the query spectra in several processes
#'
#' Hands the chunks of query spectra to the parallel workers. A worker that
#' finishes a chunk is given the next one straight away, which matters because
#' one chunk can hold several times the work of another.
#'
#' The results are collected in this process as they come in, which is also
#' where the progress is reported from.
#'
#' @param chunks A `list` of integer vectors, as created by [match_chunks()].
#' @param arguments A named `list` with the other arguments of
#'   [match_spectra_chunk()].
#' @param db_path Character(1), path to the lipid database.
#' @param workers Integer(1), the number of parallel workers.
#' @param report A `function` that is called with the number of query spectra
#'   that are done.
#'
#' @returns A `list` with a `data.frame` of hits per query spectrum.
#'
#' @importFrom BiocParallel bpiterate
#' @importFrom DBI dbDisconnect
#' @noRd
parallel_matches <- function(chunks, arguments, db_path, workers, report) {
  sent <- 0L
  done <- 0L

  BiocParallel::bpiterate(
    ITER = function() {
      sent <<- sent + 1L

      if (sent > length(chunks)) NULL else chunks[[sent]]
    },
    FUN = function(chunk, db_path, arguments) {
      # A connection can not be shared between processes, so every worker opens
      # the database itself.
      con <- lipid_db_connect(path = db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      do.call(
        what = match_spectra_chunk,
        args = c(list(chunk = chunk, con = con), arguments)
      )
    },
    db_path = db_path,
    arguments = arguments,
    REDUCE = function(collected, chunk_hits) {
      done <<- done + length(chunk_hits)
      report(done)

      c(collected, chunk_hits)
    },
    init = list(),
    reduce.in.order = FALSE,
    BPPARAM = matching_bpparam(workers = workers)
  )
}


#' Match the extracted MS/MS spectra against the lipid database
#'
#' Looks up, for every extracted MS/MS spectrum, the reference spectra with a
#' matching precursor m/z and scores them. The reference spectra are read from
#' the database one query spectrum at a time, so that the memory use does not
#' depend on the size of the database.
#'
#' The search is spread over `workers` processes, since every query spectrum is
#' scored on its own. A `SQLiteConnection` can not be shared between processes,
#' so every worker opens the database itself, which is why this function takes
#' the path of the database and not a connection.
#'
#' @param sps A `Spectra` object with the extracted MS/MS spectra.
#' @param spectra_info A `data.frame` as created by [ms2_spectra_table()], with
#'   one row per spectrum of `sps`.
#' @param db_path Character(1), path to the lipid database.
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
#' @param workers Integer(1), the number of parallel workers.
#' @param progress A `function` that is called with the fraction of the work
#'   that is done, or `NULL` when no progress should be reported.
#'
#' @returns A `data.frame` with one row per hit.
#'
#' @importFrom Spectra peaksData polarity precursorMz
#' @noRd
match_ms2_spectra <- function(sps,
                              spectra_info,
                              db_path,
                              precursor_ppm = 10,
                              tolerance = 0.01,
                              ppm = 20,
                              m = 3,
                              n = 0.6,
                              min_matched = 2,
                              top_n = 5,
                              rank_by = "weighted_dot",
                              workers = 1,
                              progress = NULL) {
  if (length(sps) == 0) {
    return(empty_match_table())
  }

  peaks <- Spectra::peaksData(sps)
  precursor <- Spectra::precursorMz(sps)
  ion_mode <- polarity_to_ion_mode(Spectra::polarity(sps))

  workers <- max(1L, as.integer(workers))
  chunks <- lapply(
    X = match_chunks(n = length(sps), workers = workers),
    FUN = spectra_chunk,
    peaks = peaks,
    precursor = precursor,
    ion_mode = ion_mode,
    spectra_info = spectra_info
  )
  arguments <- list(
    precursor_ppm = precursor_ppm,
    tolerance = tolerance,
    ppm = ppm,
    m = m,
    n = n,
    min_matched = min_matched,
    top_n = top_n,
    rank_by = rank_by
  )

  report <- function(done) {
    if (!is.null(progress)) {
      progress(done / length(sps))
    }
  }

  hits <- if (workers > 1L) {
    parallel_matches(
      chunks = chunks,
      arguments = arguments,
      db_path = db_path,
      workers = workers,
      report = report
    )
  } else {
    serial_matches(
      chunks = chunks,
      arguments = arguments,
      db_path = db_path,
      report = report
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
#'
#' @returns A `data.frame` with the columns `mz`, `intensity`, `spectrum` and
#'   `matched`, with one row per peak of both spectra.
#'
#' @noRd
mirror_spectrum_data <- function(query,
                                 reference,
                                 tolerance = 0.01,
                                 ppm = 20) {
  query <- prepare_peaks(query)
  reference <- prepare_peaks(reference)

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
                                 title = "",
                                 subtitle = "") {
  plot_data <- mirror_spectrum_data(
    query = query,
    reference = reference,
    tolerance = tolerance,
    ppm = ppm
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
