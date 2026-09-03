#' Connect to the lipid database
#'
#' Opens a read only connection to the SQLite database with the reference MS/MS
#' spectra. The connection must be closed by the caller with
#' [DBI::dbDisconnect()].
#'
#' @param path Character(1), path to the SQLite database file.
#'
#' @returns A `SQLiteConnection` object.
#'
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite SQLITE_RO
#' @noRd
lipid_db_connect <- function(path) {
  if (!file.exists(path)) {
    stop(
      "The lipid database was not found at '", path, "'.",
      call. = FALSE
    )
  }

  DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
}


#' Summarise the lipid database
#'
#' Collects the number of reference spectra, lipid classes and ion modes of the
#' lipid database, so that the user can see which database is being used.
#'
#' @param con A `SQLiteConnection` object.
#'
#' @returns A named `list` with the elements `n_spectra`, `n_classes` and
#'   `ion_modes`.
#'
#' @importFrom DBI dbGetQuery
#' @noRd
lipid_db_info <- function(con) {
  list(
    n_spectra = DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM lipid")$n,
    n_classes = DBI::dbGetQuery(
      con,
      "SELECT COUNT(*) AS n FROM lipid_class_summary"
    )$n,
    ion_modes = DBI::dbGetQuery(
      con,
      "SELECT DISTINCT ion_mode FROM lipid WHERE ion_mode IS NOT NULL AND ion_mode != ''"
    )$ion_mode
  )
}


#' Decode the peak list of a reference spectrum
#'
#' The peak lists of the lipid database are stored as a binary blob with, for
#' every peak, the m/z and the intensity as a 4 byte floating point number.
#'
#' @param blob A `raw` vector with the encoded peak list.
#' @param n_peaks Integer(1), the number of peaks in the blob.
#'
#' @returns A `matrix` with the columns `mz` and `intensity`.
#'
#' @noRd
decode_library_peaks <- function(blob, n_peaks) {
  n_peaks <- as.integer(n_peaks)

  if (is.na(n_peaks) || n_peaks < 1L || length(blob) < 8L * n_peaks) {
    return(
      matrix(
        numeric(0),
        ncol = 2L,
        dimnames = list(NULL, c("mz", "intensity"))
      )
    )
  }

  matrix(
    readBin(
      con = blob,
      what = "double",
      n = 2L * n_peaks,
      size = 4L,
      endian = "little"
    ),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mz", "intensity"))
  )
}


#' Reference spectra with a matching precursor
#'
#' Looks up the reference spectra whose precursor m/z lies within a tolerance
#' of the precursor m/z of a query spectrum. The query uses the index on the
#' precursor m/z of the database.
#'
#' @param con A `SQLiteConnection` object.
#' @param precursor_mz Numeric(1), the precursor m/z of the query spectrum.
#' @param ion_mode Character(1), `"Positive"` or `"Negative"`.
#' @param ppm Numeric(1), the precursor m/z tolerance in ppm.
#'
#' @returns A `data.frame` with one row per candidate reference spectrum.
#'
#' @importFrom DBI dbGetQuery
#' @noRd
lipid_db_candidates <- function(con, precursor_mz, ion_mode, ppm = 10) {
  delta <- precursor_mz * ppm / 1e6

  DBI::dbGetQuery(
    con,
    paste(
      "SELECT id, name, lipid_class, precursor_mz, adduct, retention_time,",
      "n_peaks, peaks FROM lipid",
      "WHERE ion_mode = ? AND precursor_mz BETWEEN ? AND ?"
    ),
    params = list(ion_mode, precursor_mz - delta, precursor_mz + delta)
  )
}


#' A single reference spectrum
#'
#' Reads one reference spectrum from the lipid database, by its identifier.
#'
#' @param con A `SQLiteConnection` object.
#' @param id Integer(1), the identifier of the reference spectrum.
#'
#' @returns A `list` with the meta data of the reference spectrum and its peak
#'   list in the element `peaks`, or `NULL` when the identifier is unknown.
#'
#' @importFrom DBI dbGetQuery
#' @noRd
lipid_db_spectrum <- function(con, id) {
  hit <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT id, name, lipid_class, precursor_mz, adduct, ion_mode,",
      "retention_time, formula, inchikey, n_peaks, peaks FROM lipid",
      "WHERE id = ?"
    ),
    params = list(id)
  )

  if (nrow(hit) == 0) {
    return(NULL)
  }

  list(
    id = hit$id[1],
    name = hit$name[1],
    lipid_class = hit$lipid_class[1],
    precursor_mz = hit$precursor_mz[1],
    adduct = hit$adduct[1],
    ion_mode = hit$ion_mode[1],
    retention_time = hit$retention_time[1],
    formula = hit$formula[1],
    inchikey = hit$inchikey[1],
    peaks = decode_library_peaks(hit$peaks[[1]], hit$n_peaks[1])
  )
}


#' Translate a polarity into an ion mode
#'
#' Translates the polarity of a spectrum, as stored in an mzML file, into the
#' ion mode used by the lipid database.
#'
#' @param polarity Integer vector with the polarity, `1` for positive and `0`
#'   for negative.
#'
#' @returns A character vector with `"Positive"`, `"Negative"` or `NA` for every
#'   element of `polarity`.
#'
#' @noRd
polarity_to_ion_mode <- function(polarity) {
  ion_mode <- rep(NA_character_, length(polarity))
  ion_mode[!is.na(polarity) & polarity == 1L] <- "Positive"
  ion_mode[!is.na(polarity) & polarity == 0L] <- "Negative"

  ion_mode
}
