#' Determine the type of a meta data file
#'
#' Determines how a meta data file should be read, based on its file
#' extension.
#'
#' @param path Character(1), path to the meta data file.
#'
#' @returns Character(1), one of `"excel"`, `"csv"`, `"tsv"` or `"unknown"`.
#'
#' @importFrom tools file_ext
#' @noRd
metadata_file_type <- function(path) {
  switch(
    tolower(tools::file_ext(path)),
    "xlsx" = ,
    "xlsm" = ,
    "xls" = "excel",
    "csv" = "csv",
    "tsv" = ,
    "tab" = ,
    "txt" = "tsv",
    "unknown"
  )
}


#' List the sheets of a meta data file
#'
#' Lists the worksheet names of an Excel meta data file. Delimited text files
#' do not have sheets and therefore return an empty vector.
#'
#' @param path Character(1), path to the meta data file.
#'
#' @returns A character vector with the sheet names.
#'
#' @importFrom readxl excel_sheets
#' @noRd
metadata_sheets <- function(path) {
  if (metadata_file_type(path) != "excel") {
    return(character(0))
  }

  readxl::excel_sheets(path = path)
}


#' Read a meta data file
#'
#' Reads the sample meta data from an Excel worksheet or from a delimited text
#' file. For text files the delimiter is guessed from the first line.
#'
#' @param path Character(1), path to the meta data file.
#' @param sheet Character(1), name of the worksheet to read. Only used for
#'   Excel files, ignored otherwise.
#'
#' @returns A `data.frame` with the sample meta data.
#'
#' @importFrom readxl read_excel
#' @importFrom utils read.delim
#' @noRd
read_metadata <- function(path, sheet = NULL) {
  type <- metadata_file_type(path)

  meta_data <- switch(
    type,
    "excel" = readxl::read_excel(
      path = path,
      sheet = sheet,
      .name_repair = "unique_quiet"
    ),
    "csv" = utils::read.delim(
      file = path,
      sep = guess_delimiter(path = path),
      check.names = TRUE,
      stringsAsFactors = FALSE
    ),
    "tsv" = utils::read.delim(
      file = path,
      sep = guess_delimiter(path = path),
      check.names = TRUE,
      stringsAsFactors = FALSE
    ),
    stop(
      "Unsupported meta data file type. Use an Excel (.xlsx, .xls) or a ",
      "delimited text (.csv, .tsv, .txt) file.",
      call. = FALSE
    )
  )

  meta_data <- as.data.frame(meta_data, stringsAsFactors = FALSE)

  if (nrow(meta_data) == 0 || ncol(meta_data) == 0) {
    stop("The meta data file does not contain any data.", call. = FALSE)
  }

  meta_data
}


#' Guess the delimiter of a text file
#'
#' Guesses the field delimiter of a delimited text file by counting the
#' candidate delimiters in the first line of the file.
#'
#' @param path Character(1), path to the text file.
#'
#' @returns Character(1), the delimiter with the highest count in the header
#'   line. Defaults to `","` when no candidate is found.
#'
#' @noRd
guess_delimiter <- function(path) {
  first_line <- readLines(con = path, n = 1L, warn = FALSE)

  if (length(first_line) == 0) {
    return(",")
  }

  candidates <- c("\t", ";", ",")
  counts <- vapply(
    X = candidates,
    FUN = function(x) lengths(regmatches(first_line, gregexpr(x, first_line))),
    FUN.VALUE = integer(1)
  )

  if (max(counts) == 0) {
    return(",")
  }

  candidates[which.max(counts)]
}


#' Guess a meta data column
#'
#' Guesses which column of the meta data holds a particular kind of
#' information, by matching the column names against a set of patterns.
#'
#' @param column_names Character vector with the column names of the meta data.
#' @param patterns Character vector with regular expressions, tried in order.
#'
#' @returns Character(1), the name of the first matching column, or the empty
#'   string when nothing matches.
#'
#' @noRd
guess_metadata_column <- function(column_names, patterns) {
  for (pattern in patterns) {
    hit <- grep(pattern = pattern, x = column_names, ignore.case = TRUE)

    if (length(hit) > 0) {
      return(column_names[hit[1]])
    }
  }

  ""
}


#' Validate the meta data
#'
#' Checks whether the selected meta data columns can be used to describe the
#' samples, i.e. whether they exist and contain unique, non missing values.
#'
#' @param meta_data A `data.frame` with the sample meta data.
#' @param sample_column Character(1), name of the column with the sample names.
#' @param file_column Character(1), name of the column with the file names.
#'
#' @returns A character vector with the problems that were found. An empty
#'   vector means that the meta data is valid.
#'
#' @noRd
check_metadata <- function(meta_data, sample_column, file_column) {
  problems <- character(0)

  if (is.null(meta_data)) {
    return("No meta data has been read yet.")
  }

  for (column in c(sample_column, file_column)) {
    if (!nzchar(column) || !column %in% colnames(meta_data)) {
      problems <- c(problems, sprintf("Column '%s' is not available.", column))
    }
  }

  if (length(problems) > 0) {
    return(problems)
  }

  if (anyDuplicated(meta_data[[sample_column]]) > 0) {
    problems <- c(
      problems,
      sprintf("Column '%s' contains duplicated sample names.", sample_column)
    )
  }

  if (anyDuplicated(meta_data[[file_column]]) > 0) {
    problems <- c(
      problems,
      sprintf("Column '%s' contains duplicated file names.", file_column)
    )
  }

  if (anyNA(meta_data[[file_column]]) || any(!nzchar(meta_data[[file_column]]))) {
    problems <- c(
      problems,
      sprintf("Column '%s' contains empty file names.", file_column)
    )
  }

  problems
}


#' Normalise a file name for comparison
#'
#' Reduces a file name to lower case letters and digits separated by single
#' underscores, so that file names that only differ in the characters that
#' conversion tools and file systems tend to replace still compare as equal.
#'
#' @param x Character vector with the file names, without their extension.
#'
#' @returns A character vector with the normalised file names.
#'
#' @noRd
normalise_file_key <- function(x) {
  key <- gsub(pattern = "[^a-z0-9]+", replacement = "_", x = tolower(x))

  gsub(pattern = "^_|_$", replacement = "", x = key)
}


#' Link mzML files to the meta data
#'
#' Matches the names of the mzML files to the file names in the meta data. The
#' comparison ignores the directory and the file extension, so that both
#' `"QC_pool_006"` and `"QC_pool_006.mzML"` are accepted in the meta data.
#'
#' Files that have no exact match are matched a second time on their
#' normalised name, see [normalise_file_key()]. Conversion software and file
#' systems regularly replace spaces and other characters in a file name, which
#' would otherwise break the link with a meta data file that still holds the
#' original name.
#'
#' @param meta_data A `data.frame` with the sample meta data.
#' @param file_column Character(1), name of the column with the file names.
#' @param file_names Character vector with the names of the mzML files.
#'
#' @returns A `data.frame` with one row per mzML file, with the columns
#'   `file_name` and `meta_row`. `meta_row` is the matching row of the meta
#'   data, or `NA` when the file is not described in the meta data.
#'
#' @importFrom tools file_path_sans_ext
#' @noRd
link_metadata_files <- function(meta_data, file_column, file_names) {
  meta_keys <- tools::file_path_sans_ext(
    basename(as.character(meta_data[[file_column]]))
  )
  file_keys <- tools::file_path_sans_ext(basename(file_names))

  meta_row <- match(x = file_keys, table = meta_keys)
  unmatched <- is.na(meta_row)

  if (any(unmatched)) {
    meta_row[unmatched] <- match(
      x = normalise_file_key(file_keys[unmatched]),
      table = normalise_file_key(meta_keys)
    )
  }

  data.frame(
    file_name = file_names,
    meta_row = meta_row,
    stringsAsFactors = FALSE
  )
}
