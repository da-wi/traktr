#' Read GPS CSV files into a single tibble
#'
#' Reads subject-level GPS CSV files from a directory, matching filenames
#' against a regex pattern and concatenating them into one tibble.
#'
#' @param pattern Regular expression used to match filenames containing GPS data.
#'   Filenames must contain a subject ID like "Z001".
#' @param dir Optional path to a directory containing the GPS files. If supplied,
#'   this overrides the `GPS_SMB_ROOT` environment variable. If not supplied,
#'   the function attempts to use the `GPS_SMB_ROOT` environment variable.
#'
#' @details
#' Each CSV must contain at least columns `timestamp` and `lat`, and either
#' `lon` or `lng`. If a `lng` column is present, it is renamed to `lon`.
#' A `subject` column is extracted from the filename using a regex pattern.
#'
#' If both `dir` and the `GPS_SMB_ROOT` environment variable are missing or empty,
#' the function will abort with an informative error.
#'
#' Progress is displayed using the \pkg{cli} package. Files are read with
#' [readr::read_csv()] using `show_col_types = FALSE`.
#'
#' @return A tibble with one row per GPS observation. Includes at least the columns
#' `subject`, `timestamp`, `lat`, and (if present) `lon`. Additional columns are preserved.
#'
#' @examples
#' \dontrun{
#' # Using GPS_SMB_ROOT set in ~/.Renviron
#' gps <- read_gps_csv(pattern = "^Z001.*\\.csv$")
#'
#' # Or explicitly specify a directory
#' gps <- read_gps_csv(pattern = "^Z001.*\\.csv$", dir = "/path/to/gps_data")
#'
#' dplyr::glimpse(gps)
#' }
#'
#' @export
read_gps_csv <- function(pattern, dir = NULL) {
  rlang::check_installed(c("readr", "dplyr", "stringr", "cli", "tibble"))
  
  root_env <- Sys.getenv("GPS_SMB_ROOT", "")
  src_dir <- if (!is.null(dir)) {
    dir
  } else if (nzchar(root_env)) {
    root_env
  } else {
    cli::cli_abort(
      "Neither 'dir' is provided nor 'GPS_SMB_ROOT' is set. Please set one of them."
    )
  }
  
  # All matching files
  files <- list.files(src_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    cli::cli_abort("No files matching '{pattern}' found in '{src_dir}'.")
  }
  
  cli::cli_alert_info("Data directory: {src_dir}")
  
  pb <- cli::cli_progress_bar(
    format = paste0(
      "{cli::pb_spin} Reading {basename(files[i])} ",
      "[{cli::pb_current}/{cli::pb_total}]   ETA:{cli::pb_eta}"
    ),
    format_done = paste0(
      "{cli::col_green(cli::symbol$tick)} Read {cli::pb_total} files ",
      "in {cli::pb_elapsed}."
    ),
    total = length(files), clear = F)
  on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  
  result <- tibble::tibble()
  
  for (i in seq_along(files)) {
    #cli::cli_progress_update(id = pb, status = paste0("Reading ", basename(files[i])))
    
    subj <- stringr::str_extract(basename(files[i]), "^[A-Za-z]\\d{3}")
    if (is.na(subj)) {
      cli::cli_warn("Couldn't extract subject ID from: {basename(files[i])}")
    }
    
    res <- readr::read_csv(files[i], show_col_types = FALSE, progress = FALSE) |>
      dplyr::mutate(subject = subj)
    
    # Required columns check (adjust to your schema as needed)
    required_cols <- c("timestamp", "lat")
    if (!all(required_cols %in% names(res))) {
      cli::cli_warn("Skipping file with missing required columns: {basename(files[i])}")
      cli::cli_progress_update(id = pb, inc = 1)
      next
    }
    
    # Normalize lng -> lon if needed
    if ("lng" %in% names(res) && !"lon" %in% names(res)) {
      res$lon <- res$lng
      res <- res[ , !(names(res) %in% "lng")]
    }
    
    # Move 'subject' before 'id' without bare symbols
    if ("id" %in% names(res)) {
      res <- dplyr::relocate(res, dplyr::all_of("subject"),
                             .before = dplyr::all_of("id"))
    } else {
      res <- dplyr::relocate(res, dplyr::all_of("subject"))
    }
    
    result <- dplyr::bind_rows(result, res)
    cli::cli_progress_update(id = pb, inc = 1)
    
  }
  
  result
}


