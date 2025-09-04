#' Read GPS CSV files into a single tibble
#'
#' Reads subject-level GPS CSV files from a source directory. Matches filenames
#' against a regex pattern and concatenates them into one tibble.
#'
#' @param pattern Regular expression used to match filenames in the GPS data directory.
#'   Must match files whose names contain a subject ID like "Z001".
#'
#' @details
#' Each CSV must contain at least columns `timestamp`, `lat`, and
#' either `lon` or `lng`. If a `lng` column is present, it is renamed
#' to `lon`. A `subject` column is created from the filename.
#'
#' Progress is displayed with \pkg{cli}. All files are read with
#' \pkg{readr::read_csv()} with `show_col_types = FALSE`.
#'
#' @return A tibble with one row per GPS observation, including at least columns
#' `subject`, `timestamp`, `lat`, and `lon` (if present). Additional columns are preserved.
#' 
#' @section SMB source:
#' This function requires an environment variable `GPS_SMB_ROOT` to be set,
#' pointing to the base directory containing GPS CSV files (e.g., a mounted SMB path).
#'
#' @examples
#' \dontrun{
#' # Requires access to SMB share configured in ~/.Renviron
#' gps <- read_gps_csv(pattern = "^Z001.*\\.csv$")
#' dplyr::glimpse(gps)
#' }
#'
#' @export
read_gps_csv <- function(pattern) {
  rlang::check_installed(c("readr","dplyr","purrr","stringr"))
  
  root <- Sys.getenv("GPS_SMB_ROOT", "")
  if (root == "") {
    cli::cli_abort("Env var 'GPS_SMB_ROOT' not set. Add it to your ~/.Renviron (e.g. //myserver/GPS).")
  }
  src_dir <- file.path(root)
  
  # All matching files
  files <- list.files(src_dir, pattern = pattern, full.names = TRUE)
  
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
  for (i in 1:length(files)) {
    subj <- stringr::str_extract(basename(files[i]), "^[A-Za-z]\\d{3}")
    if (is.na(subj)) cli::cli_warn("Couldn't extract subject ID from: {basename(files[i])}")
    
    res <- readr::read_csv(files[i], show_col_types = FALSE, progress = F) |>
      dplyr::mutate(subject = subj) |>
      relocate(subject, .before = id) 
    
    required_cols <- c("timestamp", "lat")
    if (!all(required_cols %in% names(res))) {
      cli::cli_warn("Skipping file with missing required columns: {basename(files[i])}")
      next
    }
    
    if ("lng" %in% names(res) && !"lon" %in% names(res)) {
      res <- dplyr::rename(res, lon = lng)
    }
    
    result <- result |> dplyr::bind_rows(res)
    cli::cli_progress_update(id = pb, inc = 1)
    
  }
  
  result
}
