#' Attach Human Footprint statistics to an AOI
#'
#' Extracts Human Footprint (HFP) values for each AOI, computes area-weighted
#' summary statistics (mean, sd, min, max, median) and joins them back to the
#' input data. Either one year is used for all rows, or each row gets the year
#' given in a column (e.g. the year a participant lived at a postcode).
#'
#' @param df A data object carrying the AOI (produced by \code{mutate_aoi_*()}).
#' @param out_dir Directory with \code{hfp_<year>.tif} files. Used if no files
#'   are recorded in \code{attr(df, "hfp_tiles")}.
#' @param year Single year used for all rows. Defaults to the latest available
#'   year. Ignored if \code{year_col} is given.
#' @param year_col Optional name of a column in \code{df} holding the year per
#'   row. Rows are matched to the nearest available year (a message lists any
#'   such substitutions); rows with missing year get \code{NA}.
#' @param deduplicate Logical. If \code{TRUE} (default), AOIs with identical
#'   geometry are extracted only once, then mapped back.
#' @param crs_planar EPSG code of the projected CRS the AOI is transformed to
#'   before extraction. Defaults to Web Mercator (EPSG:3857).
#' @param prefix Prefix for output column names. Default is \code{"hfp_"}.
#' @param progress Logical. If \code{TRUE} (default), show progress messages.
#'
#' @details
#' - Output: \code{hfp_mean}, \code{hfp_sd}, \code{hfp_min}, \code{hfp_max},
#'   \code{hfp_median}, plus \code{hfp_year} (the year actually used).
#' - Resolution is 1 km: small buffers (a few hundred meters) overlap only
#'   one to four cells. Values are still area-weighted by cell coverage.
#'
#' @section Attributes:
#' Attributes of \code{df} are preserved.
#'
#' @return The input \code{df}, augmented with HFP summary columns.
#'
#' @seealso \code{\link{fetch_hfp}}, \code{\link{attach_ndvi}}
#'
#' @examples
#' \dontrun{
#' # one row per participant x childhood postcode, with the year lived there
#' df <- mutate_aoi_buffer(addresses, dist = 1000)
#' df <- fetch_hfp(df, out_dir = "data/hfp", year = unique(addresses$year))
#' df <- attach_hfp(df, year_col = "year")
#' }
#'
#' @export
attach_hfp <- function(df,
                       out_dir     = NULL,
                       year        = NULL,
                       year_col    = NULL,
                       deduplicate = TRUE,
                       crs_planar  = 3857,
                       prefix      = "hfp_",
                       progress    = TRUE) {
  rlang::check_installed(c("sf", "terra", "cli", "exactextractr"))

  if (isTRUE(progress)) {
    cli::cli_progress_step("Summarizing Human Footprint (mean, sd, min, max, median)", spinner = TRUE)
  }

  info   <- .aoi_detect(df)
  id_col <- info$id_col
  avail  <- .hfp_available_files(df, out_dir)

  # Year per row (nearest available year)
  want <- if (!is.null(year_col)) {
    if (!year_col %in% names(df)) cli::cli_abort("Column {.field {year_col}} not found in `df`.")
    as.integer(df[[year_col]])
  } else {
    rep(if (is.null(year)) max(avail$year) else as.integer(year[1]), nrow(df))
  }
  used <- .nearest_year(want, avail$year)

  subst <- unique(stats::na.omit(cbind(want, used)[want != used, , drop = FALSE]))
  if (nrow(subst)) {
    cli::cli_inform(c("i" = "Using nearest available HFP year: {paste(subst[, 1], subst[, 2], sep = ' -> ')}."))
  }

  # Extract per year, only for the AOIs needed in that year
  ids <- df[[id_col]]
  stats <- lapply(sort(unique(stats::na.omit(used))), function(y) {
    need <- unique(stats::na.omit(ids[used %in% y]))
    if (!length(need)) return(NULL)
    s <- .extract_continuous(
      info$aoi, id_col, avail$file[avail$year == y],
      deduplicate = deduplicate,
      crs_planar  = crs_planar,
      ids         = need
    )
    s$year <- y
    s
  })
  stats <- do.call(rbind, stats)

  # Join by id + year actually used
  df$.__hfp_year <- used
  names(stats)[names(stats) == "year"] <- ".__hfp_year"
  out <- .join_stats(df, stats, by = c(id_col, ".__hfp_year"), prefix = prefix)
  names(out)[names(out) == ".__hfp_year"] <- paste0(prefix, "year")

  if (isTRUE(progress)) cli::cli_progress_done()
  invisible(out)
}


#' Available HFP GeoTIFFs, recorded or found in \code{out_dir} (internal)
#'
#' @return A data frame with \code{year} and \code{file}.
#' @keywords internal
#' @noRd
.hfp_available_files <- function(df, out_dir) {
  idx <- attr(df, "hfp_tiles")
  if (!is.null(idx)) {
    idx <- idx[!is.na(idx$file), , drop = FALSE]
    if (!nrow(idx)) cli::cli_abort("No files recorded in attr(df, 'hfp_tiles').")
    return(data.frame(year = idx$year, file = idx$file, stringsAsFactors = FALSE))
  }
  if (is.null(out_dir)) cli::cli_abort("No recorded HFP files; please supply 'out_dir'.")
  files <- list.files(out_dir, "^hfp_\\d{4}\\.tif$", full.names = TRUE)
  if (!length(files)) cli::cli_abort("No hfp_<year>.tif files found in 'out_dir'.")
  data.frame(
    year = as.integer(sub("^hfp_(\\d{4})\\.tif$", "\\1", basename(files))),
    file = files,
    stringsAsFactors = FALSE
  )
}


#' Map years to the nearest available year (internal)
#' @keywords internal
#' @noRd
.nearest_year <- function(want, avail) {
  avail <- sort(unique(avail))
  vapply(want, function(y) {
    if (is.na(y)) return(NA_integer_)
    as.integer(avail[which.min(abs(avail - y))])
  }, integer(1))
}
