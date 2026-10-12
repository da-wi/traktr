#' Fetch annual Human Footprint rasters
#'
#' Downloads the global annual Human Footprint (HFP) maps of Mu et al. (2022)
#' from figshare (or copies them from a local source) and records the files on
#' the input object for \code{\link{attach_hfp}}. The HFP combines built
#' environment, population density, nighttime lights, croplands, pastures,
#' roads, railways and navigable waterways into one index (0 = no human
#' pressure, 50 = maximum).
#'
#' @param df A data object carrying the AOI (produced by \code{mutate_aoi_*()}).
#'   The AOI itself is not needed for downloading (each year is one global
#'   file) but the files are recorded on \code{df}.
#' @param out_dir Directory where the GeoTIFFs are written
#'   (\code{hfp_<year>.tif}). Will be created if it does not exist.
#' @param year Integer vector of years to fetch. Available years are listed by
#'   the figshare record (2000-2024 at the time of writing).
#' @param src Optional local source: a directory holding \code{hfp<year>.zip}
#'   or \code{hfp<year>.tif} / \code{hfp_<year>.tif} files, as obtained from
#'   figshare. If supplied, nothing is downloaded.
#' @param overwrite Logical. If \code{FALSE} (default), existing GeoTIFFs in
#'   \code{out_dir} are kept.
#' @param progress Logical. If \code{TRUE} (default), show progress messages.
#'
#' @details
#' - Source: figshare record \url{https://doi.org/10.6084/m9.figshare.16571064},
#'   one zip archive (~0.5 GB) per year containing a single GeoTIFF in World
#'   Mollweide projection at 1 km resolution.
#' - The file list is queried from the figshare API, so newly added years are
#'   picked up automatically.
#' - Archives are downloaded to a temporary file, the GeoTIFF is extracted to
#'   \code{out_dir/hfp_<year>.tif} and the archive is deleted.
#' - Reading the zip remotely via GDAL is possible but very slow (one full
#'   decompression per read), hence the download.
#'
#' @section Attributes:
#' - \code{attr(df, "hfp_tiles")}: a tibble with \code{year}, \code{file},
#'   \code{status}.
#' - \code{attr(df, "traktr_layers")$hfp}: same index plus fetch parameters.
#'
#' @return The input \code{df} (invisibly) with the attributes above.
#'
#' @references
#' Mu, H., Li, X., Wen, Y. et al. (2022). A global record of annual terrestrial
#' Human Footprint dataset from 2000 to 2018. \emph{Scientific Data}, 9, 176.
#' \doi{10.1038/s41597-022-01284-8}
#'
#' @seealso \code{\link{attach_hfp}}, \code{\link{fetch_landuse}}
#'
#' @examples
#' \dontrun{
#' df <- mutate_aoi_buffer(my_points, dist = 1000)
#' df <- fetch_hfp(df, out_dir = "data/hfp", year = c(2005, 2010, 2020))
#' df <- attach_hfp(df, year = 2020)
#' }
#'
#' @export
fetch_hfp <- function(df,
                      out_dir,
                      year      = 2020,
                      src       = NULL,
                      overwrite = FALSE,
                      progress  = TRUE) {
  rlang::check_installed(c("curl", "cli", "tibble"))

  year <- sort(unique(as.integer(year)))
  if (!length(year) || anyNA(year)) cli::cli_abort("`year` must contain valid years.")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  status <- character(length(year))
  dests  <- file.path(out_dir, sprintf("hfp_%d.tif", year))

  # ask figshare only for years that actually need downloading
  need <- year[isTRUE(overwrite) | !file.exists(dests)]
  remote <- if (is.null(src) && length(need)) .hfp_remote_files() else NULL
  if (!is.null(remote)) {
    missing <- setdiff(need, remote$year)
    if (length(missing)) {
      cli::cli_abort(c(
        "Year{?s} {missing} not available on figshare.",
        "i" = "Available: {min(remote$year)}-{max(remote$year)}."
      ))
    }
  }

  for (i in seq_along(year)) {
    y <- year[i]
    d <- dests[i]

    if (file.exists(d) && !isTRUE(overwrite)) {
      status[i] <- "exists"
      next
    }
    if (isTRUE(progress)) cli::cli_progress_step("Human Footprint {y}")

    status[i] <- tryCatch({
      if (!is.null(src)) {
        .hfp_from_local(src, y, d)
      } else {
        zip <- tempfile(fileext = ".zip")
        on.exit(unlink(zip), add = TRUE)
        curl::curl_download(remote$url[remote$year == y], zip, mode = "wb", quiet = TRUE)
        .hfp_unzip_tif(zip, d)
        unlink(zip)
        "downloaded"
      }
    }, error = function(e) {
      cli::cli_alert_warning("Human Footprint {y} failed: {conditionMessage(e)}")
      "error"
    })
  }
  if (isTRUE(progress)) cli::cli_progress_done()

  ok <- status != "error"
  idx <- tibble::tibble(
    year   = year,
    file   = ifelse(ok, dests, NA_character_),
    status = status
  )

  attr(df, "hfp_tiles") <- idx
  df <- .set_traktr_layer(
    df, "hfp",
    index  = idx,
    params = list(year = year, out_dir = out_dir, src = src, overwrite = overwrite)
  )
  invisible(df)
}


#' Human Footprint files on figshare (internal)
#'
#' @return A data frame with \code{year} and download \code{url}.
#' @keywords internal
#' @noRd
.hfp_remote_files <- function(article = 16571064) {
  rlang::check_installed("jsonlite")
  api <- sprintf("https://api.figshare.com/v2/articles/%d/files?page_size=100", article)
  resp <- curl::curl_fetch_memory(api)
  if (resp$status_code != 200L) {
    cli::cli_abort("figshare API returned HTTP {resp$status_code}.")
  }
  files <- jsonlite::fromJSON(rawToChar(resp$content))
  hit <- grepl("^hfp\\d{4}\\.zip$", files$name)
  data.frame(
    year = as.integer(sub("^hfp(\\d{4})\\.zip$", "\\1", files$name[hit])),
    url  = files$download_url[hit],
    stringsAsFactors = FALSE
  )
}


#' Extract the single GeoTIFF of an HFP archive to \code{dest} (internal)
#' @keywords internal
#' @noRd
.hfp_unzip_tif <- function(zip, dest) {
  inner <- utils::unzip(zip, list = TRUE)$Name
  tif <- inner[grepl("\\.tif{1,2}$", inner, ignore.case = TRUE)]
  if (length(tif) != 1L) {
    cli::cli_abort("Expected one GeoTIFF in {.file {basename(zip)}}, found {length(tif)}.")
  }
  exdir <- tempfile("hfp_")
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  utils::unzip(zip, files = tif, exdir = exdir, junkpaths = TRUE)
  if (!file.copy(file.path(exdir, basename(tif)), dest, overwrite = TRUE)) {
    cli::cli_abort("Could not write {.file {dest}}.")
  }
  invisible(dest)
}


#' Copy one HFP year from a local directory (internal)
#' @keywords internal
#' @noRd
.hfp_from_local <- function(src, year, dest) {
  cand <- file.path(src, c(sprintf("hfp%d.tif", year), sprintf("hfp_%d.tif", year),
                           sprintf("hfp%d.zip", year)))
  cand <- cand[file.exists(cand)]
  if (!length(cand)) cli::cli_abort("No HFP file for {year} in {.path {src}}.")
  if (grepl("\\.zip$", cand[1])) {
    .hfp_unzip_tif(cand[1], dest)
  } else if (normalizePath(cand[1]) != normalizePath(dest, mustWork = FALSE)) {
    file.copy(cand[1], dest, overwrite = TRUE)
  }
  "copied"
}
