#' Fetch ESA WorldCover 10 m tiles covering an AOI
#'
#' Downloads ESA WorldCover 10 m tiles from AWS S3 for all 3 deg grid cells
#' that intersect the AOI. The function optionally HEAD-checks URLs, honors
#' local caching, and records output metadata on the input object.
#'
#' @param df A data object carrying the AOI. If \code{aoi} is not supplied,
#'   the function reads \code{attr(df, "traktr_aoi")}. The AOI must be an
#'   \code{sf} polygon/multipolygon with a valid CRS.
#' @param out_dir Directory where output GeoTIFF files are written. Will be
#'   created if it does not exist.
#' @param year WorldCover reference year. For example, 2020 (v100) or
#'   2021 (v200).
#' @param version WorldCover product version string, e.g. \code{"v100"} or
#'   \code{"v200"}.
#' @param layer Product layer name. Use \code{"map"} for the 10 m land-cover
#'   classification layer.
#' @param check_remote Logical. If \code{TRUE}, issue HTTP HEAD requests to
#'   verify remote availability before downloading.
#' @param overwrite Logical. If \code{FALSE} (default), existing files are
#'   skipped. If \code{TRUE}, files are re-downloaded.
#' @param progress Logical. If \code{TRUE} (default), show progress messages
#'   and a progress bar.
#' @param aoi Optional \code{sf} polygon/multipolygon. If supplied, overrides
#'   \code{attr(df, "traktr_aoi")}.
#'
#' @param wc_year Deprecated. Use \code{year}.
#' @param wc_version Deprecated. Use \code{version}.
#' @param check_exists_remote Deprecated. Use \code{check_remote}.
#'
#' @details
#' - AOI handling: The AOI is validated and transformed to EPSG:4326.
#' - Tiling: A regular 3 deg grid (aligned to integer multiples of 3 deg) is
#'   generated to cover the AOI bounding box. Only cells that intersect the AOI
#'   are requested. Tile codes follow the lower-left corner convention, e.g.
#'   N48E006 or S03W072.
#' - Caching: Output filenames are the official WorldCover basenames, which
#'   remain stable across runs.
#' - Attributes:
#'   * \code{attr(df, "worldcover_tiles")} stores a tibble with per-tile URL,
#'     file path, status, size, and geometry.
#'   * \code{attr(df, "worldcover_aoi")} stores the AOI used (WGS84).
#'   * \code{attr(df, "traktr_layers")$worldcover} stores the same information
#'     in a harmonized container for downstream discovery.
#'
#' @return The input \code{df}, returned invisibly, augmented with attributes
#'   describing outputs and call parameters (see Details).
#'
#' @seealso
#' \code{\link{fetch_ndvi}}
#'
#' @examples
#' \dontrun{
#' df <- some_traktr_object_with_aoi()
#' df <- fetch_landuse(
#'   df,
#'   out_dir      = "data/worldcover",
#'   year         = 2021,
#'   version      = "v200",
#'   layer        = "map",
#'   check_remote = FALSE,
#'   overwrite    = FALSE,
#'   progress     = TRUE
#' )
#' # Inspect results table
#' attr(df, "traktr_layers")$worldcover$index
#' }
#'
#' @export
fetch_landuse <- function(
    df,
    out_dir,
    year          = 2021,      # temporal selector (2020 "v100" or 2021 "v200")
    version       = "v200",    # "v100" for 2020, "v200" for 2021
    layer         = "map",     # "map" is the 10m land-cover classification
    check_remote  = FALSE,     # optional HEAD checks before download
    overwrite     = FALSE,
    progress      = TRUE,
    aoi           = NULL,
    # ---- deprecated compat args (will be removed in a later release) ----
    wc_year = NULL, wc_version = NULL, check_exists_remote = NULL
) {
 
  
  # -- package prerequisites -------------------------------------------
  rlang::check_installed(c("curl", "cli", "tibble", "sf"))
  
  # -- deprecated args shims -------------------------------------------
  if (!is.null(wc_year)) {
    warning("`wc_year` is deprecated; use `year`.", call. = FALSE)
    year <- wc_year
  }
  if (!is.null(wc_version)) {
    warning("`wc_version` is deprecated; use `version`.", call. = FALSE)
    version <- wc_version
  }
  if (!is.null(check_exists_remote)) {
    warning("`check_exists_remote` is deprecated; use `check_remote`.", call. = FALSE)
    check_remote <- isTRUE(check_exists_remote)
  }
  
  # -- AOI, dirs, and progress -----------------------------------------
  aoi <- .resolve_aoi(df, aoi)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(progress)) {
    cli::cli_progress_step(
      "Fetching WorldCover data from AWS",
      spinner   = TRUE,
      msg_done  = "Fetching WorldCover completed."
    )
  }
  
  # Always work in EPSG:4326 for 3 deg tiling and ensure validity
  aoi_ll <- sf::st_make_valid(sf::st_transform(aoi, 4326))
  
  # Build a 3 degree grid aligned to multiples of 3 deg that covers the AOI bbox
  bb  <- sf::st_bbox(aoi_ll)
  llx <- as.numeric(floor(bb["xmin"] / 3) * 3)
  lly <- as.numeric(floor(bb["ymin"] / 3) * 3)
  urx <- as.numeric(ceiling(bb["xmax"] / 3) * 3)
  ury <- as.numeric(ceiling(bb["ymax"] / 3) * 3)
  
  cover_sfc <- sf::st_as_sfc(sf::st_bbox(c(xmin = llx, ymin = lly, xmax = urx, ymax = ury), crs = 4326))
  grid3 <- sf::st_make_grid(cover_sfc, cellsize = c(3, 3), what = "polygons", square = TRUE) |>
    sf::st_as_sf()
  
  # Keep grid cells intersecting the AOI (union to speed intersection)
  hits <- grid3[sf::st_intersects(grid3, sf::st_union(aoi_ll), sparse = FALSE)[, 1], , drop = FALSE]
  
  if (nrow(hits) == 0) {
    cli::cli_abort("No tiles to request (AOI too small or outside the 3 deg grid).")
  }
  
  # Helper: tile code from lower-left corner (e.g., N48E006, S03W072)
  tile_code_from_poly <- function(poly) {
    bb <- sf::st_bbox(poly)
    lat_ll <- as.integer(bb["ymin"])
    lon_ll <- as.integer(bb["xmin"])
    lat_tag <- if (lat_ll < 0) paste0("S", sprintf("%02d", abs(lat_ll))) else paste0("N", sprintf("%02d", lat_ll))
    lon_tag <- if (lon_ll < 0) paste0("W", sprintf("%03d", abs(lon_ll))) else paste0("E", sprintf("%03d", lon_ll))
    paste0(lat_tag, lon_tag)
  }
  
  # Be robust to geometry column name ("x" vs "geometry")
  geom_col <- attr(hits, "sf_column")
  codes <- unique(vapply(hits[[geom_col]], tile_code_from_poly, character(1)))
  
  # S3 URLs
  s3_base <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com"
  urls <- sprintf(
    "%s/%s/%d/%s/ESA_WorldCover_10m_%d_%s_%s_%s.tif",
    s3_base, version, year, tolower(layer),
    year, version, codes, tools::toTitleCase(layer) # e.g., "Map"
  )
  dests <- file.path(out_dir, basename(urls))
  
  # Optional HEAD check (skip obviously missing URLs)
  good <- rep(TRUE, length(urls))
  if (isTRUE(check_remote)) {
    if (isTRUE(progress)) cli::cli_alert_info("Checking remote availability with HEAD requests ...")
    for (i in seq_along(urls)) {
      u <- urls[i]
      h <- curl::new_handle(nobody = TRUE, timeout = 10)
      res <- tryCatch(curl::curl_fetch_memory(u, handle = h), error = function(e) NULL)
      good[i] <- !is.null(res) && identical(res$status_code, 200L)
    }
  }
  
  # Download loop with optional progress
  n <- length(codes)
  status <- character(n)
  size_bytes <- rep(NA_real_, n)
  
  it <- if (isTRUE(progress)) {
    cli::cli_progress_along(
      seq_len(n),
      total        = n,
      clear        = FALSE,
      format       = "{cli::pb_spin} Processing {cli::pb_current}/{cli::pb_total}",
      format_done  = "{cli::col_green(cli::symbol$tick)} Done: {cli::pb_total} items."
    )
  } else {
    seq_len(n)
  }
  
  for (i in it) {
    u <- urls[i]
    d <- dests[i]
    ok_remote <- good[i]
    
    # Skip if remote is flagged missing
    if (!ok_remote) {
      status[i] <- "missing_remote"
      next
    }
    
    # Respect overwrite flag
    if (file.exists(d) && !isTRUE(overwrite)) {
      status[i] <- "exists"
      size_bytes[i] <- if (file.exists(d)) file.info(d)$size else NA_real_
      next
    }
    
    # Ensure parent dir
    parent <- dirname(d)
    if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    
    # Attempt download
    dl_ok <- tryCatch(
      {
        curl::curl_download(u, d, mode = "wb", quiet = TRUE)
        TRUE
      },
      error = function(e) {
        cli::cli_alert_warning("Download failed: {basename(d)} - {conditionMessage(e)}")
        FALSE
      }
    )
    
    status[i]    <- if (isTRUE(dl_ok) && file.exists(d)) "downloaded" else "error"
    size_bytes[i] <- if (file.exists(d)) file.info(d)$size else NA_real_
  }
  
  # Assemble results (tibble for nice printing; keep geometry list-column)
  results <- tibble::tibble(
    tile_id   = codes,
    file      = dests,
    geometry  = hits[[geom_col]],
    url       = urls,
    to_get    = good,
    status    = status,
    size_bytes = size_bytes
  )
  
  # Summary (shown even if progress = FALSE; comment out if you want silence)
  cli::cli_alert_info(sprintf(
    "WorldCover tiles: %d total (%d downloaded, %d existed, %d missing, %d errors)",
    nrow(results),
    sum(results$status == "downloaded"),
    sum(results$status == "exists"),
    sum(results$status == "missing_remote"),
    sum(results$status == "error")
  ))
  
  # Preserve legacy attributes for downstream compatibility
  attr(df, "worldcover_tiles") <- results
  attr(df, "worldcover_aoi")   <- aoi_ll
  
  # Also store in unified layer container
  df <- .set_traktr_layer(
    df, "worldcover",
    index  = results,
    params = list(
      year = year, version = version, layer = layer,
      out_dir = out_dir, check_remote = check_remote, overwrite = overwrite
    )
  )
  
  if (isTRUE(progress)) cli::cli_progress_done()
  invisible(df)
}