#' Fetch ESA WorldCover tiles for a given AOI
#'
#' Downloads 10-meter land-cover classification tiles from the
#' [ESA WorldCover](https://esa-worldcover.org) AWS Open Data repository
#' for a given AOI. The AOI is split into 3 deg x 3 deg grid tiles, and matching
#' WorldCover TIFFs are downloaded to a local directory.
#'
#' @param df An `sf` polygon or multipolygon representing the AOI. Must have a valid CRS.
#' @param out_dir Output directory where downloaded `.tif` files will be saved.
#' @param wc_year WorldCover year, either `2020` or `2021` (default: `2021`).
#' @param wc_version Version tag; `"v100"` for 2020 or `"v200"` for 2021 (default: `"v200"`).
#' @param layer Character; layer type (e.g., `"map"` for land-cover classification).
#' @param overwrite Logical; if `FALSE` (default), skip existing files.
#' @param check_exists_remote Logical; if `TRUE`, perform HTTP HEAD check on all tile URLs before download.
#'
#' @return The input AOI with two attributes attached:
#' \describe{
#'   \item{`worldcover_tiles`}{A tibble with one row per tile containing `tile_id`, `file`, `url`, `geometry`, `status`, and `size_bytes`.}
#'   \item{`worldcover_aoi`}{The reprojected AOI in EPSG:4326 used for tile intersection.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' poly <- st_polygon(list(rbind(c(8,47), c(8.5,47), c(8.5,47.5), c(8,47.5), c(8,47))))
#' aoi <- st_as_sf(st_sfc(poly, crs = 4326))
#' result <- fetch_worldcover(aoi, out_dir = "worldcover_data")
#' attr(result, "worldcover_tiles")
#' }
#'
#' @seealso \href{https://esa-worldcover.org}{ESA WorldCover}, [curl::curl_download()]
#' @export
fetch_landuse <- function(
    df,
    out_dir,
    wc_year = 2021,                 # 2020 (v100) or 2021 (v200)
    wc_version = "v200",            # "v100" for 2020, "v200" for 2021
    layer = "map",                  # "map" is the 10m land-cover classification
    overwrite = FALSE,
    check_exists_remote = FALSE     # set TRUE to HEAD-check URL before downloading
){
  cli::cli_progress_step("Fetching worldcover data from AWS", spinner = F, msg_done = "Fetching worldcover completed.")
  
  aoi <- attr(df, "traktr_aoi")
  
  stopifnot(inherits(aoi, "sf"))
  if (is.na(sf::st_crs(aoi))) cli::cli_abort("`aoi` has no CRS; assign one first.")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Always work in EPSG:4326 for the 3 deg. tiling
  aoi_ll <- sf::st_make_valid(sf::st_transform(aoi, 4326))
  # If your AOI could cross the dateline, consider:
  # aoi_ll <- st_wrap_dateline(aoi_ll, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  
  # Build a 3 deg. grid aligned to multiples of 3 deg. that covers the AOI bbox
  bb  <- sf::st_bbox(aoi_ll)
  llx <- as.numeric(floor(bb["xmin"]/3) * 3)
  lly <- as.numeric(floor(bb["ymin"]/3) * 3)
  urx <- as.numeric(ceiling(bb["xmax"]/3) * 3)
  ury <- as.numeric(ceiling(bb["ymax"]/3) * 3)
  
  cover_sfc <- sf::st_as_sfc(sf::st_bbox(c(xmin = llx, ymin = lly, xmax = urx, ymax = ury), crs = 4326))
  grid3 <- sf::st_make_grid(cover_sfc, cellsize = c(3,3), what = "polygons", square = TRUE) |>
    sf::st_as_sf()
  
    # Keep grid cells that intersect your AOI (union to speed it up)
  hits <- grid3[ sf::st_intersects(grid3, sf::st_union(aoi_ll), sparse = FALSE)[,1], , drop = FALSE]
  
  # Helper: tile code from lower-left corner (e.g., N48E006, S03W072)
  tile_code_from_poly <- function(poly){
    bb <- st_bbox(poly)
    lat_ll <- as.integer(bb["ymin"])
    lon_ll <- as.integer(bb["xmin"])
    lat_tag <- if (lat_ll < 0) {
      paste0("S", sprintf("%02d", abs(lat_ll)))
    } else {
      paste0("N", sprintf("%02d", lat_ll))
    }
    
    lon_tag <- if (lon_ll < 0) {
      paste0("W", sprintf("%03d", abs(lon_ll)))
    } else {
      paste0("E", sprintf("%03d", lon_ll))
    }
    #lat_tag <- if (lat_ll < 0) paste0("S", stringr::str_pad(abs(lat_ll), 2, pad = "0"))
    #else            paste0("N", stringr::str_pad(lat_ll,      2, pad = "0"))
    #lon_tag <- if (lon_ll < 0) paste0("W", stringr::str_pad(abs(lon_ll), 3, pad = "0"))
    #else            paste0("E", stringr::str_pad(lon_ll,      3, pad = "0"))
    paste0(lat_tag, lon_tag)
  }
  

  #codes <- unique(purrr::map_chr(hits$x, tile_code_from_poly))
  codes <- unique(vapply(hits$x, tile_code_from_poly, character(1)))
  
  
  s3_base <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com"
  urls <- sprintf(
    "%s/%s/%d/%s/ESA_WorldCover_10m_%d_%s_%s_%s.tif",
    s3_base, wc_version, wc_year, tolower(layer),
    wc_year, wc_version, codes, tools::toTitleCase(layer)  # Map -> 'Map'
  )
  
  dests <- file.path(out_dir, basename(urls))
  
  # Optional HEAD check (skip obviously missing URLs)
  #good <- rep(TRUE, length(urls))
  #if (check_exists_remote) {
  #  cli::cli_alert_info("Checking remote availability with HEAD requests ...")
  #  good <- map_lgl(urls, function(u){
  #    h <- new_handle(nobody = TRUE)
  #    res <- tryCatch(curl_fetch_memory(u, handle = h), error = function(e) NULL)
  #    !is.null(res) && identical(res$status_code, 200L)
  #  })
  #  }
  
  good <- rep(TRUE, length(urls))
  if (check_exists_remote) {
    cli::cli_alert_info("Checking remote availability with HEAD requests ...")
    
    for (i in seq_along(urls)) {
      u <- urls[i]
      h <- curl::new_handle(nobody = TRUE, timeout = 10)
      # optional: curl::handle_setheaders(h, "User-Agent" = "yourpkg/1.0")
      res <- tryCatch(
        curl::curl_fetch_memory(u, handle = h),
        error = function(e) NULL
      )
      good[i] <- !is.null(res) && identical(res$status_code, 200L)
      
      # optional: update progress
      # cli::cli_progress_update(status = paste0("Checking ", basename(u)))
    }
  }
  
  # Download
# 
#   results <- tibble::tibble(
#     tile_id = codes,
#     file      = dests,
#     geometry  = hits,
#     url       = urls,
#     to_get    = good
#   ) |>
#     mutate(
#       status = purrr::pmap_chr(list(url, file, to_get), function(u, d, ok){
#         if (!ok) return("missing_remote")
#         if (file.exists(d) && !overwrite) return("exists")
#         # Try download
#         ok_dl <- tryCatch({ curl::curl_download(u, d, mode = "wb", quiet = F); TRUE },
#                           error = function(e) { 
#                             cli::cli_alert_warning("Download failed: {basename(d)} - {conditionMessage(e)}"); FALSE }
#                          )
#         if (!ok_dl) return("error")
#         "downloaded"
#       }),
#       size_bytes = ifelse(file.exists(file), file.info(file)$size, NA_real_)
#     )
  
  rlang::check_installed(c("curl", "cli", "tibble"))
  
  n <- length(codes)
  status <- character(n)
  size_bytes <- rep(NA_real_, n)
  
  # progress iterator (optional)
  it <- cli::cli_progress_along(
    seq_len(n),
    total = n,
    clear = FALSE,
    format = "{cli::pb_spin} Processing {cli::pb_current}/{cli::pb_total}",
    format_done = "{cli::col_green(cli::symbol$tick)} Done: {cli::pb_total} items."
  )
  
  for (i in it) {
    u <- urls[i]
    d <- dests[i]
    ok_remote <- good[i]
    
    # skip if remote is flagged missing
    if (!ok_remote) {
      status[i] <- "missing_remote"
      next
    }
    
    # if file exists and we should not overwrite
    if (file.exists(d) && !isTRUE(overwrite)) {
      status[i] <- "exists"
      # record size if present
      size_bytes[i] <- if (file.exists(d)) file.info(d)$size else NA_real_
      next
    }
    
    # ensure parent dir exists
    parent <- dirname(d)
    if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    
    # attempt download
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
    
    status[i] <- if (isTRUE(dl_ok) && file.exists(d)) "downloaded" else "error"
    size_bytes[i] <- if (file.exists(d)) file.info(d)$size else NA_real_
  }
  
  # assemble result (as tibble if tibble is available)
  results <- data.frame(
    tile_id = codes,
    file = dests,
    geometry = I(hits),   # I() to keep list-columns if geometry/hits are lists
    url = urls,
    to_get = good,
    status = status,
    size_bytes = size_bytes,
    stringsAsFactors = FALSE
  )
  
  # optional: return as tibble for nicer printing
  results <- tibble::as_tibble(results)
  
  
  cli::cli_alert_info(sprintf(
    "WorldCover tiles: %d total (%d downloaded, %d existed, %d missing, %d errors)",
    nrow(results),
    sum(results$status == "downloaded"),
    sum(results$status == "exists"),
    sum(results$status == "missing_remote"),
    sum(results$status == "error")
  ) )
  
  attr(df, "worldcover_tiles") <- results
  attr(df, "worldcover_aoi") <- aoi_ll
  
  
  cli::cli_progress_done()
  
  df
}
