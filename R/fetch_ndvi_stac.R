#' Fetch historic NDVI composites from a STAC catalogue (Landsat / MODIS)
#'
#' Tiles an AOI into chips, searches a STAC API (by default Microsoft
#' Planetary Computer) for Landsat Collection 2 Level-2 or MODIS MOD13Q1
#' scenes in a date range, masks clouds, computes NDVI and writes a per-pixel
#' median composite per chip as GeoTIFF. The output layout matches
#' \code{\link{fetch_ndvi}}, so \code{\link{attach_ndvi}} works unchanged.
#'
#' Use this for years before Sentinel-2 (CDSE coverage starts 2015-2017):
#' Landsat 4-9 reaches back to 1982 at 30 m, MODIS to 2000 at 250 m.
#'
#' @param df A data object carrying the AOI. If \code{aoi} is not supplied,
#'   the function reads \code{attr(df, "traktr_aoi")}.
#' @param out_dir Directory where output GeoTIFF files are written.
#' @param from,to Date range (inclusive) of scenes to composite. Use the same
#'   season every year (e.g. June-August) for comparable values.
#' @param collection STAC collection: \code{"landsat-c2-l2"} (Landsat 4-9,
#'   30 m) or \code{"modis-13Q1-061"} (MODIS Terra 16-day NDVI, 250 m).
#' @param resolution_m Output pixel size in meters (approximate; chips are
#'   written in EPSG:4326 on a global grid). Defaults to the native
#'   resolution of \code{collection}.
#' @param max_cloud Maximum scene cloud cover in percent (Landsat only).
#' @param max_items Maximum number of scenes per chip (least cloudy first).
#' @param platforms Optional character vector to restrict Landsat platforms,
#'   e.g. \code{c("landsat-5", "landsat-8")}. Landsat 7 scenes after May 2003
#'   have scan-line gaps; these are masked as missing and usually filled by
#'   the composite, but can be excluded here.
#' @param tile_max_px Maximum tile width/height in pixels.
#' @param stac_url STAC API endpoint.
#' @param overwrite Logical. If \code{FALSE} (default), existing files are
#'   skipped.
#' @param progress Logical. If \code{TRUE} (default), show a progress bar.
#' @param aoi Optional \code{sf} polygon overriding \code{attr(df, "traktr_aoi")}.
#'
#' @details
#' - Landsat: surface reflectance is scaled (\code{DN * 2.75e-5 - 0.2}) before
#'   computing \code{(nir08 - red) / (nir08 + red)}; pixels flagged in
#'   \code{qa_pixel} as fill, dilated cloud, cirrus, cloud, cloud shadow or
#'   snow are removed.
#' - MODIS: the delivered NDVI (\code{* 1e-4}) is used; only pixels with
#'   reliability 0 (good) or 1 (marginal) are kept.
#' - Each chip is a 2-band GeoTIFF: band 1 median NDVI, band 2 number of
#'   valid observations. \code{attach_ndvi()} uses band 1.
#' - Planetary Computer asset URLs are signed per chip (signatures expire after
#'   about an hour), no account is needed.
#'
#' @section Attributes:
#' - \code{attr(df, "ndvi_aoi")}: \code{sf} tile index with \code{tile_id},
#'   \code{file}, \code{status}, \code{n_items}, (\code{msg}).
#' - \code{attr(df, "traktr_layers")$ndvi}: same index plus fetch parameters.
#'
#' @return The input \code{df} (invisibly) with the attributes above.
#'
#' @seealso \code{\link{fetch_ndvi}}, \code{\link{attach_ndvi}}
#'
#' @examples
#' \dontrun{
#' # Summer NDVI 2005 around childhood addresses
#' df <- mutate_aoi_buffer(addresses, dist = 500)
#' df <- fetch_ndvi_stac(df, out_dir = "data/ndvi_2005",
#'                       from = "2005-06-01", to = "2005-08-31")
#' df <- attach_ndvi(df)
#'
#' # One year per row: run per year and bind
#' res <- lapply(split(addresses, addresses$year), function(d) {
#'   y <- d$year[1]
#'   d |>
#'     mutate_aoi_buffer(dist = 500) |>
#'     fetch_ndvi_stac(out_dir = file.path("data/ndvi", y),
#'                     from = paste0(y, "-06-01"), to = paste0(y, "-08-31")) |>
#'     attach_ndvi()
#' })
#' res <- do.call(rbind, res)
#' }
#'
#' @export
fetch_ndvi_stac <- function(df,
                            out_dir,
                            from, to,
                            collection   = c("landsat-c2-l2", "modis-13Q1-061"),
                            resolution_m = NULL,
                            max_cloud    = 60,
                            max_items    = 30,
                            platforms    = NULL,
                            tile_max_px  = 1000,
                            stac_url     = "https://planetarycomputer.microsoft.com/api/stac/v1",
                            overwrite    = FALSE,
                            progress     = TRUE,
                            aoi          = NULL) {
  rlang::check_installed(c("rstac", "terra", "sf", "cli", "tibble"))

  collection <- match.arg(collection)
  spec <- .stac_ndvi_spec(collection)
  if (is.null(resolution_m)) resolution_m <- spec$res_m
  from <- as.Date(from); to <- as.Date(to)

  aoi <- .resolve_aoi(df, aoi)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  tiles <- .make_request_tiles(aoi, res_m = resolution_m, max_px = tile_max_px, overlap_px = 0)
  n_tiles <- nrow(tiles$chips_wgs84)
  if (!n_tiles) cli::cli_abort("No tiles to request (AOI too small or filtered).")

  # Faster remote COG access
  old_env <- Sys.getenv(c("GDAL_DISABLE_READDIR_ON_OPEN", "CPL_VSIL_CURL_ALLOWED_EXTENSIONS"),
                        unset = NA)
  Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
             CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.TIF")
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(old_env[nm]))
    }
  }, add = TRUE)

  skipped <- 0L
  pb <- if (isTRUE(progress)) {
    cli::cli_progress_bar(
      total  = n_tiles,
      format = "Compositing {collection} {cli::pb_bar} {cli::pb_percent} (skipped {skipped})",
      clear  = FALSE
    )
  } else NULL

  res_list <- vector("list", n_tiles)
  for (i in seq_len(n_tiles)) {
    chip <- tiles$chips_wgs84[i, ]
    bbox <- sf::st_bbox(chip)
    out_file <- file.path(out_dir, sprintf(
      "ndvi_%s_%s_%s_%s_%s_%s_%s.tif",
      collection,
      .sanitize_num(bbox$xmin), .sanitize_num(bbox$ymin),
      .sanitize_num(bbox$xmax), .sanitize_num(bbox$ymax),
      as.character(from), as.character(to)
    ))

    if (!overwrite && file.exists(out_file)) {
      skipped <- skipped + 1L
      res_list[[i]] <- tibble::tibble(tile_id = chip$tile_id, file = out_file,
                                      status = "skipped", n_items = NA_integer_,
                                      msg = NA_character_)
    } else {
      res_list[[i]] <- tryCatch({
        n <- .stac_ndvi_chip(chip, from, to, spec, resolution_m, max_cloud,
                             max_items, platforms, stac_url, out_file)
        tibble::tibble(tile_id = chip$tile_id,
                       file    = if (n > 0) out_file else NA_character_,
                       status  = if (n > 0) "ok" else "no_data",
                       n_items = n, msg = NA_character_)
      }, error = function(e) {
        tibble::tibble(tile_id = chip$tile_id, file = NA_character_, status = "error",
                       n_items = NA_integer_, msg = conditionMessage(e))
      })
    }
    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }
  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  results <- do.call(rbind, res_list)
  n_bad <- sum(results$status %in% c("error", "no_data"))
  if (n_bad) {
    cli::cli_warn("{n_bad} of {n_tiles} chip{?s} without NDVI (see {.code attr(df, \"ndvi_aoi\")}).")
  }

  left <- tiles$chips_wgs84[, intersect(names(tiles$chips_wgs84), c("tile_id", "geometry")), drop = FALSE]
  tile_index <- sf::st_as_sf(merge(left, results, by = "tile_id", all.x = TRUE))

  params <- list(
    source = "stac", collection = collection, stac_url = stac_url,
    from = from, to = to, out_dir = out_dir,
    resolution_m = resolution_m, max_cloud = max_cloud, max_items = max_items,
    platforms = platforms, tile_max_px = tile_max_px, overwrite = overwrite
  )
  attr(df, "ndvi_aoi") <- tile_index
  attr(df, "ndvi_fetch_params") <- params
  df <- .set_traktr_layer(df, "ndvi", index = tile_index, params = params)

  invisible(df)
}


#' Collection-specific settings for STAC NDVI (internal)
#' @keywords internal
#' @noRd
.stac_ndvi_spec <- function(collection) {
  switch(collection,
    "landsat-c2-l2" = list(
      collection = collection, res_m = 30, cloud_filter = TRUE,
      assets = c(red = "red", nir = "nir08", qa = "qa_pixel"),
      ndvi = function(b) {
        red <- b$red * 2.75e-5 - 0.2
        nir <- b$nir * 2.75e-5 - 0.2
        # QA_PIXEL bits 0-5: fill, dilated cloud, cirrus, cloud, shadow, snow
        clear <- terra::app(b$qa, function(q) bitwAnd(as.integer(q), 63L) == 0L)
        v <- (nir - red) / (nir + red)
        v <- terra::mask(v, clear, maskvalues = c(0, NA))
        terra::clamp(v, -1, 1, values = FALSE)
      }
    ),
    "modis-13Q1-061" = list(
      collection = collection, res_m = 250, cloud_filter = FALSE,
      assets = c(ndvi = "250m_16_days_NDVI", qa = "250m_16_days_pixel_reliability"),
      ndvi = function(b) {
        v <- b$ndvi * 1e-4
        ok <- b$qa == 0 | b$qa == 1
        terra::mask(v, ok, maskvalues = c(0, NA))
      }
    )
  )
}


#' Composite NDVI for one chip and write it to disk (internal)
#'
#' @return Number of scenes that contributed (0 = nothing written).
#' @keywords internal
#' @noRd
.stac_ndvi_chip <- function(chip, from, to, spec, resolution_m, max_cloud,
                            max_items, platforms, stac_url, out_file) {
  bb <- as.numeric(sf::st_bbox(chip))

  items <- rstac::stac(stac_url) |>
    rstac::stac_search(
      collections = spec$collection,
      bbox        = bb,
      datetime    = paste0(format(from), "T00:00:00Z/", format(to), "T23:59:59Z"),
      limit       = 1000
    ) |>
    rstac::post_request()

  feats <- items$features
  props <- function(f, k) { v <- f$properties[[k]]; if (is.null(v)) NA else v }
  if (spec$cloud_filter && length(feats)) {
    cc <- vapply(feats, props, numeric(1), k = "eo:cloud_cover")
    keep <- is.na(cc) | cc <= max_cloud
    feats <- feats[keep]; cc <- cc[keep]
    feats <- feats[order(cc)]
  }
  if (!is.null(platforms) && length(feats)) {
    pf <- vapply(feats, props, character(1), k = "platform")
    feats <- feats[pf %in% platforms]
  }
  if (!length(feats)) return(0L)
  feats <- utils::head(feats, max_items)
  items$features <- feats
  items <- rstac::items_sign(items, sign_fn = rstac::sign_planetary_computer())

  template <- .stac_template(chip, resolution_m)

  layers <- list()
  for (f in items$features) {
    v <- tryCatch(.stac_item_ndvi(f, spec, template), error = function(e) NULL)
    if (!is.null(v) && !all(is.na(terra::minmax(v)))) layers[[length(layers) + 1L]] <- v
  }
  if (!length(layers)) return(0L)

  s <- terra::rast(layers)
  med <- terra::app(s, stats::median, na.rm = TRUE)
  cnt <- terra::app(!is.na(s), sum)
  out <- c(med, cnt)
  names(out) <- c("ndvi", "n_obs")
  terra::writeRaster(out, out_file, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE"))
  length(layers)
}


#' Output grid for a chip: EPSG:4326, globally aligned (internal)
#'
#' All chips share origin and resolution, so they can be combined into one
#' VRT by \code{attach_ndvi()}.
#'
#' @keywords internal
#' @noRd
.stac_template <- function(chip, resolution_m) {
  res <- resolution_m / 111320
  bb  <- sf::st_bbox(chip)
  e <- terra::ext(
    floor(bb$xmin / res) * res, ceiling(bb$xmax / res) * res,
    floor(bb$ymin / res) * res, ceiling(bb$ymax / res) * res
  )
  terra::rast(e, resolution = res, crs = "EPSG:4326")
}


#' NDVI of one STAC item resampled to the chip template (internal)
#' @keywords internal
#' @noRd
.stac_item_ndvi <- function(feature, spec, template) {
  first <- NULL
  bands <- lapply(spec$assets, function(key) {
    href <- feature$assets[[key]]$href
    if (is.null(href)) cli::cli_abort("Asset {key} missing.")
    r <- terra::rast(paste0("/vsicurl/", href))
    # use raw DNs; scaling is defined in `spec` (MODIS GDAL metadata uses
    # the inverse HDF convention, scale = 10000)
    terra::scoff(r) <- cbind(1, 0)
    if (is.null(first)) {
      # crop window in the item's native CRS (+ one native pixel margin)
      win <- terra::ext(terra::project(terra::as.polygons(terra::ext(template), crs = "EPSG:4326"),
                                       terra::crs(r)))
      first <<- terra::extend(win, terra::res(r)[1])
    }
    terra::crop(r, first, snap = "out")
  })
  v <- spec$ndvi(bands)
  terra::project(v, template, method = "bilinear")
}
