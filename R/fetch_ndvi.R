#' Fetch NDVI chips from CDSE over an AOI
#'
#' Tiles an area of interest (AOI) into WGS84 chips at a target resolution,
#' then downloads NDVI GeoTIFFs from the Copernicus Data Space Ecosystem (CDSE)
#' using a provided evalscript. The function is rate-limited, honors local
#' caching, and records output metadata on the input object.
#'
#' @param df A data object carrying the AOI. If \code{aoi} is not supplied,
#'   the function reads \code{attr(df, "traktr_aoi")}. The AOI must be an
#'   \code{sf} polygon/multipolygon with a valid CRS.
#' @param out_dir Directory where output GeoTIFF files are written. Will be
#'   created if it does not exist.
#' @param from,to Date range (inclusive) used to query NDVI imagery. Coerced
#'   with \code{as.Date()}.
#' @param resolution_m Target pixel size in meters used when requesting chips.
#' @param tile_max_px Maximum tile width/height in pixels used by the tiler.
#'   Controls chip size and request count.
#' @param tile_overlap_px Optional overlap between adjacent tiles in pixels.
#' @param client OAuth client for CDSE. Typically an \code{httr2} OAuth client
#'   or equivalent object.
#' @param evalscript Path or character string containing the Sentinel Hub
#'   evalscript to compute NDVI.
#' @param rate A \code{ratelimitr::rate()} object controlling request rate,
#'   e.g. \code{ratelimitr::rate(n = 2, period = 1)}.
#' @param overwrite Logical. If \code{FALSE} (default), existing files are
#'   skipped. If \code{TRUE}, files are re-downloaded.
#' @param progress Logical. If \code{TRUE} (default), show progress messages
#'   and a progress bar.
#' @param aoi Optional \code{sf} polygon/multipolygon. If supplied, overrides
#'   \code{attr(df, "traktr_aoi")}.
#'
#' @param res_m Deprecated. Use \code{resolution_m}.
#' @param max_px Deprecated. Use \code{tile_max_px}.
#' @param overlap_px Deprecated. Use \code{tile_overlap_px}.
#' @param rate_n,rate_period Deprecated. Use \code{rate = ratelimitr::rate(n, period)}.
#'
#' @details
#' - AOI handling: The AOI is validated, transformed to EPSG:4326, and used to
#'   build request chips. If the AOI has no CRS, the call aborts.
#' - Caching: Filenames encode the chip bounding box (xmin, ymin, xmax, ymax)
#'   and the date range, which provides a stable cache key across runs.
#' - Rate limiting: All remote requests are wrapped with \code{ratelimitr}.
#' - Attributes:
#'   * \code{attr(df, "ndvi_aoi")} stores an \code{sf} index of chips with
#'     status and file paths.
#'   * \code{attr(df, "ndvi_fetch_params")} stores call parameters.
#'   * \code{attr(df, "traktr_layers")$ndvi} stores the same information in a
#'     harmonized container for downstream discovery.
#'
#' @return The input \code{df}, returned invisibly, augmented with attributes
#'   describing outputs and call parameters (see Details).
#'
#' @seealso
#' \code{\link{fetch_landuse}}, \code{ratelimitr::rate}
#'
#' @examples
#' \dontrun{
#' df <- some_traktr_object_with_aoi()
#' df <- fetch_ndvi(
#'   df,
#'   out_dir       = "data/ndvi",
#'   from          = as.Date("2024-06-01"),
#'   to            = as.Date("2024-06-30"),
#'   client        = my_cdse_client,
#'   resolution_m  = 10,
#'   overwrite     = FALSE,
#'   progress      = TRUE
#' )
#' # Inspect chip index
#' attr(df, "traktr_layers")$ndvi$index
#' }
#'
#' @export
fetch_ndvi <- function(
    df,
    out_dir,
    from, to,                         # temporal selector
    resolution_m     = 10,
    tile_max_px      = 2500,
    tile_overlap_px  = 0,
    client,                           # required: OAuth client for CDSE
    evalscript       = .ndvi_evalscript_path(),
    rate             = ratelimitr::rate(n = 2, period = 1),
    overwrite        = FALSE,
    progress         = TRUE,
    aoi              = NULL,
    # ---- deprecated compat args (kept for one release cycle) ----
    res_m = NULL, max_px = NULL, overlap_px = NULL, rate_n = NULL, rate_period = NULL
) {
  
  # -- package prerequisites --------------------------------------------
  rlang::check_installed(c("cli", "ratelimitr", "tibble", "sf"))
  
  # -- deprecated args shims --------------------------------------------
  if (!is.null(res_m)) {
    warning("`res_m` is deprecated; use `resolution_m`.", call. = FALSE)
    resolution_m <- res_m
  }
  if (!is.null(max_px)) {
    warning("`max_px` is deprecated; use `tile_max_px`.", call. = FALSE)
    tile_max_px <- max_px
  }
  if (!is.null(overlap_px)) {
    warning("`overlap_px` is deprecated; use `tile_overlap_px`.", call. = FALSE)
    tile_overlap_px <- overlap_px
  }
  if (!is.null(rate_n) || !is.null(rate_period)) {
    warning("`rate_n`/`rate_period` are deprecated; use `rate = ratelimitr::rate(n, period)`.", call. = FALSE)
    n <- if (!is.null(rate_n)) rate_n else 2
    p <- if (!is.null(rate_period)) rate_period else 1
    rate <- ratelimitr::rate(n = n, period = p)
  }
  
  # -- AOI, dirs, and progress ------------------------------------------
  aoi <- .resolve_aoi(df, aoi)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(progress)) {
    cli::cli_progress_step(
      "Fetching NDVI data from CDSE",
      spinner  = TRUE,
      msg_done = "Fetching NDVI data completed."
    )
  }
  
  # 1) Build request tiles (keep original API of helper)
  tiles <- .make_request_tiles(
    aoi,
    res_m      = resolution_m,      # pass with legacy arg names expected by helper
    max_px     = tile_max_px,
    overlap_px = tile_overlap_px
  )
  
  n_tiles <- nrow(tiles$chips_wgs84)
  if (!n_tiles) cli::cli_abort("No tiles to request (AOI too small or filtered).")
  idx <- seq_len(n_tiles)
  
  # 2) Safe + rate-limited chip fetch
  safe_one <- function(i) {
    tryCatch({
      bbox <- sf::st_bbox(tiles$chips_wgs84[i, ])
      
      filename <- sprintf(
        "ndvi_%s_%s_%s_%s_%s_%s.tif",
        .sanitize_num(bbox$xmin),
        .sanitize_num(bbox$ymin),
        .sanitize_num(bbox$xmax),
        .sanitize_num(bbox$ymax),
        as.character(as.Date(from)),
        as.character(as.Date(to))
      )
      
      out_file <- file.path(out_dir, filename)
      
      result <- .cdse_fetch_chip(
        tiles$chips_wgs84[i, ],
        from, to,
        evalscript,
        out_file,
        resolution = resolution_m,
        client     = client
      )
      
      list(result = result, error = NULL)
    }, error = function(e) {
      list(result = NULL, error = e)
    })
  }
  
  limited <- ratelimitr::limit_rate(safe_one, rate)
  
  # 3) Progress + loop
  skipped <- 0L
  pb <- if (isTRUE(progress)) {
    cli::cli_progress_bar(
      total  = n_tiles,
      format = "Downloading {cli::pb_bar} {cli::pb_percent} (skipped {skipped})",
      clear  = FALSE
    )
  } else NULL
  
  res_list <- vector("list", n_tiles)
  
  for (j in seq_along(idx)) {
    i <- idx[j]
    
    bbox <- sf::st_bbox(tiles$chips_wgs84[i, ])
    filename <- sprintf(
      "ndvi_%s_%s_%s_%s_%s_%s.tif",
      .sanitize_num(bbox$xmin),
      .sanitize_num(bbox$ymin),
      .sanitize_num(bbox$xmax),
      .sanitize_num(bbox$ymax),
      as.character(as.Date(from)),
      as.character(as.Date(to))
    )
    out_file <- file.path(out_dir, filename)
    
    if (!overwrite && file.exists(out_file)) {
      skipped <- skipped + 1L
      res_list[[j]] <- tibble::tibble(
        tile_id = tiles$chips_wgs84$tile_id[i],
        file    = out_file,
        status  = "skipped"
      )
    } else {
      if (is.null(client)) {
        cli::cli_abort("No OAuth client set. Pass a valid `client`.")
      }
      ans <- limited(i)
      
      if (!is.null(ans$error)) {
        res_list[[j]] <- tibble::tibble(
          tile_id = tiles$chips_wgs84$tile_id[i],
          file    = NA_character_,
          status  = "error",
          msg     = paste(ans$error)
        )
      } else {
        res_list[[j]] <- ans$result
      }
    }
    
    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }
  
  if (!is.null(pb)) cli::cli_progress_done(id = pb)
  
  # 4) Bind results and normalize columns
  results <- do.call(rbind, res_list)
  
  if ("id" %in% names(results) && !"tile_id" %in% names(results)) {
    results$tile_id <- results$id
    results <- results[, !(names(results) %in% "id")]
  }
  
  # 5) Left-join chip geometry with results
  left <- tiles$chips_wgs84[, intersect(names(tiles$chips_wgs84), c("tile_id", "geometry")), drop = FALSE]
  right <- results[, intersect(names(results), c("tile_id", "file", "status", "msg")), drop = FALSE]
  
  tile_index <- merge(left, right, by = "tile_id", all.x = TRUE)
  tile_index <- sf::st_as_sf(tile_index)
  
  # 6) Legacy attrs + unified layer container
  attr(df, "ndvi_aoi") <- tile_index
  attr(df, "ndvi_fetch_params") <- list(
    from = as.Date(from), to = as.Date(to),
    out_dir = out_dir,
    resolution_m = resolution_m,
    tile_max_px = tile_max_px,
    tile_overlap_px = tile_overlap_px,
    evalscript = evalscript,
    rate = rate,
    overwrite = overwrite
  )
  
  df <- .set_traktr_layer(
    df, "ndvi",
    index  = tile_index,
    params = list(
      from = as.Date(from), to = as.Date(to),
      out_dir = out_dir,
      resolution_m = resolution_m,
      tile_max_px = tile_max_px,
      tile_overlap_px = tile_overlap_px,
      evalscript = evalscript,
      rate = rate,
      overwrite = overwrite
    )
  )
  
  invisible(df)
}




#' Split an AOI into request tiles under a pixel-size limit
#'
#' Splits an \pkg{sf}/\pkg{sfc} AOI into square tiles in a metric CRS so that the
#' projected width/height in pixels at a given ground resolution does not exceed
#' a maximum (e.g., API limit). Optionally adds overlap (to reduce seam artifacts)
#' and drops tiny slivers. Returns both metric (e.g., EPSG:3857) and WGS84 versions.
#'
#' @param aoi An \code{sf} or \code{sfc} polygon/multipolygon representing the AOI.
#'   Can be in any CRS; it will be reprojected to \code{crs_metric}.
#' @param res_m Numeric. Ground resolution in meters per pixel to plan for
#'   (e.g., \code{10} or \code{20}).
#' @param max_px Integer. Maximum allowed pixel width/height per request
#'   (default \code{2500}).
#' @param overlap_px Integer. Overlap in pixels to add around each tile
#'   (default \code{0}).
#' @param min_area_m2 Numeric. Minimum tile area (m^2) to keep after clipping
#'   (default \code{0}).
#' @param crs_metric Integer or CRS. Metric CRS used to build tiles
#'   (default \code{3857}, Web Mercator).
#'
#' @return A list with two \code{sf} objects:
#' \describe{
#'   \item{\code{chips_metric}}{Tiles in \code{crs_metric} with columns
#'     \code{tile_id}, \code{area_m2}, \code{width_m}, \code{height_m},
#'     \code{width_px}, \code{height_px}, and geometry.}
#'   \item{\code{chips_wgs84}}{Same tiles transformed to EPSG:4326.}
#' }
#'
#' @details
#' The grid cell size is computed as \code{(max_px - 2 * overlap_px) * res_m} meters.
#' Each grid cell is optionally buffered by \code{overlap_px * res_m}, intersected
#' with the AOI, tiny slivers are removed, and per-tile pixel dimensions are derived
#' from the tile bounding box.  
#' Note: this implementation references the geometry column name \code{"x"} created
#' by \code{sf::st_as_sf()} when converting an \code{sfc}. If your environment uses
#' a different geometry column name, update references (or replace \code{x} with
#' \code{sf::st_geometry(.)} calls).
#'
#'
#' @keywords internal
#' @noRd
.make_request_tiles <- function(aoi, res_m = 20, max_px = 2500,
                               overlap_px = 2, min_area_m2 = 0,
                               crs_metric = 3857) {
  
  cli::cli_progress_step("Creating AOIs for CSDE", spinner = TRUE)
  
  
  stopifnot(inherits(aoi, "sf") || inherits(aoi, "sfc"))
  
  aoi_m <- aoi |>
    sf::st_make_valid() |>
    sf::st_transform(crs_metric)
  
  # Ziel-Kantenlaenge je Kachel (Meter)
  tile_size_m <- (max_px - 2*overlap_px) * res_m
  ov_m <- overlap_px * res_m
  
  # Grid nur ueber die AOI-BBox
  grid <- sf::st_make_grid(aoi_m,
                       cellsize = c(tile_size_m, tile_size_m),
                       what = "polygons", square = TRUE)
  
  #grid_sf <- sf::st_as_sf(grid) |>
    #dplyr::mutate(tile_id = row_number())
  grid_sf <- sf::st_as_sf(grid)
  grid_sf$tile_id <- seq_len(nrow(grid_sf))
  
  # nur Zellen behalten, die AOI schneiden
  has_int <- lengths(sf::st_intersects(grid_sf, aoi_m)) > 0
  grid_sel <- grid_sf[has_int, , drop = FALSE]
  
  grid_sel_buf_geom <- if (ov_m > 0) {
    sf::st_buffer(sf::st_geometry(grid_sel), ov_m)
  } else {
    sf::st_geometry(grid_sel)
  }
  grid_sel_buf <- sf::st_sf(tile_id = grid_sel$tile_id,
                            geometry = grid_sel_buf_geom,
                            crs = sf::st_crs(aoi_m))
  
  suppressWarnings({
    aoi_union <- sf::st_union(sf::st_geometry(aoi_m))
    # auf AOI clippen -> "leere" Bereiche entfallen
    # Perform spatial intersection
    chips_m <- sf::st_intersection(sf::st_as_sf(grid_sel_buf), aoi_union)
    chips_m <- sf::st_make_valid(chips_m)
    chips_m$tile_id <- grid_sel$tile_id
  })
  
  # Tiny slivers entfernen
  #chips_m <- suppressWarnings(
  #  chips_m |>
  #    dplyr::mutate(area_m2 = as.numeric(sf::st_area(sf::st_geometry(chips_m)))) |>
  #    dplyr::filter(dplyr::all_of("area_m2") >= min_area_m2) |>
  #    sf::st_as_sf()
  #)
  suppressWarnings(
    {
    chips_m$area_m2 <- as.numeric(sf::st_area(sf::st_geometry(chips_m)))
    
    # subset rows safely — preserves sf class
    keep <- !is.na(chips_m$area_m2) & (chips_m$area_m2 >= min_area_m2)
    chips_m <- chips_m[keep, , drop = FALSE]
    chips_m <- sf::st_as_sf(chips_m)
    }
  )  
  
  
  # BBox-Dimensionen (Meter) & erwartete Pixelmasse je Chip
  # chips_m <- chips_m |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(
  #     bb = list(sf::st_bbox(sf::st_geometry(chips_m))),
  #     width_m  = bb$xmax - bb$xmin,
  #     height_m = bb$ymax - bb$ymin,
  #     width_px  = ceiling(width_m  / res_m),
  #     height_px = ceiling(height_m / res_m)
  #   ) |>
  #   # rename(geometry = x) |>
  #   dplyr::ungroup() |>
  #   dplyr::select(dplyr::all_of(c("tile_id", "area_m2", "width_m", "height_m", "width_px", "height_px", "geometry")))
  # 
  # compute as above
  geoms <- sf::st_geometry(chips_m)
  bboxes <- lapply(geoms, sf::st_bbox)
  width_m  <- vapply(bboxes, function(bb) as.numeric(bb["xmax"] - bb["xmin"]), numeric(1))
  height_m <- vapply(bboxes, function(bb) as.numeric(bb["ymax"] - bb["ymin"]), numeric(1))
  width_px  <- ceiling(width_m / res_m)
  height_px <- ceiling(height_m / res_m)
  
  # chips_m <- chips_m |>
  #   dplyr::mutate(
  #     width_m  = width_m,
  #     height_m = height_m,
  #     width_px  = width_px,
  #     height_px = height_px,
  #     bb = bboxes
  #   ) |>
  #   dplyr::select(dplyr::all_of(c("tile_id", "area_m2", "width_m", "height_m", "width_px", "height_px", attr(chips_m, "sf_column"))))
  # # Add/overwrite columns
  chips_m$width_m   <- width_m
  chips_m$height_m  <- height_m
  chips_m$width_px  <- width_px
  chips_m$height_px <- height_px
  chips_m$bb        <- bboxes
  
  # Dynamically get the geometry column name
  geom_col <- attr(chips_m, "sf_column")
  
  # Select only specified columns
  cols_to_keep <- c("tile_id", "area_m2", "width_m", "height_m", "width_px", "height_px", geom_col)
  chips_m <- chips_m[, intersect(cols_to_keep, names(chips_m)), drop = FALSE]
  
  
  # Hinweis, falls etwas noch zu gross ist
  if (any(chips_m$width_px > max_px | chips_m$height_px > max_px)) {
    cli::cli_alert_warning("Some chips exceed max_px - reduce tile_size_m (e.g. increase overlap_px or reduce max_px)")
  }
  
  # Fuer die Anfrage ggf. wieder nach 4326 oder direkt im metrischen CRS lassen,
  # je nachdem ob dein Client resx/resy in m (metrisches CRS) oder Grad (4326) erwartet.
  chips_4326 <- sf::st_transform(chips_m, 4326)
  
  cli::cli_progress_done()
  
  list(
    chips_metric = chips_m,   # EPSG:3857 (Pixelmasse stimmen zu res_m in m)
    chips_wgs84  = chips_4326 # EPSG:4326 (falls dein Client 4326-Geometrie braucht)
  )
}


#' Path to the js file
#'
#' @keywords internal
#' @noRd
.ndvi_evalscript_path <- function() {
  return("inst/ext/ndvi_default_2band.js")
}



#' Fetch a single chip from CDSE and write to disk
#'
#' @keywords internal
#' @noRd
.cdse_fetch_chip <- function(aoi_wgs84,
                            from,
                            to, 
                            script_file, 
                            out_file,
                            collection = "sentinel-2-l2a",
                            resolution = 20,
                            mosaicking_order = "leastCC",
                            mask = TRUE,
                            client = NULL,
                            format = "image/tiff",
                            buffer = 0) {
  if (!requireNamespace("CDSE", quietly = TRUE)) {
    cli::cli_abort("Package 'CDSE' is required but not installed.")
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  # coerce AOI to sf (EPSG:4326)
  if (inherits(aoi_wgs84, "sfc")) aoi_wgs84 <- sf::st_as_sf(aoi_wgs84)
  stopifnot(inherits(aoi_wgs84, "sf"))
  if (is.na(sf::st_crs(aoi_wgs84))) {
    aoi_wgs84 <- sf::st_set_crs(aoi_wgs84, 4326)
  } else if (sf::st_crs(aoi_wgs84)$epsg != 4326) {
    aoi_wgs84 <- sf::st_transform(aoi_wgs84, 4326)
  }
  
  # call API
  res <- CDSE::GetImage(
    aoi              = sf::st_geometry(aoi_wgs84),
    time_range       = c(as.Date(from), as.Date(to)),
    collection       = collection,
    script           = script_file,
    format           = format,
    mosaicking_order = mosaicking_order,
    resolution       = resolution,
    mask             = mask,
    client           = client
  )
  
  #print(out_file)
  # write result
  if (is.raw(res)) {
    writeBin(res, out_file)
  } else if (inherits(res, "SpatRaster")) {
    terra::writeRaster(res, out_file, overwrite = TRUE)
  } else {
    # try to coerce to SpatRaster
    maybe <- try(terra::rast(res), silent = TRUE)
    if (inherits(maybe, "try-error")) {
      cli::cli_abort("CDSE::GetImage() returned an unsupported object; cannot write to file.")
    }
    terra::writeRaster(maybe, out_file, overwrite = TRUE)
  }
  
  id <- tryCatch(aoi_wgs84$tile_id[1], error = function(e) NA_integer_)
  tibble::tibble(tile_id = id, file = out_file, status = "ok")
}
