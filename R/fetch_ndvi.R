#' Fetch NDVI chips for an AOI
#'
#' Splits the AOI from `df` (in `attr(df,"ndvi_aoi")`) into request tiles,
#' downloads Sentinel-2 NDVI chips from the Copernicus Data Space Ecosystem (CDSE),
#' and stores a tile index as `attr(df,"ndvi_tiles")`. A record of fetch parameters
#' is stored in `attr(df,"ndvi_fetch_params")`.  
#' 
#' This function is AOI-agnostic: it works for H3-derived AOIs
#' (via [mutate_aoi_h3()]) or buffer-based AOIs (via [mutate_aoi_buffer()]).
#'
#' @param df A data frame previously passed through [mutate_aoi_h3()] or
#'   [mutate_aoi_buffer()], i.e. with an AOI stored in `attr(df,"ndvi_aoi")`.
#' @param from,to Date range for Sentinel-2 search (inclusive). Must be coercible
#'   to `Date`.
#' @param out_dir Directory to which NDVI GeoTIFF chips will be written.
#' @param res_m Numeric ground resolution in meters per pixel (default `20`).
#'   This is passed both to tile generation and to CDSE `resolution`.
#' @param max_px Maximum allowed pixel width/height per request (default `2500`).
#'   Used when splitting AOIs into tiles.
#' @param overlap_px Overlap (in pixels) to add around each tile (default `0`).
#' @param evalscript Path to SentinelHub evalscript file.
#' @param client OAuth client object for CDSE, e.g. from
#'   `CDSE::GetOAuthClient(id=..., secret=...)`.
#' @param rate_n,rate_period Passed to \pkg{ratelimitr} to throttle requests
#'   (default 2 requests per second).
#' @param overwrite Logical; if `FALSE` (default), existing `.tif` files in
#'   `out_dir` are skipped. If `TRUE`, files are re-fetched.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Creates request tiles on the AOI.
#'   \item Iterates through each tile and fetches the satellite data.
#'   \item Applies a rate limit to avoid hitting API caps.
#'   \item Builds an `sf` index of all tiles with columns:
#'     `tile_id`, `geometry`, `file`, `status`, `msg`.
#' }
#'
#' The tile index is attached as `attr(df,"ndvi_tiles")`.
#' Parameters used are stored as a list in `attr(df,"ndvi_fetch_params")`.
#'
#' Progress is displayed using \pkg{cli}.
#'
#' @return The input `df`, with additional attributes:
#' * `ndvi_tiles` - an `sf` object of tiles with download status and file paths.
#' * `ndvi_fetch_params` - list of parameters used.
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(lat = c(47.38,47.39), lon = c(8.53,8.54)) |>
#'   mutate_aoi_h3(lat, lon, res = 10)
#'
#' oauth <- CDSE::GetOAuthClient(
#'   id     = Sys.getenv("CDSE_CLIENT_ID"),
#'   secret = Sys.getenv("CDSE_CLIENT_SECRET")
#' )
#'
#' df <- fetch_ndvi(df,
#'                  from = as.Date("2025-04-01"),
#'                  to   = as.Date("2025-04-30"),
#'                  out_dir = "output/ndvi",
#'                  res_m   = 20,
#'                  client  = oauth)
#'
#' attr(df, "ndvi_tiles")
#' }
#'
#' @seealso [mutate_aoi_h3()], [mutate_aoi_buffer()], [attach_ndvi()],
#' @export
fetch_ndvi <- function(df, from, to, out_dir,
                       res_m = 10, max_px = 2500, overlap_px = 0,
                       evalscript = ndvi_evalscript_path(),
                       client, rate_n = 2, rate_period = 1,
                       overwrite = FALSE) {
  rlang::check_installed(c("cli","ratelimitr","purrr","tibble","dplyr","sf","glue"))
  
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    rlang::abort("No AOI found: please call a mutate_aoi_*() function first.")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1) Tiles erzeugen
  tiles <- make_request_tiles(aoi, res_m = res_m, max_px = max_px, overlap_px = overlap_px)
  n_tiles <- nrow(tiles$chips_wgs84)
  if (!n_tiles) rlang::abort("No tiles to request (AOI too small or filtered).")
  idx <- seq_len(n_tiles)
  
  sanitize <- function(x) gsub("[-.]", "_", format(x, scientific = FALSE))
  
  # 2) Safe + Rate-limit
  safe_one <- purrr::safely(function(i) {
    bbox <- sf::st_bbox(tiles$chips_wgs84[i, ])
    
    filename <- glue::glue(
      "ndvi_{sanitize(bbox$xmin)}_{sanitize(bbox$ymin)}_{sanitize(bbox$xmax)}_{sanitize(bbox$ymax)}_{as.character(as.Date(from))}_{as.character(as.Date(to))}.tif"
    )
    
    out_file <- file.path(out_dir, filename)
    
    #out_file <- file.path(out_dir, sprintf("ndvi_%06d.tif", tiles$chips_wgs84$tile_id[i]))
    
    cdse_fetch_chip(tiles$chips_wgs84[i, ], from, to, evalscript, out_file,  resolution = res_m, client = client)
  })
  
  
  limited <- ratelimitr::limit_rate(safe_one, ratelimitr::rate(n = rate_n, period = rate_period))
  
  # 3) Fortschritt
  skipped <- 0
  cli::cli_progress_step("Fetching data from CDSE",
                         msg_done = "Fetching NDVI data completed", spinner = TRUE)
  
  pb <- cli::cli_progress_bar(total = n_tiles,
                              format = "Downloading {cli::pb_bar} {cli::pb_percent} (skipped {skipped})",
                              clear = FALSE)
  
  res_list <- vector("list", n_tiles)
  for (j in seq_along(idx)) {
    i <- idx[j]
    
    bbox <- sf::st_bbox(tiles$chips_wgs84[i, ])
    
    filename <- glue::glue(
      "ndvi_{sanitize(bbox$xmin)}_{sanitize(bbox$ymin)}_{sanitize(bbox$xmax)}_{sanitize(bbox$ymax)}_{as.character(as.Date(from))}_{as.character(as.Date(to))}.tif"
    )
    
    out_file <- file.path(out_dir, filename)
    
    #out_file <- file.path(out_dir, sprintf("ndvi_%06d.tif", tiles$chips_wgs84$tile_id[i]))
    
    if (!overwrite && file.exists(out_file)) {
      skipped <- skipped + 1
      res_list[[j]] <- tibble::tibble(tile_id = tiles$chips_wgs84$tile_id[i],
                                      file = out_file, status = "skipped")
    } else {
      if (is.null(client )) {
        cli::cli_abort("No OAuth client set. ")
      } else {
        ans <- limited(i)
      }
      
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
    cli::cli_progress_update(id = pb)
  }
  cli::cli_progress_done(id = pb)
  
  results <- dplyr::bind_rows(res_list)
  
  if ("id" %in% names(results) && !"tile_id" %in% names(results)) {
    #results <- dplyr::rename(results, tile_id = id)
    results$tile_id <- results$id
    #results <- subset(results, select = -id)
    results <- results[ , !(names(results) %in% "id")]
    
  }
  
  tile_index <- dplyr::left_join(
    tiles$chips_wgs84 |> dplyr::select(dplyr::any_of(c("tile_id", "geometry"))),
    results |> dplyr::select(dplyr::any_of(c("tile_id","file","status","msg"))),
    by = "tile_id"
  )
  tile_index <- sf::st_as_sf(tile_index)
  attr(df, "ndvi_aoi") <- tile_index
  attr(df, "ndvi_fetch_params") <- list(
    from = as.Date(from), to = as.Date(to),
    out_dir = out_dir, res_m = res_m, max_px = max_px, overlap_px = overlap_px,
    evalscript = evalscript, rate_n = rate_n, rate_period = rate_period,
    overwrite = overwrite
  )
  class(df) <- union("traktr_df", class(df))
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
#' @examples
#' \dontrun{
#' library(sf)
#' poly <- st_polygon(list(rbind(c(8,47), c(8.5,47), c(8.5,47.5), c(8,47.5), c(8,47))))
#' aoi  <- st_as_sf(st_sfc(poly, crs = 4326))
#' tiles <- make_request_tiles(aoi, res_m = 10, max_px = 2500, overlap_px = 2)
#' nrow(tiles$chips_metric)
#' }
#'
#' @importFrom sf st_make_valid st_transform st_make_grid st_as_sf st_intersects
#' @importFrom sf st_buffer st_intersection st_area st_bbox
#' @importFrom dplyr mutate row_number filter select rowwise ungroup
#' @keywords internal
#' @noRd
make_request_tiles <- function(aoi, res_m = 20, max_px = 2500,
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
  
  grid_sf <- sf::st_as_sf(grid) |>
    dplyr::mutate(tile_id = row_number())
  
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
    chips_m <- sf::st_intersection(st_as_sf(grid_sel_buf), aoi_union) |>
      sf::st_make_valid() |>
      dplyr::mutate(tile_id = grid_sel$tile_id)
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
  
  chips_m <- chips_m |>
    dplyr::mutate(
      width_m  = width_m,
      height_m = height_m,
      width_px  = width_px,
      height_px = height_px,
      bb = bboxes
    ) |>
    dplyr::select(dplyr::all_of(c("tile_id", "area_m2", "width_m", "height_m", "width_px", "height_px", attr(chips_m, "sf_column"))))
  
  
  
  # Hinweis, falls etwas noch zu gross ist
  if (any(chips_m$width_px > max_px | chips_m$height_px > max_px)) {
    cli::cli_alert_warning("Some chips exceed max_px - reduce tile_size_m (e.g. increase overlap_px or reduce max_px) and call `make_request_tiles` again.")
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
ndvi_evalscript_path <- function() {
  return("inst/ext/ndvi_default_2band.js")
}



#' Fetch a single chip from CDSE and write to disk
#'
#' @keywords internal
#' @noRd
cdse_fetch_chip <- function(aoi_wgs84,
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
