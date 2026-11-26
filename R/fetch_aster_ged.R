#' Fetch and convert ASTER Global Emissivity Dataset tiles
#'
#' Retrieves one-degree AG100 v003 HDF5 tiles covering an area of interest
#' and converts them to seven-band GeoTIFFs in WGS 84. A local folder of AG100
#' files can be used instead of downloading.
#'
#' @param df An object carrying a traktr AOI, for example from [aoi_from_sf()].
#' @param out_dir Output directory for `ag100_v003_<north>_<west>.tif` files.
#'   Created if needed.
#' @param src Optional directory searched recursively for local AG100 HDF5 files.
#' @param keep_h5 Whether to retain downloaded HDF5 files in `out_dir`,
#'   including when GeoTIFF conversion fails.
#' @param overwrite Whether to replace all existing GeoTIFFs. Files with an
#'   incompatible grid or band layout are recreated automatically.
#' @param use_cmr Whether to look up remote tiles via NASA's Common Metadata
#'   Repository before falling back to constructed URLs.
#' @param workers Maximum number of tiles to process in parallel, capped at four.
#'   Parallel execution requires future.apply and a parallel future plan set
#'   by the caller.
#' @param progress Whether to show a tile progress bar.
#' @param aoi Optional `sf` polygon or multipolygon AOI overriding the AOI on `df`.
#'
#' @details
#' The output preserves the product's integer values: elevation in metres
#' (`elev`), NDVI mean and standard deviation as NDVI times 100 (`ndvi_mean`,
#' `ndvi_sd`), temperature mean and standard deviation in centikelvin
#' (`temp_mean`, `temp_sd`), the land/water flag (`lwm`), and observation
#' count (`nobs`). The fill value of -9999 is recorded as NoData.
#'
#' Downloads require Earthdata Login credentials. See [earthdata_auth_check()]
#' for the supported methods.
#'
#' @return Invisibly, `df` with an `aster_ged_tiles` index and an
#'   `aster_ged_map` mapping AOI IDs to tile IDs. The AOI used for fetch and
#'   subsequent attachment is recorded in `attr(df, "aster_ged_aoi")`.
#'   The index and fetch parameters are also recorded in
#'   `attr(df, "traktr_layers")$aster_ged`. Tile statuses
#'   are `downloaded`, `converted` (from a local HDF5 file), `exists`,
#'   `missing_remote`, `missing_local`, and `error`.
#' @examples
#' \dontrun{
#' df <- aoi_from_sf(my_polygons, id_col = "id")
#' df <- fetch_aster_ged(df, out_dir = "data/ag100")
#' }
#' @export
fetch_aster_ged <- function(df, out_dir, src = NULL, keep_h5 = FALSE, overwrite = FALSE,
                            use_cmr = TRUE, workers = 1, progress = TRUE, aoi = NULL) {
  aoi <- .resolve_aoi(df, aoi)
  id_col <- attr(df, "traktr_aoi_params")$id_col
  if (!length(id_col) || !id_col %in% names(aoi)) {
    id_col <- .aoi_detect(structure(aoi, traktr_aoi = aoi))$id_col
  }
  aoi <- sf::st_transform(aoi[, id_col, drop = FALSE], 4326)
  tl <- .ag100_tiles_for_aoi(aoi, id_col)
  tiles <- tl$tiles
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tif <- file.path(out_dir, sprintf("ag100_v003_%d_%d.tif", tiles$lat_n, tiles$lon_w))
  available <- file.exists(tif) & !dir.exists(tif)
  need <- which(isTRUE(overwrite) | !available)
  if (!isTRUE(overwrite)) {
    candidates <- which(available)
    outdated <- candidates[!vapply(candidates, function(i) {
      .ag100_tif_current(tif[i], tiles$lat_n[i], tiles$lon_w[i])
    }, logical(1))]
    if (length(outdated)) {
      tile_ids <- tiles$tile_id[outdated]
      cli::cli_inform("ASTER GED outdated grid, re-creating tile{?s} {.val {tile_ids}}.")
      need <- sort(c(need, outdated))
    }
  }
  index <- tibble::tibble(
    tile_id = tiles$tile_id, lat_n = tiles$lat_n, lon_w = tiles$lon_w,
    file = tif, status = "exists", url_source = NA_character_,
    size_bytes = NA_real_, message = NA_character_
  )
  if (isTRUE(progress) && nrow(index)) {
    progress_id <- cli::cli_progress_bar("ASTER GED", total = nrow(index), auto_terminate = FALSE)
    on.exit(cli::cli_progress_done(id = progress_id), add = TRUE)
    cli::cli_progress_update(set = nrow(index) - length(need), id = progress_id)
  }

  if (length(need)) {
    pending <- tiles[need, ]
    pending$status <- pending$url <- pending$url_source <- NA_character_
    pending$source <- NA_character_
    if (is.null(src)) {
      pending <- .ag100_resolve(pending, use_cmr = use_cmr)
    } else {
      pending$source <- .ag100_local_files(src, pending$tile_id)
    }
    pending$tif <- tif[need]

    workers <- min(4L, as.integer(workers))
    parallel <- workers > 1L && .has_future_apply()
    if (workers > 1L && !parallel) {
      cli::cli_warn("Package {.pkg future.apply} is unavailable. Processing ASTER GED sequentially.")
    }
    handle <- NULL
    if (!parallel && is.null(src) && anyNA(pending$status)) handle <- .earthdata_handle()

    # Limit work in flight so an authentication error stops subsequent batches.
    batch_size <- if (parallel) workers else 1L
    for (start in seq.int(1L, nrow(pending), by = batch_size)) {
      rows <- seq.int(start, min(start + batch_size - 1L, nrow(pending)))
      if (parallel) {
        jobs <- lapply(rows, function(i) pending[i, ])
        results <- future.apply::future_lapply(
          jobs, .ag100_fetch_tile, local = !is.null(src), keep_h5 = keep_h5,
          future.chunk.size = 1L, future.seed = TRUE
        )
      } else {
        results <- list(.ag100_fetch_tile(pending[rows, ], !is.null(src), keep_h5, handle))
      }
      for (j in seq_along(rows)) {
        i <- need[rows[j]]
        result <- results[[j]]
        index$status[i] <- result$status
        index$url_source[i] <- pending$url_source[rows[j]]
        index$message[i] <- result$message
        if (result$status == "error") {
          cli::cli_warn("ASTER GED tile {index$tile_id[i]} failed: {result$message}")
        }
      }
      if (isTRUE(progress)) cli::cli_progress_update(inc = length(rows), id = progress_id)
    }
  }

  ok <- index$status %in% c("exists", "downloaded", "converted")
  index$file[!ok] <- NA_character_
  index$size_bytes[ok] <- file.size(index$file[ok])
  attr(df, "aster_ged_tiles") <- index
  attr(df, "aster_ged_map") <- tl$map
  attr(df, "aster_ged_aoi") <- aoi
  df <- .set_traktr_layer(df, "aster_ged", index, params = list(
    product = "AG100", version = "003", concept_id = "C2763266348-LPCLOUD",
    out_dir = out_dir, src = src, overwrite = overwrite, retrieved = Sys.time()
  ))
  invisible(df)
}

.ag100_local_files <- function(src, tile_ids) {
  pattern <- "^AG100\\.v003\\.0*(-?\\d+)\\.0*(-?\\d+)\\.0001.*\\.h5$"
  paths <- list.files(src, pattern = pattern, recursive = TRUE, full.names = TRUE)
  # Prefer a file at the root over a duplicate in a subdirectory.
  paths <- paths[order(nchar(gsub("[^/]", "", paths)))]
  coords <- utils::strcapture(pattern, basename(paths), list(lat = integer(), lon = integer()))
  paths[match(tile_ids, paste(coords$lat, coords$lon, sep = "."))]
}

.has_future_apply <- function() requireNamespace("future.apply", quietly = TRUE)

.sleep <- function(seconds) Sys.sleep(seconds)

.ag100_fetch_tile <- function(tile, local, keep_h5, handle = NULL) {
  if (identical(tile$status, "missing_remote")) {
    return(list(status = "missing_remote", message = "Tile unavailable in the remote catalogue."))
  }
  if (identical(tile$status, "error")) {
    return(list(status = "error", message = "Could not resolve an AG100 download URL."))
  }
  if (local && is.na(tile$source)) {
    return(list(status = "missing_local", message = "No matching local HDF5 file."))
  }

  if (!local && is.null(handle)) handle <- .earthdata_handle()
  tryCatch({
    h5 <- tile$source
    if (!local) {
      staged_h5 <- tempfile("ag100-", tmpdir = dirname(tile$tif), fileext = ".h5")
      on.exit(unlink(staged_h5), add = TRUE)
      .ag100_download_retry(tile$url, staged_h5, handle)
      h5 <- staged_h5
    }
    conversion_error <- tryCatch({
      .ag100_h5_to_tif(h5, tile$tif, tile$lat_n, tile$lon_w)
      NULL
    }, error = identity)
    msg <- if (is.null(conversion_error)) NA_character_ else .ag100_error_message(conversion_error)
    if (!local && isTRUE(keep_h5)) {
      kept <- sub("\\.tif$", ".h5", tile$tif)
      if (!is.null(conversion_error) && file.exists(kept)) {
        kept <- tempfile(sub("\\.h5$", "-failed-", basename(kept)),
                         tmpdir = dirname(kept), fileext = ".h5")
      }
      if (!suppressWarnings(file.rename(staged_h5, kept))) {
        retention_error <- paste0("The HDF5 file could not be retained at ", kept, ".")
        msg <- if (is.null(conversion_error)) retention_error else paste(msg, retention_error)
        if (is.null(conversion_error)) cli::cli_warn("{msg}")
      } else if (!is.null(conversion_error)) {
        msg <- paste0(msg, " Downloaded HDF5 retained as ", basename(kept), ".")
      }
    }
    if (!is.null(conversion_error)) return(list(status = "error", message = msg))
    list(status = if (local) "converted" else "downloaded", message = msg)
  }, error = function(e) {
    if (inherits(e, "ag100_auth_error")) stop(e)
    list(status = "error", message = .ag100_error_message(e))
  })
}

.ag100_download <- function(url, dest_h5, handle) {
  tmp <- tempfile("ag100-download-", tmpdir = dirname(dest_h5), fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  curl::curl_download(url, tmp, handle = handle, mode = "wb", quiet = TRUE)
  if (!file.rename(tmp, dest_h5)) cli::cli_abort("Could not stage the downloaded HDF5 file.")
}

.ag100_download_retry <- function(url, dest_h5, handle) {
  signature <- as.raw(c(137, 72, 68, 70, 13, 10, 26, 10))
  for (attempt in seq_len(3L)) {
    error <- tryCatch({
      .ag100_download(url, dest_h5, handle)
      NULL
    }, error = identity)
    if (is.null(error)) {
      if (!identical(readBin(dest_h5, "raw", n = 8L), signature)) .ag100_auth_abort()
      return(invisible(dest_h5))
    }
    auth_status <- "(HTTP(/[0-9.]+)?( status( code)?)?[: ]+|returned error: )(401|403)\\b"
    if (any(error$status_code %in% c(401L, 403L)) ||
        grepl(auth_status, conditionMessage(error), ignore.case = TRUE)) .ag100_auth_abort()
    if (attempt == 3L) stop(error)
    unlink(dest_h5)
    .sleep(2^attempt)
  }
}

.ag100_auth_abort <- function() {
  cli::cli_abort(.earthdata_auth_help(), class = "ag100_auth_error")
}

.ag100_error_message <- function(error) {
  message <- conditionMessage(error)
  token <- Sys.getenv("EARTHDATA_TOKEN")
  if (nzchar(token)) message <- gsub(token, "[redacted]", message, fixed = TRUE)
  message
}

.ag100_read_array <- function(address) {
  if (file.exists(address)) address <- normalizePath(address, winslash = "/")
  info <- sf::gdal_utils("info", address, options = "-nomd", quiet = TRUE)
  lines <- strsplit(info, "\n", fixed = TRUE)[[1L]]
  size <- utils::strcapture("^Size is ([0-9]+), ([0-9]+)$",
                           lines[startsWith(lines, "Size is ")],
                           list(ncol = integer(), nrow = integer()))
  if (nrow(size) != 1L || anyNA(size) || any(size <= 0L)) {
    cli::cli_abort("Could not determine AG100 raster dimensions.")
  }
  bands <- grep("^Band [0-9]+ ", lines, value = TRUE)
  if (length(bands) != 1L) cli::cli_abort("Expected one raster layer.")
  type <- utils::strcapture(" Type=([[:alnum:]]+)", bands, list(type = character()))$type
  if (is.na(type)) cli::cli_abort("Could not determine the AG100 raster data type.")
  nodata <- grep("^NoData Value=", trimws(lines), value = TRUE)
  nodata <- sub("^NoData Value=", "", nodata)

  vrt <- tempfile("ag100-", fileext = ".vrt")
  on.exit(unlink(vrt), add = TRUE)
  # Keep the source's NoData marker. Read the remaining pixels without
  # borrowing the source transform, projection, scale or offset.
  source <- gsub("&", "&amp;", address, fixed = TRUE)
  source <- gsub("<", "&lt;", source, fixed = TRUE)
  source <- gsub(">", "&gt;", source, fixed = TRUE)
  writeLines(c(
    sprintf('<VRTDataset rasterXSize="%d" rasterYSize="%d">', size$ncol, size$nrow),
    sprintf("  <GeoTransform>0,1,0,%d,0,-1</GeoTransform>", size$nrow),
    sprintf('  <VRTRasterBand dataType="%s" band="1">', type),
    if (length(nodata)) sprintf("    <NoDataValue>%s</NoDataValue>", nodata),
    "    <SimpleSource>",
    sprintf('      <SourceFilename relativeToVRT="0">%s</SourceFilename>', source),
    "      <SourceBand>1</SourceBand>",
    "    </SimpleSource>",
    "  </VRTRasterBand>",
    "</VRTDataset>"
  ), vrt)
  r <- withCallingHandlers(terra::rast(vrt), warning = function(w) {
    if (grepl("unknown extent|the data in this file are rotated", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
  matrix(terra::values(r, mat = FALSE), size$nrow, size$ncol, byrow = TRUE)
}

.ag100_band_paths <- c(
  elev = "ASTER GDEM/ASTGDEM", ndvi_mean = "NDVI/Mean", ndvi_sd = "NDVI/SDev",
  temp_mean = "Temperature/Mean", temp_sd = "Temperature/SDev",
  lwm = "Land Water Map/LWmap", nobs = "Observations/NumObs"
)

.ag100_tif_current <- function(tif, lat_n, lon_w) {
  tryCatch({
    r <- suppressWarnings(terra::rast(tif, guessCRS = FALSE))
    nr <- nrow(r)
    nc <- ncol(r)
    if (min(nr, nc) < 2L || isTRUE(terra::is.rotated(r)) ||
        !identical(names(r), names(.ag100_band_paths)) ||
        !isTRUE(sf::st_crs(terra::crs(r)) == sf::st_crs(4326))) return(FALSE)

    res <- 1 / c(nc - 1L, nr - 1L)
    bounds <- .ag100_grid_bounds(lat_n, lon_w, res)
    isTRUE(all(abs(terra::res(r) - res) <= 1e-9) &&
             all(abs(as.vector(terra::ext(r)) - bounds) <= 1e-9))
  }, error = function(e) FALSE)
}

.ag100_h5_to_tif <- function(h5, tif, lat_n, lon_w, tol = 1e-6) {
  bands <- .ag100_band_paths
  h5 <- normalizePath(h5, winslash = "/")
  subdatasets <- terra::describe(h5, sds = TRUE)
  read_dataset <- function(path, geolocation = FALSE) {
    pattern <- paste0("/", gsub(" ", "[ _]", path, fixed = TRUE), "$")
    address <- subdatasets$name[grepl(pattern, subdatasets$name)][1L]
    if (is.na(address) && geolocation) {
      # GDAL can open rank-one arrays directly when the listing omits them.
      address <- sprintf('HDF5:"%s"://%s', h5, path)
    }
    if (is.na(address)) cli::cli_abort("Missing AG100 dataset {.val {path}}.")
    tryCatch(
      .ag100_read_array(address),
      error = function(e) cli::cli_abort("Could not read AG100 dataset {.val {path}}.", parent = e)
    )
  }

  first <- read_dataset(bands[[1L]])
  nr <- nrow(first)
  nc <- ncol(first)
  lat <- read_dataset("Geolocation/Latitude", TRUE)
  lon <- read_dataset("Geolocation/Longitude", TRUE)
  if (any(!is.finite(lat))) cli::cli_abort("Invalid AG100 Geolocation/Latitude values.")
  if (any(!is.finite(lon))) cli::cli_abort("Invalid AG100 Geolocation/Longitude values.")

  vectors <- min(dim(lat)) == 1L && min(dim(lon)) == 1L
  transposed <- if (vectors) {
    nr != nc && length(lat) == nc && length(lon) == nr
  } else {
    identical(dim(lat), dim(first)) && identical(dim(lon), dim(first)) &&
      any(abs(lat[1L, ] - lat[1L, 1L]) > tol) && any(abs(lon[, 1L] - lon[1L, 1L]) > tol) &&
      all(abs(t(lat) - lat[1L, ]) <= tol) && all(abs(lon - lon[, 1L]) <= tol)
  }
  if (transposed) {
    cli::cli_inform("AG100 tile {lat_n}.{lon_w}: transposed layout.")
    lat <- t(lat)
    lon <- t(lon)
  }
  output_dim <- if (transposed) c(nc, nr) else c(nr, nc)
  lat <- .ag100_coordinate(lat, output_dim, "Geolocation/Latitude", tol)
  lon <- .ag100_coordinate(lon, output_dim, "Geolocation/Longitude", tol)
  bounds <- .ag100_extent(lat, lon, lat_n, lon_w, check_swap = vectors && nr == nc)

  raw <- matrix(NA_real_, nr * nc, length(bands))
  for (i in seq_along(bands)) {
    path <- bands[[i]]
    band <- if (i == 1L) first else read_dataset(path)
    if (!identical(dim(band), c(nr, nc))) {
      cli::cli_abort("Unexpected dimensions in AG100 dataset {.val {path}}.")
    }
    band[band == -9999] <- NA
    if (any(!is.na(band) & (!is.finite(band) | band != trunc(band)))) {
      cli::cli_abort("Values in AG100 dataset {.val {path}} must be integers.")
    }
    if (any(band < -2147483648 | band > 2147483647, na.rm = TRUE)) {
      cli::cli_abort("Values in AG100 dataset {.val {path}} must fit in signed 32-bit integers.")
    }
    if (transposed) band <- t(band)
    raw[, i] <- as.vector(t(band))
  }
  r <- terra::rast(nrows = output_dim[1L], ncols = output_dim[2L], nlyrs = length(bands),
                   extent = terra::ext(bounds), crs = "EPSG:4326")
  r <- terra::setValues(r, raw)
  names(r) <- names(bands)

  tmp_tif <- tempfile("ag100-", tmpdir = dirname(tif), fileext = ".tif")
  on.exit(unlink(tmp_tif), add = TRUE)
  terra::writeRaster(r, tmp_tif, datatype = "INT4S", NAflag = -9999,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"), overwrite = TRUE)
  if (!file.rename(tmp_tif, tif)) cli::cli_abort("Could not write {.file {tif}}.")
  invisible(tif)
}

.ag100_extent <- function(lat, lon, lat_n, lon_w, check_swap = FALSE) {
  endpoint_tol <- 1e-4
  observed <- c(lat[1L], lat[length(lat)], lon[1L], lon[length(lon)])
  expected <- c(lat_n, lat_n - 1, lon_w, lon_w + 1)
  endpoints <- function(x) {
    sprintf("Geolocation/Latitude %.9g to %.9g, Geolocation/Longitude %.9g to %.9g",
            x[1L], x[2L], x[3L], x[4L])
  }
  lat_step <- -diff(lat)
  lon_step <- diff(lon)
  valid_direction <- length(lat_step) && length(lon_step) &&
    all(lat_step > 0) && all(lon_step > 0)
  if (!valid_direction || any(abs(observed - expected) > endpoint_tol)) {
    hint <- if (check_swap) " Check direction and extent: lat/lon axes swapped?" else ""
    cli::cli_abort(paste0(
      "Invalid AG100 extent or orientation. Observed endpoints: {endpoints(observed)}. ",
      "Expected endpoints: {endpoints(expected)}.", hint
    ))
  }
  lat_res <- 1 / length(lat_step)
  lon_res <- 1 / length(lon_step)
  # Float32 coordinates can round successive steps differently.
  if (any(abs(lat_step - lat_res) > lat_res / 100) ||
      any(abs(lon_step - lon_res) > lon_res / 100) ||
      abs(lat_res - lon_res) > min(lat_res, lon_res) / 100) {
    cli::cli_abort("AG100 Geolocation/Latitude and Geolocation/Longitude must have a constant step on both axes.")
  }

  .ag100_grid_bounds(lat_n, lon_w, c(lon_res, lat_res))
}

.ag100_grid_bounds <- function(lat_n, lon_w, res) {
  half <- res / 2
  c(lon_w - half[1L], lon_w + 1 + half[1L],
    lat_n - 1 - half[2L], lat_n + half[2L])
}

.ag100_coordinate <- function(grid, dimensions, path, tol) {
  latitude <- path == "Geolocation/Latitude"
  expected <- dimensions[if (latitude) 1L else 2L]
  if (min(dim(grid)) == 1L && length(grid) == expected) return(as.vector(grid))
  if (!identical(dim(grid), dimensions)) {
    cli::cli_abort("Unexpected dimensions for AG100 dataset {.val {path}}.")
  }
  if (!latitude) grid <- t(grid)
  if (any(abs(grid - grid[, 1L]) > tol)) {
    cli::cli_abort("AG100 dataset {.val {path}} must be constant across the other axis.")
  }
  grid[, 1L]
}
