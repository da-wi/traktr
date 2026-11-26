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
#' @param keep_h5 Whether to retain downloaded HDF5 files in `out_dir`.
#' @param overwrite Whether to replace existing GeoTIFFs.
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
  need <- which(isTRUE(overwrite) | !file.exists(tif) | dir.exists(tif))
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
    .ag100_h5_to_tif(h5, tile$tif, tile$lat_n, tile$lon_w)
    # The GeoTIFF is the product; a failed keep_h5 rename only warns, so the
    # tile status stays consistent with the published (cached) GeoTIFF.
    msg <- NA_character_
    if (!local && isTRUE(keep_h5)) {
      kept <- sub("\\.tif$", ".h5", tile$tif)
      if (!suppressWarnings(file.rename(staged_h5, kept))) {
        msg <- paste0("GeoTIFF written, but the HDF5 file could not be retained at ", kept, ".")
        cli::cli_warn(msg)
      }
    }
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

.ag100_h5_to_tif <- function(h5, tif, lat_n, lon_w, tol = 1e-6) {
  bands <- c(elev = "ASTER GDEM/ASTGDEM", ndvi_mean = "NDVI/Mean", ndvi_sd = "NDVI/SDev",
             temp_mean = "Temperature/Mean", temp_sd = "Temperature/SDev",
             lwm = "Land Water Map/LWmap", nobs = "Observations/NumObs")
  subdatasets <- terra::describe(h5, sds = TRUE)
  read_dataset <- function(path, geolocation = FALSE) {
    pattern <- paste0("/", gsub(" ", "[ _]", path, fixed = TRUE), "$")
    address <- subdatasets$name[grepl(pattern, subdatasets$name)][1L]
    if (is.na(address) && geolocation) {
      # GDAL can open rank-one arrays directly when the listing omits them.
      address <- sprintf('HDF5:"%s"://%s', normalizePath(h5, winslash = "/"), path)
    }
    if (is.na(address)) cli::cli_abort("Missing AG100 dataset {.val {path}}.")
    r <- tryCatch(
      withCallingHandlers(terra::rast(address), warning = function(w) {
        if (grepl("unknown extent", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }),
      error = function(e) cli::cli_abort("Could not read AG100 dataset {.val {path}}.", parent = e)
    )
    if (terra::nlyr(r) != 1L) cli::cli_abort("Expected one layer in AG100 dataset {.val {path}}.")
    if (!geolocation) terra::scoff(r) <- NULL
    r
  }

  rasters <- lapply(bands, read_dataset)
  r <- do.call(c, unname(rasters))
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  lat <- .ag100_coordinate(read_dataset("Geolocation/Latitude", TRUE), nr, nc, "latitude", tol)
  lon <- .ag100_coordinate(read_dataset("Geolocation/Longitude", TRUE), nr, nc, "longitude", tol)
  bounds <- .ag100_extent(lat, lon, lat_n, lon_w)

  raw <- terra::values(r, mat = TRUE)
  if (any(!is.na(raw) & (!is.finite(raw) | raw != trunc(raw)))) {
    cli::cli_abort("AG100 band values must be integers.")
  }
  if (any(raw < -2147483648 | raw > 2147483647, na.rm = TRUE)) {
    cli::cli_abort("AG100 band values must fit in signed 32-bit integers.")
  }
  raw[raw == -9999] <- NA
  r <- terra::setValues(r, raw)
  terra::ext(r) <- bounds
  terra::crs(r) <- "EPSG:4326"
  names(r) <- names(bands)

  tmp_tif <- tempfile("ag100-", tmpdir = dirname(tif), fileext = ".tif")
  on.exit(unlink(tmp_tif), add = TRUE)
  terra::writeRaster(r, tmp_tif, datatype = "INT4S", NAflag = -9999,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"), overwrite = TRUE)
  if (!file.rename(tmp_tif, tif)) cli::cli_abort("Could not write {.file {tif}}.")
  invisible(tif)
}

.ag100_extent <- function(lat, lon, lat_n, lon_w) {
  lat_step <- -diff(lat)
  lon_step <- diff(lon)
  if (!length(lat_step) || !length(lon_step) || any(lat_step <= 0) || any(lon_step <= 0)) {
    cli::cli_abort("Invalid AG100 orientation: latitude must be decreasing and longitude increasing.")
  }
  lat_res <- mean(lat_step)
  lon_res <- mean(lon_step)
  # Float32 coordinates can round successive steps differently.
  if (any(abs(lat_step - lat_res) > lat_res / 100) ||
      any(abs(lon_step - lon_res) > lon_res / 100) ||
      abs(lat_res - lon_res) > min(lat_res, lon_res) / 100) {
    cli::cli_abort("AG100 coordinates must have a constant step on both axes.")
  }

  west <- min(lon) - lon_res / 2
  east <- max(lon) + lon_res / 2
  north <- max(lat) + lat_res / 2
  south <- min(lat) - lat_res / 2
  if (abs(north - lat_n) > 1e-6 * length(lat) + 1e-9 || abs(west - lon_w) > lon_res / 10 ||
      abs(east - west - 1) > lon_res / 10 || abs(north - south - 1) > lat_res / 10) {
    cli::cli_abort("AG100 extent does not match the expected one-degree tile at {lat_n}, {lon_w}.")
  }
  c(west, east, south, north)
}

.ag100_coordinate <- function(r, nr, nc, axis, tol) {
  coord <- terra::values(r, mat = FALSE)
  latitude <- axis == "latitude"
  expected <- if (latitude) nr else nc
  if (any(!is.finite(coord))) cli::cli_abort("Invalid AG100 {axis} values.")
  if (min(dim(r)[1:2]) == 1L && length(coord) == expected) return(coord)
  if (!all(dim(r)[1:2] == c(nr, nc))) {
    cli::cli_abort("Unexpected dimensions for AG100 {axis}.")
  }
  grid <- matrix(coord, nr, nc, byrow = TRUE)
  if (!latitude) grid <- t(grid)
  if (any(abs(grid - grid[, 1L]) > tol)) {
    cli::cli_abort("AG100 {axis} must be constant across the other axis.")
  }
  grid[, 1L]
}
