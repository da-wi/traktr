#' Attach area-weighted ASTER GED statistics
#'
#' Summarizes raster cells from [fetch_aster_ged()] within each AOI polygon,
#' weighting valid values by the area covered. Joins the results by AOI ID
#' and preserves input row order and AOI metadata.
#'
#' @param df Output of [fetch_aster_ged()].
#' @param bands Character vector selecting one or more ASTER GED bands.
#' @param naming Column naming scheme. `"county"` uses the ecolo-zip names
#'   with `ASTER_county_elevation_*` for elevation. `"ecolozip"` uses
#'   `ASTER_zip_elevation_*`. `"traktr"` uses lower-case names such as
#'   `aster_elev_mean`.
#' @param extras Character vector selecting area-weighted medians (except for
#'   the land/water band), the number of valid cells in elevation or the first
#'   requested band, and geodesic AOI area in square kilometres. Use
#'   `character()` to omit all extras.
#' @param convert Convert NDVI from NDVI times 100 to NDVI, and temperature
#'   from centikelvin to degrees Celsius for temperature means and kelvin
#'   for temperature standard deviations and spatial SDs of temperature means.
#' @param ecolozip_na_quirk Report zero missing percentages as `NA`, matching
#'   the ecolo-zip convention.
#' @param allow_missing Warn and skip failed or missing local tiles if
#'   `TRUE`. Otherwise stop and report their tile IDs. Tiles recorded as
#'   missing remotely are skipped silently.
#' @param workers Maximum number of tiles processed in parallel. Parallel
#'   execution requires future.apply and a parallel future plan set by the caller.
#' @param progress Whether to show a tile extraction progress bar.
#' @param cache Whether to cache per-tile moments in the `_moments`
#'   subdirectory of the tile directory.
#'
#' @details
#' Extraction and geodesic area use the AOI recorded by [fetch_aster_ged()],
#' including any AOI override. Older fetch results use the original traktr AOI.
#' Rows whose ID is absent from the AOI receive `NA` for all summary columns.
#'
#' Standard deviations use the area-weighted population formula. Missing
#' percentages refer to the area covered by available tiles and are `NA`
#' where there is no coverage. Elevation `cin` measures the percentage of
#' valid area at exactly zero metres. Elevation `sum` is always `NA` in
#' polygon mode. The observation flag `lt_4` marks a mean below four.
#' Land/water summaries report the percentage of valid area classified as
#' water and, as `cin`, the percentage of covered area with valid values.
#' Existing columns with the selected statistic names are refreshed.
#'
#' @return Invisibly, `df` with the selected statistics and an
#'   `aster_ged_summary` entry in its `traktr_layers` attribute containing
#'   the tile index and summary parameters.
#' @examples
#' \dontrun{
#' df <- aoi_from_sf(my_polygons, id_col = "id")
#' df <- fetch_aster_ged(df, out_dir = "data/ag100")
#' df <- attach_aster_ged(df, naming = "traktr", convert = TRUE)
#' }
#' @export
attach_aster_ged <- function(df,
                             bands = c("elev", "ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd", "lwm", "nobs"),
                             naming = c("county", "ecolozip", "traktr"),
                             extras = c("median", "n_cells", "area"),
                             convert = FALSE,
                             ecolozip_na_quirk = FALSE,
                             allow_missing = FALSE,
                             workers = 1,
                             progress = TRUE,
                             cache = TRUE) {
  bands <- unique(rlang::arg_match(bands, multiple = TRUE))
  naming <- rlang::arg_match(naming)
  extras <- unique(rlang::arg_match(extras, multiple = TRUE))
  if (!length(bands)) cli::cli_abort("{.arg bands} must select at least one ASTER GED band.")

  tiles <- attr(df, "aster_ged_tiles")
  map <- attr(df, "aster_ged_map")
  aoi <- attr(df, "aster_ged_aoi")
  if (is.null(aoi)) {
    aoi <- attr(df, "traktr_aoi")
    id_col <- attr(df, "traktr_aoi_params")$id_col
  } else {
    id_col <- setdiff(names(aoi), attr(aoi, "sf_column"))
  }
  if (is.null(tiles) || is.null(map) || !inherits(aoi, "sf") ||
      length(id_col) != 1L || !id_col %in% intersect(names(df), names(aoi))) {
    cli::cli_abort("Missing ASTER GED or AOI metadata: run fetch_aster_ged() first.")
  }

  aoi <- aoi[aoi[[id_col]] %in% df[[id_col]], , drop = FALSE]
  map <- map[map$aoi_id %in% aoi[[id_col]], , drop = FALSE]
  tiles <- tiles[tiles$tile_id %in% map$tile_id & !(tiles$status %in% "missing_remote"),
                 , drop = FALSE]
  missing <- tiles$status %in% c("error", "missing_local") |
    !file.exists(tiles$file) | dir.exists(tiles$file)
  if (any(missing)) {
    missing_ids <- tiles$tile_id[missing]
    msg <- "Unavailable ASTER GED tiles: {.val {missing_ids}}."
    if (!isTRUE(allow_missing)) cli::cli_abort(msg)
    cli::cli_warn(msg)
  }
  tiles <- tiles[!missing, , drop = FALSE]

  cache_dir <- NULL
  if (isTRUE(cache) && nrow(tiles)) {
    tile_dir <- attr(df, "traktr_layers")$aster_ged$params$out_dir
    if (is.null(tile_dir)) tile_dir <- dirname(tiles$file[1L])
    cache_dir <- file.path(tile_dir, "_moments")
  }

  dict <- .aster_ged_dictionary(convert)
  core <- dict[!dict$extra & dict$band %in% bands, , drop = FALSE]
  medians <- dict$stat == "median" & dict$band %in% bands & "median" %in% extras
  other_extras <- is.na(dict$band) & dict$stat %in% extras
  columns <- rbind(core, dict[medians | other_extras, , drop = FALSE])
  output_names <- columns[[paste0("name_", naming)]]
  if (id_col %in% output_names) {
    cli::cli_abort("AOI ID column {.val {id_col}} conflicts with an ASTER GED summary name.")
  }

  layers <- unique(core$band)
  extraction_tiles <- tiles
  if (all(c("lat_n", "lon_w") %in% names(tiles))) {
    extraction_tiles$xmin <- tiles$lon_w
    extraction_tiles$xmax <- tiles$lon_w + 1
    extraction_tiles$ymin <- tiles$lat_n - 1
    extraction_tiles$ymax <- tiles$lat_n
  }
  moments <- .extract_tiled(aoi, id_col, extraction_tiles, map, layers = layers, workers = workers,
                            cache_dir = cache_dir, progress = progress)
  ids <- aoi[[id_col]]
  stats <- data.frame(ids = ids[!duplicated(ids)], check.names = FALSE)
  names(stats) <- id_col
  summaries <- lapply(layers, function(band) {
    band_stats <- columns$stat[which(columns$band == band)]
    .aster_ged_band_summary(lapply(moments, `[[`, band), band, band_stats,
                            convert, ecolozip_na_quirk)
  })
  names(summaries) <- layers
  count_band <- if ("elev" %in% bands) "elev" else bands[1L]
  for (i in seq_len(nrow(columns))) {
    band <- columns$band[i]
    stat <- columns$stat[i]
    stats[[output_names[i]]] <- if (!is.na(band)) {
      summaries[[band]][[stat]]
    } else if (stat == "n_cells") {
      vapply(moments, function(m) m[[count_band]]$n, numeric(1))
    } else {
      # Sum over all rows of an ID, matching how the moments combine them
      area <- rowsum(.aster_ged_area(aoi), as.character(ids), reorder = FALSE)
      unname(area[as.character(ids[!duplicated(ids)]), 1L])
    }
  }

  df[intersect(output_names, names(df))] <- NULL
  out <- .join_stats(df, stats, by = id_col, prefix = "")
  out <- .set_traktr_layer(out, "aster_ged_summary", tiles, params = list(
    bands = bands, naming = naming, extras = extras, convert = convert,
    ecolozip_na_quirk = ecolozip_na_quirk
  ))
  invisible(out)
}

# Confirmed on tile 40.-105 (smoke run 2, 2026-11-26).
.ag100_water_value <- 1L

.aster_ged_band_summary <- function(moments, band, stats, convert, ecolozip_na_quirk) {
  base <- c("mean", "min", "max", "sd", "median", "missing_perc")
  summary <- lapply(moments, .finalize_moments, stats = base)
  if (!length(summary)) {
    summary <- list(.finalize_moments(.moments_fn(numeric(), numeric()), base)[FALSE, ])
  }
  out <- do.call(rbind, unname(summary))

  if (band == "elev") {
    out$cin <- vapply(moments, .hist_share, numeric(1), pred = function(value) value == 0)
    out$sum <- rep(NA_real_, nrow(out))
  } else if (band == "nobs") {
    out$lt_4 <- as.integer(out$mean < 4)
  } else if (band == "lwm") {
    out$perc <- vapply(moments, .hist_share, numeric(1),
                       pred = function(value) value == .ag100_water_value)
    out$cin <- 100 - out$missing_perc
  }

  if (isTRUE(ecolozip_na_quirk)) out$missing_perc[which(out$missing_perc == 0)] <- NA_real_
  if (isTRUE(convert) && band %in% c("ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd")) {
    out[base[1:5]] <- out[base[1:5]] / 100
    if (band == "temp_mean") {
      location <- c("mean", "min", "max", "median")
      out[location] <- out[location] - 273.15
    }
  }
  out[, stats, drop = FALSE]
}

.aster_ged_area <- function(aoi) {
  previous_s2 <- getOption("sf_use_s2")
  on.exit(options(sf_use_s2 = previous_s2), add = TRUE)
  suppressMessages(sf::sf_use_s2(TRUE))
  as.numeric(sf::st_area(sf::st_transform(aoi, 4326))) / 1e6
}
