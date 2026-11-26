#' @noRd
.spiderweb_ecolozip <- function(radii_km = c(8.0467, 16.0934, 24.1401, 32.1869),
                                n_points = 36, drop_duplicate = TRUE) {
  list(radii_km = radii_km, n_points = n_points, drop_duplicate = drop_duplicate)
}

#' @noRd
.spiderweb_points <- function(lat, lon, spec, lat_ref = NULL) {
  if (length(lat) != length(lon)) cli::cli_abort("Supply one longitude per latitude.")
  if (is.null(lat_ref)) lat_ref <- lat
  if (length(lat_ref) != 1L && length(lat_ref) != length(lat)) {
    cli::cli_abort("Supply one reference latitude or one per centroid.")
  }
  lat_ref <- rep_len(lat_ref, length(lat))

  angle <- seq(0, 2 * pi, length.out = spec$n_points)
  if (isTRUE(spec$drop_duplicate) && length(angle) > 1L) angle <- angle[-length(angle)]
  offsets <- expand.grid(k = seq_along(angle), ring = seq_along(spec$radii_km),
                         KEEP.OUT.ATTRS = FALSE)
  id <- rep(seq_along(lat), each = nrow(offsets))
  k <- rep(offsets$k, times = length(lat))
  ring <- rep(offsets$ring, times = length(lat))
  radius_deg <- spec$radii_km[ring] / 111

  sample_lon <- lon[id] + radius_deg / cos(lat_ref[id] / 57.3) * cos(angle[k])
  wrap <- which(sample_lon < -180 | sample_lon >= 180)
  sample_lon[wrap] <- (sample_lon[wrap] + 180) %% 360 - 180
  data.frame(id = id, ring = ring, k = k,
             lat = lat[id] + radius_deg * sin(angle[k]), lon = sample_lon)
}

#' @noRd
.spiderweb_tile_id <- function(lat, lon) {
  if (length(lat) != length(lon)) cli::cli_abort("Supply one longitude per latitude.")
  id <- rep(NA_character_, length(lat))
  located <- which(lat > -90 & lat <= 90 & lon >= -180 & lon < 180)
  # Match asterlookup: the southern tile owns the latitude seam and the
  # eastern tile owns the longitude seam.
  id[located] <- sprintf("%d.%d", ceiling(lat[located]), floor(lon[located]))
  id
}

#' @noRd
.aster_ged_sample_points <- function(points, tile_dir, layers = names(.ag100_band_paths)) {
  sampled <- matrix(NA_real_, nrow(points), length(layers))
  tile_missing <- rep(TRUE, nrow(points))

  tiles <- .spiderweb_tile_id(points$lat, points$lon)
  groups <- split(seq_along(tiles), tiles)
  for (tile in names(groups)) {
    path <- file.path(tile_dir, paste0("ag100_v003_", sub(".", "_", tile, fixed = TRUE), ".tif"))
    if (!file.exists(path) || dir.exists(path)) next
    rows <- groups[[tile]]
    raster <- terra::rast(path)[[layers]]
    xy <- as.matrix(points[rows, c("lon", "lat")])
    sampled[rows, ] <- as.matrix(terra::extract(raster, xy, method = "simple"))
    tile_missing[rows] <- FALSE
  }
  points[layers] <- as.data.frame(sampled)
  points$tile_missing <- tile_missing
  points
}

#' @noRd
.spss_aggregate <- function(values, id) {
  if (nrow(values) != length(id)) cli::cli_abort("Supply one ID per sampled point.")
  result <- data.frame(id = unique(id))
  groups <- split(seq_along(id), factor(match(id, result$id), levels = seq_len(nrow(result))))
  dict <- .aster_ged_dictionary()
  dict <- dict[!dict$extra & dict$band %in% names(values), , drop = FALSE]

  for (band in unique(dict$band)) {
    columns <- dict[dict$band == band, , drop = FALSE]
    x <- values[[band]]
    summaries <- vapply(groups, function(rows) {
      samples <- x[rows]
      valid <- samples[!is.na(samples)]
      n <- length(valid)
      total <- length(rows)
      missing <- total - n
      missing_perc <- if (missing) 100 * missing / total else NA_real_
      if (band == "lwm") {
        return(c(perc = if (n) 100 * sum(valid == .ag100_water_value) / total else NA_real_,
                 cin = if (n) 100 * n / total else NA_real_, missing_perc = missing_perc))
      }
      mean_value <- if (n) mean(valid) else NA_real_
      stats <- c(mean = mean_value, min = if (n) min(valid) else NA_real_,
                 max = if (n) max(valid) else NA_real_, sd = stats::sd(valid),
                 missing_perc = missing_perc)
      if (band == "elev") {
        stats <- c(stats, cin = if (n) 100 * sum(valid == 0) / total else NA_real_,
                   sum = if (n) sum(valid) else NA_real_)
      }
      if (band == "nobs") stats <- c(stats, lt_4 = mean_value < 4)
      stats[columns$stat]
    }, numeric(nrow(columns)))
    for (i in seq_len(nrow(columns))) {
      summary <- unname(summaries[i, ])
      if (columns$stat[i] == "lt_4") summary <- as.integer(summary)
      result[[columns$name_ecolozip[i]]] <- summary
    }
  }
  result
}
