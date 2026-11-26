#' Area-weighted summary of a continuous raster (internal)
#'
#' Summary callback for \code{exactextractr::exact_extract()}. Computes the
#' coverage-weighted mean, sd, min, max and median of the cell values that
#' fall into one polygon.
#'
#' @param values Cell values (vector, or data frame for multi-layer rasters,
#'   in which case layers are averaged per cell).
#' @param coverage_fraction Fraction of each cell covered by the polygon.
#'
#' @return A one-row data frame with columns \code{mean}, \code{sd},
#'   \code{min}, \code{max}, \code{median}.
#'
#' @keywords internal
#' @noRd
.continuous_summary_fn <- function(values, coverage_fraction) {
  if (is.data.frame(values)) {
    v <- rowMeans(as.matrix(values), na.rm = TRUE)
  } else {
    v <- values
  }
  w <- coverage_fraction
  ok <- is.finite(v) & is.finite(w) & (w > 0)
  if (!any(ok)) {
    return(data.frame(
      mean   = NA_real_,
      sd     = NA_real_,
      min    = NA_real_,
      max    = NA_real_,
      median = NA_real_
    ))
  }
  v <- v[ok]; w <- w[ok]
  w <- w / sum(w)  # normalize weights

  m     <- sum(v * w)
  var_w <- sum(w * (v - m)^2)

  ord   <- order(v)
  cw    <- cumsum(w[ord])
  med   <- v[ord][which(cw >= 0.5)[1]]

  data.frame(
    mean   = m,
    sd     = sqrt(var_w),
    min    = min(v),
    max    = max(v),
    median = med
  )
}


#' Extract continuous raster statistics per AOI (internal)
#'
#' Builds a raster (VRT if several files) from \code{files}, selects
#' \code{band}, and computes area-weighted statistics for every AOI id.
#'
#' @param aoi AOI \code{sf} (id column + geometry).
#' @param id_col Name of the AOI id column.
#' @param files Character vector of GeoTIFF paths sharing one CRS.
#' @param band Band (layer) index to summarise.
#' @param count_band,min_count Optional band holding a per-pixel observation
#'   count; pixels with a count below \code{min_count} are masked out.
#' @param extra_bands Optional named list of further bands to summarise over
#'   the pixels that have a valid value in \code{band}, each as
#'   \code{list(band = <index>, fun = <exactextractr operation>)}, e.g.
#'   \code{list(n_obs = list(band = 2, fun = "mean"))}. Bands the raster does
#'   not have are skipped.
#' @param deduplicate Extract identical geometries only once.
#' @param crs_planar Projected CRS the AOI is transformed to before extraction,
#'   or `NULL` to transform directly to the raster CRS.
#' @param ids Optional subset of AOI ids to extract (all if \code{NULL}).
#' @param weights Use cell coverage fractions or covered cell areas.
#'
#' @return A data frame with \code{id_col} and the unprefixed statistic
#'   columns (\code{mean}, \code{sd}, \code{min}, \code{max}, \code{median},
#'   plus one column per entry of \code{extra_bands}), one row per AOI id.
#'
#' @keywords internal
#' @noRd
.extract_continuous <- function(aoi, id_col, files,
                                band = 1,
                                count_band = NULL,
                                min_count = NULL,
                                extra_bands = NULL,
                                deduplicate = TRUE,
                                crs_planar = 3857,
                                ids = NULL,
                                weights = c("fraction", "area")) {
  weights <- match.arg(weights)
  r_all <- if (length(files) == 1L) terra::rast(files) else terra::vrt(terra::sprc(files))
  nl <- terra::nlyr(r_all)
  if (band > nl) {
    cli::cli_abort("Raster has {nl} band{?s}; cannot select band {band}.")
  }
  if (!is.null(count_band) && count_band > nl) {
    cli::cli_abort("Raster has no observation-count band {count_band}; cannot apply a minimum count.")
  }
  extra_bands <- Filter(function(e) e$band <= nl, extra_bands)

  # one extraction pass over the needed bands; masking happens per polygon
  # inside the summary function, never on the whole (possibly huge) mosaic
  layers <- unique(c(band, count_band, vapply(extra_bands, `[[`, numeric(1), "band")))
  r <- r_all[[layers]]
  col_of <- function(b) match(b, layers)

  summary_fn <- function(values, coverage_fraction) {
    values <- as.data.frame(values)
    v <- values[[col_of(band)]]
    if (!is.null(count_band)) v[!(values[[col_of(count_band)]] >= min_count)] <- NA
    out <- .continuous_summary_fn(v, coverage_fraction)
    ok <- is.finite(v) & is.finite(coverage_fraction) & coverage_fraction > 0
    for (nm in names(extra_bands)) {
      e <- extra_bands[[nm]]
      x <- values[[col_of(e$band)]][ok]
      w <- coverage_fraction[ok]
      out[[nm]] <- if (!length(x) || all(is.na(x))) NA_real_ else switch(e$fun,
        mean = sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)]),
        max  = max(x, na.rm = TRUE),
        cli::cli_abort("Unsupported summary {.val {e$fun}}.")
      )
    }
    out
  }

  if (!is.null(ids)) aoi <- aoi[aoi[[id_col]] %in% ids, , drop = FALSE]

  aoi_planar <- if (is.null(crs_planar)) aoi else sf::st_transform(aoi, crs_planar)
  dd <- .aoi_deduplicate(aoi_planar, id_col = id_col, deduplicate = deduplicate,
                         to_wgs84 = !is.null(crs_planar))

  # match the raster CRS (e.g. Mollweide for HFP) to avoid on-the-fly transforms
  au <- sf::st_transform(dd$aoi_unique, terra::crs(r))

  res <- exactextractr::exact_extract(
    r, au, fun = summary_fn, coverage_area = weights == "area", progress = FALSE
  )
  if (is.list(res) && !is.data.frame(res)) {
    res <- do.call(rbind, res)
    rownames(res) <- NULL
  }
  if (!nrow(res)) {
    cli::cli_abort("Raster extraction returned no results.")
  }
  res$.aoi_row <- seq_len(nrow(res))

  per_aoi <- merge(dd$map_back, res, by = ".aoi_row", all.x = TRUE, sort = FALSE)
  per_aoi$.aoi_row <- NULL
  per_aoi[!duplicated(per_aoi[[id_col]]), , drop = FALSE]
}

#' @noRd
.weighted_hist <- function(values, weights) {
  distinct <- sort(unique(values))
  sums <- if (length(distinct)) {
    as.numeric(rowsum(weights, match(values, distinct)))
  } else {
    numeric()
  }
  data.frame(value = distinct, w = sums)
}

#' @noRd
.moments_fn <- function(values, weights) {
  if (is.data.frame(values)) {
    return(lapply(values, .moments_fn, weights = weights))
  }

  covered <- is.finite(weights) & weights > 0
  valid <- covered & is.finite(values)
  v <- as.double(values[valid])
  wt <- weights[valid]
  w <- sum(wt)
  mean <- if (w > 0) sum(v * (wt / w)) else NA_real_

  list(
    w = w, mean = mean, m2 = sum(wt * (v - mean)^2),
    min = if (w > 0) min(v) else NA_real_,
    max = if (w > 0) max(v) else NA_real_,
    w_na = sum(weights[covered & is.na(values)]), n = length(v),
    hist = .weighted_hist(v, wt)
  )
}

#' @noRd
.combine_moments <- function(a, b) {
  if (a$w == 0 || b$w == 0) {
    out <- if (a$w == 0) b else a
    out$w_na <- a$w_na + b$w_na
    return(out)
  }

  w <- a$w + b$w
  delta <- b$mean - a$mean
  a_edges <- attr(a, "edge_cells")
  b_edges <- attr(b, "edge_cells")
  out <- list(
    w = w, mean = a$mean + delta * (b$w / w),
    m2 = a$m2 + b$m2 + delta^2 * a$w * (b$w / w),
    min = min(a$min, b$min), max = max(a$max, b$max),
    w_na = a$w_na + b$w_na, n = sum(a$n, b$n) - length(intersect(a_edges, b_edges)),
    hist = .weighted_hist(c(a$hist$value, b$hist$value), c(a$hist$w, b$hist$w))
  )
  if (length(a_edges) || length(b_edges)) attr(out, "edge_cells") <- union(a_edges, b_edges)
  out
}

#' @noRd
.finalize_moments <- function(m, stats) {
  total <- m$w + m$w_na
  median <- NA_real_
  if (m$w > 0) {
    # Resolve near-half cumulative weights to the lower value. This can differ
    # from the continuous summary within 1e-9 relative to the half-weight.
    # Tile order cannot change the median unless the cumulative weight lies
    # within a few ULP of this cutoff itself.
    cutoff <- (m$w / 2) * (1 - 1e-9)
    median <- m$hist$value[which(cumsum(m$hist$w) >= cutoff)[1L]]
  }
  data.frame(
    mean = m$mean,
    sd = if (m$w > 0) sqrt(m$m2 / m$w) else NA_real_,
    min = m$min,
    max = m$max,
    median = median,
    missing_perc = if (total > 0) 100 * (m$w_na / total) else NA_real_,
    n_cells = m$n,
    area_valid_km2 = m$w / 1e6
  )[, stats, drop = FALSE]
}

#' @noRd
.hist_share <- function(m, pred) {
  if (m$w == 0) return(NA_real_)
  100 * (sum(m$hist$w[pred(m$hist$value)]) / m$w)
}

#' @noRd
.tile_moments <- function(r, polys, box = NULL) {
  layers <- names(r)
  # Keep arbitrary product layer names clear of exact_extract's column names.
  raster_names <- paste0("layer_", seq_along(layers))
  names(r) <- raster_names
  value_columns <- if (length(layers) == 1L) "value" else raster_names
  track_edges <- !is.null(box) && any(abs(as.vector(terra::ext(r)) - box) > 1e-9)
  res <- terra::res(r)
  summarize <- function(cells) {
    weights <- cells$coverage_area
    values <- cells[value_columns]
    moments <- .moments_fn(values, weights)
    if (track_edges) {
      # Only cells centred on a box boundary can occur in adjacent AG100 tiles.
      near <- function(x, edge) abs(x - edge) <= 1e-9
      edge <- which(is.finite(weights) & weights > 0 &
                      (near(cells$x, box[1L]) | near(cells$x, box[2L]) |
                         near(cells$y, box[3L]) | near(cells$y, box[4L])))
      x <- cells$x[edge]
      x[near(x, 180)] <- -180
      keys <- sprintf("%.10f,%.10f,%.10f,%.10f",
                      round(x, 10L) + 0, round(cells$y[edge], 10L) + 0, res[1L], res[2L])
      for (i in seq_along(moments)) {
        valid <- is.finite(values[[i]][edge])
        if (any(valid)) attr(moments[[i]], "edge_cells") <- keys[valid]
      }
    }
    # A length-one list preserves the nested result through exact_extract().
    list(stats::setNames(moments, layers))
  }
  exactextractr::exact_extract(
    r, polys, fun = summarize, coverage_area = TRUE, summarize_df = TRUE,
    include_xy = track_edges, progress = FALSE
  )
}

.moments_cache_version <- 2L

#' @noRd
.is_cached_moment <- function(m) {
  fields <- c("w", "mean", "m2", "min", "max", "w_na", "n")
  if (!is.list(m) || !identical(names(m), c(fields, "hist"))) return(FALSE)
  if (!all(vapply(m[fields], function(x) {
    is.vector(x, mode = "numeric") && length(x) == 1L
  }, logical(1)))) return(FALSE)
  if (!all(is.finite(c(m$w, m$m2, m$w_na, m$n))) ||
      min(m$w, m$m2, m$w_na, m$n) < 0 || m$n != floor(m$n)) return(FALSE)
  location <- c(m$mean, m$min, m$max)
  if (m$w > 0 && !all(is.finite(location))) return(FALSE)
  if (m$w == 0 && !all(is.na(location))) return(FALSE)

  edges <- attr(m, "edge_cells")
  if (!is.null(edges) && (!is.character(edges) || anyNA(edges) ||
                          anyDuplicated(edges) || length(edges) > m$n)) return(FALSE)
  hist <- m$hist
  is.data.frame(hist) && identical(names(hist), c("value", "w")) &&
    is.vector(hist$value, mode = "double") && is.vector(hist$w, mode = "numeric") &&
    (nrow(hist) > 0L) == (m$w > 0) &&
    all(is.finite(hist$value)) && !is.unsorted(hist$value, strictly = TRUE) &&
    all(is.finite(hist$w) & hist$w > 0)
}

#' @noRd
.is_tile_cache <- function(cached, ids, layers) {
  if (!is.list(cached) || !identical(cached[["version"]], .moments_cache_version)) return(FALSE)
  moments <- cached[["moments"]]
  is.list(moments) && identical(names(moments), ids) &&
    all(vapply(moments, function(m) {
      is.list(m) && identical(names(m), layers) &&
        all(vapply(m, .is_cached_moment, logical(1)))
    }, logical(1)))
}

#' @noRd
.read_or_extract_tile <- function(job, layers, cache_dir) {
  path <- job$file
  tile_id <- job$tile_id
  polys <- job$polys
  poly_ids <- job$ids
  if (!is.null(cache_dir)) {
    info <- file.info(path)
    fingerprint <- unname(tools::md5sum(path))
    if (is.na(fingerprint)) {
      cli::cli_abort("Could not fingerprint the raster for tile {.val {tile_id}}.")
    }
    key <- rlang::hash(list(.moments_cache_version, fingerprint,
                            info$size, info$mtime, layers, poly_ids, job$box, job$rows,
                            sf::st_as_binary(sf::st_geometry(polys))))
    cache_file <- file.path(cache_dir, paste0(tile_id, "_", key, ".rds"))
    if (file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL,
                         warning = function(w) NULL)
      if (.is_tile_cache(cached, poly_ids, layers)) {
        return(cached$moments)
      }
      cli::cli_warn("Moments cache for tile {.val {tile_id}} is unreadable or has an invalid schema. Recomputing.")
    }
  }

  r <- terra::rast(path, guessCRS = FALSE)[[layers]]
  out <- stats::setNames(.tile_moments(r, polys, job$box), poly_ids)
  if (!is.null(job$box)) {
    for (j in seq_along(out)) {
      for (i in seq_along(out[[j]])) {
        edges <- attr(out[[j]][[i]], "edge_cells")
        if (length(edges)) {
          attr(out[[j]][[i]], "edge_cells") <- paste(job$rows[j], edges, sep = ":")
        }
      }
    }
  }
  if (!is.null(cache_dir) && !identical(unname(tools::md5sum(path)), fingerprint)) {
    # The raster was replaced during extraction; do not cache under the old key.
    cli::cli_warn("Raster for tile {.val {tile_id}} changed during extraction. Result not cached.")
  } else if (!is.null(cache_dir)) {
    tmp <- tempfile(tmpdir = cache_dir)
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(list(version = .moments_cache_version, moments = out), tmp)
    if (!suppressWarnings(file.rename(tmp, cache_file))) {
      cli::cli_abort("Could not publish the moments cache for tile {.val {tile_id}}.")
    }
  }
  out
}

#' @noRd
.validate_tile_extents <- function(extents, tile_ids) {
  if (anyDuplicated(tile_ids)) {
    cli::cli_abort("Tile ID {.val {tile_ids[anyDuplicated(tile_ids)]}} occurs more than once.")
  }

  rectangles <- sf::st_sfc(lapply(extents, function(e) {
    x <- e[1:2]
    y <- e[3:4]
    sf::st_polygon(list(cbind(x[c(1, 2, 2, 1, 1)], y[c(1, 1, 2, 2, 1)])))
  }))
  hits <- sf::st_intersects(rectangles)
  for (i in seq_along(hits)) {
    a <- extents[[i]]
    for (j in hits[[i]][hits[[i]] > i]) {
      b <- extents[[j]]
      # Shared boundaries are allowed. Any intersection with area is an overlap.
      if (max(a[1L], b[1L]) < min(a[2L], b[2L]) &&
          max(a[3L], b[3L]) < min(a[4L], b[4L])) {
        cli::cli_abort("Tile {.val {tile_ids[i]}} overlaps tile {.val {tile_ids[j]}}.")
      }
    }
  }
}

#' @noRd
.clip_tile_polygons <- function(polys, box) {
  previous_s2 <- getOption("sf_use_s2")
  on.exit(options(sf_use_s2 = previous_s2), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  # Most AOIs can be clipped in one pass without splitting multipart features.
  if (!any(vapply(polys, .ag100_crosses_antimeridian, logical(1)))) {
    return(.intersect_tile_box(polys, box))
  }

  # Classify each polygon independently so an already split multipart AOI
  # can also contain a ring that still needs its longitudes shifted.
  parts <- sf::st_cast(sf::st_sf(row = seq_along(polys), geometry = polys),
                       "POLYGON", warn = FALSE)
  rings <- sf::st_geometry(parts)
  shifted <- vapply(rings, .ag100_crosses_antimeridian, logical(1))
  pieces <- polys[0L]
  rows <- integer()
  for (shift in c(FALSE, TRUE)) {
    i <- which(shifted == shift)
    if (!length(i)) next
    tile_rings <- rings[i]
    if (shift) tile_rings <- sf::st_shift_longitude(tile_rings)
    offset <- if (shift && box[1L] < 0) 360 else 0
    clipped <- .intersect_tile_box(tile_rings, box + c(offset, offset, 0, 0))
    if (!length(clipped$index)) next
    if (offset != 0) {
      clipped$polys <- sf::st_set_crs(clipped$polys - c(offset, 0), 4326)
    }
    pieces <- c(pieces, clipped$polys)
    rows <- c(rows, parts$row[i[clipped$index]])
  }
  # Reassemble each original AOI row before counting its raster cells.
  by_row <- split(seq_along(rows), rows)
  merged <- lapply(by_row, function(i) sf::st_combine(pieces[i])[[1L]])
  list(polys = sf::st_sfc(merged, crs = 4326), index = as.integer(names(by_row)))
}

#' @noRd
.intersect_tile_box <- function(polys, box) {
  bounds <- stats::setNames(box, c("xmin", "xmax", "ymin", "ymax"))
  rectangle <- sf::st_as_sfc(sf::st_bbox(bounds, crs = 4326))
  clipped <- suppressMessages(sf::st_intersection(polys, rectangle))
  keep <- which(!sf::st_is_empty(clipped) & sf::st_dimension(clipped) == 2L)
  index <- attr(clipped, "idx")[keep, 1L]
  clipped <- clipped[keep]
  if (any(sf::st_geometry_type(clipped) == "GEOMETRYCOLLECTION")) {
    clipped <- suppressWarnings(sf::st_collection_extract(clipped, "POLYGON"))
  }
  list(polys = clipped, index = index)
}

#' @noRd
.extract_tiled <- function(aoi, id_col, tiles, map, layers = NULL, workers = 1,
                           cache_dir = NULL, progress = TRUE) {
  tiles <- tiles[!is.na(tiles$file), , drop = FALSE]
  ids <- as.character(aoi[[id_col]])
  unique_ids <- unique(ids)
  extents <- vector("list", nrow(tiles))
  box_cols <- c("xmin", "xmax", "ymin", "ymax")
  has_box <- all(box_cols %in% names(tiles))
  if (!has_box && any(box_cols %in% names(tiles))) {
    cli::cli_abort("Nominal tile boxes require xmin, xmax, ymin and ymax.")
  }
  for (i in seq_len(nrow(tiles))) {
    tile_id <- tiles$tile_id[i]
    r <- terra::rast(tiles$file[i], guessCRS = FALSE)
    crs <- terra::crs(r)
    if (!nzchar(crs) || sf::st_crs(crs) != sf::st_crs(4326)) {
      cli::cli_abort("Tile {.val {tile_id}} must use EPSG:4326.")
    }
    extent <- as.vector(terra::ext(r))
    box <- if (has_box) unlist(tiles[i, box_cols], use.names = FALSE) else extent
    if (!is.numeric(box) || any(!is.finite(box))) {
      cli::cli_abort("Tile {.val {tile_id}} must have finite numeric bounds.")
    }
    edges <- round(box)
    if (any(abs(box - edges) > 1e-9) || any(edges[c(1L, 3L)] >= edges[c(2L, 4L)])) {
      cli::cli_abort("Tile {.val {tile_id}} must have whole-degree bounds.")
    }
    if (has_box && any((extent - box) * c(1, -1, 1, -1) > 1e-9)) {
      cli::cli_abort("Raster for tile {.val {tile_id}} does not cover its nominal box.")
    }
    extents[[i]] <- box
    if (is.null(layers)) layers <- names(r)
    unknown <- setdiff(layers, names(r))
    if (length(unknown)) {
      cli::cli_abort("Unknown layer selection {.val {unknown}} in tile {.val {tile_id}}.")
    }
  }

  empty <- .moments_fn(numeric(), numeric())
  empty_layers <- stats::setNames(rep(list(empty), length(layers)), layers)
  out <- stats::setNames(rep(list(empty_layers), length(unique_ids)), unique_ids)
  if (!nrow(tiles)) return(out)
  .validate_tile_extents(extents, tiles$tile_id)
  if (!length(layers)) return(out)

  aoi_rows <- split(seq_along(ids), match(ids, unique_ids))
  tile_rows <- lapply(split(match(map$aoi_id, unique_ids), map$tile_id), function(mapped) {
    sort(unlist(aoi_rows[unique(stats::na.omit(mapped))], use.names = FALSE))
  })
  tile_rows <- unname(tile_rows[as.character(tiles$tile_id)])
  # Pass only a tile's AOI geometries to its worker.
  geometries <- sf::st_geometry(aoi)
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  workers <- as.integer(workers)
  parallel <- workers > 1L && .has_future_apply()
  if (workers > 1L && !parallel) {
    cli::cli_warn("Package {.pkg future.apply} is unavailable. Processing tiles sequentially.")
  }
  if (isTRUE(progress)) {
    progress_id <- cli::cli_progress_bar("Tile extraction", total = nrow(tiles), auto_terminate = FALSE)
    on.exit(cli::cli_progress_done(id = progress_id), add = TRUE)
  }

  # Keep only a batch of tile summaries in memory while reducing into the AOIs.
  batch_size <- if (parallel) workers else 1L
  for (start in seq.int(1L, nrow(tiles), by = batch_size)) {
    batch <- seq.int(start, min(start + batch_size - 1L, nrow(tiles)))
    jobs <- lapply(batch[lengths(tile_rows[batch]) > 0L], function(i) {
      rows <- tile_rows[[i]]
      polys <- geometries[rows]
      box <- if (has_box) extents[[i]]
      if (has_box) {
        clipped <- .clip_tile_polygons(polys, box)
        polys <- clipped$polys
        rows <- rows[clipped$index]
      }
      list(file = tiles$file[i], tile_id = tiles$tile_id[i],
           polys = polys, ids = ids[rows], box = box, rows = if (has_box) rows)
    })
    jobs <- Filter(function(job) length(job$ids) > 0L, jobs)
    tile_results <- if (parallel && length(jobs)) {
      future.apply::future_lapply(
        jobs, .read_or_extract_tile, layers = layers, cache_dir = cache_dir,
        future.seed = TRUE, future.chunk.size = 1L
      )
    } else {
      lapply(jobs, .read_or_extract_tile, layers = layers, cache_dir = cache_dir)
    }
    for (tile in tile_results) {
      for (j in seq_along(tile)) {
        id <- names(tile)[j]
        out[[id]] <- Map(.combine_moments, out[[id]], tile[[j]])
      }
    }
    if (isTRUE(progress)) cli::cli_progress_update(inc = length(batch), id = progress_id)
  }
  out
}


#' Join per-AOI statistics back to the input data (internal)
#'
#' Prefixes the statistic columns, merges them into \code{df} by \code{by},
#' restores row order, puts geometry last and restores attributes/class.
#'
#' @param df The input data (carrying the AOI attributes).
#' @param stats Data frame with the \code{by} columns plus statistics.
#' @param by Join columns.
#' @param prefix Prefix for the statistic columns.
#'
#' @return \code{df} augmented with the prefixed statistic columns.
#'
#' @keywords internal
#' @noRd
.join_stats <- function(df, stats, by, prefix) {
  stat_cols <- setdiff(names(stats), by)
  names(stats)[match(stat_cols, names(stats))] <- paste0(prefix, stat_cols)

  existing_names <- c(names(df), names(stats))
  row_order <- ".__row_order"
  while (row_order %in% existing_names) row_order <- paste0(row_order, "_")

  df[[row_order]] <- seq_len(nrow(df))
  out <- merge(df, stats, by = by, all.x = TRUE, sort = FALSE)
  out <- out[order(out[[row_order]]), , drop = FALSE]
  out[[row_order]] <- NULL

  if ("geometry" %in% names(out)) {
    geom_col <- "geometry"
    out <- out[, c(setdiff(names(out), geom_col), geom_col), drop = FALSE]
  }

  .restore_attrs_class(df, out)
}
