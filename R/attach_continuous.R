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
#' @param deduplicate Extract identical geometries only once.
#' @param crs_planar Projected CRS the AOI is transformed to before extraction.
#' @param ids Optional subset of AOI ids to extract (all if \code{NULL}).
#'
#' @return A data frame with \code{id_col} and the unprefixed statistic
#'   columns (\code{mean}, \code{sd}, \code{min}, \code{max}, \code{median}),
#'   one row per AOI id.
#'
#' @keywords internal
#' @noRd
.extract_continuous <- function(aoi, id_col, files,
                                band = 1,
                                count_band = NULL,
                                min_count = NULL,
                                deduplicate = TRUE,
                                crs_planar = 3857,
                                ids = NULL) {
  r <- if (length(files) == 1L) terra::rast(files) else terra::vrt(terra::sprc(files))
  if (band > terra::nlyr(r)) {
    cli::cli_abort("Raster has {terra::nlyr(r)} band{?s}; cannot select band {band}.")
  }
  if (!is.null(count_band)) {
    if (count_band > terra::nlyr(r)) {
      cli::cli_abort("Raster has no observation-count band {count_band}; cannot apply a minimum count.")
    }
    enough <- r[[count_band]] >= min_count
    r <- terra::mask(r[[band]], enough, maskvalues = c(0, NA))
  } else {
    r <- r[[band]]
  }

  if (!is.null(ids)) aoi <- aoi[aoi[[id_col]] %in% ids, , drop = FALSE]

  aoi_planar <- sf::st_transform(aoi, crs_planar)
  dd <- .aoi_deduplicate(aoi_planar, id_col = id_col, deduplicate = deduplicate)

  # match the raster CRS (e.g. Mollweide for HFP) to avoid on-the-fly transforms
  au <- sf::st_transform(dd$aoi_unique, terra::crs(r))

  res <- exactextractr::exact_extract(r, au, fun = .continuous_summary_fn,
                                      progress = FALSE)
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

  df$.__row_order <- seq_len(nrow(df))
  out <- merge(df, stats, by = by, all.x = TRUE, sort = FALSE)
  out <- out[order(out$.__row_order), , drop = FALSE]
  out$.__row_order <- NULL

  if ("geometry" %in% names(out)) {
    geom_col <- "geometry"
    out <- out[, c(setdiff(names(out), geom_col), geom_col), drop = FALSE]
  }

  .restore_attrs_class(df, out)
}
