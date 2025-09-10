#' Attach NDVI statistics to an AOI
#'
#' Extracts NDVI values from GeoTIFF chips covering an AOI, computes
#' area-weighted summary statistics (mean, sd, min, max, median), and joins
#' them back to the input data. All statistics are weighted by covered surface
#' area in square meters, ensuring comparability across latitudes and AOI sizes.
#'
#' @param df A data object carrying the AOI (produced by \code{mutate_aoi_*()}).
#'   Must contain an AOI identifier column and have an attached AOI geometry.
#' @param out_dir Directory where NDVI GeoTIFF files are located. Used if
#'   no recorded tiles are found in \code{attr(df, "ndvi_aoi")}.
#' @param deduplicate Logical. If \code{TRUE} (default), AOIs with identical
#'   geometry are extracted only once, then mapped back.
#' @param crs_planar EPSG code of the projected CRS used to compute surface
#'   areas in meters. Defaults to Web Mercator (EPSG:3857).
#' @param prefix Prefix for output column names. Default is \code{"ndvi_"}.
#' @param progress Logical. If \code{TRUE} (default), show progress messages
#'   during extraction.
#'
#' @details
#' - AOI handling: The AOI is validated and transformed to \code{crs_planar}
#'   to ensure that surface area is measured in square meters.
#' - Summaries: Each NDVI value is weighted by its overlap area (m²) when
#'   computing the statistics.
#' - Output: Five new columns are created per AOI:
#'   \code{ndvi_mean}, \code{ndvi_sd}, \code{ndvi_min}, \code{ndvi_max},
#'   \code{ndvi_median}.
#'
#' @section Attributes:
#' Attributes of \code{df} are preserved. AOI and fetch metadata (e.g.,
#' \code{attr(df, "traktr_aoi")}, \code{attr(df, "ndvi_aoi")}) remain available.
#'
#' @return The input \code{df}, augmented with NDVI summary columns.
#'
#' @seealso
#' \code{\link{attach_landuse}},
#' \code{\link{mutate_aoi_h3}}, \code{\link{mutate_aoi_buffer}},
#' \code{\link{fetch_ndvi}}
#'
#' @examples
#' \dontrun{
#' df <- mutate_aoi_h3(my_points, res = 10)
#' df <- fetch_ndvi(df, out_dir = "data/ndvi", from = Sys.Date()-30, to = Sys.Date())
#' df <- attach_ndvi(df, out_dir = "data/ndvi", crs_planar = 3857)
#'
#' # Inspect statistics
#' head(df[, c("h3", "ndvi_mean", "ndvi_sd")])
#' }
#'
#' @export
attach_ndvi <- function(df,
                        out_dir = NULL,
                        deduplicate = TRUE,
                        crs_planar = 3857,
                        prefix = "ndvi_",
                        progress = TRUE) {
  # check dependencies
  rlang::check_installed(c("sf", "terra", "cli", "exactextractr"))
  
  if (isTRUE(progress)) {
    cli::cli_progress_step("Summarizing NDVI (mean, sd, min, max, median)", spinner = TRUE)
  }
  
  # 1. Detect AOI and ID column
  info   <- .aoi_detect(df)
  aoi    <- info$aoi
  id_col <- info$id_col
  
  # 2. Collect NDVI raster files
  files <- .collect_raster_files(df, "ndvi_aoi", out_dir)
  vrt   <- terra::vrt(terra::sprc(files))
  
  # 3. Transform AOI to planar CRS for area-weighted stats
  aoi_planar <- sf::st_transform(aoi, crs_planar)
  
  # 4. Deduplicate AOIs for faster extraction
  dd   <- .aoi_deduplicate(aoi_planar, id_col = id_col, deduplicate = deduplicate)
  au   <- dd$aoi_unique
  map_back <- dd$map_back
  
  # 5. Custom summary function (area-weighted)
  ndvi_summary_fn <- function(values, coverage_fraction) {
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
    s     <- sqrt(var_w)
    vmin  <- min(v)
    vmax  <- max(v)
    
    ord   <- order(v)
    cw    <- cumsum(w[ord])
    med   <- v[ord][which(cw >= 0.5)[1]]
    
    data.frame(
      mean   = m,
      sd     = s,
      min    = vmin,
      max    = vmax,
      median = med
    )
  }
  
  # 6. Extract values
  res <- exactextractr::exact_extract(vrt, au, fun = ndvi_summary_fn, progress = FALSE)
  
  if (is.list(res) && !is.data.frame(res)) {
    res <- do.call(rbind, res)
    rownames(res) <- NULL
  }
  if (!nrow(res)) {
    cli::cli_abort("NDVI extraction returned no results.")
  }
  
  res$.aoi_row <- seq_len(nrow(res))
  names(res)[names(res) != ".aoi_row"] <- paste0(prefix, names(res)[names(res) != ".aoi_row"])
  
  # 7. Map back to all AOIs
  per_aoi <- merge(map_back, res, by = ".aoi_row", all.x = TRUE, sort = FALSE)
  per_aoi$.aoi_row <- NULL
  
  # 8. Join back to df; preserve order and attrs
  df$.__row_order <- seq_len(nrow(df))
  out <- merge(df, per_aoi, by = id_col, all.x = TRUE, sort = FALSE)
  out <- out[order(out$.__row_order), , drop = FALSE]
  out$.__row_order <- NULL
  
  if ("geometry" %in% names(out)) {
    geom_col <- "geometry"
    out <- out[, c(setdiff(names(out), geom_col), geom_col), drop = FALSE]
  }
  
  out <- .restore_attrs_class(df, out)
  
  if (isTRUE(progress)) cli::cli_progress_done()
  invisible(out)
}