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
#' @param band Band (layer) of the GeoTIFFs holding NDVI. Defaults to 1;
#'   further bands (e.g. a data mask or observation count) are not averaged in.
#' @param min_obs Optional minimum number of observations per pixel. Pixels
#'   whose observation count (band 2 of chips from
#'   \code{fetch_ndvi(composite = "median")} or \code{\link{fetch_ndvi_stac}})
#'   is below \code{min_obs} are dropped before summarising. With
#'   \code{fallback = TRUE} the count refers to the tier used for that pixel.
#'   \code{NULL} (default) keeps every pixel with a valid NDVI.
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
#'   \code{ndvi_median}. For median composites (\code{fetch_ndvi()} default,
#'   \code{\link{fetch_ndvi_stac}}) also \code{ndvi_n_obs}, the area-weighted
#'   mean number of observations per pixel, and for
#'   \code{fetch_ndvi(fallback = TRUE)} \code{ndvi_tier}, the worst tier used
#'   in the AOI (1 clear, 2 relaxed cloud mask, 3 maximum value composite).
#'   Use \code{ndvi_tier} for sensitivity analyses, e.g. \code{ndvi_tier == 1}.
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
                        band = 1,
                        min_obs = NULL,
                        progress = TRUE) {
  # check dependencies
  rlang::check_installed(c("sf", "terra", "cli", "exactextractr"))
  
  if (isTRUE(progress)) {
    cli::cli_progress_step("Summarizing NDVI (mean, sd, min, max, median)", spinner = TRUE)
  }
  
  info  <- .aoi_detect(df)
  files <- .collect_raster_files(df, "ndvi_aoi", out_dir)

  # band 2 is an observation count for median/STAC chips; for least-cloudy
  # chips (historical names without method tag) it is only a data mask
  has_count <- !all(grepl("^ndvi_[_0-9]", basename(files)))
  extra <- if (has_count) {
    list(n_obs = list(band = 2L, fun = "mean"),   # coverage-weighted mean
         tier  = list(band = 3L, fun = "max"))    # worst tier in the AOI
  }

  stats <- .extract_continuous(
    info$aoi, info$id_col, files,
    band        = band,
    count_band  = if (is.null(min_obs)) NULL else 2L,
    min_count   = min_obs,
    extra_bands = extra,
    deduplicate = deduplicate,
    crs_planar  = crs_planar
  )
  out <- .join_stats(df, stats, by = info$id_col, prefix = prefix)
  
  if (isTRUE(progress)) cli::cli_progress_done()
  invisible(out)
}
