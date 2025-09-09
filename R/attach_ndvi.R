#' Attach NDVI summaries back to the input data frame
#'
#' Reads fetched NDVI chip files (GeoTIFFs), mosaics them into a virtual raster,
#' extracts NDVI statistics per AOI (from `attr(df, "traktr_aoi")`), and joins the
#' results back to the input tibble.
#'
#' @param df A data frame of GPS or other observations that already has an AOI id
#'   column (e.g. created via [mutate_aoi_h3()] or [mutate_aoi_buffer()]) and
#'   AOI metadata in `attr(df, "traktr_aoi")` and `attr(df, "traktr_aoi_params")`.
#' @param out_dir Directory containing previously fetched NDVI chip `.tif` files
#'   (as created by [fetch_ndvi()]).
#' @param summary Summary function to apply in \pkg{exactextractr}
#'   (e.g. `"mean"`, `"median"`, `"stdev"`). Default `"mean"`.
#' @param drop_geometry Logical; if `TRUE`, the AOI geometry column is
#'   dropped after joining. If `FALSE` (default), an `sf` object is returned with geometries
#'   retained.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Collects NDVI chip files either from `attr(df,"ndvi_aoi")` or by
#'   listing `out_dir`.
#'   \item Builds a lazy mosaic (VRT) using \pkg{terra}.
#'   \item Reprojects the AOIs to match the raster CRS.
#'   \item Uses \pkg{exactextractr::exact_extract()} to compute per-AOI summaries
#'   on the first band (NDVI), ignoring masked areas.
#'   \item Joins the summary values back into `df` by the AOI id column
#'   (`attr(df,"ndvi_aoi_params")$id_col`).
#' }
#'
#' For H3-based AOIs, duplicate geometries are deduplicated before extraction.
#'
#' The created NDVI column is named `ndvi_<summary>` (e.g. `ndvi_mean`).
#'
#' @return A data frame (class `"ndvi_df"`) with an additional NDVI summary column.
#' If `drop_geometry = FALSE`, the result is an `sf` object with geometries.
#'
#' @examples
#' \dontrun{
#' gps <- tibble::tibble(lat = c(47.38, 47.39), lon = c(8.53, 8.54))
#' df <- gps |>
#'   mutate_aoi_h3(lat, lon, res = 10) |>
#'   fetch_ndvi(from = as.Date("2025-04-01"),
#'              to   = as.Date("2025-07-31"),
#'              out_dir = "output/ndvi") |>
#'   attach_ndvi(out_dir = "output/ndvi", summary = "mean")
#'
#' head(df$ndvi_mean)
#' }
#'
#' @seealso [mutate_aoi_h3()], [mutate_aoi_buffer()], [fetch_ndvi()]
#' @export
attach_ndvi <- function(df, out_dir, summary = "mean", drop_geometry = FALSE) {
  cli::cli_progress_step("Summarizing NDVI data using {cli::col_blue(summary)}",
                         msg_done = "Summarized NDVI data using {cli::col_blue(summary)}",
                         spinner = TRUE)
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    cli::cli_abort("No AOI found: call a mutate_aoi_*() function first.")
  }
  params <- attr(df, "traktr_aoi_params")
  id_col <- if (!is.null(params) && !is.null(params$id_col)) params$id_col else {
    # fallback: first non-geometry column in AOI
    setdiff(names(aoi), attr(aoi, "sf_column"))[1]
  }
  if (is.null(id_col) || !id_col %in% names(aoi))
    cli::cli_abort("Cannot determine AOI id column; set attr(df, 'traktr_aoi_params')$id_col.")
  
  # ensure df has the id column to join on
  if (!id_col %in% names(df)) {
    cli::cli_abort(paste0("Input df does not contain AOI id column '", id_col, "'."))
  }
  
  # Collect chip files
  tiles <- attr(df, "ndvi_aoi")
  if (is.null(tiles)) {
    files <- list.files(out_dir, "\\.tif$", full.names = TRUE)
    if (!length(files)) cli::cli_abort("No NDVI files found (empty out_dir).")
  } else {
    files <- dplyr::filter(tiles, !is.na(tiles$file))$file
    if (!length(files)) cli::cli_abort("No NDVI files recorded in attr(df,'ndvi_aoi').")
  }
  
  if ("geometry" %in% names(df)) {
    #df <- subset(df, select = -geometry)
    df <- df[ , !(names(df) %in% "geometry")]
  }
  
  # create virtual raster from files
  virtual_raster <- terra::vrt(terra::sprc(files))
  aoi_r <- sf::st_transform(aoi, terra::crs(virtual_raster))
  
  # special case h3 
  if (params$id_col == "h3") {
    # extract for each distinct h3
    # keep first occurrence of each h3 (preserves all columns, including geometry)
    aoi_r <- aoi_r[!duplicated(aoi_r$h3), , drop = FALSE]
  }
  

  new_col <- sprintf("ndvi_%s", summary)
  
  vals <- exactextractr::exact_extract(virtual_raster, aoi_r, summary, progress = FALSE)[, 1]
  
  if (length(vals) != nrow(aoi_r)) {
    cli::cli_abort("Length mismatch: exact_extract returned ", length(vals),
         " rows but aoi_r has ", nrow(aoi_r))
  }
  
  # assign (works for sf objects too)
  aoi_r[[new_col]] <- vals
  
  # then join back
  df <- df |> sf::st_drop_geometry() |> dplyr::inner_join(aoi_r, by = id_col)
  
  cli::cli_process_done()
  
  # drop geometry
  if (drop_geometry == TRUE) {
    df |> sf::st_drop_geometry()
  } else {
    df |> dplyr::relocate(dplyr::all_of("geometry"), .after = dplyr::last_col()) |> sf::st_as_sf()
  }
}
