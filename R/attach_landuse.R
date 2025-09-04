#' Attach ESA WorldCover land-use composition per AOI (area-weighted)
#'
#' Computes per-class fractions from WorldCover (10 m, categorical) for each AOI
#' stored in `attr(df, "ndvi_aoi")` and joins results back to `df` by the AOI id.
#'
#' @param df A data frame previously augmented with AOIs
#'           (class "ndvi_df", attributes "ndvi_aoi" and "ndvi_aoi_params").
#' @param out_dir Out dir
#'           If NULL, the function tries `attr(df, "worldcover")` or `attr(df, "landcover")`.
#' @param exact Use area-weighted extraction (`terra::extract(..., exact=TRUE)`).
#' @param deduplicate If TRUE, AOIs with identical geometry are computed once.
#' @param crs_planar CRS used for area computations (default 3857).
#' @param include_cropland If TRUE, cropland (40) is counted in `frac_green_plus`.
#' @param prefix Column name prefix for per-class fractions (default "lc_").
#'
#' @return `df` with new columns:
#'   - per-class fractions `lc_<code>` (0..1),
#'   - rollups: `frac_green_core`, `frac_green_plus`, `frac_blue`, `frac_grey`,
#'   - `shannon` (land-cover diversity), `wc_cover` (covered fraction), `aoi_area_m2`.
#'   Also adds attribute `attr(df, "landuse_table")` (per-AOI table).
#'
#' @examples
#' # gps |> mutate_aoi_buffer(lat, lon, dist=75) |> fetch_worldcover() |> attach_landuse()
#' @export
attach_landuse <- function(df,
                           out_dir = NULL,
                           exact = TRUE,
                           deduplicate = TRUE,
                           crs_planar = 3857,
                           include_cropland = FALSE,
                           prefix = "lc_") {
  rlang::check_installed(c("sf","terra","dplyr","tidyr","purrr","cli"))
  
  cli::cli_progress_step(paste0("Summarizing landuse data"), spinner = TRUE)
  aoi <- attr(df, "ndvi_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    rlang::abort("No AOI found – call a mutate_aoi_*() function first.")
  }
  params <- attr(df, "ndvi_aoi_params")
  id_col <- if (!is.null(params) && !is.null(params$id_col)) params$id_col else {
    # fallback: first non-geometry column in AOI
    setdiff(names(aoi), attr(aoi, "sf_column"))[1]
  }
  if (is.null(id_col) || !id_col %in% names(aoi))
    rlang::abort("Cannot determine AOI id column; set attr(df, 'ndvi_aoi_params')$id_col.")
  
  # ensure df has the id column to join on
  if (!id_col %in% names(df)) {
    rlang::abort(paste0("Input df does not contain AOI id column '", id_col, "'."))
  }
  
  tiles <- attr(df, "worldcover_tiles")
  if (is.null(tiles)) {
    files <- list.files(out_dir, "\\.tif$", full.names = TRUE)
    if (!length(files)) rlang::abort("No worldcover files found (empty out_dir).")
  } else {
    files <- dplyr::filter(tiles, !is.na(.data$file))$file
    if (!length(files)) rlang::abort("No worldcover files recorded in attr(df,'worldcover').")
  }
  
  # ---- 1) Optional deduplication of AOIs to save work ----
  aoi <- sf::st_make_valid(aoi)
  aoi4326 <- sf::st_transform(aoi, 4326)
  
  
  # special case h3 
  if (params$id_col == "h3") {
    # extract for each distinct h3
    aoi4326 <- aoi4326 |>
      dplyr::distinct(h3, .keep_all = T)
  }
  
  if (deduplicate) {
    aoi4326$.__wkb <- vapply(sf::st_as_binary(aoi4326$geometry), paste, "", collapse = "")
    aoi_unique <- aoi4326 |>
      dplyr::distinct(.__wkb, .keep_all = TRUE) |>
      dplyr::mutate(.aoi_row = dplyr::row_number())
    map_back <- aoi4326 |>
      sf::st_drop_geometry() |>
      dplyr::select(!!id_col, .__wkb) |>
      dplyr::left_join(sf::st_drop_geometry(aoi_unique)[, c(".__wkb", ".aoi_row")], by=".__wkb")
  } else {
    aoi_unique <- dplyr::mutate(aoi4326, .aoi_row = dplyr::row_number())
    map_back <- sf::st_drop_geometry(aoi_unique)[, c(id_col, ".aoi_row")]
  }
  
  # ---- 2) Area (for reference) ----
  aoi_area <- aoi_unique |>
    sf::st_transform(crs_planar) |>
    dplyr::transmute(.aoi_row, aoi_area_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry()
  
  # ---- 3) Extract categorical values with area weights ----
  #cli::cli_progress_step("Extracting WorldCover per AOI")
  
  files <- list.files(out_dir, "\\.tif$", full.names = TRUE)

  wc <- terra::vrt(terra::sprc(files))
  
  # Keep name of the categorical band
  wc_band <- names(wc)[1]
  ext <- terra::extract(wc, terra::vect(aoi_unique), exact = exact, cells = TRUE)
  # ext: columns: ID (row idx of aoi_unique), cell, <band>, fraction
  names(ext)[names(ext) == wc_band] <- "lc"
  ext <- dplyr::as_tibble(ext) |>
    dplyr::filter(!is.na(lc), fraction > 0)
  
  # total covered fraction per AOI (handle NODATA / partial coverage)
  cov_tbl <- ext |>
    dplyr::group_by(ID) |>
    dplyr::summarise(weight_sum = sum(fraction), .groups = "drop") |>
    dplyr::rename(.aoi_row = ID)
  
  # ---- 4) Per-class fractions (normalized within AOI) ----
  comp_long <- ext |>
    dplyr::group_by(ID, lc) |>
    dplyr::summarise(w = sum(fraction), .groups = "drop") |>
    dplyr::rename(.aoi_row = ID) |>
    dplyr::left_join(cov_tbl, by = ".aoi_row") |>
    dplyr::mutate(p = dplyr::if_else(weight_sum > 0, w / weight_sum, NA_real_)) |>
    dplyr::select(.aoi_row, lc, p) |>
    dplyr::distinct()
  
  
  comp_wide <- comp_long |>
    dplyr::mutate(lc_name = paste0(prefix, lc)) |>
    tidyr::pivot_wider(
      id_cols    = .aoi_row,
      names_from = lc_name,
      values_from= p,
      values_fill= 0
    )
  
  # ---- 5) Rollups + Shannon diversity ----
  # ESA WorldCover v200 codes
  veg_core <- c(10, 20, 30, 90, 95, 100)
  veg_plus <- if (include_cropland) c(veg_core, 40) else veg_core
  blue <- c(80)
  grey <- c(50, 60)
  
  # Build helper to sum a set of class columns if present
  sum_cols <- function(df, codes) {
    cols <- paste0(prefix, codes)
    cols <- cols[cols %in% names(df)]
    if (!length(cols)) return(rep(0, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  stats_tbl <- comp_wide |>
    dplyr::rowwise() |>
    dplyr::mutate(
      frac_green_core = sum_cols(cur_data(), veg_core)[1],
      frac_green_plus = sum_cols(cur_data(), veg_plus)[1],
      frac_blue       = sum_cols(cur_data(), blue)[1],
      frac_grey       = sum_cols(cur_data(), grey)[1],
      shannon = {
        probs <- as.numeric(dplyr::c_across(tidyselect::starts_with(prefix)))
        probs <- probs[probs > 0]
        if (length(probs)) -sum(probs * log(probs)) else 0
      }
    ) |>
    dplyr::ungroup()
  
  # ---- 6) Bring everything together per AOI row ----
  per_aoi <- aoi_unique |>
    sf::st_drop_geometry() |>
    dplyr::select(.aoi_row, !!id_col) |>
    dplyr::left_join(cov_tbl,  by = ".aoi_row") |>
    dplyr::left_join(aoi_area, by = ".aoi_row") |>
    dplyr::left_join(stats_tbl, by = ".aoi_row") |>
    dplyr::select(-.aoi_row)
  
  # Map back to all AOI ids (if deduplicated)
  per_aoi_full <- map_back |>
    dplyr::left_join(per_aoi, by = id_col) |>
    dplyr::select(-dplyr::any_of(c(".__wkb")))
  
  # ---- 7) Join to df ----
  out <- df |>
    dplyr::left_join(per_aoi_full, by = id_col)
  
  cli::cli_progress_done()
  
  # keep a tidy per-AOI table as attribute
  # attr(out, "landuse_table") <- per_aoi
  out
  
}
