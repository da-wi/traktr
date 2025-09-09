#' Attach ESA WorldCover land-use composition per AOI (area-weighted)
#'
#' Computes per-class fractions from WorldCover (10 m, categorical) for each AOI
#' stored in `attr(df, "traktr_aoi")` and joins results back to `df` by the AOI id.
#'
#' @param df A data frame of GPS or other observations that already has an AOI id
#'   column (e.g. created via [mutate_aoi_h3()] or [mutate_aoi_buffer()]) and
#'   AOI metadata in `attr(df, "traktr_aoi")` and `attr(df, "traktr_aoi_params")`.
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
  rlang::check_installed(c("sf", "terra", "dplyr", "tidyr", "purrr", "cli", "exactextractr"))
  
  pb <- cli::cli_progress_bar(
    "attach_landuse()",
    total = 3,
    type = "iterator",
    clear = FALSE,
    status = "attach_landuse()",
    format = "{cli::pb_spin} {cli::pb_name} | {cli::pb_status} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}  {cli::pb_eta}",
    format_done = "{cli::col_green(cli::symbol$tick)} Summarized landuse data {cli::style_dim(paste0('[', cli::pb_elapsed, ']'))}"
  )
  on.exit(cli::cli_progress_done(pb), add = TRUE)
  
  cli::cli_progress_update(pb, set = 1, status = "Validating AOIs", force = TRUE)
  
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    cli::cli_abort("No AOI found: call a mutate_aoi_*() function first.")
  }
  params <- attr(df, "traktr_aoi_params")
  id_col <- if (!is.null(params) && !is.null(params$id_col)) params$id_col else {
    setdiff(names(aoi), attr(aoi, "sf_column"))[1]
  }
  if (is.null(id_col) || !id_col %in% names(aoi)) {
    cli::cli_abort("Cannot determine AOI id column; set attr(df, 'traktr_aoi_params')$id_col.")
  }
  if (!id_col %in% names(df)) {
    cli::cli_abort(paste0("Input df does not contain AOI id column '", id_col, "'."))
  }
  
  tiles <- attr(df, "worldcover_tiles")
  if (is.null(tiles)) {
    files <- list.files(out_dir, "\\.tif$", full.names = TRUE)
    if (!length(files)) cli::cli_abort("No worldcover files found (empty out_dir).")
  } else {
    # avoid NSE/.data: use base subsetting
    files <- tiles$file[!is.na(tiles$file)]
    if (!length(files)) cli::cli_abort("No worldcover files recorded in attr(df, 'worldcover_tiles').")
  }
  
  if (".aoi_row" %in% names(df)) {
    df <- df[ , !(names(df) %in% ".aoi_row")]
  }
  
  # AOIs (dedup + validity)
  aoi <- sf::st_make_valid(aoi)
  aoi4326 <- sf::st_transform(aoi, 4326)
  
  # h3 special case: drop duplicates by h3 id (base to avoid NSE)
  if (identical(id_col, "h3") && "h3" %in% names(aoi4326)) {
    aoi4326 <- aoi4326[!duplicated(aoi4326$h3), , drop = FALSE]
  }
  
  if (deduplicate) {
    aoi4326$.__wkb <- vapply(sf::st_as_binary(aoi4326$geometry), paste, "", collapse = "")
    # base unique to avoid NSE in distinct()
    aoi_unique <- aoi4326[!duplicated(aoi4326$.__wkb), , drop = FALSE] |>
      dplyr::mutate(.aoi_row = dplyr::row_number())
    
    map_back <- aoi4326 |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(c(id_col, ".__wkb"))) |>
      dplyr::left_join(
        sf::st_drop_geometry(aoi_unique)[, c(".__wkb", ".aoi_row")],
        by = ".__wkb"
      )
  } else {
    aoi_unique <- dplyr::mutate(aoi4326, .aoi_row = dplyr::row_number())
    map_back <- sf::st_drop_geometry(aoi_unique)[, c(id_col, ".aoi_row")]
  }
  
  # Area per AOI (use explicit geometry to avoid NSE)
  #aoi_area <- aoi_unique |>
  # sf::st_transform(crs_planar) |>
  #  dplyr::transmute(
  #    .aoi_row,
  #    aoi_area_m2 = as.numeric(sf::st_area(aoi_unique$geometry))
  #  ) |>
  #  sf::st_drop_geometry()
  
  # ensure .aoi_row exists
  if (!".aoi_row" %in% names(aoi_unique)) {
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
  }
  
  # 1) transform to planar CRS
  aoi_transformed <- sf::st_transform(aoi_unique, crs_planar)
  
  # 2) compute areas (sf::st_area returns a 'units' object)
  areas <- as.numeric(sf::st_area(aoi_transformed))  # numeric in same units as CRS (m if projected in metres)
  
  # 3) build a plain data.frame with .aoi_row and area
  aoi_area <- sf::st_drop_geometry(aoi_transformed)   # keep non-geometry columns
  aoi_area <- data.frame(
    .aoi_row = aoi_area$.aoi_row,
    aoi_area_m2 = areas,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Build VRT and extract
  wc <- terra::vrt(terra::sprc(files))
  
  # callback for exactextractr
  cat_summary_fn <- function(values, coverage_fraction) {
    # vectorized validity mask
    valid <- !is.na(values) & !is.na(coverage_fraction) & coverage_fraction > 0
    if (!any(valid)) return(data.frame())
    
    values  <- values[valid]
    weights <- coverage_fraction[valid]
    
    tab <- tapply(weights, values, sum)
    if (length(tab) == 0L) return(data.frame())
    
    nm <- paste0(prefix, names(tab))
    out <- as.list(as.numeric(tab))
    names(out) <- nm
    as.data.frame(out)
  }
  
  cli::cli_progress_update(pb, set = 2, status = "Building VRT / reading tiles", force = TRUE)
  
  per_aoi_class_sums <- exactextractr::exact_extract(
    wc,
    aoi_unique,
    fun = cat_summary_fn,
    progress = FALSE
  )
  
  if (is.list(per_aoi_class_sums) && !is.data.frame(per_aoi_class_sums)) {
    per_aoi_class_sums <- do.call(
      what = function(...) dplyr::bind_rows(list(...)),
      args = per_aoi_class_sums
    )
  }
  
  # Ensure 1 row per AOI
  n_aoi <- nrow(aoi_unique)
  if (nrow(per_aoi_class_sums) < n_aoi) {
    missing_rows <- n_aoi - nrow(per_aoi_class_sums)
    per_aoi_class_sums <- dplyr::bind_rows(
      per_aoi_class_sums,
      as.data.frame(matrix(NA_real_, nrow = missing_rows, ncol = ncol(per_aoi_class_sums)))
    )
    names(per_aoi_class_sums) <- names(per_aoi_class_sums)
  }
  
  # Replace NA class weights with 0; add AOI row id
  if (ncol(per_aoi_class_sums) > 0) per_aoi_class_sums[is.na(per_aoi_class_sums)] <- 0
  per_aoi_class_sums$.aoi_row <- seq_len(n_aoi)
  
  # Totals per AOI
  lc_cols <- grep(paste0("^", prefix), names(per_aoi_class_sums), value = TRUE)
  if (length(lc_cols) == 0L) {
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row, weight_sum = numeric(n_aoi))
    comp_long <- dplyr::tibble(.aoi_row = integer(0), lc = integer(0), p = numeric(0))
  } else {
    weight_sum <- rowSums(per_aoi_class_sums[, lc_cols, drop = FALSE], na.rm = TRUE)
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row, weight_sum = weight_sum)
    
    # Long format without NSE
    tmp <- dplyr::select(per_aoi_class_sums, dplyr::all_of(c(".aoi_row", lc_cols)))
    comp_long <- tidyr::pivot_longer(
      tmp,
      cols = dplyr::all_of(lc_cols),
      names_to = "lc_name",
      values_to = "w"
    )
    comp_long <- comp_long[!is.na(comp_long$w) & comp_long$w > 0, , drop = FALSE]
    comp_long$lc <- as.integer(sub(prefix, "", comp_long$lc_name))
    
    comp_long <- dplyr::select(comp_long, dplyr::all_of(c(".aoi_row", "lc", "w"))) |>
      dplyr::left_join(cov_tbl, by = ".aoi_row")
    
    # proportions (leave NA when weight_sum is 0 to match previous logic)
    comp_long$p <- ifelse(comp_long$weight_sum > 0, comp_long$w / comp_long$weight_sum, NA_real_)
    
    comp_long <- dplyr::select(comp_long, dplyr::all_of(c(".aoi_row", "lc", "p"))) |>
      dplyr::distinct()
  }
  
  # Wide composition (strings in pivot args to avoid NSE)
  if (nrow(comp_long) > 0) {
    comp_wide <- comp_long
    comp_wide$lc_name <- paste0(prefix, comp_wide$lc)
    comp_wide <- tidyr::pivot_wider(
      comp_wide,
      id_cols     = ".aoi_row",
      names_from  = "lc_name",
      values_from = "p",
      values_fill = 0
    )
  } else {
    comp_wide <- dplyr::tibble(.aoi_row = integer(0))
  }
  
  # Rollups + Shannon diversity
  veg_core <- c(10, 20, 30, 90, 95, 100)
  veg_plus <- if (include_cropland) c(veg_core, 40) else veg_core
  blue <- c(80)
  grey <- c(50, 60)
  
  sum_cols <- function(df, codes) {
    cols <- paste0(prefix, codes)
    cols <- cols[cols %in% names(df)]
    if (!length(cols)) return(rep(0, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  cli::cli_progress_update(id = pb, set = 3, status = "Rollups & Shannon / join", force = TRUE)
  
  stats_tbl <- comp_wide |>
    dplyr::rowwise() |>
    dplyr::mutate(
      frac_green_core = sum_cols(dplyr::cur_data(), veg_core)[1],
      frac_green_plus = sum_cols(dplyr::cur_data(), veg_plus)[1],
      frac_blue       = sum_cols(dplyr::cur_data(), blue)[1],
      frac_grey       = sum_cols(dplyr::cur_data(), grey)[1],
      shannon = {
        probs <- as.numeric(dplyr::c_across(tidyselect::starts_with(prefix)))
        probs <- probs[probs > 0]
        if (length(probs)) -sum(probs * log(probs)) else 0
      }
    ) |>
    dplyr::ungroup()
  
  # Per-AOI table
  per_aoi <- aoi_unique |>
    sf::st_drop_geometry() |>
    dplyr::select(dplyr::all_of(c(".aoi_row", id_col))) |>
    dplyr::left_join(cov_tbl,  by = ".aoi_row") |>
    dplyr::left_join(aoi_area, by = ".aoi_row") |>
    dplyr::left_join(stats_tbl, by = ".aoi_row") |>
    dplyr::select(-dplyr::all_of(".aoi_row"))
  
  # Map back (if deduplicated)
  per_aoi_full <- map_back |>
    dplyr::left_join(per_aoi, by = id_col) |>
    dplyr::select(-dplyr::any_of(c(".__wkb")))
  
  # Join to df; move geometry to end if present
  out <- df |>
    dplyr::left_join(per_aoi_full, by = id_col) |>
    dplyr::select(-dplyr::any_of(c(".aoi_row", "weight_sum", "aoi_area_m2")))
  
  if ("geometry" %in% names(out)) {
    out <- dplyr::relocate(out, dplyr::all_of("geometry"), .after = dplyr::last_col())
  }
  
  out
}