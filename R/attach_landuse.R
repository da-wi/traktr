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
  rlang::check_installed(c("sf","terra","dplyr","tidyr","purrr","cli"))
  
  #pb <- cli::cli_progress_bar(
  #  name   = "attach_landuse()", 
  #  total  = 3,                                   # number of major steps
  #  clear  = FALSE,                               # keep final line; set TRUE to remove
  #  format = "{cli::pb_spin} {name} | {msg} {pb} {n}/{total}  {eta}"
  #)
  pb <- cli::cli_progress_bar("attach_landuse()",
                              total = 3,
                              type = "iterator", 
                              clear = F,
                              status = "attach_landuse()",
                              format = "{cli::pb_spin} {cli::pb_name} | {cli::pb_status} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}  {cli::pb_eta}",
                              format_done = "{cli::col_green(cli::symbol$tick)} Summarized landuse data {cli::style_dim(paste0('[', cli::pb_elapsed, ']'))}"
  )
  on.exit(cli::cli_progress_done(pb), add = TRUE) # always close, even on error
  
  #cli::cli_progress_step("Summarizing landuse data", spinner = F, msg_done = "Landuse data written. ")
  cli::cli_progress_update(pb, set = 1, status = "Validating AOIs", force = TRUE)
  #cli::cli_progress_message(msg = "Validating AOIs1" )
  
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    cli::cli_abort("No AOI found – call a mutate_aoi_*() function first.")
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
  
  tiles <- attr(df, "worldcover_tiles")
  if (is.null(tiles)) {
    files <- list.files(out_dir, "\\.tif$", full.names = TRUE)
    if (!length(files)) cli::cli_abort("No worldcover files found (empty out_dir).")
  } else {
    files <- dplyr::filter(tiles, !is.na(.data$file))$file
    if (!length(files)) cli::cli_abort("No worldcover files recorded in attr(df,'worldcover_tiles').")
  }
  
  if(".aoi_row" %in% names(df)) { df <- dplyr::select(-.aoi_row) }
  

  # Optional deduplication of AOIs to save work
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
  
  # Area (for reference)
  aoi_area <- aoi_unique |>
    sf::st_transform(crs_planar) |>
    dplyr::transmute(.aoi_row, aoi_area_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry()
  
  # Extract categorical values with area weights

  # files <- list.files(out_dir, "\\.tif$", full.names = TRUE)

  # Extract categorical values with area weights (fast path via exactextractr)
  wc <- terra::vrt(terra::sprc(files))
  
  # require exactextractr for the fast callback approach
  rlang::check_installed("exactextractr")
  
  # callback: for one polygon, receives (values, coverage_fraction)
  # returns a one-row data.frame with columns lc_<code> containing summed coverage fractions
  cat_summary_fn <- function(values, coverage_fraction) {
    valid <- !is.na(values) & coverage_fraction > 0
    values <- values[valid]
    weights <- coverage_fraction[valid]
    
    if (length(values) == 0L) {
      return(data.frame()) # no data for this AOI -> exact_extract will produce an empty row
    }
    
    # weighted sum of coverage per category (values are integer class codes)
    tab <- tapply(weights, values, sum)
    if (length(tab) == 0L) return(data.frame())
    
    nm <- paste0(prefix, names(tab))
    df <- as.data.frame(as.list(setNames(as.numeric(tab), nm)))
    df
  }
  
  cli::cli_progress_update(pb, set = 2, status = "Building VRT / reading tiles", force = T)
  
  # run exact_extract with our callback -> returns one-row summary per AOI (columns vary by classes present)
  per_aoi_class_sums <- exactextractr::exact_extract(
    wc,
    aoi_unique,
    fun = cat_summary_fn,
    progress = FALSE
  )
  
  # exact_extract may return a list or data.frame depending on callback; normalize to data.frame
  if (is.list(per_aoi_class_sums) && !is.data.frame(per_aoi_class_sums)) {
    # bind rows; missing columns will become NA (we'll set to 0)
    per_aoi_class_sums <- do.call(
      what = function(...) dplyr::bind_rows(list(...)),
      args = per_aoi_class_sums
    )
  }
  
  # Ensure we have one row per AOI in same order as aoi_unique
  # If exact_extract produced fewer rows (rare), pad with empty rows
  n_aoi <- nrow(aoi_unique)
  if (nrow(per_aoi_class_sums) < n_aoi) {
    missing_rows <- n_aoi - nrow(per_aoi_class_sums)
    per_aoi_class_sums <- dplyr::bind_rows(
      per_aoi_class_sums,
      as.data.frame(matrix(NA, nrow = missing_rows, ncol = ncol(per_aoi_class_sums)))
    )
    names(per_aoi_class_sums) <- names(per_aoi_class_sums)
  }
  
  # replace NA with 0 (no coverage for that class in that AOI)
  if (ncol(per_aoi_class_sums) > 0) per_aoi_class_sums[is.na(per_aoi_class_sums)] <- 0
  
  # add .aoi_row so downstream joins work exactly like before
  per_aoi_class_sums$.aoi_row <- seq_len(n_aoi)
  
  # build cov_tbl: total covered fraction per AOI (sum over all lc_* columns)
  lc_cols <- grep(paste0("^", prefix), names(per_aoi_class_sums), value = TRUE)
  if (length(lc_cols) == 0L) {
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row,
                          weight_sum = numeric(n_aoi))
    comp_long <- dplyr::tibble(.aoi_row = integer(0), lc = integer(0), p = numeric(0))
  } else {
    weight_sum <- rowSums(per_aoi_class_sums[, lc_cols, drop = FALSE], na.rm = TRUE)
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row,
                          weight_sum = weight_sum)
    
    # melt per_aoi_class_sums into long (.aoi_row, lc, w)
    comp_long <- per_aoi_class_sums |>
      dplyr::select(.aoi_row, dplyr::all_of(lc_cols)) |>
      tidyr::pivot_longer(
        cols = dplyr::starts_with(prefix),
        names_to = "lc_name",
        values_to = "w"
      ) |>
      dplyr::mutate(
        lc = as.integer(sub(prefix, "", lc_name))
      ) |>
      dplyr::filter(!is.na(w) & w > 0) |>
      dplyr::select(.aoi_row, lc, w) |>
      dplyr::left_join(cov_tbl, by = ".aoi_row") |>
      dplyr::mutate(p = dplyr::if_else(weight_sum > 0, w / weight_sum, NA_real_)) |>
      dplyr::select(.aoi_row, lc, p) |>
      dplyr::distinct()
  }
  

  comp_wide <- comp_long |>
    dplyr::mutate(lc_name = paste0(prefix, lc)) |>
    tidyr::pivot_wider(
      id_cols    = .aoi_row,
      names_from = lc_name,
      values_from= p,
      values_fill= 0
    )

  # Rollups + Shannon diversity 
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
  
  cli::cli_progress_update(id=pb, set = 3, status = "Rollups & Shannon / join", force= T)
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
  
  # Bring everything together per AOI row 
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
  
  # Join to df 
  out <- df |>
    dplyr::left_join(per_aoi_full, by = id_col) |>
    dplyr::select(-c(.aoi_row,weight_sum,aoi_area_m2)) |>
    relocate(geometry, .after = dplyr::last_col())
  
  out
  
}
