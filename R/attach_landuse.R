#' Attach land-cover composition to an AOI
#'
#' Extracts ESA WorldCover land-cover classes from GeoTIFF tiles covering an
#' AOI, computes per-class proportions, rollups (green, blue, grey fractions),
#' Shannon diversity, and joins them back to the input data. All fractions are
#' area-weighted in square meters, ensuring comparability across latitudes and
#' AOI sizes.
#'
#' @param df A data object carrying the AOI (produced by \code{mutate_aoi_*()}).
#'   Must contain an AOI identifier column and have an attached AOI geometry.
#' @param out_dir Directory where WorldCover GeoTIFF files are located. Used if
#'   no recorded tiles are found in \code{attr(df, "worldcover_tiles")}.
#' @param deduplicate Logical. If \code{TRUE} (default), AOIs with identical
#'   geometry are extracted only once, then mapped back.
#' @param crs_planar EPSG code of the projected CRS used to compute surface
#'   areas in meters. Defaults to Web Mercator (EPSG:3857).
#' @param include_cropland Logical. If \code{TRUE}, cropland (class 40) is
#'   included in the "green+" rollup.
#' @param prefix Prefix for output column names. Default is \code{"lc_"}.
#' @param progress Logical. If \code{TRUE} (default), show progress messages
#'   during extraction.
#'
#' @details
#' - AOI handling: The AOI is validated and transformed to \code{crs_planar}
#'   to ensure that surface area is measured in square meters.
#' - Summaries: Each land-cover class is weighted by its overlap area (m²).
#' - Output: Per-class proportion columns (e.g. \code{lc_10}, \code{lc_20}),
#'   rollups (\code{frac_green_core}, \code{frac_green_plus},
#'   \code{frac_blue}, \code{frac_grey}), and \code{shannon} diversity.
#'
#' @section Attributes:
#' Attributes of \code{df} are preserved. AOI and fetch metadata (e.g.,
#' \code{attr(df, "traktr_aoi")}, \code{attr(df, "worldcover_tiles")})
#' remain available.
#'
#' @return The input \code{df}, augmented with land-cover composition columns.
#'
#' @seealso
#' \code{\link{attach_ndvi}},
#' \code{\link{mutate_aoi_h3}}, \code{\link{mutate_aoi_buffer}},
#' \code{\link{fetch_landuse}}
#'
#' @examples
#' \dontrun{
#' df <- mutate_aoi_buffer(my_points, dist = 100)
#' df <- fetch_landuse(df, out_dir = "data/worldcover", year = 2021, version = "v200")
#' df <- attach_landuse(df, out_dir = "data/worldcover", crs_planar = 3857)
#'
#' # Inspect composition
#' head(df[, c("aoi_id", "frac_green_core", "shannon")])
#' }
#'
#' @export
attach_landuse <- function(df,
                           out_dir = NULL,
                           deduplicate = TRUE,
                           crs_planar = 3857,
                           include_cropland = FALSE,
                           prefix = "lc_",
                           progress = TRUE) {
  rlang::check_installed(c("sf", "terra", "cli", "exactextractr"))
  
  pb <- if (isTRUE(progress)) {
    cli::cli_progress_bar(
      "attach_landuse()",
      total = 3,
      type = "iterator",
      clear = FALSE,
      status = "attach_landuse()",
      format = "{cli::pb_spin} {cli::pb_name} | {cli::pb_status} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}  {cli::pb_eta}",
      format_done = "{cli::col_green(cli::symbol$tick)} Summarized landuse data {cli::style_dim(paste0('[', cli::pb_elapsed, ']'))}"
    )
  } else NULL
  if (!is.null(pb)) on.exit(cli::cli_progress_done(pb), add = TRUE)
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 1, status = "Validating AOIs", force = TRUE)
  
  # AOI, id
  info <- .aoi_detect(df)
  aoi  <- info$aoi
  id_col <- info$id_col
  
  # files -> VRT
  files <- .collect_raster_files(df, "worldcover_tiles", out_dir)
  wc    <- terra::vrt(terra::sprc(files))
  
  # de-dup AOIs for faster extraction
  dd   <- .aoi_deduplicate(sf::st_transform(aoi, 4326), id_col = id_col, deduplicate = deduplicate)
  au   <- dd$aoi_unique
  map_back <- dd$map_back
  
  # areas (for reference; projected meters)
  aoi_proj <- sf::st_transform(au, crs_planar)
  aoi_area <- data.frame(
    .aoi_row = seq_len(nrow(au)),
    aoi_area_m2 = as.numeric(sf::st_area(aoi_proj)),
    row.names = NULL
  )
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 2, status = "Building VRT / reading tiles", force = TRUE)
  
  # category summary callback -> area-weighted sums per class
  cat_summary_fn <- function(values, coverage_fraction) {
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
  
  per_aoi_class_sums <- exactextractr::exact_extract(wc, au, fun = cat_summary_fn, progress = FALSE)
  if (is.list(per_aoi_class_sums) && !is.data.frame(per_aoi_class_sums)) {
    per_aoi_class_sums <- do.call(rbind, per_aoi_class_sums)
    rownames(per_aoi_class_sums) <- NULL
  }
  
  n_aoi <- nrow(au)
  if (nrow(per_aoi_class_sums) < n_aoi) {
    per_aoi_class_sums[n_aoi, ] <- NA
    rownames(per_aoi_class_sums) <- NULL
  }
  if (ncol(per_aoi_class_sums) > 0) per_aoi_class_sums[is.na(per_aoi_class_sums)] <- 0
  per_aoi_class_sums$.aoi_row <- seq_len(n_aoi)
  
  # totals and proportions
  lc_cols <- grep(paste0("^", prefix), names(per_aoi_class_sums), value = TRUE)
  if (length(lc_cols) == 0L) {
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row, weight_sum = numeric(n_aoi))
    comp_long <- tibble::tibble(.aoi_row = integer(0), lc = integer(0), p = numeric(0))
  } else {
    weight_sum <- rowSums(per_aoi_class_sums[, lc_cols, drop = FALSE], na.rm = TRUE)
    cov_tbl <- data.frame(.aoi_row = per_aoi_class_sums$.aoi_row, weight_sum = weight_sum)
    
    tmp <- per_aoi_class_sums[, c(".aoi_row", lc_cols), drop = FALSE]
    id_cols <- setdiff(names(tmp), lc_cols)
    comp_long <- cbind(tmp[id_cols], utils::stack(tmp[lc_cols]))
    names(comp_long)[names(comp_long) == "ind"] <- "lc_name"
    names(comp_long)[names(comp_long) == "values"] <- "w"
    
    comp_long <- comp_long[!is.na(comp_long$w) & comp_long$w > 0, , drop = FALSE]
    comp_long$lc <- as.integer(sub(prefix, "", comp_long$lc_name))
    
    # join totals for proportions
    x <- comp_long[, c(".aoi_row", "lc", "w"), drop = FALSE]
    idx <- match(x$.aoi_row, cov_tbl$.aoi_row)
    comp_long <- cbind(x, cov_tbl[idx, "weight_sum", drop = FALSE])
    comp_long$p <- ifelse(comp_long$weight_sum > 0, comp_long$w / comp_long$weight_sum, NA_real_)
    comp_long <- unique(comp_long[, c(".aoi_row", "lc", "p"), drop = FALSE])
  }
  
  # wide composition (proportions per class)
  if (nrow(comp_long) > 0) {
    comp_wide <- comp_long
    comp_wide$lc_name <- paste0(prefix, comp_wide$lc)
    tab <- stats::xtabs(p ~ .aoi_row + lc_name, comp_wide)
    comp_wide <- cbind(.aoi_row = utils::type.convert(rownames(tab), as.is = TRUE),
                       as.data.frame.matrix(tab))
  } else {
    comp_wide <- tibble::tibble(.aoi_row = integer(0))
  }
  
  # rollups (+ optional cropland)
  veg_core <- c(10, 20, 30, 90, 95, 100)
  veg_plus <- if (include_cropland) c(veg_core, 40) else veg_core
  blue     <- 80
  grey     <- c(50, 60)
  
  sum_cols <- function(df, codes) {
    cols <- paste0(prefix, codes)
    cols <- cols[cols %in% names(df)]
    if (!length(cols)) return(rep(0, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  stats_tbl <- comp_wide
  stats_tbl$frac_green_core <- sum_cols(stats_tbl, veg_core)
  stats_tbl$frac_green_plus <- sum_cols(stats_tbl, veg_plus)
  stats_tbl$frac_blue       <- sum_cols(stats_tbl, blue)
  stats_tbl$frac_grey       <- sum_cols(stats_tbl, grey)
  
  # Shannon entropy across all land-cover proportion columns
  prob_cols <- names(stats_tbl)[startsWith(names(stats_tbl), prefix)]
  if (length(prob_cols)) {
    P <- data.matrix(stats_tbl[, prob_cols, drop = FALSE])
    P[P <= 0] <- NA_real_
    stats_tbl$shannon <- -rowSums(P * log(P), na.rm = TRUE)
  } else {
    stats_tbl$shannon <- numeric(nrow(stats_tbl))
  }
  
  if (!is.null(pb)) cli::cli_progress_update(id = pb, set = 3, status = "Rollups & join", force = TRUE)
  
  # Assemble per-AOI table and map back (if deduped)
  per_aoi <- data.frame(.aoi_row = seq_len(nrow(au)), check.names = FALSE)
  per_aoi <- merge(per_aoi, cov_tbl,   by = ".aoi_row", all.x = TRUE, sort = FALSE)
  per_aoi <- merge(per_aoi, aoi_area,  by = ".aoi_row", all.x = TRUE, sort = FALSE)
  per_aoi <- merge(per_aoi, stats_tbl, by = ".aoi_row", all.x = TRUE, sort = FALSE)
  
  # add id_col for join back
  per_aoi$id_col_tmp__ <- sf::st_drop_geometry(au)[[id_col]]
  per_aoi <- per_aoi[, c("id_col_tmp__", setdiff(names(per_aoi), "id_col_tmp__")), drop = FALSE]
  names(per_aoi)[1] <- id_col
  
  per_aoi_full <- merge(map_back, per_aoi, by = c(id_col, ".aoi_row"), all.x = TRUE, sort = FALSE)
  per_aoi_full$.aoi_row <- NULL
  
  # Join to df; preserve order & class/attrs
  df$.__row_order <- seq_len(nrow(df))
  out <- merge(df, per_aoi_full, by = id_col, all.x = TRUE, sort = FALSE)
  out <- out[order(out$.__row_order), , drop = FALSE]
  out$.__row_order <- NULL
  
  # geometry last, if present
  if ("geometry" %in% names(out)) {
    geom_col <- "geometry"
    out <- out[, c(setdiff(names(out), geom_col), geom_col), drop = FALSE]
  }
  
  out <- .restore_attrs_class(df, out)
  
  invisible(out)
}