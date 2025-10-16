#' Resolve AOI from input or object attribute
#'
#' Returns an AOI \code{sf} object from the explicit \code{aoi_obj} argument
#' or from \code{attr(df, "traktr_aoi")}. Validates that the AOI is an \code{sf}
#' polygon/multipolygon with a defined CRS. Aborts with a clear message if the
#' AOI is missing or has no CRS.
#'
#' @param df An object that may carry \code{attr(df, "traktr_aoi")}.
#' @param aoi_obj Optional \code{sf} polygon/multipolygon to use instead of the
#'   attribute.
#' @param error_call Passed to \code{cli::cli_abort()} to attribute the error
#'   to the caller; typically \code{rlang::caller_env()}.
#'
#' @return An \code{sf} object representing the AOI.
#'
#' @details
#' If \code{aoi_obj} is \code{NULL}, the function looks up
#' \code{attr(df, "traktr_aoi")}. The AOI must have a non-NA CRS.
#'
#' @seealso \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @keywords internal
#' @noRd
.resolve_aoi <- function(df, aoi_obj, error_call = rlang::caller_env()) {
  x <- if (!is.null(aoi_obj)) aoi_obj else attr(df, "traktr_aoi")
  if (is.null(x) || !inherits(x, "sf")) {
    cli::cli_abort(
      "No AOI found: pass `aoi` or call a mutate_aoi_*() function first.",
      call = error_call
    )
  }
  if (is.na(sf::st_crs(x))) {
    cli::cli_abort("`aoi` has no CRS; assign one first.", call = error_call)
  }
  x
}

#' Store layer outputs and parameters on a traktr object
#'
#' Inserts or updates a named entry in \code{attr(df, "traktr_layers")}, storing
#' the layer index (e.g., an \code{sf} tibble) and the parameters used to create
#' it. Also ensures class \code{"traktr_df"} is present on \code{df}.
#'
#' @param df A traktr-like object.
#' @param name Character scalar; the layer name (e.g., "ndvi", "worldcover").
#' @param index The per-layer index object to store (e.g., \code{sf} or tibble).
#' @param params A named list of parameters used to build the layer.
#'
#' @return \code{df} with updated \code{traktr_layers} attribute (and class).
#'
#' @details
#' The \code{traktr_layers} attribute is a list. Each element is a list with
#' fields \code{index} and \code{params}. Existing entries with the same
#' \code{name} are replaced.
#'
#' @seealso \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @keywords internal
#' @noRd
.set_traktr_layer <- function(df, name, index, params) {
  layers <- attr(df, "traktr_layers")
  if (is.null(layers)) layers <- list()
  layers[[name]] <- list(index = index, params = params)
  attr(df, "traktr_layers") <- layers
  class(df) <- union("traktr_df", class(df))
  df
}

#' Sanitize numeric values for stable filenames
#'
#' Formats numeric-like values without scientific notation and replaces hyphens
#' and dots with underscores, producing cache-friendly filename components.
#'
#' @param x A numeric or character vector to sanitize.
#'
#' @return A character vector with \code{"-"} and \code{"."} replaced by
#'   \code{"_"} after non-scientific formatting.
#'
#' @examples
#' sanitize_num(c(-7.5, 12.0, 1e-6))
#' # " _7_5", "12_0", "0_000001"
#'
#' @keywords internal
#' @noRd
.sanitize_num <- function(x) {
  gsub("[-.]", "_", format(x, scientific = FALSE))
}


#' AOI helpers (internal)
#'
#' Internal utilities used by AOI mutate/attach functions for validation,
#' raster file collection, deduplication, and attribute restoration.
#'
#' These functions are not exported and are intended for internal use only.
#'
#' @keywords internal
#' @rdname aoi-helpers
#' @noRd
.aoi_detect <- function(df) {
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    cli::cli_abort("No AOI found: call a mutate_aoi_*() function first.")
  }
  params <- attr(df, "traktr_aoi_params")
  id_col <- if (!is.null(params) && !is.null(params$id_col)) {
    params$id_col
  } else {
    setdiff(names(aoi), attr(aoi, "sf_column"))[1]
  }
  if (is.null(id_col) || !id_col %in% names(aoi)) {
    cli::cli_abort("Cannot determine AOI id column; set attr(df, 'traktr_aoi_params')$id_col.")
  }
  if (!id_col %in% names(df)) {
    cli::cli_abort(paste0("Input df does not contain AOI id column '", id_col, "'."))
  }
  list(aoi = aoi, id_col = id_col)
}

#' @description
#' \code{.collect_raster_files()} finds raster files either from a recorded
#' attribute (e.g., \code{"ndvi_aoi"}, \code{"worldcover_tiles"}) or from a
#' directory on disk.
#'
#' @param df A data object carrying tile metadata as attributes.
#' @param attr_name Name of the attribute containing tile metadata.
#' @param out_dir Optional directory path where raster files are stored.
#'
#' @return A character vector of raster file paths.
#'
#' @rdname aoi-helpers
#' @noRd
.collect_raster_files <- function(df, attr_name, out_dir) {
  tiles <- attr(df, attr_name)
  if (is.null(tiles)) {
    if (is.null(out_dir)) {
      cli::cli_abort("No recorded tiles; please supply 'out_dir'.")
    }
    files <- list.files(out_dir, "\\.tif(f)?$", full.names = TRUE)
    if (!length(files)) cli::cli_abort("No raster files found in 'out_dir'.")
  } else {
    files <- tiles$file[!is.na(tiles$file)]
    if (!length(files)) cli::cli_abort(paste0("No files recorded in attr(df, '", attr_name, "')."))
  }
  files
}

#' @description
#' \code{.aoi_deduplicate()} reduces AOIs to unique geometries for efficient
#' extraction, with special handling for H3 cells (unique by ID). It also
#' creates a mapping table back to the full AOI set so results can be joined
#' correctly.
#'
#' @param aoi An \code{sf} object of AOI polygons.
#' @param id_col Column name containing AOI identifiers.
#' @param deduplicate Logical. If \code{TRUE}, geometrically identical AOIs
#'   are deduplicated before extraction.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{aoi_unique}: deduplicated AOIs (with \code{.aoi_row}).
#'   \item \code{map_back}: table mapping original AOIs to \code{.aoi_row}.
#' }
#'
#' @rdname aoi-helpers
#' @noRd
.aoi_deduplicate <- function(aoi, id_col, deduplicate = TRUE) {
  aoi <- sf::st_make_valid(aoi)
  aoi4326 <- sf::st_transform(aoi, 4326)
  
  if (identical(id_col, "h3") && "h3" %in% names(aoi4326)) {
    # H3 cells are intrinsically unique by id
    aoi_unique <- aoi4326[!duplicated(aoi4326$h3), , drop = FALSE]
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    map_back <- sf::st_drop_geometry(aoi4326)[, c(id_col), drop = FALSE]
    map_back$.aoi_row <- match(map_back[[id_col]], aoi_unique[[id_col]])
    return(list(aoi_unique = aoi_unique, map_back = map_back))
  }
  
  if (deduplicate) {
    aoi4326$.__wkb <- vapply(sf::st_as_binary(aoi4326$geometry), paste, "", collapse = "")
    aoi_unique <- aoi4326[!duplicated(aoi4326$.__wkb), , drop = FALSE]
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    
    df1 <- sf::st_drop_geometry(aoi4326)[, c(id_col, ".__wkb"), drop = FALSE]
    df2 <- sf::st_drop_geometry(aoi_unique)[, c(".__wkb", ".aoi_row"), drop = FALSE]
    map_back <- merge(df1, df2, by = ".__wkb", all.x = TRUE, sort = FALSE)
    map_back$.__wkb <- NULL
  } else {
    aoi_unique <- aoi4326
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    map_back <- sf::st_drop_geometry(aoi_unique)[, c(id_col, ".aoi_row"), drop = FALSE]
  }
  
  list(aoi_unique = aoi_unique, map_back = map_back)
}

#' @description
#' \code{.restore_attrs_class()} copies attributes and classes from the
#' original \code{df} to a modified version after joins or transformations.
#'
#' @param df The original input object.
#' @param out The modified output object.
#'
#' @return The output object with attributes and class restored.
#'
#' @rdname aoi-helpers
#' @noRd
.restore_attrs_class <- function(df, out) {
  attrs <- attributes(df)
  attrs <- attrs[!names(attrs) %in% c("names", "row.names", "class")]
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  class(out) <- class(df)
  out
}


# Great-circle distance (meters), Haversine; vectorized, paired points
.haversine_vec <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  R <- 6371000
  dlon <- (lon2 - lon1) * rad
  dlat <- (lat2 - lat1) * rad
  a <- sin(dlat/2)^2 + cos(lat1*rad) * cos(lat2*rad) * sin(dlon/2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

# Convert POSIXct to another tz (base R)
.with_tz <- function(x, tz) {
  # keep the instant, display in tz
  as.POSIXct(format(x, tz = "UTC", usetz = TRUE), tz = tz)
}

# Fast overlap (seconds) of [s,e) with daily window [day_start_hr, day_end_hr) in tz
.fast_day_overlap <- function(s, e, tz = "UTC", day_start_hr = 7, day_end_hr = 22) {
  stopifnot(inherits(s, "POSIXct"), inherits(e, "POSIXct"))
  n <- length(s)
  out <- numeric(n)
  if (!n) return(out)
  
  s_loc <- .with_tz(s, tz)
  e_loc <- .with_tz(e, tz)
  
  valid <- is.finite(s_loc) & is.finite(e_loc) & (e_loc > s_loc)
  if (!any(valid)) return(out)
  
  s_loc <- s_loc[valid]; e_loc <- e_loc[valid]
  
  # floor to local day
  sd <- as.POSIXct(format(s_loc, "%Y-%m-%d 00:00:00", tz = tz), tz = tz)
  ed <- as.POSIXct(format(e_loc, "%Y-%m-%d 00:00:00", tz = tz), tz = tz)
  same_day <- sd == ed
  
  to_num <- function(x) as.numeric(x)
  
  # first day
  f_start <- pmax(to_num(s_loc), to_num(sd + day_start_hr*3600))
  f_end   <- pmin(to_num(e_loc), to_num(sd + day_end_hr*3600))
  first   <- pmax(f_end - f_start, 0)
  
  # last day (if span > 1 day)
  l_start <- pmax(to_num(ed + day_start_hr*3600), to_num(s_loc))
  l_end   <- pmin(to_num(e_loc), to_num(ed + day_end_hr*3600))
  last    <- pmax(l_end - l_start, 0)
  last[same_day] <- 0
  
  # middle full days
  n_days  <- as.integer((to_num(ed) - to_num(sd)) / 86400)
  mid_days <- pmax(n_days - 1L, 0L)
  window_len <- (day_end_hr - day_start_hr) * 3600
  mid <- mid_days * window_len
  
  res <- first + last + mid
  out[valid] <- res
  out
}

.wtd_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

.wtd_sd_unbiased <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  W  <- sum(w); if (W <= 0) return(NA_real_)
  mu <- sum(w * x) / W
  num <- sum(w * (x - mu)^2)
  W2  <- sum(w^2)
  denom <- W - (W2 / W)
  if (denom <= 0) return(NA_real_)
  sqrt(num / denom)
}

.safe_sum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  sum(x)
}

#' Summarize greenness exposure with (auto) time-weighting
#'
#' Computes per-group exposure summaries from GPS-derived intervals, with
#' optional **auto-computation of daytime dwell weights** based on a local
#' daytime window and gap handling rules. Returns **time-weighted means**,
#' **unbiased weighted SDs**, optional **unweighted unique-AOI** statistics,
#' and **coverage** diagnostics.
#'
#' @param data A \code{data.frame} containing at least the grouping variables
#'   (see \code{group}), the exposure value columns (see \code{value_cols}),
#'   and either a weight column (see \code{weight_col}) or the fields required
#'   to auto-compute weights: \code{subject}, \code{timestamp} (POSIXct),
#'   \code{lat}, \code{lon}.
#' @param group Character vector of grouping columns. Defaults to
#'   \code{c("subject","date")}. If \code{"date"} is included but not present
#'   in \code{data}, it is derived from \code{timestamp} in \code{tz_local}
#'   when \code{add_date_if_missing = TRUE}.
#' @param value_cols Named character vector mapping metric names to columns,
#'   e.g. \code{c(ndvi = "ndvi_mean", green = "frac_green_core",
#'   grey = "frac_grey")}. For each name \code{m}, the function produces
#'   \code{m_time_weighted}, \code{m_sd_weighted}, and, if enabled,
#'   \code{m_mean_unweighted}, \code{m_sd_unweighted}.
#' @param weight_col Name of the weight column (seconds), default
#'   \code{"day_dwell_sec"}. If absent and \code{auto_compute_weights = TRUE},
#'   weights are derived (see Details).
#' @param weight_fun Optional function applied to the weights before
#'   summarization; signature \code{function(w) w}. Examples: cap weights
#'   \code{function(w) pmin(w, 300)}, presence-only
#'   \code{function(w) as.numeric(w > 0)}, or \code{function(w) sqrt(w)}.
#' @param aoi_id Optional AOI identifier column (e.g., \code{"h3"} or
#'   \code{"buffer_id"}). If supplied and \code{include_unweighted = TRUE},
#'   unweighted unique-AOI summaries are computed per group.
#' @param include_unweighted Logical; if \code{TRUE} (default), compute
#'   unweighted statistics across unique AOIs per group (one row per AOI).
#' @param coverage_cols Named character vector with names
#'   \code{c(valid = "day_dwell_sec", missing = "missing_day_sec")}. If both
#'   columns exist in the group, coverage metrics are reported (see Details).
#' @param tz_local Time zone used to derive \code{date} from \code{timestamp}
#'   and to evaluate daytime windows. Default \code{"UTC"}.
#' @param add_date_if_missing Logical; if \code{TRUE} (default) and
#'   \code{"date"} is in \code{group} but missing in \code{data}, derive it
#'   from \code{timestamp} in \code{tz_local}.
#' @param auto_compute_weights Logical; if \code{TRUE} (default) and
#'   \code{weight_col} is absent, compute \code{day_dwell_sec} and
#'   \code{missing_day_sec} from successive fixes (see Details).
#' @param nominal_fix_sec Intended sampling interval in seconds; used to detect
#'   gaps. Default \code{10}.
#' @param delta_t_sec Standby-gap threshold in seconds (no movement for at
#'   least this long). Default \code{3*60}.
#' @param delta_d_m Distance threshold in meters below which we consider
#'   "no movement" across a gap. Default \code{150}.
#' @param gap_big_sec Threshold (seconds) for classifying a gap as "big/uncertain".
#'   Default \code{15*60}.
#' @param gap_strategy Strategy for big/uncertain gaps: \code{"exclude"} (mark
#'   as missing) or \code{"carry_forward"} (cap carry-forward dwell). Default
#'   \code{c("exclude","carry_forward")}, first value is used.
#' @param cap_carry_forward Cap in seconds when \code{gap_strategy = "carry_forward"}.
#'   Default \code{5*60}.
#' @param day_start_hr,day_end_hr Local hours defining the daytime window
#'   \code{[start, end)}. Defaults \code{7} and \code{22}.
#' @param progress Logical; if \code{TRUE} (default), display a progress bar
#'   via \pkg{cli}.
#'
#' @details
#' **Auto-computation of weights (when `weight_col` is missing).**  
#' For each subject, consecutive GPS fixes define intervals
#' \eqn{[t_i, t_{i+1})} with duration \eqn{d_i = t_{i+1} - t_i} (seconds) and
#' displacement \eqn{\Delta s_i} (meters; Haversine).
#'
#' Gaps are identified by \eqn{d_i > \mathrm{nominal\_fix\_sec}}.  
#' A gap is *standby* if
#' \eqn{d_i \ge \mathrm{delta\_t\_sec}} and
#' \eqn{\Delta s_i \le \mathrm{delta\_d\_m}}.  
#' A gap is *big/uncertain* if
#' \eqn{\Delta s_i > \mathrm{delta\_d\_m}} or
#' \eqn{d_i \ge \mathrm{gap\_big\_sec}}.
#'
#' Effective dwell for interval \eqn{i} is:
#' \deqn{
#'   \tilde{d}_i =
#'   \begin{cases}
#'     \mathrm{NA}, & \text{if big/uncertain and strategy = exclude} \\
#'     \min(d_i, C), & \text{if big/uncertain and strategy = carry\_forward} \\
#'     d_i, & \text{otherwise,}
#'   \end{cases}
#' }
#' where \eqn{C} is the carry-forward cap (seconds).
#'
#' Let \eqn{O_i} be the overlap (seconds) of \eqn{[t_i, t_{i+1})} with the
#' local daytime window \eqn{[H_{\mathrm{start}}, H_{\mathrm{end}})}
#' (evaluated in `tz_local`). The daytime dwell weight and missing daytime
#' seconds are:
#' \deqn{w_i = \min(\tilde{d}_i, O_i),}
#' \deqn{
#'   m_i =
#'   \begin{cases}
#'     O_i, & \text{if } w_i \text{ is NA and } O_i > 0 \\
#'     0,   & \text{otherwise.}
#'   \end{cases}
#' }
#' The function stores \eqn{w_i} in `day_dwell_sec` (or `weight_col`) and
#' \eqn{m_i} in `missing_day_sec`. Intervals with \eqn{O_i = 0} are dropped.
#' Weights may be transformed via `weight_fun(w)` before use.
#'
#' **Time-weighted summaries.**  
#' For any attribute \eqn{x_i} (e.g., NDVI) with weights \eqn{w_i}, the
#' time-weighted mean is
#' \deqn{\bar{x}_w = \frac{\sum_i w_i x_i}{\sum_i w_i}.}
#'
#' The unbiased weighted standard deviation for frequency-like weights is
#' \deqn{
#'   s_w = \sqrt{
#'     \frac{\sum_i w_i (x_i - \bar{x}_w)^2}{
#'       \sum_i w_i - \frac{\sum_i w_i^2}{\sum_i w_i}
#'     }
#'   }.
#' }
#' If the denominator is non-positive, `NA` is returned.
#'
#' **Coverage.**  
#' If both `coverage_cols$valid` (e.g., `day_dwell_sec`) and
#' `coverage_cols$missing` (e.g., `missing_day_sec`) exist, coverage is
#' \deqn{
#'   \mathrm{coverage} =
#'   \frac{\sum_i \mathrm{valid}_i}{\sum_i \mathrm{valid}_i + \sum_i \mathrm{missing}_i}.
#' }
#' 
#'
#' @section Output columns:
#' For each name \code{m} in \code{value_cols}, the output includes:
#' \itemize{
#'   \item \code{m_time_weighted}: time-weighted mean.
#'   \item \code{m_sd_weighted}: unbiased weighted SD.
#'   \item \code{m_mean_unweighted}, \code{m_sd_unweighted} (if requested):
#'         unweighted mean/SD across unique AOIs in the group.
#' }
#' If coverage inputs are present:
#' \itemize{
#'   \item \code{valid_day_sec}, \code{missing_day_sec}: sums across the group.
#'   \item \code{coverage_valid_share}: \eqn{\in [0,1]}.
#' }
#' Grouping columns (as in \code{group}) are preserved and appear first.
#'
#' @return A \code{data.frame} with one row per group containing the summary
#'   columns described above. Returned \emph{invisibly}.
#'
#' @seealso
#' \code{\link{attach_ndvi}},
#' \code{\link{attach_landuse}},
#' \code{\link{mutate_aoi_h3}},
#' \code{\link{mutate_aoi_buffer}}
#'
#' @examples
#' \dontrun{
#' # 1) Using precomputed weights (day_dwell_sec) and H3 AOIs
#' expo_daily <- gps_data |>
#'   summarize_exposure(
#'     group      = c("subject","date"),
#'     value_cols = c(ndvi="ndvi_mean", green="frac_green_core", grey="frac_grey"),
#'     weight_col = "day_dwell_sec",
#'     aoi_id     = "h3",
#'     tz_local   = "Europe/Zurich"
#'   )
#'
#' # 2) Auto-compute weights from GPS fixes, exclude big/uncertain gaps
#' expo_daily2 <- gps_data |>
#'   summarize_exposure(
#'     group            = c("subject","date"),
#'     aoi_id           = "h3",
#'     auto_compute_weights = TRUE,
#'     gap_strategy     = "exclude",
#'     day_start_hr     = 7, day_end_hr = 22,
#'     tz_local         = "Europe/Zurich"
#'   )
#'
#' # 3) Presence-only weighting and unique-buffer AOI stats
#' expo_buf <- gps_data |>
#'   summarize_exposure(
#'     group      = c("subject","date"),
#'     aoi_id     = "buffer_id",
#'     weight_fun = function(w) as.numeric(w > 0)
#'   )
#' }
#'
#' @export
summarize_exposure <- function(
    data,
    group = c("subject", "date"),
    value_cols = c(ndvi = "ndvi_mean",
                   green = "frac_green_core",
                   grey  = "frac_grey"),
    weight_col = "day_dwell_sec",
    weight_fun = function(w) w,                # optional transform on weights
    aoi_id = NULL,
    include_unweighted = TRUE,
    coverage_cols = c(valid = "day_dwell_sec", missing = "missing_day_sec"),
    tz_local = "UTC",
    add_date_if_missing = TRUE,
    
    # --- auto-computation controls (used if weight_col is absent) ---
    auto_compute_weights = TRUE,
    nominal_fix_sec = 10,
    delta_t_sec     = 3*60,     # standby threshold
    delta_d_m       = 150,
    gap_big_sec     = 15*60,
    gap_strategy    = c("exclude","carry_forward"),
    cap_carry_forward = 5*60,
    day_start_hr = 7,
    day_end_hr   = 22,
    progress = TRUE
) {
  stopifnot(is.data.frame(data))
  gap_strategy <- match.arg(gap_strategy)
  
  # Add date from timestamp if requested and needed
  if (add_date_if_missing && ("date" %in% group) && !("date" %in% names(data))) {
    if (!("timestamp" %in% names(data))) {
      cli::cli_abort("`date` requested in `group` but `timestamp` is missing for derivation.")
    }
    if (!inherits(data$timestamp, "POSIXt")) {
      cli::cli_abort("`timestamp` must be POSIXct/POSIXlt to derive `date`.")
    }
    data$date <- as.Date(format(as.POSIXct(data$timestamp, tz = tz_local), tz = tz_local))
  }
  
  pb <- if (isTRUE(progress)) {
    cli::cli_progress_bar(
      "summarize_exposure()",
      total = 4,
      type = "iterator",
      clear = FALSE,
      status = "summarize_exposure()",
      format = "{cli::pb_spin} {cli::pb_name} | {cli::pb_status} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}  {cli::pb_eta}",
      format_done = "{cli::col_green(cli::symbol$tick)} Summarized exposure {cli::style_dim(paste0('[', cli::pb_elapsed, ']'))}"
    )
  } else NULL
  if (!is.null(pb)) on.exit(cli::cli_progress_done(pb), add = TRUE)
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 1, status = "Calculating time overlaps", force = TRUE)
  
  
  # If weights missing and allowed, derive them
  if (!(weight_col %in% names(data)) && isTRUE(auto_compute_weights)) {
    # Require minimal columns
    req <- c("subject","timestamp","lat","lon")
    miss <- setdiff(req, names(data))
    if (length(miss)) cli::cli_abort("Cannot auto-compute weights; missing: ", paste(miss, collapse = ", "))
    
    # Order by subject,timestamp
    ord <- order(data$subject, data$timestamp)
    data <- data[ord, , drop = FALSE]
    
    # Lead timestamp/lat/lon within subject
    sid <- data$subject
    ts  <- data$timestamp
    lat <- data$lat; lon <- data$lon
    
    # lead indices per subject
    lead_idx <- ave(seq_along(sid), sid, FUN = function(ix) c(ix[-1], NA))
    next_time <- ts[lead_idx]
    next_lat  <- lat[lead_idx]
    next_lon  <- lon[lead_idx]
    
    dt_sec <- as.numeric(difftime(next_time, ts, units = "secs"))
    dist_m <- .haversine_vec(lon, lat, next_lon, next_lat)
    
    is_gap       <- dt_sec > nominal_fix_sec
    standby_gap  <- is_gap & (dt_sec >= delta_t_sec) & (dist_m <= delta_d_m)
    big_uncertain<- is_gap & ((dist_m > delta_d_m) | (dt_sec >= gap_big_sec))
    
    dwell_eff_sec <- dt_sec
    if (gap_strategy == "exclude") {
      dwell_eff_sec[big_uncertain] <- NA_real_
    } else {
      # carry_forward, but cap
      dwell_eff_sec[big_uncertain] <- pmin(dt_sec[big_uncertain], cap_carry_forward)
    }
    
    # Daytime overlap
    daytime_overlap_sec <- .fast_day_overlap(ts, next_time, tz = tz_local,
                                             day_start_hr = day_start_hr, day_end_hr = day_end_hr)
    
    # Final weights and missingness
    day_dwell_sec   <- pmin(daytime_overlap_sec, dwell_eff_sec)
    missing_day_sec <- ifelse(is.na(day_dwell_sec) & (daytime_overlap_sec > 0), daytime_overlap_sec, 0)
    
    data[[weight_col]] <- day_dwell_sec
    data[["missing_day_sec"]] <- missing_day_sec
    
    # Keep only intervals that intersect daytime at all (optional; mirrors your note)
    keep <- (daytime_overlap_sec > 0)
    data <- data[keep, , drop = FALSE]
  }
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 2, status = "Validating columns", force = TRUE)
  
  
  # Validate required columns now
  miss_g <- setdiff(group, names(data))
  if (length(miss_g)) cli::cli_abort("Missing grouping columns: ", paste(miss_g, collapse = ", "))
  miss_v <- setdiff(unname(value_cols), names(data))
  if (length(miss_v)) cli::cli_abort("Missing value columns: ", paste(miss_v, collapse = ", "))
  if (!(weight_col %in% names(data))) cli::cli_abort("Missing weight column '", weight_col, "' (and auto-compute disabled).")
  
  # Prepare grouping key
  key <- if (length(group) == 1L) data[[group]] else do.call(interaction, c(data[group], drop = TRUE))
  split_idx <- split(seq_len(nrow(data)), key)
  
  out_rows <- vector("list", length(split_idx))
  i <- 0L
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 3, status = "Summarizing per group", force = TRUE)
  
  
  for (idx in split_idx) {
    i <- i + 1L
    df <- data[idx, , drop = FALSE]
    
    # Extract group values
    group_vals <- lapply(group, function(g) {
      v <- df[[g]]
      v[which(!is.na(v))[1L]]
    })
    names(group_vals) <- group
    
    # Weights (allow custom transform)
    w <- tryCatch(weight_fun(df[[weight_col]]), error = function(e) {
      cli::cli_abort("`weight_fun` failed: ", conditionMessage(e))
    })
    
    res <- list()
    for (nm in names(value_cols)) {
      col <- value_cols[[nm]]
      x <- df[[col]]
      res[[paste0(nm, "_time_weighted")]] <- .wtd_mean(x, w)
      res[[paste0(nm, "_sd_weighted")]]   <- .wtd_sd_unbiased(x, w)
    }
    
    # Coverage if present
    if (all(unname(coverage_cols) %in% names(df))) {
      vcol <- coverage_cols[["valid"]]
      mcol <- coverage_cols[["missing"]]
      valid_sec   <- .safe_sum(df[[vcol]])
      missing_sec <- .safe_sum(df[[mcol]])
      denom <- valid_sec + missing_sec
      res[["valid_day_sec"]]        <- valid_sec
      res[["missing_day_sec"]]      <- missing_sec
      res[["coverage_valid_share"]] <- if (denom > 0) valid_sec / denom else NA_real_
    }
    
    # Unweighted unique-AOI stats (optional)
    if (include_unweighted && !is.null(aoi_id) && (aoi_id %in% names(df))) {
      keep <- !duplicated(df[[aoi_id]])
      dfu  <- df[keep, , drop = FALSE]
      for (nm in names(value_cols)) {
        col <- value_cols[[nm]]
        x <- dfu[[col]]
        x <- x[is.finite(x)]
        res[[paste0(nm, "_mean_unweighted")]] <- if (length(x)) mean(x) else NA_real_
        res[[paste0(nm, "_sd_unweighted")]]   <- if (length(x) > 1L) stats::sd(x) else NA_real_
      }
    }
    
    out_rows[[i]] <- c(group_vals, res)
  }
  
  if (!is.null(pb)) cli::cli_progress_update(pb, set = 4, status = "Cleaning", force = TRUE)
  
  # Assemble, group columns first
  col_names <- unique(unlist(lapply(out_rows, names), use.names = FALSE))
  col_names <- c(group, setdiff(col_names, group))
  
  out <- do.call(rbind, lapply(out_rows, function(x) {
    row <- setNames(vector("list", length(col_names)), col_names)
    for (nm in names(x)) row[[nm]] <- x[[nm]]
    as.data.frame(row, stringsAsFactors = FALSE)
  }))
  
  # Coerce non-group columns to numeric
  num_cols <- setdiff(col_names, group)
  for (nm in num_cols) {
    if (!is.null(out[[nm]]) && !inherits(out[[nm]], "numeric")) {
      suppressWarnings(out[[nm]] <- as.numeric(out[[nm]]))
    }
  }
  rownames(out) <- NULL
  if (isTRUE(progress)) cli::cli_progress_done()
  
  invisible(out)
}