#' Coordinate validation and extraction (internal)
#'
#' Validates latitude/longitude columns and extracts only the rows with
#' finite coordinates for downstream AOI computation.
#'
#' @param df A data frame (or tibble) containing coordinate columns.
#' @param lat,lon Column names (strings or symbols) pointing to latitude and
#'   longitude in WGS84 degrees.
#'
#' @details
#' - Columns are captured with \code{rlang::ensym()} and looked up in \code{df}.
#' - Rows with non-finite coordinates are marked as invalid.
#' - If no rows are valid, the function aborts with a clear message.
#' - Returned coordinates use column names \code{.lon} and \code{.lat}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{ok}: logical vector (length \code{nrow(df)}) indicating valid rows.
#'   \item \code{coords}: tibble with \code{.lon}, \code{.lat} for valid rows.
#'   \item \code{lat_sym}, \code{lon_sym}: captured symbols for the input columns.
#' }
#'
#' @seealso \code{\link{mutate_aoi_h3}}, \code{\link{mutate_aoi_buffer}}
#' @keywords internal
#' @rdname aoi-internals
#' @noRd
.aoi_validate_coords <- function(df, lat, lon) {
  lat_sym <- rlang::ensym(lat)
  lon_sym <- rlang::ensym(lon)
  
  latv <- df[[rlang::as_name(lat_sym)]]
  lonv <- df[[rlang::as_name(lon_sym)]]
  
  ok <- is.finite(latv) & is.finite(lonv)
  if (!any(ok)) {
    cli::cli_abort("No valid coordinates (lat/lon) found.")
  }
  
  coords <- tibble::tibble(.lon = lonv[ok], .lat = latv[ok])
  list(ok = ok, coords = coords, lat_sym = lat_sym, lon_sym = lon_sym)
}


#' Move AOI id column to the front (internal)
#'
#' Reorders columns so the AOI identifier appears first, preserving the
#' geometry column and CRS of the input \code{sf} object.
#'
#' @param aoi An \code{sf} object.
#' @param id_col Name of the identifier column to move to the front.
#'
#' @details
#' - Only column order is changed; geometry column and data are preserved.
#' - Used to produce a consistent \code{sf} layout for attributes and outputs.
#'
#' @return The input \code{aoi} with columns reordered.
#'
#' @seealso \code{\link{mutate_aoi_h3}}, \code{\link{mutate_aoi_buffer}}
#' @keywords internal
#' @rdname aoi-internals
#' @noRd
.aoi_relocate_id_first <- function(aoi, id_col) {
  cols <- names(aoi)
  aoi[, c(id_col, setdiff(cols, id_col)), drop = FALSE]
}

#' Attach AOI attributes and id column (internal)
#'
#' Writes an AOI identifier column back to \code{df}, stores a minimal AOI
#' \code{sf} (id + geometry) in \code{attr(df, "traktr_aoi")}, records
#' parameters in \code{attr(df, "traktr_aoi_params")}, and tags the object
#' with class \code{"traktr_df"}.
#'
#' @param df A data frame (or tibble) to receive the id column and attributes.
#' @param id_vec Vector of AOI identifiers aligned to \code{df} rows.
#' @param aoi An \code{sf} object containing AOI geometries (already valid and
#'   in EPSG:4326) with an \code{id_col} column.
#' @param id_col Name of the identifier column to write and keep in attributes.
#' @param params Named list of method-specific parameters to record.
#'
#' @details
#' - Overwrites \code{df[[id_col]]} if it already exists.
#' - The AOI attribute keeps only \code{id_col} and the geometry column (id first).
#' - The \code{"traktr_df"} class is prepended to \code{class(df)}.
#' - Callers are responsible for validating \code{aoi} (e.g., \code{st_make_valid()}).
#'
#' @return \code{df}, augmented with the id column and attributes.
#'
#' @seealso \code{\link{mutate_aoi_h3}}, \code{\link{mutate_aoi_buffer}}
#' @keywords internal
#' @rdname aoi-internals
#' @noRd
.aoi_attach_attrs <- function(df, id_vec, aoi, id_col, params) {
  # Write id column back to df (row order preserved)
  df[[id_col]] <- id_vec
  
  # Keep only id + geometry on the AOI attribute, id first
  geom_col <- attr(aoi, "sf_column")
  aoi_attr <- aoi[, c(id_col, geom_col), drop = FALSE]
  
  attr(df, "traktr_aoi") <- aoi_attr
  attr(df, "traktr_aoi_params") <- c(params, list(id_col = id_col))
  
  # Mark class for downstream verbs
  class(df) <- union("traktr_df", class(df))
  df
}

#' Derive AOI polygons from point coordinates using H3
#'
#' Computes H3 cell IDs from latitude/longitude columns and builds one AOI
#' polygon per unique cell. The H3 ID is written back to the input data as a
#' new column, and an AOI \code{sf} is attached as an attribute for downstream
#' fetching functions.
#'
#' @param df A data frame (or tibble) with point coordinates.
#' @param lat,lon Names (strings or symbols) of the latitude and longitude
#'   columns in \code{df}. Must be in WGS84 (EPSG:4326).
#' @param res H3 resolution (integer). Higher values yield smaller cells.
#' @param .name Name of the output column in \code{df} that stores the H3 ID.
#'   If a column of the same name exists, it will be overwritten.
#'
#' @details
#' - Coordinates are validated; rows with non-finite \code{lat}/\code{lon} get
#'   \code{NA} in the ID column and do not contribute polygons.
#' - The function preserves row order and does not reshuffle \code{df}.
#' - One AOI polygon is produced per **unique** H3 ID.
#' - The AOI is created as an \code{sf} object in WGS84 and geometries are
#'   run through \code{sf::st_make_valid()}.
#' - A standard set of attributes is written to \code{df} for sibling functions.
#'
#' @section Attributes:
#' - \code{attr(df, "traktr_aoi")}: an \code{sf} with two columns:
#'   \code{.name} (H3 ID) and \code{geometry}.
#' - \code{attr(df, "traktr_aoi_params")}: a list with
#'   \code{list(method = "h3", res = res, id_col = .name)}.
#' - The input is also marked with class \code{"traktr_df"} for downstream verbs.
#'
#' @return The input \code{df}, with a new \code{.name} column (character H3 ID)
#'   and AOI-related attributes (see Attributes).
#'
#' @seealso
#' \code{\link{mutate_aoi_buffer}},
#' \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(
#'   id  = 1:3,
#'   lon = c(8.54, 8.55, NA_real_),
#'   lat = c(47.37, 47.38, 47.39)
#' )
#'
#' df <- mutate_aoi_h3(df, lat = lat, lon = lon, res = 10, .name = "h3")
#'
#' # Inspect H3 column and AOI polygons
#' df$h3
#' attr(df, "traktr_aoi")
#'
#' # Use AOI with a fetching sibling afterwards
#' # df <- fetch_ndvi(df, out_dir = "data/ndvi", from = Sys.Date()-30, to = Sys.Date())
#' }
#'
#' @export
mutate_aoi_h3 <- function(df, lat = "lat", lon = "lon", res = 10, .name = "h3") {
  rlang::check_installed(c("h3jsr", "sf", "tibble", "rlang", "cli"))
  
  cli::cli_progress_step("Creating H3 tiles", spinner = TRUE)
  
  prep <- .aoi_validate_coords(df, lat, lon)
  
  # Compute H3 cell ids for valid rows
  pts <- sf::st_as_sf(prep$coords, coords = c(".lon", ".lat"), crs = 4326)
  h3_ok <- h3jsr::point_to_cell(sf::st_geometry(pts), res = res)
  
  # Map back to full vector
  id_vec <- rep(NA_character_, nrow(df))
  id_vec[prep$ok] <- h3_ok
  
  # Build AOI polygons (one row per unique H3 id)
  uniq_ids <- unique(stats::na.omit(id_vec))
  if (!length(uniq_ids)) {
    cli::cli_abort("No H3 cells could be derived from the input.")
  }
  
  poly <- h3jsr::cell_to_polygon(uniq_ids)       # sfc
  aoi  <- sf::st_sf(geometry = poly, crs = 4326) # sf
  aoi[[.name]] <- uniq_ids
  aoi <- .aoi_relocate_id_first(aoi, .name)
  aoi <- sf::st_make_valid(aoi)
  
  out <- .aoi_attach_attrs(
    df = df,
    id_vec = id_vec,
    aoi = aoi,
    id_col = .name,
    params = list(method = "h3", res = res)
  )
  
  cli::cli_progress_done()
  out
}


#' Derive AOI polygons by buffering point coordinates
#'
#' Buffers each valid point row by a given metric distance (after projecting to
#' a planar CRS) and writes a sequential AOI ID back to the data. The per-row
#' buffer polygons are attached as an AOI \code{sf} for downstream fetching
#' functions.
#'
#' @param df A data frame (or tibble) with point coordinates.
#' @param lat,lon Names (strings or symbols) of the latitude and longitude
#'   columns in \code{df}. Must be in WGS84 (EPSG:4326).
#' @param dist Buffer radius in meters (applied in a planar CRS).
#' @param .name Name of the output column in \code{df} that stores the AOI ID
#'   (e.g., \code{"aoi_id"}). If a column of the same name exists, it will be
#'   overwritten.
#' @param crs_planar EPSG code for the planar (projected) CRS used during
#'   buffering. Defaults to Web Mercator (EPSG:3857).
#'
#' @details
#' - Coordinates are validated; rows with non-finite \code{lat}/\code{lon} get
#'   \code{NA} in the ID column and no polygon is produced for them.
#' - The function preserves row order and does not reshuffle \code{df}.
#' - One AOI polygon is produced **per valid input row** (not deduplicated).
#' - Buffering is performed in \code{crs_planar}, then geometries are
#'   transformed back to WGS84 and run through \code{sf::st_make_valid()}.
#' - AOI IDs are sequential strings (e.g., \code{"aoi_001"}).
#' - A standard set of attributes is written to \code{df} for sibling functions.
#'
#' @section Attributes:
#' - \code{attr(df, "traktr_aoi")}: an \code{sf} with two columns:
#'   \code{.name} (sequential AOI ID) and \code{geometry}.
#' - \code{attr(df, "traktr_aoi_params")}: a list with
#'   \code{list(method = "buffer", dist_m = dist, crs_planar = crs_planar, id_col = .name)}.
#' - The input is also marked with class \code{"traktr_df"} for downstream verbs.
#'
#' @return The input \code{df}, with a new \code{.name} column (character AOI ID)
#'   and AOI-related attributes (see Attributes).
#'
#' @seealso
#' \code{\link{mutate_aoi_h3}},
#' \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(
#'   id  = 1:3,
#'   lon = c(8.54, 8.55, NA_real_),
#'   lat = c(47.37, 47.38, 47.39)
#' )
#'
#' df <- mutate_aoi_buffer(
#'   df,
#'   lat = lat, lon = lon,
#'   dist = 100,
#'   .name = "aoi_id",
#'   crs_planar = 3857
#' )
#'
#' # Inspect AOI IDs and polygons
#' df$aoi_id
#' attr(df, "traktr_aoi")
#' }
#'
#' @export
mutate_aoi_buffer <- function(df,
                              lat = "lat",
                              lon = "lon",
                              dist = 100,
                              .name = "aoi_id",
                              crs_planar = 3857) {
  rlang::check_installed(c("sf", "tibble", "rlang", "cli"))
  
  cli::cli_progress_step("Creating buffers", spinner = TRUE)
  
  prep <- .aoi_validate_coords(df, lat, lon)
  
  # Points for valid rows, buffer in planar, back to WGS84
  pts_planar <- sf::st_as_sf(prep$coords, coords = c(".lon", ".lat"), crs = 4326) |>
    sf::st_transform(crs_planar)
  
  bufs <- sf::st_buffer(pts_planar, dist) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  
  # Sequential AOI ids only for valid rows; map back to all rows
  ids_ok <- sprintf("aoi_%03d", seq_len(nrow(bufs)))
  id_vec <- rep(NA_character_, nrow(df))
  id_vec[prep$ok] <- ids_ok
  
  # AOI sf: one polygon per valid row, id first
  aoi <- bufs
  aoi[[.name]] <- ids_ok
  aoi <- .aoi_relocate_id_first(aoi, .name)
  
  out <- .aoi_attach_attrs(
    df = df,
    id_vec = id_vec,
    aoi = aoi,
    id_col = .name,
    params = list(method = "buffer", dist_m = dist, crs_planar = crs_planar)
  )
  
  cli::cli_progress_done()
  out
}
