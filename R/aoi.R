#' Create H3-based AOIs and attach an AOI id column
#'
#' Adds an H3 index column to `df` from latitude/longitude and stores the
#' corresponding H3 polygons as an `sf` object in `attr(df, "ndvi_aoi")`.
#' Also records parameters in `attr(df, "ndvi_aoi_params")` and adds class
#' `"ndvi_df"` to the returned data. Rows with missing/invalid coordinates
#' get `NA` in the id column.
#'
#' @param df A data frame / tibble containing coordinates.
#' @param lat,lon Columns in `df` holding latitude and longitude (tidy-eval).
#'   Coordinates are interpreted as WGS84 (EPSG:4326), decimal degrees.
#' @param res Integer H3 resolution (default `10`).
#' @param .name Character scalar; name of the AOI id column to create
#'   in `df` (default `"h3"`).
#'
#' @details
#' The created attribute `attr(df, "ndvi_aoi")` is an `sf` object in EPSG:4326
#' with one row per unique H3 cell and columns `(.name, geometry)`. No union is
#' performed across cells. The attribute `attr(df, "ndvi_aoi_params")` is a list
#' with elements:
#' \itemize{
#'   \item `method = "h3"`
#'   \item `res` (the H3 resolution)
#'   \item `id_col` (the chosen column name in `df`)
#' }
#'
#' This function emits progress messages via \pkg{cli}. It requires packages
#' \pkg{h3jsr}, \pkg{sf}, \pkg{dplyr}, \pkg{tibble}, and \pkg{rlang}.
#'
#' @return `df` with a new id column named `.name` (character) and the
#' attributes `ndvi_aoi` and `ndvi_aoi_params` set; class `"ndvi_df"` is added.
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(lat = c(47.38, 47.39), lon = c(8.53, 8.54))
#' df2 <- mutate_aoi_h3(df, lat, lon, res = 10)
#' attr(df2, "ndvi_aoi")            # sf with per-cell polygons
#' unique(df2$h3)[1:3]              # id column in df
#' }
#'
#' @seealso [mutate_aoi_buffer()], [fetch_ndvi()], [attach_ndvi()]
#' @family AOI functions
#' @rdname mutate_aoi_h3
#' @export
mutate_aoi_h3 <- function(df, lat = "lat", lon = "lon", res = 10, .name = "h3") {
  rlang::check_installed(c("h3jsr", "sf", "tibble", "dplyr", "rlang"))
  
  cli::cli_progress_step("Creating H3 tiles", spinner = TRUE)
  
  
  lat <- rlang::ensym(lat)
  lon <- rlang::ensym(lon)
  id_sym <- rlang::sym(.name)
  
  # validate coords (avoid failing on NA/Inf)
  latv <- df[[rlang::as_name(lat)]]
  lonv <- df[[rlang::as_name(lon)]]
  ok <- is.finite(latv) & is.finite(lonv)
  if (!any(ok)) cli::cli_abort("No valid coordinates (lat/lon) found.")
  
  # compute H3 for valid rows using a temporary sf of just the coordinates
  coords <- tibble::tibble(.lon = lonv[ok], .lat = latv[ok])
  pts <- sf::st_as_sf(coords, coords = c(".lon", ".lat"), crs = 4326)
  h3_ok <- h3jsr::point_to_cell(sf::st_geometry(pts), res = res)
  
  # attach id column (character) to df, preserving row order; NA where not ok
  out <- rep(NA_character_, nrow(df))
  out[ok] <- h3_ok
  df[[.name]] <- out
  
  # build AOI polygons: one row per unique H3 id
  uniq <- unique(out[!is.na(out)])
  if (!length(uniq)) cli::cli_abort("No H3 cells could be derived from the input.")
  
  poly <- h3jsr::cell_to_polygon(uniq)  # sfc
  aoi  <- sf::st_sf(geometry = poly, crs = 4326)       # sf
  aoi[[.name]] <- uniq                                 # Spalte nachträglich setzen
  aoi  <- dplyr::relocate(aoi, !!id_sym) |> sf::st_make_valid()
  #aoi  <- sf::st_sf(!!id_sym := uniq, geometry = poly, crs = 4326) |>
  #  sf::st_make_valid()
  
  # store AOIs and parameters on df

  attr(df, "traktr_aoi") <- dplyr::select(aoi, dplyr::all_of(c(.name, attr(aoi, "sf_column")   )))
  
  
  #attr(df, "traktr_aoi") <- aoi |> dplyr::select(dplyr::all_of(c("id_sym", "geometry")))
  attr(df, "traktr_aoi_params") <- list(
    method = "h3",
    res = res,
    id_col = .name
  )
  
  class(df) <- union("traktr_df", class(df))
  cli::cli_progress_done()
  
  
  df
}

#' Create buffer-based AOIs (meters) and attach an AOI id column
#'
#' Builds circular AOIs (metric buffers) around each row's point
#' and stores the resulting polygons as an `sf` object in
#' `attr(df, "ndvi_aoi")`. Also records parameters in
#' `attr(df, "ndvi_aoi_params")` and adds class `"ndvi_df"`.
#' Rows with missing/invalid coordinates get `NA` in the id column.
#'
#' @param df A data frame / tibble containing coordinates.
#' @param lat,lon Columns in `df` holding latitude and longitude (tidy-eval).
#'   Coordinates are interpreted as WGS84 (EPSG:4326), decimal degrees.
#' @param dist Numeric buffer radius in **meters** (default `100`).
#' @param .name Character scalar; name of the AOI id column to create
#'   in `df` (default `"aoi_id"`).
#' @param crs_planar Integer or CRS used for buffering in meters
#'   (default `3857`, Web Mercator). Buffers are computed in `crs_planar`
#'   and reprojected back to EPSG:4326 in the output attribute.
#'
#' @details
#' The created attribute `attr(df, "ndvi_aoi")` is an `sf` object in EPSG:4326
#' with one row per input row that had valid coordinates, containing columns
#' `(.name, geometry)`. The attribute `attr(df, "ndvi_aoi_params")` is a list
#' with elements:
#' \itemize{
#'   \item `method = "buffer"`
#'   \item `dist_m` (buffer radius in meters)
#'   \item `crs_planar` (CRS used for buffering)
#'   \item `id_col` (the chosen column name in `df`)
#' }
#'
#' This function emits progress messages via \pkg{cli}. It requires packages
#' \pkg{sf}, \pkg{dplyr}, \pkg{tibble}, and \pkg{rlang}.
#'
#' @return `df` with a new id column named `.name` (character IDs for rows with 
#' valid coordinates) and the attributes `ndvi_aoi` and `ndvi_aoi_params` set.
#' Class `"ndvi_df"` is added.
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(lat = c(47.38, 47.39, NA), lon = c(8.53, 8.54, 8.55))
#' df2 <- mutate_aoi_buffer(df, lat, lon, dist = 150)
#' attr(df2, "ndvi_aoi")            # sf with per-row buffers
#' df2$aoi_id[is.na(df$lat)]        # NA where coords were missing
#' }
#'
#' @seealso [mutate_aoi_h3()], [fetch_ndvi()], [attach_ndvi()]
#' @family AOI functions
#' @rdname mutate_aoi_buffer
#' @export
mutate_aoi_buffer <- function(df, lat = "lat", lon = "lon",
                              dist = 100,
                              .name = "aoi_id",
                              crs_planar = 3857) {
  rlang::check_installed(c("sf","dplyr","tibble","rlang"))
  
  cli::cli_progress_step("Creating buffers", spinner = TRUE)
  
  
  lat <- rlang::ensym(lat)
  lon <- rlang::ensym(lon)
  
  # Valid rows (coords present)
  latv <- df[[rlang::as_name(lat)]]
  lonv <- df[[rlang::as_name(lon)]]
  ok <- is.finite(latv) & is.finite(lonv)
  if (!any(ok)) cli::cli_abort("No valid coordinates (lat/lon) found.")
  
  # Keep a stable row id for joining back
  df$.row_id__ <- seq_len(nrow(df))
  
  # Points (only valid rows), buffer in metric CRS
  pts <- sf::st_as_sf(df[ok, , drop = FALSE],
                      coords = c(rlang::as_name(lon), rlang::as_name(lat)),
                      crs = 4326, remove = FALSE) |>
    sf::st_transform(crs_planar)
  
  # One AOI per row (only where coords are valid)
  bufs <- sf::st_buffer(pts, dist) |> sf::st_make_valid()
  #aoi <- bufs |>
  #  dplyr::transmute(.row_id__ = .data$.row_id__) |>
  #  sf::st_transform(4326)
  
  aoi <- bufs[, c(".row_id__", attr(bufs, "sf_column")), drop = FALSE] |>
    sf::st_transform(4326)
  
  # Assign sequential AOI IDs following the valid-row order
  # aoi[[.name]] <- seq_len(nrow(aoi))
  aoi[[.name]] <- sprintf("aoi_%03d", seq_len(nrow(aoi)))
  #aoi <- aoi |> dplyr::select(dplyr::all_of(.name), .row_id__, geometry)
  #print(colnames(aoi))
  
  aoi <- aoi[ , c(.name, ".row_id__", attr(bufs, "sf_column")), drop = FALSE]
  
  # Map AOI ID back to df via row id (rows with NA coords get NA aoi_id)
  df <- df |>
    dplyr::left_join(sf::st_drop_geometry(aoi), by = ".row_id__")
  
  
  # Clean helper column
  df$.row_id__ <- NULL
  
  # Store AOIs (only the necessary columns) as attribute
  aoi_attr <- aoi[ , c(.name, attr(bufs, "sf_column")), drop = FALSE]
  
  #aoi |> dplyr::select(dplyr::all_of(c(.name, "geometry")))
  
  attr(df, "traktr_aoi") <- aoi_attr
  attr(df, "traktr_aoi_params") <- list(
    method = "buffer",
    dist_m = dist,
    crs_planar = crs_planar,
    id_col = .name
  )
  
  # Mark class for downstream verbs
  class(df) <- union("traktr_df", class(df))
  
  cli::cli_progress_done()
  
  invisible(df)
}
