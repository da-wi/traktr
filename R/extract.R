#' Mittel-NDVI in Pufferflächen
mean_ndvi_in_buffers <- function(r, pts_ll, radius_m = 300) {
  pts <- sf::st_as_sf(pts_ll, coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
    sf::st_transform(terra::crs(r))
  buf <- sf::st_buffer(pts, radius_m)
  d  <- exactextractr::exact_extract(r, buf, "mean")
  sf::st_drop_geometry(buf) |>
    dplyr::mutate(ndvi_mean = d)
}