plot_ndvi_points <- function(r, df_ll) {
  pts <- sf::st_as_sf(df_ll, coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
    sf::st_transform(terra::crs(r))
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_sf(data = pts, ggplot2::aes(color = ndvi_mean)) +
    ggplot2::scale_color_viridis_c()
}