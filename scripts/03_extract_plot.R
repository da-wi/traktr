source("R/extract.R")
source("R/plot.R")


r   <- terra::rast("output/ndvi_2025_AprJul.tif")
gps <- traj_sf |> dplyr::filter(id == "P01")

res <- mean_ndvi_in_buffers(r, gps, radius_m = 300)

readr::write_csv(res, "output/ndvi_home_summary.csv")
p <- plot_ndvi_points(r, res)
ggplot2::ggsave("output/ndvi_map.png", p, width=8, height=6, dpi=300)

p
