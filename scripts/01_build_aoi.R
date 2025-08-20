source("R/aoi.R")
library(readr)

gps <- traj_sf |> dplyr::filter(id == "P01") #readr::read_csv("data/traj.csv") 

aoi <- aoi_from_gps(gps)

sf::st_write(sf::st_sf(geometry = aoi), "output/aoi.geojson", delete_dsn = TRUE)
