source("R/ndvi_composite.R")
library(CDSE)

aoi <- sf::st_read("output/aoi.geojson", quiet = TRUE)$geom

id <- "sh-09a6a314-1c14-4612-be75-6ea7354309eb"
secret <- "wKZwZZUtWzYCfwhgYYxz0aTBLeEvZPdL"
oauth <- GetOAuthClient(id = id, secret = secret)

script_file <- system.file("scripts", "NDVI_float32.js", package = "CDSE")
r <- make_seasonal_ndvi(aoi, 2025, oauth, script_file)
terra::writeRaster(r, "output/ndvi_2025_AprJul.tif", overwrite = TRUE)
