library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(lubridate)

gps <- traj_sf |> filter(id == "P01")
# --- Daten ---
# gps: id, ts (UTC), lon, lat
# ndvi_month.tif: monthly, masked, same CRS
crs_m <- "EPSG:2056"  # LV95

pts <- st_as_sf(gps, coords = c("lon","lat"), crs = 4326) |>
  st_transform(crs_m) |>
  mutate(ts = with_tz(time, "Europe/Zurich"),
         month = floor_date(ts, "month"))
#####################################################################
id <- "sh-09a6a314-1c14-4612-be75-6ea7354309eb"
secret <- "wKZwZZUtWzYCfwhgYYxz0aTBLeEvZPdL"
oauth <- GetOAuthClient(id = id, secret = secret)

aoi <- zurich_osm$bbox |> st_as_sfc() |> st_transform(4326)

images <- SearchCatalog(aoi = aoi,
                        from = "2022-01-01",
                        to = "2025-07-07", 
                        collection = "sentinel-2-l2a",
                        with_geometry = TRUE, 
                        client = OAuthClient) 

# Pick the least-cloudy day
day <- images |> 
  arrange(tileCloudCover) |> 
  head(1) |> 
  pull(acquisitionDate)


#script_file <- system.file("scripts", "TrueColor.js", package = "CDSE")
# day <- "2025-04-04"
script_file <- system.file("scripts", "NDVI_float32.js", package = "CDSE")

fp <- GetImage(
  aoi = aoi, time_range = day, script = script_file,
  collection = "sentinel-2-l2a", format = "image/tiff",
  mosaicking_order = "leastCC", resolution = 10, mask = TRUE,
  buffer = 0, client = oauth
)

ndvi <- fp |> project(crs_m)
# ndvi <- rast("ndvi_2024-06.tif") |> project(crs_m)

# Optional glätten (3x3 Median)
ndvi_smooth <- focal(ndvi, w = matrix(1,3,3), fun = median, na.policy = "omit", na.rm = TRUE)

pts_r <- st_transform(pts, crs(ndvi))      # transform to raster CRS (meter-based)
polys <- st_buffer(pts_r, 50)              # 50 m buffer

ndvi_mean <- exactextractr::exact_extract(ndvi_smooth, polys, "mean")


# Optional: green-share (NDVI > 0.4), area-weighted
green_share <- exactextractr::exact_extract(
  ndvi_smooth, polys,
  function(values, coverage_fraction) {
    ok <- is.finite(values) & coverage_fraction > 0
    if (!any(ok)) return(NA_real_)
    sum(coverage_fraction[ok] * (values[ok] > 0.4)) / sum(coverage_fraction[ok])
  }
)

# Bind back to your (non-geometry) points table
out <- st_drop_geometry(pts) |>
  mutate(ndvi_mean = ndvi_mean,
         green_share = green_share)


# Zusammenführen
out <- bind_cols(st_drop_geometry(pts), data.frame(percell=percell)) |>
  mutate(is_valid = percell >= 0.2)
