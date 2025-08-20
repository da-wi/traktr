library(cli)
make_ndvi_composite <- function(aoi, start_date, end_date, oauth, script_file) {
  # 1) list acquisitions in window
  days <- SearchCatalog(
    aoi = aoi,
    from = as.character(start_date),
    to   = as.character(end_date),
    collection = "sentinel-2-l2a",
    with_geometry = TRUE,
    client = oauth
  ) %>%
    arrange(tileCloudCover) %>%
    pull(acquisitionDate) %>%
    unique()
  
  if (length(imgs) == 0) return(NULL)
  
  # 2) fetch NDVI tif for each acquisition day (cloud-masked by your NDVI script)
  fetch_ndvi_day <- function(day, aoi, oauth, script_file) {
    x <- GetImage(
      aoi = aoi,
      time_range = day,
      script = script_file,
      collection = "sentinel-2-l2a",
      format = "image/tiff",   # could still return SpatRaster in this client
      mosaicking_order = "leastCC",
      resolution = 10, mask = TRUE,
      buffer = 0, client = oauth
    )
    if (inherits(x, "SpatRaster")) return(x)
    terra::rast(x)  # if a file path is returned instead
  }
  
  # Fetch rasters (list of SpatRaster)
  rs_list <- map(days, ~ fetch_ndvi_day(.x, aoi, oauth, script_file), .progress = list(
    type = "iterator", 
    format = "Downloading {cli::pb_bar} {cli::pb_percent}",
    clear = TRUE))
  
  # Optional: drop NULL or failed downloads
  rs_list <- keep(rs_list, ~ inherits(.x, "SpatRaster"))
  
  # Stack (assumes same grid since same AOI/resolution)
  rs <- rast(rs_list)
  
  # Pixel-wise median across dates (robust to residual clouds)
  ndvi_med <- app(rs, median, na.rm = TRUE)
  
  
  # Optional: gentle 3×3 median filter
  ndvi_med_smooth <- focal(ndvi_med, w = matrix(1,3,3), fun = median, na.rm = TRUE, na.policy = "omit")
}

make_seasonal_ndvi <- function(aoi, year, oauth, script_file) {
  start <- as.Date(paste0(year,"-04-01"))
  end   <- as.Date(paste0(year,"-07-31"))
  make_ndvi_composite(aoi, start, end, oauth, script_file)
}

aoi_from_gps <- function(gps, buffer_point_m = 50, margin_m = 500, metric_crs = 2056) {
  stopifnot(all(c("lon","lat") %in% names(gps)))
  pts_ll <- st_as_sf(gps, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pts_m  <- st_transform(pts_ll, metric_crs)
  
  # per-point buffer representing your exposure radius
  bufs   <- st_buffer(pts_m, buffer_point_m)
  
  # dissolve to one polygon, then get bbox, then expand by margin
  u      <- st_union(st_make_valid(bufs))
  bbm    <- st_bbox(u)
  
  bbm_expanded <- bbm
  bbm_expanded["xmin"] <- bbm["xmin"] - margin_m
  bbm_expanded["ymin"] <- bbm["ymin"] - margin_m
  bbm_expanded["xmax"] <- bbm["xmax"] + margin_m
  bbm_expanded["ymax"] <- bbm["ymax"] + margin_m
  
  aoi_m  <- st_as_sfc(bbm_expanded, crs = metric_crs)
  aoi_ll <- st_transform(aoi_m, 4326)  # CDSE/Sentinel APIs expect lon/lat
  aoi_ll
}

aoi <- aoi_from_gps(gps %>% filter(id == "P01"))

ndvi_season <- make_seasonal_ndvi(aoi, 2025, oauth,
                                  system.file("scripts", "NDVI_float32.js", package = "CDSE"))

plot(ndvi_season)


gps <- traj_sf |> filter(id == "P01")

pts <- st_as_sf(gps, coords = c("lon","lat"), crs = 4326) |>
  st_transform(crs(ndvi_season))        # match raster CRS
home_buf <- st_buffer(pts, 300)          # 300 m home buffer

ndvi_fin <- home_buf |> mutate(ndvi_home = exact_extract(ndvi_season, home_buf, "mean"))
ggplot(ndvi_fin)+
  tidyterra::geom_spatraster(data=ndvi_season)+
  ggnewscale::new_scale_color()+
  geom_point(aes(x = lon, y=lat,color=ndvi_home))+
  scale_color_viridis_c()


ndvi_fin