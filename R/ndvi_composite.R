#' Tagesliste & NDVI-Komposit bauen
make_ndvi_composite <- function(aoi, start_date, end_date, oauth, script_file) {
  days <- SearchCatalog(
    aoi = aoi, from = as.character(start_date), to = as.character(end_date),
    collection = "sentinel-2-l2a", with_geometry = TRUE, client = oauth
  ) |>
    dplyr::arrange(tileCloudCover) |>
    dplyr::pull(acquisitionDate) |>
    unique()
  
  if (length(days) == 0) return(NULL)  # <— Bugfix: nicht 'imgs'
  
  fetch_ndvi_day <- function(day) {
    x <- GetImage(
      aoi = aoi, time_range = c(day, day),  # 1-Tages-Fenster
      script = script_file, collection = "sentinel-2-l2a",
      format = "image/tiff", mosaicking_order = "leastCC",
      resolution = 10, mask = TRUE, buffer = 0, client = oauth
    )
    if (inherits(x, "SpatRaster")) return(x)
    terra::rast(x)
  }
  
  rs_list <- purrr::map(days, fetch_ndvi_day, .progress = list(
    type = "iterator", 
    format = "Downloading {cli::pb_bar} {cli::pb_percent}",
    clear = TRUE))
  rs_list <- purrr::keep(rs_list, ~ inherits(.x, "SpatRaster"))
  stopifnot(length(rs_list) > 0)
  
  rs <- terra::rast(rs_list)
  ndvi_med <- terra::app(rs, median, na.rm = TRUE)
  terra::focal(ndvi_med, w = matrix(1,3,3), fun = median, na.rm = TRUE, na.policy = "omit")
}

make_seasonal_ndvi <- function(aoi, year, oauth, script_file) {
  start <- as.Date(paste0(year,"-04-01"))
  end   <- as.Date(paste0(year,"-07-31"))
  make_ndvi_composite(aoi, start, end, oauth, script_file)
}