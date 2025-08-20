#' AOI aus GPS-Punkten erzeugen
#' @param gps data.frame/sf mit Spalten lon, lat (WGS84)
#' @param buffer_point_m numerisch, Punktpuffer in Metern
#' @param margin_m numerisch, BBox-Marge in Metern
#' @param metric_crs EPSG-Code für metrisches CRS (z.B. 2056/3857)
#' @return sfc POLYGON im EPSG:4326
aoi_from_gps <- function(gps, buffer_point_m = 50, margin_m = 500, metric_crs = 2056) {
  stopifnot(all(c("lon","lat") %in% names(gps)))
  pts_ll <- sf::st_as_sf(gps, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pts_m  <- sf::st_transform(pts_ll, metric_crs)
  bufs   <- sf::st_buffer(pts_m, buffer_point_m)
  u      <- sf::st_union(sf::st_make_valid(bufs))
  bbm    <- sf::st_bbox(u)
  bbm[c("xmin","ymin","xmax","ymax")] <- c(bbm["xmin"]-margin_m, bbm["ymin"]-margin_m,
                                           bbm["xmax"]+margin_m, bbm["ymax"]+margin_m)
  aoi_ll <- sf::st_transform(sf::st_as_sfc(bbm, crs = metric_crs), 4326)
  aoi_ll
}