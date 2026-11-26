# Run from the package root to regenerate the synthetic test boundaries.
rectangle <- function(xmin, xmax, ymin, ymax) {
  matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
         ncol = 2L, byrow = TRUE)
}

polygon <- function(xmin, xmax, ymin, ymax) {
  sf::st_multipolygon(list(list(rectangle(xmin, xmax, ymin, ymax))))
}

geoid <- c("08031", "06071", "02016", "02185", "15001", "72127")
name <- c("Denver", "San Bernardino", "Aleutians West", "North Slope",
          "Hawaii", "San Juan")
lsad <- c("06", "06", "07", "04", "06", "13")
counties <- sf::st_sf(
  STATEFP = substr(geoid, 1L, 2L),
  COUNTYFP = substr(geoid, 3L, 5L),
  COUNTYNS = sprintf("%08d", seq_along(geoid)),
  AFFGEOID = paste0("0500000US", geoid),
  GEOID = geoid,
  NAME = name,
  NAMELSAD = paste(name, c("County", "County", "Census Area", "Borough",
                          "County", "Municipio")),
  STUSPS = c("CO", "CA", "AK", "AK", "HI", "PR"),
  STATE_NAME = c("Colorado", "California", "Alaska", "Alaska", "Hawaii",
                 "Puerto Rico"),
  LSAD = lsad,
  ALAND = rep(1000000, length(geoid)),
  AWATER = rep(1000, length(geoid)),
  geometry = sf::st_sfc(
    polygon(-105.1, -104.9, 39.6, 39.8),
    polygon(-117.9, -114.1, 33.9, 35.8),
    sf::st_multipolygon(list(
      list(rectangle(172, 179.5, 51.8, 52.2)),
      list(rectangle(-179.5, -176, 51.8, 52.2))
    )),
    polygon(-162, -141, 68, 71),
    polygon(-155.6, -155.4, 19.5, 19.7),
    polygon(-66.2, -66.0, 18.3, 18.5),
    crs = 4269
  )
)

sf::st_write(
  counties, "tests/testthat/fixtures/counties5.gpkg", layer = "counties",
  delete_dsn = TRUE, quiet = TRUE, layer_options = "SPATIAL_INDEX=NO",
  config_options = c(OGR_SQLITE_PRAGMA = "page_size=1024")
)
