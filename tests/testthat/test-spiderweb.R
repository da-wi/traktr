test_that("ecolo-zip spiderwebs use the released spacing and coordinates", {
  spec <- .spiderweb_ecolozip()
  points <- .spiderweb_points(c(40, 30), c(-105, -95), spec)
  expect_equal(as.integer(table(points$id)), c(140L, 140L))
  expect_true(all(table(points$id, points$ring) == 35L))
  expect_false(anyDuplicated(points[c("id", "lat", "lon")]) > 0L)

  lat <- c(40, 30)[points$id]
  lon <- c(-105, -95)[points$id]
  radius <- 111 * sqrt((points$lat - lat)^2 +
                         ((points$lon - lon) * cos(lat / 57.3))^2)
  expect_equal(radius, spec$radii_km[points$ring], tolerance = 1e-10)
  expect_equal(points$lat[1], 40)
  expect_equal(points$lon[1], -105 + 8.0467 / 111 / cos(40 / 57.3))

  quarter_turns <- .spiderweb_ecolozip(8.0467, n_points = 5)
  batch <- .spiderweb_points(c(40, 30), c(-105, -95), quarter_turns, lat_ref = 35)
  expect_equal(batch$lat[2], 40 + 8.0467 / 111)
  expect_equal(batch$lon[2], -105)
  expect_equal(batch$lon[batch$k == 1], c(-105, -95) + 8.0467 / 111 / cos(35 / 57.3))
  per_point <- .spiderweb_points(c(40, 30), c(-105, -95), quarter_turns, lat_ref = c(35, 25))
  expect_equal(per_point$lon[per_point$k == 1],
               c(-105, -95) + 8.0467 / 111 / cos(c(35, 25) / 57.3))
  closed <- .spiderweb_points(40, -105, .spiderweb_ecolozip(drop_duplicate = FALSE))
  expect_equal(nrow(closed), 144L)
  expect_equal(closed[closed$k == 1, "lon"], closed[closed$k == 36, "lon"])
})

write_spiderweb_tile <- function(tile_dir, lat_n, lon_w, elev) {
  bounds <- .ag100_grid_bounds(lat_n, lon_w, c(0.5, 0.5))
  tile <- terra::rast(nrows = 3, ncols = 3, nlyrs = length(.ag100_band_paths),
                      xmin = bounds[1], xmax = bounds[2], ymin = bounds[3], ymax = bounds[4],
                      crs = "EPSG:4326")
  terra::values(tile) <- cbind(rep_len(elev, 9), matrix(1, 9, length(.ag100_band_paths) - 1))
  names(tile) <- names(.ag100_band_paths)
  path <- file.path(tile_dir, sprintf("ag100_v003_%d_%d.tif", lat_n, lon_w))
  terra::writeRaster(tile, path)
  path
}

test_that("AG100 points read nearest cells and flag absent tiles", {
  tile_dir <- withr::local_tempdir()
  write_spiderweb_tile(tile_dir, 40, -105, c(11:18, NA))
  write_spiderweb_tile(tile_dir, 41, -105, 21:29)
  points <- data.frame(id = 1:4, lat = c(40.01, 39.62, 39.8, 39.01),
                       lon = c(-104.5, -104.8, -103.999, -104.01))

  sampled <- .aster_ged_sample_points(points, tile_dir)
  expect_identical(sampled[names(points)], points)
  expect_equal(sampled$elev, c(28, 14, NA, NA))
  expect_equal(sampled$ndvi_mean, c(1, 1, NA, 1))
  expect_identical(sampled$tile_missing, c(FALSE, FALSE, TRUE, FALSE))
  expect_true(all(is.na(sampled[3, names(.ag100_band_paths)])))
})

test_that("latitude and longitude seams use asterlookup tile ownership", {
  tile_dir <- withr::local_tempdir()
  south <- write_spiderweb_tile(tile_dir, 40, -105, 10)
  write_spiderweb_tile(tile_dir, 41, -105, 20)
  east <- write_spiderweb_tile(tile_dir, 40, -104, 30)
  points <- data.frame(lat = c(40, 39.5, 40.0001, 39.5, 40),
                       lon = c(-104.5, -104, -104.5, -104.0001, -104))
  expect_identical(.spiderweb_tile_id(points$lat, points$lon),
                   c("40.-105", "40.-104", "41.-105", "40.-105", "40.-104"))
  expect_equal(.aster_ged_sample_points(points, tile_dir, "elev")$elev, c(10, 30, 20, 10, 30))

  unlink(east)
  seam <- .aster_ged_sample_points(points[2:4, ], tile_dir, "elev")
  expect_equal(seam$elev, c(NA, 20, 10))
  expect_identical(seam$tile_missing, c(TRUE, FALSE, FALSE))

  unlink(south)
  seam <- .aster_ged_sample_points(points[c(1, 3), ], tile_dir, "elev")
  expect_equal(seam$elev, c(NA, 20))
  expect_identical(seam$tile_missing, c(TRUE, FALSE))
})

test_that("spiderweb longitudes wrap and sample both sides of the antimeridian", {
  tile_dir <- withr::local_tempdir()
  write_spiderweb_tile(tile_dir, 41, 179, 10)
  write_spiderweb_tile(tile_dir, 41, -180, 20)
  spec <- .spiderweb_ecolozip(32.1869, n_points = 5)
  points <- .spiderweb_points(40.5, 179.9, spec)
  expect_true(all(points$lon >= -180 & points$lon < 180))
  expect_equal(points$lon[1], 179.9 + 32.1869 / 111 / cos(40.5 / 57.3) - 360)
  expect_equal(.aster_ged_sample_points(points, tile_dir, "elev")$elev, c(20, 10, 10, 10))
  expect_equal(.aster_ged_sample_points(
    .spiderweb_points(40.5, -179.9, spec), tile_dir, "elev"
  )$elev, c(20, 20, 10, 20))
})

test_that("SPSS statistics use valid samples and the published column names", {
  values <- data.frame(
    elev = c(0, 9, 2, NA, NA, NA),
    ndvi_mean = c(NA, 2, 8, 4, 6, NA),
    nobs = c(2, 4, 4, 4, NA, NA),
    lwm = c(1, 0, 0, 1, NA, NA)
  )
  id <- c("a", "b", "a", "b", "a", "c")
  result <- .spss_aggregate(values, id)
  dict <- .aster_ged_dictionary()
  expect_identical(names(result), c("id", dict$name_ecolozip[!dict$extra & dict$band %in% names(values)]))
  expect_identical(result$id, c("a", "b", "c"))
  expect_equal(result$ASTER_zip_elevation_mean, c(1, 9, NA))
  expect_equal(result$ASTER_zip_elevation_sd, c(sqrt(2), NA, NA))
  expect_equal(result$ASTER_zip_elevation_min, c(0, 9, NA))
  expect_equal(result$ASTER_zip_elevation_max, c(2, 9, NA))
  expect_equal(result$ASTER_zip_elevation_cin, c(100 / 3, 0, NA))
  expect_equal(result$ASTER_zip_elevation_sum, c(2, 9, NA))
  expect_equal(result$ASTER_zip_elevation_missing_perc, c(100 / 3, 50, 100))
  expect_equal(result$ASTER_ndvi_mean_missing_perc, c(100 / 3, NA, 100))
  expect_identical(result$ASTER_observation_number_lt_4, c(1L, 0L, NA_integer_))
  expect_equal(result$ASTER_lwm_water_perc, c(100 / 3, 50, NA))
  expect_equal(result$ASTER_lwm_water_cin, c(200 / 3, 100, NA))
  expect_equal(result$ASTER_lwm_water_missing_perc, c(100 / 3, NA, 100))
  expect_identical(.spss_aggregate(values[FALSE, ], character()), result[FALSE, ])
})

test_that("SPSS percentages reproduce the published 140-point cases", {
  values <- data.frame(
    elev = c(rep(0, 84), rep(1, 34), rep(NA_real_, 22 + 140)),
    lwm = c(.ag100_water_value, rep(NA_real_, 139 + 140))
  )
  result <- .spss_aggregate(values, rep(1:2, each = 140))
  expect_equal(result$ASTER_zip_elevation_cin, c(60, NA))
  expect_equal(result$ASTER_zip_elevation_missing_perc, c(100 * 22 / 140, 100))
  expect_equal(result$ASTER_lwm_water_perc, c(100 / 140, NA))
  expect_equal(result$ASTER_lwm_water_cin, c(100 / 140, NA))
  expect_equal(result$ASTER_lwm_water_missing_perc, c(100 * 139 / 140, 100))
})
