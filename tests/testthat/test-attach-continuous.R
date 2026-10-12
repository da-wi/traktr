make_df <- function() {
  pts <- data.frame(pid = 1:3, year = c(2005L, 2012L, NA),
                    lat = c(48.400, 48.401, 48.420), lon = c(15.600, 15.601, 15.620))
  mutate_aoi_buffer(pts, dist = 300)
}

write_const_raster <- function(path, value, crs = "EPSG:4326", nlyr = 1) {
  r <- terra::rast(terra::ext(15.5, 15.7, 48.3, 48.5), resolution = 0.001, crs = crs, nlyrs = nlyr)
  terra::values(r) <- rep(value, each = terra::ncell(r))[seq_len(terra::ncell(r) * nlyr)]
  terra::writeRaster(r, path, overwrite = TRUE)
  path
}

test_that("weighted summary handles empty and weighted input", {
  s <- .continuous_summary_fn(c(NA, NA), c(1, 1))
  expect_true(all(is.na(s)))
  s <- .continuous_summary_fn(c(0, 1), c(1, 3))
  expect_equal(s$mean, 0.75)
  expect_equal(s$median, 1)
  expect_equal(s$min, 0)
  expect_equal(s$max, 1)
})

test_that("attach_ndvi uses band 1 of multi-band chips", {
  skip_if_not_installed("exactextractr")
  dir <- withr::local_tempdir()
  # band 1 = 0.5 (NDVI), band 2 = 7 (observation count)
  write_const_raster(file.path(dir, "ndvi_chip.tif"), c(0.5, 7), nlyr = 2)
  out <- attach_ndvi(make_df(), out_dir = dir, progress = FALSE)
  expect_equal(out$ndvi_mean, rep(0.5, 3), tolerance = 1e-6)
  expect_equal(out$pid, 1:3)
  expect_s3_class(out, "traktr_df")
  expect_false(is.null(attr(out, "traktr_aoi")))
})

test_that("attach_hfp picks the nearest year per row", {
  skip_if_not_installed("exactextractr")
  dir <- withr::local_tempdir()
  write_const_raster(file.path(dir, "hfp_2005.tif"), 10)
  write_const_raster(file.path(dir, "hfp_2010.tif"), 20)

  out <- attach_hfp(make_df(), out_dir = dir, year_col = "year", progress = FALSE)
  expect_equal(out$hfp_year, c(2005L, 2010L, NA))
  expect_equal(out$hfp_mean, c(10, 20, NA), tolerance = 1e-6)

  out <- attach_hfp(make_df(), out_dir = dir, progress = FALSE)
  expect_equal(out$hfp_year, rep(2010L, 3))
})

test_that("nearest year mapping", {
  expect_equal(.nearest_year(c(1999L, 2007L, NA, 2030L), c(2000L, 2005L, 2010L)),
               c(2000L, 2005L, NA, 2010L))
})

test_that("attach_ndvi min_obs drops pixels with too few observations", {
  skip_if_not_installed("exactextractr")
  dir <- withr::local_tempdir()
  # band 1 NDVI: left half 0.2, right half 0.8; band 2 count: left 1, right 10
  r <- terra::rast(terra::ext(15.5, 15.7, 48.3, 48.5), resolution = 0.001, crs = "EPSG:4326", nlyrs = 2)
  xy <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  left <- xy < 15.6005
  terra::values(r) <- cbind(ifelse(left, 0.2, 0.8), ifelse(left, 1, 10))
  terra::writeRaster(r, file.path(dir, "ndvi_chip.tif"))
  
  all <- attach_ndvi(make_df(), out_dir = dir, progress = FALSE)
  few <- attach_ndvi(make_df(), out_dir = dir, min_obs = 5, progress = FALSE)
  expect_gt(all$ndvi_max[1], all$ndvi_min[1])         # AOI 1 straddles both halves
  expect_equal(few$ndvi_mean[1], 0.8, tolerance = 1e-6) # only well-observed pixels left
})

test_that("chip file names tag the composite method, leastcc keeps old names", {
  chip <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = 1, ymin = 2, xmax = 3, ymax = 4), crs = 4326)))
  expect_equal(basename(.ndvi_chip_file("d", chip, "2024-06-01", "2024-08-31", "leastcc")),
               "ndvi_1_2_3_4_2024-06-01_2024-08-31.tif")
  expect_equal(basename(.ndvi_chip_file("d", chip, "2024-06-01", "2024-08-31", "median")),
               "ndvi_median_1_2_3_4_2024-06-01_2024-08-31.tif")
  expect_equal(basename(.ndvi_chip_file("d", chip, "2024-06-01", "2024-08-31", "medianfb")),
               "ndvi_medianfb_1_2_3_4_2024-06-01_2024-08-31.tif")
})

test_that("attach_ndvi reports observation count and worst tier", {
  skip_if_not_installed("exactextractr")
  dir <- withr::local_tempdir()
  # left half: tier 1 with 3 obs; right half: tier 3 with 1 obs
  r <- terra::rast(terra::ext(15.5, 15.7, 48.3, 48.5), resolution = 0.001, crs = "EPSG:4326", nlyrs = 3)
  left <- terra::xFromCell(r, seq_len(terra::ncell(r))) < 15.6005
  terra::values(r) <- cbind(ifelse(left, 0.5, 0.3), ifelse(left, 3, 1), ifelse(left, 1, 3))
  terra::writeRaster(r, file.path(dir, "ndvi_medianfb_chip.tif"))

  out <- attach_ndvi(make_df(), out_dir = dir, progress = FALSE)
  expect_equal(out$ndvi_tier, c(3, 3, 3))          # AOIs 1-2 straddle, AOI 3 right half
  expect_true(out$ndvi_n_obs[1] > 1 && out$ndvi_n_obs[1] < 3)
  expect_equal(out$ndvi_n_obs[3], 1)

  # least-cloudy chips (band 2 = data mask) get no count/tier columns
  dir2 <- withr::local_tempdir()
  write_const_raster(file.path(dir2, "ndvi_1_2_3_4_2024-06-01_2024-06-30.tif"), c(0.5, 1), nlyr = 2)
  out2 <- attach_ndvi(make_df(), out_dir = dir2, progress = FALSE)
  expect_false(any(c("ndvi_n_obs", "ndvi_tier") %in% names(out2)))
})
