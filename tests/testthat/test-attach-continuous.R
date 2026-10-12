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
