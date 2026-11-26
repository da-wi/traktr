tiled_stats <- c("mean", "sd", "min", "max", "median", "missing_perc", "n_cells", "area_valid_km2")

tiled_polygon <- function(xmin = 0, xmax = 1, ymin = 0, ymax = 1) {
  sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                            crs = 4326))
}

tiled_aoi <- function(...) sf::st_sf(region = "a", geometry = tiled_polygon(...))

tiled_raster <- function(xmin = 0, xmax = 1, ymin = 0, ymax = 1) {
  r <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    resolution = 0.01, crs = "EPSG:4326", nlyrs = 2)
  withr::local_seed(173)
  n <- terra::ncell(r)
  values <- cbind(round(stats::runif(n, 29000, 29100)), round(stats::runif(n, -10, 10)))
  values[seq.int(1L, n, by = 7L), 1L] <- NA
  values[seq.int(1L, n, by = 11L), 2L] <- NA
  terra::values(r) <- values
  names(r) <- c("temp", "elev")
  r
}

tiled_write <- function(r, dir, id) {
  path <- file.path(dir, paste0(id, ".tif"))
  terra::writeRaster(r, path, datatype = "INT4S", NAflag = -9999)
  data.frame(tile_id = id, file = path)
}

tiled_map <- function(tiles, id = "a") data.frame(aoi_id = id, tile_id = tiles$tile_id)

test_that("single-tile moments agree with continuous extraction and use square metres", {
  dir <- withr::local_tempdir()
  tile <- tiled_write(tiled_raster(ymin = 47, ymax = 48), dir, "one")
  aoi <- tiled_aoi(xmin = 0.015, xmax = 0.82, ymin = 47.12, ymax = 47.8)
  moments <- .extract_tiled(aoi, "region", tile, tiled_map(tile), progress = FALSE)
  expect_named(moments$a, c("temp", "elev"))
  stats <- tiled_stats[1:5]
  for (i in seq_along(moments$a)) {
    expected <- .extract_continuous(
      aoi, "region", tile$file, band = i, crs_planar = NULL, weights = "area"
    )
    actual <- .finalize_moments(moments$a[[i]], stats)
    expect_equal(actual, expected[, stats], tolerance = 1e-9)
  }

  single <- .extract_tiled(aoi, "region", tile, tiled_map(tile), layers = "elev", progress = FALSE)
  expect_identical(single$a, moments$a["elev"])
  r <- terra::rast(tile$file)[["elev"]]
  cells <- exactextractr::exact_extract(r, aoi, coverage_area = TRUE, progress = FALSE)[[1L]]
  valid <- !is.na(cells$value)
  expect_equal(single$a$elev$w, sum(cells$coverage_area[valid]), tolerance = 1e-9)
  expect_equal(single$a$elev$w_na, sum(cells$coverage_area[!valid]), tolerance = 1e-9)
  expect_gt(single$a$elev$w, 1e9)
})

test_that("whole and split rasters yield the same statistics regardless of tile order", {
  dir <- withr::local_tempdir()
  r <- tiled_raster(0, 2, 0, 2)
  whole <- tiled_write(r, dir, "whole")
  grid <- expand.grid(x = 0:1, y = 0:1)
  tiles <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    x <- grid$x[i]
    y <- grid$y[i]
    tiled_write(terra::crop(r, terra::ext(x, x + 1, y, y + 1)), dir, paste0("tile", i))
  }))
  aoi <- tiled_aoi(0.24, 1.71, 0.28, 1.75)
  extract <- function(index) .extract_tiled(aoi, "region", index, tiled_map(index), progress = FALSE)
  expected <- extract(whole)
  for (index in list(tiles, tiles[4:1, ])) {
    actual <- extract(index)
    for (layer in names(expected$a)) {
      result <- .finalize_moments(actual$a[[layer]], tiled_stats)
      reference <- .finalize_moments(expected$a[[layer]], tiled_stats)
      expect_equal(result, reference, tolerance = 1e-9)
      expect_identical(result$n_cells, reference$n_cells)
      expect_identical(result$median, reference$median)
    }
  }
})

test_that("centered moments combine associatively and discard invalid weights", {
  a <- .moments_fn(c(29001, 29002, NA, Inf), c(2e6, 3e6, 5e6, 1e6))
  b <- .moments_fn(c(29002, 29003, 29004, NA), c(4e6, 1e6, 2e6, Inf))
  c <- .moments_fn(c(29001, 29004, 29005, 29006), c(7e6, 6e6, 0, -1e6))
  left <- .combine_moments(.combine_moments(a, b), c)
  right <- .combine_moments(a, .combine_moments(c, b))
  expect_equal(left, right, tolerance = 1e-9)
  expect_equal(left$hist, data.frame(value = as.double(29001:29004), w = c(9, 7, 1, 8) * 1e6))
  expect_equal(left$w_na, 5e6)
  expect_identical(left$n, 7L)
  missing <- .moments_fn(c(NA, NA), c(1, 2))
  expect_equal(.combine_moments(left, missing)$w_na, 5e6 + 3)
  expect_identical(.combine_moments(left, missing), .combine_moments(missing, left))
  expect_equal(.combine_moments(missing, missing)$w_na, 6)
})

test_that("missing percentages measure covered area", {
  dir <- withr::local_tempdir()
  r <- tiled_raster()
  left <- terra::xFromCell(r, seq_len(terra::ncell(r))) < 0.5
  terra::values(r) <- cbind(ifelse(left, NA, 29000), 1)
  tile <- tiled_write(r, dir, "missing")
  moments <- .extract_tiled(tiled_aoi(), "region", tile, tiled_map(tile), progress = FALSE)$a
  expect_equal(.finalize_moments(moments$temp, "missing_perc")$missing_perc, 50, tolerance = 1e-2)
  expect_identical(.finalize_moments(moments$elev, "missing_perc")$missing_perc, 0)
  missing <- .finalize_moments(.moments_fn(NA_real_, 1), tiled_stats)
  expect_true(all(is.na(missing[, tiled_stats[1:5]])))
  expect_equal(missing$missing_perc, 100)
  expect_identical(missing$n_cells, 0L)
  expect_equal(missing$area_valid_km2, 0)
})

test_that("an exact weighted-median tie selects the smaller integer", {
  values <- c(4, 7, 4, 7)
  weights <- c(2, 3, 3, 2) * 1e6
  expected <- .continuous_summary_fn(values, weights)$median
  moments <- .combine_moments(.moments_fn(values[1:2], weights[1:2]),
                              .moments_fn(values[3:4], weights[3:4]))
  expect_equal(expected, 4)
  expect_identical(.finalize_moments(moments, "median")$median, expected)
})

test_that("the tiled median selects the lower value at a rounded half-weight", {
  dir <- withr::local_tempdir()
  r <- terra::rast(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                    nrows = 1, ncols = 12, crs = "EPSG:4326")
  terra::values(r) <- rep(c(4, 7), each = 6)
  names(r) <- "value"
  tile <- tiled_write(r, dir, "tie")
  moments <- .extract_tiled(tiled_aoi(), "region", tile, tiled_map(tile), progress = FALSE)
  expect_identical(.finalize_moments(moments$a$value, "median")$median, 4)

  parts <- lapply(list(c(0.8, 0.4), c(0.6, 0.6), c(0.4, 0.8)), function(weights) {
    .moments_fn(c(4, 7), weights)
  })
  for (order in list(1:3, c(1, 3, 2))) {
    combined <- Reduce(.combine_moments, parts[order])
    expect_identical(.finalize_moments(combined, "median")$median, 4)
  }

  below_cutoff <- .moments_fn(c(4, 7), c(0.5 - 1e-8, 0.5 + 1e-8))
  expect_identical(.finalize_moments(below_cutoff, "median")$median, 7)
})

test_that("multipart AOIs extract each side of the antimeridian once", {
  dir <- withr::local_tempdir()
  west <- tiled_write(tiled_raster(-180, -179), dir, "west")
  east <- tiled_write(tiled_raster(179, 180), dir, "east")
  tiles <- rbind(west, east)
  polygons <- c(tiled_polygon(-180, -179.6), tiled_polygon(179.6, 180))
  aoi <- sf::st_sf(region = "a", geometry = sf::st_combine(polygons))
  actual <- .extract_tiled(aoi, "region", tiles, tiled_map(tiles), progress = FALSE)$a
  parts <- lapply(1:2, function(i) {
    part <- sf::st_sf(region = "a", geometry = polygons[i])
    .extract_tiled(part, "region", tiles[i, ], tiled_map(tiles[i, ]), progress = FALSE)$a
  })
  expected <- Map(.combine_moments, parts[[1L]], parts[[2L]])
  expect_equal(actual, expected, tolerance = 1e-9)
  expect_identical(actual$elev$n, parts[[1L]]$elev$n + parts[[2L]]$elev$n)
})

test_that("missing tiles leave mapped AOIs with empty moments", {
  dir <- withr::local_tempdir()
  tiles <- rbind(tiled_write(tiled_raster(), dir, "one"),
                 data.frame(tile_id = "absent", file = NA_character_))
  aoi <- sf::st_sf(region = c(2L, 8L, 9L),
                    geometry = c(tiled_polygon(), tiled_polygon(1, 2), tiled_polygon(2, 3)))
  map <- tiled_map(tiles, id = c(2L, 8L))
  out <- .extract_tiled(aoi, "region", tiles, map, progress = FALSE)
  expect_named(out, c("2", "8", "9"))
  expect_gt(out[["2"]]$temp$w, 0)
  expect_identical(out[["8"]], out[["9"]])
  missing <- .finalize_moments(out[["8"]]$temp, tiled_stats)
  expect_true(all(is.na(missing[, tiled_stats[1:6]])))
  expect_identical(missing$n_cells, 0L)
  expect_equal(missing$area_valid_km2, 0)
  all_missing <- .extract_tiled(aoi, "region", tiles[2L, ], map,
                                layers = "temp", progress = FALSE)
  expect_named(all_missing, names(out))
  expect_identical(all_missing[["2"]], out[["8"]]["temp"])
})

test_that("overlapping, misaligned, and non-4326 tiles and unknown layers are rejected", {
  dir <- withr::local_tempdir()
  tile <- tiled_write(tiled_raster(), dir, "one")
  overlap <- tiled_write(tiled_raster(0, 2, 0, 2), dir, "overlap")
  misaligned <- tiled_write(tiled_raster(0.01, 1.01), dir, "misaligned")
  projected <- tiled_raster()
  terra::crs(projected) <- "EPSG:3857"
  projected <- tiled_write(projected, dir, "projected")
  aoi <- tiled_aoi()
  extract <- function(tiles, ...) {
    .extract_tiled(aoi, "region", tiles, tiled_map(tiles), progress = FALSE, ...)
  }
  expect_error(extract(rbind(tile, overlap)), "overlaps.*overlap")
  expect_error(extract(misaligned), "misaligned.*whole-degree")
  expect_error(extract(projected), "projected.*EPSG:4326")
  expect_error(extract(tile, layers = "unknown"), "Unknown layer.*unknown.*one")
})

test_that("nominal tile boxes count a shared AG100 row once", {
  dir <- withr::local_tempdir()
  res <- 1 / 19
  h <- res / 2
  r <- terra::rast(xmin = -h, xmax = 1 + h, ymin = 39 - h, ymax = 41 + h,
                   ncols = 20, nrows = 39, crs = "EPSG:4326")
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "value"
  tiles <- rbind(
    tiled_write(terra::crop(r, terra::ext(-h, 1 + h, 40 - h, 41 + h)), dir, "north"),
    tiled_write(terra::crop(r, terra::ext(-h, 1 + h, 39 - h, 40 + h)), dir, "south")
  )
  tiles$xmin <- 0
  tiles$xmax <- 1
  tiles$ymin <- c(40, 39)
  tiles$ymax <- c(41, 40)
  whole <- transform(tiled_write(r, dir, "whole"), xmin = 0, xmax = 1, ymin = 39, ymax = 41)
  aoi <- tiled_aoi(0.19, 0.81, 39.75, 40.375)
  extract <- function(index, ...) {
    .extract_tiled(aoi, "region", index, tiled_map(index), progress = FALSE, ...)$a$value
  }
  original_s2 <- getOption("sf_use_s2")
  cache <- file.path(dir, "moments")
  actual <- extract(tiles, cache_dir = cache)
  expect_identical(getOption("sf_use_s2"), original_s2)
  expected <- extract(whole)
  stats <- c("mean", "sd", "median", "n_cells")
  expect_equal(.finalize_moments(actual, stats), .finalize_moments(expected, stats),
               tolerance = 1e-9)
  expect_identical(actual$n, expected$n)
  expect_equal(actual$w / expected$w, 1, tolerance = 1e-6)
  expect_identical(extract(tiles, cache_dir = cache), actual)

  # Bypass the overlap guard for the countercheck and extract both rasters
  # directly with the unclipped AOI.
  unclipped <- lapply(tiles$file, function(file) {
    .tile_moments(terra::rast(file), sf::st_geometry(aoi))[[1L]]$value
  })
  expect_gt(sum(vapply(unclipped, `[[`, numeric(1), "w")), expected$w)
})

test_that("nominal boxes are validated even for the same underlying raster", {
  dir <- withr::local_tempdir()
  tile <- tiled_write(tiled_raster(0, 2), dir, "wide")
  west <- transform(tile, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  east <- transform(tile, tile_id = "east", xmin = 1, xmax = 2, ymin = 0, ymax = 1)
  aoi <- tiled_aoi(0.25, 1, 0.25, 0.75)
  extract <- function(tiles, ...) {
    .extract_tiled(aoi, "region", tiles, tiled_map(tiles), progress = FALSE, ...)
  }
  expect_identical(extract(east)$a$elev$n, 0L)
  expect_identical(extract(rbind(west, east)), extract(west))
  expect_error(extract(rbind(transform(west, xmax = 2), east)), "wide.*overlaps.*east")
  expect_error(extract(transform(west, xmin = -1)), "wide.*does not cover.*nominal box")

  cache <- file.path(dir, "moments")
  expect_identical(extract(west, cache_dir = cache),
                   extract(transform(west, xmax = 2), cache_dir = cache))
  expect_length(list.files(cache), 2L)
})

test_that("tiles without a CRS are rejected even when their bounds look geographic", {
  dir <- withr::local_tempdir()
  r <- tiled_raster()
  terra::crs(r) <- ""
  tile <- tiled_write(r, dir, "untagged")
  expect_identical(terra::crs(terra::rast(tile$file, guessCRS = FALSE)), "")
  expect_error(
    .extract_tiled(tiled_aoi(), "region", tile, tiled_map(tile), progress = FALSE),
    "untagged.*EPSG:4326"
  )
})

test_that("the cache resumes per tile, invalidates with geometry and repairs corrupt RDS", {
  dir <- withr::local_tempdir()
  cache <- file.path(dir, "moments")
  tiles <- tiled_write(tiled_raster(), dir, "one")
  # Multiple polygons in one tile exercise the callback's nested return shape.
  aoi <- sf::st_sf(region = c("a", "b"),
                    geometry = c(tiled_polygon(), tiled_polygon(0, 0.5)))
  map <- data.frame(aoi_id = aoi$region, tile_id = "one")
  extract <- function(aoi) {
    .extract_tiled(aoi, "region", tiles, map, cache_dir = cache, progress = FALSE)
  }
  first <- extract(aoi)
  files <- list.files(cache, full.names = TRUE)
  expect_length(files, 1L)
  expect_named(first$b, c("temp", "elev"))
  expect_gt(first$a$elev$n, first$b$elev$n)

  worker <- .tile_moments
  calls <- 0L
  local_mocked_bindings(.tile_moments = function(...) {
    calls <<- calls + 1L
    worker(...)
  })
  expect_identical(extract(aoi), first)
  expect_identical(calls, 0L)
  changed <- aoi
  sf::st_geometry(changed)[1L] <- tiled_polygon(0, 0.7)
  expect_false(identical(extract(changed), first))
  expect_identical(calls, 1L)
  expect_length(list.files(cache), 2L)

  writeLines("not an RDS", files)
  expect_warning(repaired <- extract(aoi), "cache.*one.*unreadable")
  expect_identical(repaired, first)
  expect_identical(calls, 2L)
  expect_identical(extract(aoi), first)
  expect_identical(calls, 2L)
})

test_that("replacing raster contents invalidates a cache with the same size and mtime", {
  dir <- withr::local_tempdir()
  cache <- file.path(dir, "moments")
  tile <- data.frame(tile_id = "one", file = file.path(dir, "one.tif"))
  replacement <- file.path(dir, "replacement.tif")
  r <- terra::rast(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                    nrows = 2, ncols = 2, crs = "EPSG:4326")
  names(r) <- "value"
  terra::values(r) <- 1
  terra::writeRaster(r, tile$file, datatype = "INT4S", gdal = "COMPRESS=NONE")
  terra::values(r) <- 2
  terra::writeRaster(r, replacement, datatype = "INT4S", gdal = "COMPRESS=NONE")
  timestamp <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  expect_true(Sys.setFileTime(tile$file, timestamp))
  stamp <- file.info(tile$file)
  extract <- function() {
    .extract_tiled(tiled_aoi(), "region", tile, tiled_map(tile),
                   cache_dir = cache, progress = FALSE)
  }
  expect_equal(extract()$a$value$mean, 1)
  expect_true(file.copy(replacement, tile$file, overwrite = TRUE))
  expect_true(Sys.setFileTime(tile$file, stamp$mtime))
  expect_identical(file.size(tile$file), stamp$size)
  expect_identical(file.mtime(tile$file), stamp$mtime)
  expect_equal(extract()$a$value$mean, 2)
  expect_length(list.files(cache), 2L)
})

test_that("readable cache records with invalid versions or schemas are recomputed", {
  dir <- withr::local_tempdir()
  cache <- file.path(dir, "moments")
  tile <- tiled_write(tiled_raster(), dir, "one")
  extract <- function() {
    .extract_tiled(tiled_aoi(), "region", tile, tiled_map(tile),
                   cache_dir = cache, progress = FALSE)
  }
  first <- extract()
  path <- list.files(cache, full.names = TRUE)
  valid <- readRDS(path)
  bad_version <- valid
  bad_version$version <- valid$version + 1L
  bad_layer <- valid
  bad_layer$moments$a["temp"] <- list(NULL)
  bad_field <- valid
  bad_field$moments$a$temp$w <- c(1, 2)
  bad_hist <- valid
  bad_hist$moments$a$temp$hist$w <- as.character(bad_hist$moments$a$temp$hist$w)

  for (bad in list(bad_version, bad_layer, bad_field, bad_hist)) {
    saveRDS(bad, path)
    expect_warning(repaired <- extract(), "cache.*one.*invalid schema")
    expect_identical(repaired, first)
    expect_identical(readRDS(path), valid)
  }
})

test_that("parallel tile extraction agrees with sequential extraction", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  skip_on_os("windows")
  skip_if_not(getExportedValue("future", "supportsMulticore")())
  plan <- getExportedValue("future", "plan")
  previous <- plan("multicore", workers = 2)
  withr::defer(plan(previous))
  dir <- withr::local_tempdir()
  tiles <- rbind(tiled_write(tiled_raster(), dir, "left"),
                 tiled_write(tiled_raster(1, 2), dir, "right"))
  aoi <- tiled_aoi(0.3, 1.6)
  map <- tiled_map(tiles)
  sequential <- .extract_tiled(aoi, "region", tiles, map, progress = FALSE)
  parallel <- .extract_tiled(aoi, "region", tiles, map, workers = 2, progress = FALSE)
  expect_identical(parallel, sequential)

  local_mocked_bindings(.has_future_apply = function() FALSE)
  expect_warning(
    fallback <- .extract_tiled(aoi, "region", tiles, map, workers = 2, progress = FALSE),
    "future.apply.*unavailable"
  )
  expect_identical(fallback, sequential)
})

test_that("continuous extraction defaults to coverage fractions", {
  dir <- withr::local_tempdir()
  tile <- tiled_write(tiled_raster(), dir, "one")
  aoi <- tiled_aoi(0.2, 0.8, 0.2, 0.8)
  expect_identical(.extract_continuous(aoi, "region", tile$file),
                   .extract_continuous(aoi, "region", tile$file, weights = "fraction"))
})

test_that("continuous extraction can transform directly from the AOI to the raster CRS", {
  dir <- withr::local_tempdir()
  tile <- tiled_write(tiled_raster(), dir, "one")
  transform <- sf::st_transform
  aoi <- transform(tiled_aoi(0.2, 0.8, 0.2, 0.8), 3857)
  transformations <- integer()
  local_mocked_bindings(.package = "sf", st_transform = function(x, crs, ...) {
    if (inherits(x, "sf")) transformations <<- c(transformations, sf::st_crs(crs)$epsg)
    transform(x, crs, ...)
  })
  out <- .extract_continuous(aoi, "region", tile$file, crs_planar = NULL)
  expect_true(is.finite(out$mean))
  expect_identical(transformations, 4326L)
})

test_that("histogram shares use only valid covered area", {
  moments <- .moments_fn(c(0, 0, 5, NA), c(1, 2, 7, 4))
  expect_equal(.hist_share(moments, function(value) value == 0), 30)
  expect_equal(.hist_share(moments, function(value) value < 4), 30)
  expect_identical(.hist_share(.moments_fn(NA_real_, 1), is.na), NA_real_)
})

test_that("a raster replaced during extraction is not cached under its old fingerprint", {
  dir <- withr::local_tempdir()
  cache <- file.path(dir, "moments")
  dir.create(cache)
  path <- file.path(dir, "one.tif")
  r <- terra::rast(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                   nrows = 2, ncols = 2, crs = "EPSG:4326")
  names(r) <- "value"
  terra::values(r) <- 1
  terra::writeRaster(r, path, datatype = "INT4S")
  replacement <- file.path(dir, "replacement.tif")
  terra::values(r) <- 2
  terra::writeRaster(r, replacement, datatype = "INT4S")
  tile <- data.frame(tile_id = "one", file = path)
  aoi <- tiled_aoi()

  worker <- .tile_moments
  local_mocked_bindings(.tile_moments = function(...) {
    out <- worker(...)
    file.copy(replacement, path, overwrite = TRUE)
    out
  })
  expect_warning(
    .extract_tiled(aoi, "region", tile, tiled_map(tile), cache_dir = cache, progress = FALSE),
    "changed during extraction"
  )
  expect_length(list.files(cache), 0L)
})
