ag100_fixture <- function(file = "AG100.v003.40.-105.0001.h5") {
  testthat::test_path("fixtures", "ag100", file)
}

ag100_skip_without_hdf5 <- function() {
  skip_if_not("HDF5" %in% terra::gdal(drivers = TRUE)$name, "GDAL has no HDF5 driver")
}

ag100_fetch_aoi <- function(xmin = -104.9, xmax = -104.1, ymin = 39.1, ymax = 39.9) {
  bounds <- sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = 4326)
  aoi_from_sf(sf::st_sf(id = "a", geometry = sf::st_as_sfc(bounds)), id_col = "id")
}

ag100_mock_urls <- function(tiles, use_cmr) {
  tiles$url <- paste0("https://example.test/", tiles$tile_id, ".h5")
  tiles$url_source <- "cmr"
  tiles$status <- NA_character_
  tiles
}

test_that("AG100 reads raw pixels from a rotated source in row order", {
  dir <- withr::local_tempdir()
  tif <- file.path(dir, "source & pixels.tif")
  vrt <- file.path(dir, "rotated & scaled.vrt")
  r <- terra::rast(nrows = 2L, ncols = 3L, vals = 1:6)
  terra::writeRaster(r, tif, datatype = "INT4S")
  writeLines(c(
    '<VRTDataset rasterXSize="3" rasterYSize="2">',
    '  <GeoTransform>0,0,0.001,40,-0.001,0</GeoTransform>',
    '  <VRTRasterBand dataType="Int32" band="1">',
    '    <Scale>10</Scale>',
    '    <Offset>100</Offset>',
    '    <SimpleSource>',
    '      <SourceFilename relativeToVRT="1">source &amp; pixels.tif</SourceFilename>',
    '      <SourceBand>1</SourceBand>',
    '    </SimpleSource>',
    '  </VRTRasterBand>',
    '</VRTDataset>'
  ), vrt)

  expect_warning(terra::rast(vrt), "rotated")
  expect_no_warning(result <- .ag100_read_array(vrt))
  expect_equal(result, matrix(1:6, nrow = 2L, byrow = TRUE))
})

test_that("AG100 conversion preserves band order, georeferencing and raw values", {
  ag100_skip_without_hdf5()
  tif <- withr::local_tempfile(fileext = ".tif")
  .ag100_h5_to_tif(ag100_fixture(), tif, 40, -105)
  r <- terra::rast(tif)

  bounds <- c(-105, -104, 39, 40) + c(-1, 1, -1, 1) / (2 * 19)
  expect_equal(unname(as.vector(terra::ext(r))), bounds, tolerance = 1e-12)
  expect_equal(terra::res(r), rep(1 / 19, 2), tolerance = 1e-12)
  expect_identical(sf::st_crs(terra::crs(r))$epsg, 4326L)
  expect_identical(names(r), c("elev", "ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd", "lwm", "nobs"))
  expect_identical(terra::datatype(r), rep("INT4S", 7))
  values <- terra::values(r)
  expect_equal(unname(values[1, ]), c(1511, 302, 21, 33001, 10, 1, 3))
  expect_equal(unname(values[400, ]), c(1720, 340, 40, 33020, 200, 0, 60))
  expect_true(all(is.na(values[c(87, 88, 107, 108), ])))
  expect_error(.ag100_h5_to_tif(ag100_fixture(), tif, 41, -105), "extent")
  expect_equal(unname(as.vector(terra::ext(terra::rast(tif)))), bounds)
})

test_that("AG100 rejects the old inner grid and reports observed and expected endpoints", {
  ag100_skip_without_hdf5()
  tif <- withr::local_tempfile(fileext = ".tif")
  h5 <- ag100_fixture(file.path("bad", "AG100.v003.40.-105.0001.h5"))
  err <- expect_error(.ag100_h5_to_tif(h5, tif, 40, -105), "Observed endpoints")
  observed <- utils::strcapture(
    paste0("Observed endpoints: Geolocation/Latitude (-?[0-9]+[.][0-9]+) to (-?[0-9]+[.][0-9]+), ",
           "Geolocation/Longitude (-?[0-9]+[.][0-9]+) to (-?[0-9]+[.][0-9]+)"),
    conditionMessage(err), list(north = double(), south = double(), west = double(), east = double())
  )
  expect_equal(unlist(observed, use.names = FALSE), c(39.975, 39.025, -104.975, -104.025),
               tolerance = 1e-7)
  expect_match(conditionMessage(err),
               "Expected endpoints: Geolocation/Latitude 40 to 39, Geolocation/Longitude -105 to -104")
  expect_false(file.exists(tif))
})

test_that("the local Denver AG100 tile has inclusive endpoints and plausible values", {
  skip_on_ci()
  skip_on_cran()
  h5 <- testthat::test_path("..", "..", "testoutput", "ag100", "ag100_v003_40_-105.h5")
  skip_if_not(file.exists(h5), "Local Denver AG100 tile is unavailable")
  ag100_skip_without_hdf5()
  tif <- withr::local_tempfile(fileext = ".tif")
  .ag100_h5_to_tif(h5, tif, 40, -105)
  r <- terra::rast(tif)
  expect_equal(unname(as.vector(terra::ext(r))),
               c(-105, -104, 39, 40) + c(-1, 1, -1, 1) / (2 * 999), tolerance = 1e-12)
  expect_equal(terra::res(r), rep(1 / 999, 2), tolerance = 1e-12)
  elev <- terra::extract(r[["elev"]], matrix(c(-104.9903, 39.7392), ncol = 2))$elev
  expect_true(elev >= 1550 && elev <= 1650)
  expect_setequal(unique(terra::values(r[["lwm"]], mat = FALSE)), c(0, 1))
})

test_that("AG100 transposes both band and coordinate arrays in a transposed layout", {
  ag100_skip_without_hdf5()
  original <- withr::local_tempfile(fileext = ".tif")
  transposed <- withr::local_tempfile(fileext = ".tif")
  .ag100_h5_to_tif(ag100_fixture(), original, 40, -105)
  read_array <- .ag100_read_array
  local_mocked_bindings(.ag100_read_array = function(address) t(read_array(address)))

  expect_message(.ag100_h5_to_tif(ag100_fixture(), transposed, 40, -105), "transposed layout")
  r <- terra::rast(original)
  s <- terra::rast(transposed)
  expect_equal(as.vector(terra::ext(s)), as.vector(terra::ext(r)))
  expect_identical(terra::values(s), terra::values(r))
})

test_that("AG100 masks the HDF5 dataset's declared fill value", {
  ag100_skip_without_hdf5()
  h5 <- ag100_fixture(file.path("nodata", "AG100.v003.40.-105.0001.h5"))
  tif <- withr::local_tempfile(fileext = ".tif")
  .ag100_h5_to_tif(h5, tif, 40, -105)
  values <- terra::values(terra::rast(tif))

  expect_true(is.na(values[1L, "elev"]))
  expect_equal(unname(values[2L, "elev"]), 1512)
  expect_equal(unname(values[1L, "ndvi_mean"]), 302)
  expect_true(all(is.na(values[c(87, 88, 107, 108), ])))
})

test_that("AG100 accepts one-dimensional coordinates and diagnoses swapped axes", {
  ag100_skip_without_hdf5()
  tif <- withr::local_tempfile(fileext = ".tif")
  h5 <- ag100_fixture("AG100.v003.05.-075.0001.h5")
  read_array <- .ag100_read_array
  swapped <- FALSE
  local_mocked_bindings(.ag100_read_array = function(address) {
    if (endsWith(address, "/Latitude")) {
      val <- if (swapped) read_array(sub("/Latitude$", "/Longitude", address))[1L, ] else
        read_array(address)[, 1L]
      return(matrix(val, nrow = 1L))
    }
    if (endsWith(address, "/Longitude")) {
      val <- if (swapped) read_array(sub("/Longitude$", "/Latitude", address))[, 1L] else
        read_array(address)[1L, ]
      return(matrix(val, nrow = 1L))
    }
    read_array(address)
  })
  expect_no_error(.ag100_h5_to_tif(h5, tif, 5, -75))
  swapped <- TRUE
  expect_error(.ag100_h5_to_tif(h5, tif, 5, -75), "lat/lon axes swapped")
})

test_that("AG100 diagnoses reversed square coordinate vectors within the tile bounds", {
  ag100_skip_without_hdf5()
  tif <- withr::local_tempfile(fileext = ".tif")
  lat <- c(39, 40)
  lon <- c(40, 39)
  local_mocked_bindings(.ag100_read_array = function(address) {
    if (endsWith(address, "/Latitude")) return(matrix(lat, nrow = 1L))
    if (endsWith(address, "/Longitude")) return(matrix(lon, nrow = 1L))
    matrix(1:4, nrow = 2L)
  })
  expect_error(.ag100_h5_to_tif(ag100_fixture(), tif, 40, 39), "lat/lon axes swapped\\?")
  lat <- rev(lat)
  expect_error(.ag100_h5_to_tif(ag100_fixture(), tif, 40, 39), "lat/lon axes swapped\\?")
  expect_false(file.exists(tif))
})

test_that("AG100 accepts float32 coordinate rounding and rejects an irregular step", {
  res <- 1 / 999
  centres <- seq(0, 1, length.out = 1000L)
  lat <- 40 - centres
  lon <- readBin(writeBin(-105 + centres, raw(), size = 4L),
                 "double", n = length(centres), size = 4L)

  expect_equal(.ag100_extent(lat, lon, 40, -105),
               c(-105, -104, 39, 40) + c(-1, 1, -1, 1) * res / 2)
  lon[500L] <- lon[500L] + res / 20
  expect_error(.ag100_extent(lat, lon, 40, -105), "constant step")
})

test_that("local AG100 conversion and the cache bypass remote lookup and authentication", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  calls <- 0L
  remote <- function(...) {
    calls <<- calls + 1L
    stop("Unexpected remote action")
  }
  local_mocked_bindings(.ag100_resolve = remote, .earthdata_handle = remote, .ag100_download = remote)

  df <- ag100_fetch_aoi()
  h5 <- ag100_fixture()
  original <- tools::md5sum(h5)
  result <- withVisible(fetch_aster_ged(df, out, src = dirname(h5), progress = FALSE))
  expect_false(result$visible)
  index <- attr(result$value, "aster_ged_tiles")
  expect_named(index, c("tile_id", "lat_n", "lon_w", "file", "status", "url_source",
                       "size_bytes", "message"))
  expect_identical(index$status, "converted")
  expect_true(file.exists(index$file))
  expect_equal(index$size_bytes, file.size(index$file))
  expect_identical(tools::md5sum(h5), original)
  expect_identical(attr(result$value, "traktr_layers")$aster_ged$params$product, "AG100")
  map <- attr(result$value, "aster_ged_map")
  expect_named(map, c("aoi_id", "tile_id"))
  expect_identical(map$tile_id, "40.-105")
  expect_identical(attr(result$value, "aster_ged_aoi"),
                   sf::st_transform(attr(df, "traktr_aoi"), 4326))

  cached <- fetch_aster_ged(df, out, progress = FALSE)
  expect_identical(attr(cached, "aster_ged_tiles")$status, "exists")
  expect_identical(calls, 0L)
})

test_that("AG100 recreates old TIFF grids and rejects incompatible cached bands", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  tif <- file.path(out, "ag100_v003_40_-105.tif")
  old <- terra::rast(nrows = 20, ncols = 20, nlyrs = 7, xmin = -105, xmax = -104,
                     ymin = 39, ymax = 40, crs = "EPSG:4326", vals = 0)
  names(old) <- c("elev", "ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd", "lwm", "nobs")
  terra::writeRaster(old, tif, datatype = "INT4S")
  fetch <- function() {
    fetch_aster_ged(ag100_fetch_aoi(), out, src = dirname(ag100_fixture()), progress = FALSE)
  }
  expect_message(fetched <- fetch(), "outdated grid, re-creating")
  expect_identical(attr(fetched, "aster_ged_tiles")$status, "converted")
  r <- terra::rast(tif)
  expect_equal(unname(as.vector(terra::ext(r))),
               c(-105, -104, 39, 40) + c(-1, 1, -1, 1) / (2 * 19))
  expect_equal(terra::res(r), rep(1 / 19, 2))
  expect_equal(unname(terra::values(r)[1L, "elev"]), 1511)
  expect_no_message(cached <- fetch())
  expect_identical(attr(cached, "aster_ged_tiles")$status, "exists")

  wrong_bands <- terra::rast(r)
  terra::values(wrong_bands) <- 0
  names(wrong_bands)[1L] <- "other"
  terra::writeRaster(wrong_bands, tif, datatype = "INT4S", overwrite = TRUE)
  expect_message(fetched <- fetch(), "outdated grid, re-creating")
  expect_identical(attr(fetched, "aster_ged_tiles")$status, "converted")
  expect_identical(names(terra::rast(tif)), names(old))
})

test_that("AG100 recreates cached TIFFs with a rotated geotransform", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  tif <- file.path(out, "ag100_v003_40_-105.tif")
  fetch <- function() {
    fetch_aster_ged(ag100_fetch_aoi(), out, src = dirname(ag100_fixture()), progress = FALSE)
  }
  fetch()
  scratch <- withr::local_tempdir()
  vrt <- file.path(scratch, "rotated.vrt")
  sf::gdal_utils("translate", tif, vrt, options = c("-of", "VRT"))
  xml <- readLines(vrt)
  line <- grep("<GeoTransform>", xml)
  gt <- as.numeric(strsplit(gsub(".*<GeoTransform>|</GeoTransform>.*", "", xml[line]), ",")[[1L]])
  gt[c(3L, 5L)] <- 0.01
  xml[line] <- paste0("  <GeoTransform>", paste(format(gt, digits = 17), collapse = ", "),
                      "</GeoTransform>")
  writeLines(xml, vrt)
  rotated <- file.path(scratch, "rotated.tif")
  sf::gdal_utils("translate", vrt, rotated)
  expect_true(terra::is.rotated(suppressWarnings(terra::rast(rotated))))
  file.copy(rotated, tif, overwrite = TRUE)

  expect_message(fetched <- fetch(), "outdated grid, re-creating")
  expect_identical(attr(fetched, "aster_ged_tiles")$status, "converted")
  expect_false(terra::is.rotated(terra::rast(tif)))
})

test_that("AG100 matches padded filenames and reads two-dimensional geolocation", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  df <- fetch_aster_ged(ag100_fetch_aoi(-74.8, -74.2, 4.2, 4.8), out,
                       src = dirname(ag100_fixture()), progress = FALSE)
  index <- attr(df, "aster_ged_tiles")
  expect_identical(index$tile_id, "5.-75")
  expect_identical(index$status, "converted")
  expect_equal(unname(as.vector(terra::ext(terra::rast(index$file)))),
               c(-75, -74, 4, 5) + c(-1, 1, -1, 1) / (2 * 19))
})

test_that("AG100 reports missing local sources in the index", {
  source <- withr::local_tempdir()
  out <- withr::local_tempdir()
  df <- fetch_aster_ged(ag100_fetch_aoi(), out, src = source, progress = FALSE)
  index <- attr(df, "aster_ged_tiles")
  expect_identical(index$status, "missing_local")
  expect_identical(index$file, NA_character_)
  expect_match(index$message, "No matching local HDF5")
})

test_that("AG100 discards or retains a downloaded HDF5 according to keep_h5", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  handles <- 0L
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() {
      handles <<- handles + 1L
      "dummy"
    },
    .ag100_download = function(url, dest_h5, handle) {
      expect_identical(handle, "dummy")
      file.copy(ag100_fixture(), dest_h5)
    }
  )
  df <- ag100_fetch_aoi()

  downloaded <- fetch_aster_ged(df, out, progress = FALSE)
  expect_identical(attr(downloaded, "aster_ged_tiles")$status, "downloaded")
  expect_length(list.files(out, "\\.h5$"), 0L)
  expect_identical(handles, 1L)

  kept <- fetch_aster_ged(df, out, keep_h5 = TRUE, overwrite = TRUE, progress = FALSE)
  expect_identical(attr(kept, "aster_ged_tiles")$status, "downloaded")
  expect_length(list.files(out, "\\.h5$"), 1L)
  expect_identical(handles, 2L)
})

test_that("a failed keep_h5 rename warns but keeps the tile consistent with its GeoTIFF", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() "dummy",
    .ag100_download = function(url, dest_h5, handle) file.copy(ag100_fixture(), dest_h5)
  )
  df <- ag100_fetch_aoi()
  # Block the retained HDF5 path with a directory so the rename fails
  dir.create(file.path(out, "ag100_v003_40_-105.h5"))

  expect_warning(
    res <- fetch_aster_ged(df, out, keep_h5 = TRUE, progress = FALSE),
    "could not be retained"
  )
  idx <- attr(res, "aster_ged_tiles")
  expect_identical(idx$status, "downloaded")
  expect_true(file.exists(idx$file))
  expect_match(idx$message, "could not be retained")
  expect_length(list.files(out, "^ag100-.*\\.h5$"), 0L)
})

test_that("AG100 keeps both HDF5 files when a replacement fails conversion", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  kept <- file.path(out, "ag100_v003_40_-105.h5")
  file.copy(ag100_fixture(), kept)
  original <- readBin(kept, "raw", n = file.size(kept))
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() "dummy",
    .ag100_download = function(url, dest_h5, handle) {
      writeBin(c(original[1:8], charToRaw("truncated HDF5")), dest_h5)
    }
  )

  expect_warning(
    df <- fetch_aster_ged(ag100_fetch_aoi(), out, keep_h5 = TRUE,
                          overwrite = TRUE, progress = FALSE),
    "ASTER GED tile .* failed"
  )
  index <- attr(df, "aster_ged_tiles")
  expect_identical(index$status, "error")
  expect_identical(readBin(kept, "raw", n = file.size(kept)), original)
  failed <- list.files(out, "^ag100_v003_40_-105-failed-.*\\.h5$")
  expect_length(failed, 1L)
  expect_match(index$message, failed, fixed = TRUE)
  expect_identical(readBin(file.path(out, failed), "raw", n = 22L),
                   c(original[1:8], charToRaw("truncated HDF5")))
})

test_that("AG100 retains a downloaded HDF5 after failed conversion only on request", {
  ag100_skip_without_hdf5()
  bad <- ag100_fixture(file.path("bad", "AG100.v003.40.-105.0001.h5"))
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() "dummy",
    .ag100_download = function(url, dest_h5, handle) file.copy(bad, dest_h5)
  )
  out <- withr::local_tempdir()
  for (keep in c(FALSE, TRUE)) {
    expect_warning(
      df <- fetch_aster_ged(ag100_fetch_aoi(), out, keep_h5 = keep, progress = FALSE),
      "Geolocation/Latitude"
    )
    index <- attr(df, "aster_ged_tiles")
    expect_identical(index$status, "error")
    expect_identical(index$file, NA_character_)
    files <- list.files(out)
    if (keep) {
      expect_identical(files, "ag100_v003_40_-105.h5")
      expect_match(index$message, files, fixed = TRUE)
      expect_identical(unname(tools::md5sum(file.path(out, files))), unname(tools::md5sum(bad)))
    } else {
      expect_length(files, 0L)
      expect_no_match(index$message, "retained")
    }
  }
})

test_that("the AG100 download wrapper stages success and removes partial files", {
  out <- withr::local_tempdir()
  dest <- file.path(out, "retained.h5")
  content <- readBin(ag100_fixture(), "raw", n = file.size(ag100_fixture()))
  fail_download <- TRUE
  url <- "https://example.test/AG100.h5"
  local_mocked_bindings(
    curl_download = function(url, destfile, quiet, mode, handle) {
      expect_identical(dirname(destfile), out)
      expect_false(file.exists(dest))
      expect_false(identical(destfile, dest))
      expect_identical(mode, "wb")
      expect_true(quiet)
      expect_identical(handle, "dummy")
      if (fail_download) {
        writeBin(content[1:100], destfile)
        stop("Interrupted transfer")
      }
      writeBin(content, destfile)
    },
    .package = "curl"
  )

  expect_error(.ag100_download(url, dest, "dummy"), "Interrupted transfer")
  expect_length(list.files(out), 0L)
  fail_download <- FALSE
  .ag100_download(url, dest, "dummy")
  expect_identical(list.files(out), basename(dest))
  expect_identical(readBin(dest, "raw", n = file.size(dest)), content)
})

test_that("AG100 reports missing remote tiles without downloading", {
  out <- withr::local_tempdir()
  local_mocked_bindings(
    .ag100_resolve = function(tiles, use_cmr) {
      tiles$status <- "missing_remote"
      tiles$url_source <- "cmr"
      tiles
    },
    .earthdata_handle = function(...) stop("Unexpected handle"),
    .ag100_download = function(...) stop("Unexpected download")
  )
  df <- fetch_aster_ged(ag100_fetch_aoi(), out, progress = FALSE)
  expect_identical(attr(df, "aster_ged_tiles")$status, "missing_remote")
  expect_length(list.files(out), 0L)
})

test_that("AG100 stops immediately on an HTML response or HTTP 401/403", {
  out <- withr::local_tempdir()
  withr::local_envvar(c(EARTHDATA_TOKEN = "mock-secret"))
  calls <- 0L
  response <- "html"
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() "dummy",
    .sleep = function(...) stop("Unexpected retry")
  )
  local_mocked_bindings(
    curl_download = function(url, destfile, quiet, mode, handle) {
      calls <<- calls + 1L
      if (response != "html") {
        writeLines("partial response", destfile)
        stop(paste("HTTP response code said error [mock-secret]: The requested URL returned error:", response))
      }
      writeLines("<html>login</html>", destfile)
    },
    .package = "curl"
  )
  for (response in c("html", "401", "403")) {
    err <- expect_error(fetch_aster_ged(ag100_fetch_aoi(), out, progress = FALSE),
                        "EARTHDATA_TOKEN", class = "ag100_auth_error")
    expect_no_match(conditionMessage(err), "mock-secret", fixed = TRUE)
  }
  expect_identical(calls, 3L)
  expect_length(list.files(out), 0L)
})

test_that("AG100 retries downloads and continues with the next tile after an error", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  calls <- 0L
  waits <- numeric()
  local_mocked_bindings(
    .ag100_resolve = ag100_mock_urls,
    .earthdata_handle = function() "dummy",
    .sleep = function(seconds) waits <<- c(waits, seconds),
    .ag100_download = function(url, dest_h5, handle) {
      calls <<- calls + 1L
      if (calls < 3L || grepl("-106", url, fixed = TRUE)) stop("Transient download failure")
      file.copy(ag100_fixture(), dest_h5)
    }
  )
  recovered <- fetch_aster_ged(ag100_fetch_aoi(), out, progress = FALSE)
  expect_identical(attr(recovered, "aster_ged_tiles")$status, "downloaded")
  expect_identical(calls, 3L)
  expect_equal(waits, c(2, 4))

  calls <- 0L
  waits <- numeric()
  expect_warning(
    mixed <- fetch_aster_ged(ag100_fetch_aoi(xmin = -105.2), out,
                             overwrite = TRUE, progress = FALSE),
    "Transient download failure"
  )
  index <- attr(mixed, "aster_ged_tiles")
  expect_identical(index$status, c("error", "downloaded"))
  expect_identical(index$message, c("Transient download failure", NA_character_))
  expect_identical(index$file[1], NA_character_)
  expect_identical(calls, 4L)
  expect_equal(waits, c(2, 4))
  expect_length(list.files(out), 1L)
})

test_that("AG100 uses sequential processing when future.apply is unavailable", {
  ag100_skip_without_hdf5()
  out <- withr::local_tempdir()
  local_mocked_bindings(.has_future_apply = function() FALSE)
  df <- ag100_fetch_aoi()
  src <- dirname(ag100_fixture())
  sequential <- fetch_aster_ged(df, out, src = src, progress = FALSE)
  expect_warning(
    fallback <- fetch_aster_ged(df, out, src = src, workers = 2, overwrite = TRUE, progress = FALSE),
    "future.apply.*unavailable"
  )
  expect_identical(attr(fallback, "aster_ged_tiles"), attr(sequential, "aster_ged_tiles"))
})
