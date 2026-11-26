aster_summary_polygon <- function(xmin, xmax, ymin, ymax) {
  sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                            crs = 4326))
}

aster_summary_data <- function(dir) {
  aoi <- sf::st_sf(
    id = c("a", "b"), label = c("west", "east"),
    geometry = c(aster_summary_polygon(0, 0.5, 0, 1), aster_summary_polygon(0.5, 1, 0, 1))
  )
  x <- aoi_from_sf(aoi, id_col = "id")
  r <- terra::rast(xmin = 0, xmax = 1, ymin = 0, ymax = 1, nrows = 1, ncols = 4,
                   crs = "EPSG:4326", nlyrs = 7)
  values <- cbind(
    elev = c(0, 20, NA, 8), ndvi_mean = c(10, 30, 40, 50), ndvi_sd = c(2, 6, 4, 8),
    temp_mean = c(27315, 27515, 27015, 28015), temp_sd = c(100, 300, 200, 400),
    lwm = c(1, 0, 1, NA), nobs = c(2, 4, 4, 6)
  )
  terra::values(r) <- values
  names(r) <- colnames(values)
  path <- file.path(dir, "one.tif")
  terra::writeRaster(r, path, datatype = "INT4S", NAflag = -9999)
  attr(x, "aster_ged_tiles") <- data.frame(tile_id = "one", file = path, status = "converted")
  attr(x, "aster_ged_map") <- data.frame(aoi_id = x$id, tile_id = "one")
  x
}

test_that("the public fetch-to-attach workflow produces the exact column names", {
  skip_if_not("HDF5" %in% terra::gdal(drivers = TRUE)$name, "GDAL has no HDF5 driver")
  aoi <- sf::st_sf(
    id = c("colorado", "colombia"),
    geometry = c(aster_summary_polygon(-104.9, -104.1, 39.1, 39.9),
                 aster_summary_polygon(-74.9, -74.1, 4.1, 4.9))
  )
  input <- fetch_aster_ged(
    aoi_from_sf(aoi, id_col = "id"), out_dir = withr::local_tempdir(),
    src = test_path("fixtures", "ag100"), use_cmr = FALSE, progress = FALSE
  )
  expect_identical(attr(input, "aster_ged_tiles")$status, rep("converted", 2))

  county <- c(
    "ASTER_county_elevation_mean", "ASTER_county_elevation_min", "ASTER_county_elevation_max",
    "ASTER_county_elevation_sd", "ASTER_county_elevation_cin", "ASTER_county_elevation_sum",
    "ASTER_county_elevation_missing_perc",
    "ASTER_ndvi_mean_mean", "ASTER_ndvi_mean_min", "ASTER_ndvi_mean_max",
    "ASTER_ndvi_mean_sd", "ASTER_ndvi_mean_missing_perc",
    "ASTER_ndvi_sd_mean", "ASTER_ndvi_sd_min", "ASTER_ndvi_sd_max",
    "ASTER_ndvi_sd_sd", "ASTER_ndvi_sd_missing_perc",
    "ASTER_temp_mean_mean", "ASTER_temp_mean_min", "ASTER_temp_mean_max",
    "ASTER_temp_mean_sd", "ASTER_temp_mean_missing_perc",
    "ASTER_temp_sd_mean", "ASTER_temp_sd_min", "ASTER_temp_sd_max",
    "ASTER_temp_sd_sd", "ASTER_temp_sd_missing_perc",
    "ASTER_observation_number_mean", "ASTER_observation_number_min", "ASTER_observation_number_max",
    "ASTER_observation_number_sd", "ASTER_observation_number_lt_4",
    "ASTER_lwm_water_perc", "ASTER_lwm_water_cin", "ASTER_lwm_water_missing_perc"
  )
  ecolozip <- c(
    "ASTER_zip_elevation_mean", "ASTER_zip_elevation_min", "ASTER_zip_elevation_max",
    "ASTER_zip_elevation_sd", "ASTER_zip_elevation_cin", "ASTER_zip_elevation_sum",
    "ASTER_zip_elevation_missing_perc",
    "ASTER_ndvi_mean_mean", "ASTER_ndvi_mean_min", "ASTER_ndvi_mean_max",
    "ASTER_ndvi_mean_sd", "ASTER_ndvi_mean_missing_perc",
    "ASTER_ndvi_sd_mean", "ASTER_ndvi_sd_min", "ASTER_ndvi_sd_max",
    "ASTER_ndvi_sd_sd", "ASTER_ndvi_sd_missing_perc",
    "ASTER_temp_mean_mean", "ASTER_temp_mean_min", "ASTER_temp_mean_max",
    "ASTER_temp_mean_sd", "ASTER_temp_mean_missing_perc",
    "ASTER_temp_sd_mean", "ASTER_temp_sd_min", "ASTER_temp_sd_max",
    "ASTER_temp_sd_sd", "ASTER_temp_sd_missing_perc",
    "ASTER_observation_number_mean", "ASTER_observation_number_min", "ASTER_observation_number_max",
    "ASTER_observation_number_sd", "ASTER_observation_number_lt_4",
    "ASTER_lwm_water_perc", "ASTER_lwm_water_cin", "ASTER_lwm_water_missing_perc"
  )
  traktr <- c(
    "aster_elev_mean", "aster_elev_min", "aster_elev_max", "aster_elev_sd",
    "aster_elev_cin", "aster_elev_sum", "aster_elev_missing_perc",
    "aster_ndvi_mean_mean", "aster_ndvi_mean_min", "aster_ndvi_mean_max",
    "aster_ndvi_mean_sd", "aster_ndvi_mean_missing_perc",
    "aster_ndvi_sd_mean", "aster_ndvi_sd_min", "aster_ndvi_sd_max",
    "aster_ndvi_sd_sd", "aster_ndvi_sd_missing_perc",
    "aster_temp_mean_mean", "aster_temp_mean_min", "aster_temp_mean_max",
    "aster_temp_mean_sd", "aster_temp_mean_missing_perc",
    "aster_temp_sd_mean", "aster_temp_sd_min", "aster_temp_sd_max",
    "aster_temp_sd_sd", "aster_temp_sd_missing_perc",
    "aster_nobs_mean", "aster_nobs_min", "aster_nobs_max", "aster_nobs_sd", "aster_nobs_lt_4",
    "aster_lwm_water_perc", "aster_lwm_cin", "aster_lwm_missing_perc"
  )
  core <- list(county = county, ecolozip = ecolozip, traktr = traktr)
  extras <- list(
    county = c(
      "ASTER_county_elevation_median", "ASTER_ndvi_mean_median", "ASTER_ndvi_sd_median",
      "ASTER_temp_mean_median", "ASTER_temp_sd_median", "ASTER_observation_number_median",
      "n_cells", "county_area_km2_geodesic"
    ),
    ecolozip = c(
      "ASTER_zip_elevation_median", "ASTER_ndvi_mean_median", "ASTER_ndvi_sd_median",
      "ASTER_temp_mean_median", "ASTER_temp_sd_median", "ASTER_observation_number_median",
      "n_cells", "county_area_km2_geodesic"
    ),
    traktr = c(
      "aster_elev_median", "aster_ndvi_mean_median", "aster_ndvi_sd_median",
      "aster_temp_mean_median", "aster_temp_sd_median", "aster_nobs_median",
      "aster_n_cells", "aoi_area_km2"
    )
  )
  dictionary <- .aster_ged_dictionary()
  for (naming in names(core)) {
    result <- attach_aster_ged(input, naming = naming, extras = character(), progress = FALSE)
    expect_length(core[[naming]], 35)
    expect_named(result, c(names(input), core[[naming]]))
    expect_true(all(is.finite(result[[core[[naming]][1L]]])))
    expect_equal(result[[core[[naming]][2L]]], rep(1533, 2))
    expect_equal(result[[core[[naming]][3L]]], rep(1698, 2))
    result <- attach_aster_ged(input, naming = naming, progress = FALSE)
    expected_names <- c(core[[naming]], extras[[naming]])
    expect_named(result, c(names(input), expected_names))
    expect_equal(result[[extras[[naming]][7L]]], rep(252, 2))
    expect_identical(dictionary[[paste0("name_", naming)]], expected_names)
    expect_false(anyDuplicated(expected_names) > 0)
  }
})

test_that("band subsets use canonical order and cell counts use the first requested band", {
  x <- aster_summary_data(withr::local_tempdir())
  out <- attach_aster_ged(x, bands = c("lwm", "nobs"), naming = "traktr",
                          extras = c("n_cells", "median"), progress = FALSE)
  expect_named(out, c(
    "id", "label", "aster_nobs_mean", "aster_nobs_min", "aster_nobs_max", "aster_nobs_sd",
    "aster_nobs_lt_4", "aster_lwm_water_perc", "aster_lwm_cin", "aster_lwm_missing_perc",
    "aster_nobs_median", "aster_n_cells"
  ))
  expect_equal(out$aster_n_cells, c(2, 1))
  expect_error(attach_aster_ged(x, bands = "elevation"), "bands.*elev")
})

test_that("ASTER statistics use population SD, area shares and observation counts", {
  x <- aster_summary_data(withr::local_tempdir())
  out <- attach_aster_ged(x, naming = "traktr", progress = FALSE)
  expect_equal(out$aster_elev_mean, c(10, 8))
  expect_equal(out$aster_elev_sd, c(10, 0))
  expect_equal(out$aster_elev_median, c(0, 8))
  expect_equal(out$aster_elev_cin, c(50, 0))
  expect_identical(out$aster_elev_sum, rep(NA_real_, 2))
  expect_equal(out$aster_elev_missing_perc, c(0, 50))
  expect_identical(out$aster_nobs_lt_4, c(1L, 0L))
  expect_equal(out$aster_lwm_water_perc, c(50, 100))
  expect_equal(out$aster_lwm_cin, c(100, 50))
  expect_equal(out$aster_lwm_missing_perc, c(0, 50))
  expect_equal(out$aster_n_cells, c(2, 1))

  quirk <- attach_aster_ged(x, naming = "traktr", ecolozip_na_quirk = TRUE, progress = FALSE)
  expect_equal(quirk$aster_elev_missing_perc, c(NA, 50))
  expect_equal(quirk$aster_lwm_missing_perc, c(NA, 50))
  expect_true(all(is.na(quirk$aster_ndvi_mean_missing_perc)))
  expect_identical(quirk$aster_lwm_cin, out$aster_lwm_cin)

  partial <- attr(x, "traktr_aoi")
  sf::st_geometry(partial)[1L] <- aster_summary_polygon(0.125, 0.5, 0, 1)
  attr(x, "traktr_aoi") <- partial
  weighted <- attach_aster_ged(x, bands = c("elev", "lwm"), progress = FALSE)
  expect_equal(weighted$ASTER_county_elevation_cin[1L], 100 / 3)
  expect_equal(weighted$ASTER_lwm_water_perc[1L], 100 / 3)
})

test_that("conversion changes the units of location and spread statistics", {
  x <- aster_summary_data(withr::local_tempdir())
  raw <- attach_aster_ged(x, naming = "traktr", progress = FALSE)
  out <- attach_aster_ged(x, naming = "traktr", convert = TRUE, progress = FALSE)
  rows_changed <- grepl("^aster_(ndvi|temp)_.*_(mean|min|max|sd|median)$", names(raw))
  expect_identical(out[!rows_changed], raw[!rows_changed], ignore_attr = TRUE)
  for (band in c("ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd")) {
    for (stat in c("mean", "min", "max", "sd", "median")) {
      name <- paste("aster", band, stat, sep = "_")
      expected <- raw[[name]] / 100
      if (band == "temp_mean" && stat != "sd") expected <- expected - 273.15
      expect_equal(out[[name]], expected)
    }
  }

  raw_units <- .aster_ged_dictionary()$unit
  converted <- .aster_ged_dictionary(convert = TRUE)
  expect_identical(raw_units[match(c("aster_ndvi_mean_sd", "aster_temp_mean_mean"),
                                  converted$name_traktr)], c("NDVI\u00d7100", "cK"))
  expect_identical(converted$unit[match(
    c("aster_ndvi_mean_sd", "aster_temp_mean_mean", "aster_temp_mean_sd", "aster_temp_sd_mean",
      "aster_temp_mean_missing_perc", "aster_nobs_lt_4", "aoi_area_km2"),
    converted$name_traktr
  )], c("NDVI", "\u00b0C", "K", "K", "%", "flag 0/1", "km\u00b2"))
})

test_that("missing tiles are reported according to status and missing coverage stays NA", {
  x <- aster_summary_data(withr::local_tempdir())
  tiles <- attr(x, "aster_ged_tiles")
  attr(x, "aster_ged_map") <- rbind(attr(x, "aster_ged_map"), data.frame(
    aoi_id = "a", tile_id = c("local", "failed", "ocean")
  ))
  attr(x, "aster_ged_tiles") <- rbind(tiles, data.frame(
    tile_id = c("local", "failed", "ocean"), file = NA_character_,
    status = c("missing_local", "error", "missing_remote")
  ))
  expect_error(attach_aster_ged(x), "local.*failed")
  expect_warning(out <- attach_aster_ged(x, allow_missing = TRUE, progress = FALSE), "local.*failed")
  expect_identical(attr(out, "traktr_layers")$aster_ged_summary$index, tiles)
  index <- attr(x, "aster_ged_tiles")
  attr(x, "aster_ged_tiles") <- index[c(1L, 4L), ]
  expect_no_warning(attach_aster_ged(x, progress = FALSE))
  attr(x, "aster_ged_tiles") <- index[4L, ]
  expect_no_warning(ocean <- attach_aster_ged(x, naming = "traktr", progress = FALSE))
  expect_true(all(is.na(ocean$aster_elev_missing_perc)))
  expect_true(all(is.na(ocean$aster_lwm_cin)))
  expect_true(all(is.na(ocean$aster_nobs_lt_4)))
  expect_equal(ocean$aster_n_cells, c(0, 0))
  expect_equal(nrow(attr(ocean, "traktr_layers")$aster_ged_summary$index), 0)
  expect_error(attach_aster_ged(data.frame(id = "a")), "run fetch_aster_ged\\(\\) first")

  attr(x, "aster_ged_tiles") <- tiles
  unlink(tiles$file)
  expect_error(attach_aster_ged(x), "Unavailable ASTER GED tiles:.*one")
})

test_that("attachment uses the AOI override recorded by fetch", {
  skip_if_not("HDF5" %in% terra::gdal(drivers = TRUE)$name, "GDAL has no HDF5 driver")
  colorado <- aster_summary_polygon(-104.9, -104.1, 39.1, 39.9)
  x <- aoi_from_sf(sf::st_sf(id = c("shared", "absent"),
                            geometry = c(colorado, colorado)), id_col = "id")
  colombia <- aster_summary_polygon(-74.9, -74.1, 4.1, 4.9)
  override <- sf::st_sf(
    id = c("shared", "unused"), label = c("colombia", "outside fixtures"),
    geometry = c(colombia, aster_summary_polygon(0.1, 0.9, 0.1, 0.9))
  )
  fetched <- fetch_aster_ged(
    x, out_dir = withr::local_tempdir(), src = test_path("fixtures", "ag100"),
    use_cmr = FALSE, progress = FALSE, aoi = sf::st_transform(override, 3857)
  )
  stored_aoi <- attr(fetched, "aster_ged_aoi")
  expect_named(stored_aoi, c("id", "geometry"))
  expect_identical(sf::st_crs(stored_aoi)$epsg, 4326L)
  expect_identical(attr(fetched, "traktr_aoi"), attr(x, "traktr_aoi"))
  expect_identical(attr(fetched, "aster_ged_tiles")$status, c("converted", "missing_local"))

  out <- attach_aster_ged(fetched, bands = "elev", progress = FALSE)
  expect_identical(out$id, x$id)
  expect_true(is.finite(out$ASTER_county_elevation_mean[1L]))
  expect_true(is.na(out$ASTER_county_elevation_mean[2L]))
  expect_equal(out$county_area_km2_geodesic, c(.aster_ged_area(override[1L, ]), NA_real_))
  expect_identical(attr(out, "traktr_layers")$aster_ged_summary$index$tile_id, "5.-75")
})

test_that("AG100 fetch and attach agree for split and unsplit antimeridian polygons", {
  dir <- withr::local_tempdir()
  src <- withr::local_tempdir()
  h <- 1 / (2 * 19)
  for (lon in c(179, -180)) {
    r <- terra::rast(nrows = 20, ncols = 20, nlyrs = 7, xmin = lon - h,
                     xmax = lon + 1 + h, ymin = -h, ymax = 1 + h,
                     crs = "EPSG:4326", vals = if (lon == 179) 10 else 20)
    names(r) <- c("elev", "ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd", "lwm", "nobs")
    path <- file.path(dir, sprintf("ag100_v003_1_%d.tif", lon))
    terra::writeRaster(r, path, datatype = "INT4S")
  }
  aoi <- function(geometry) sf::st_sf(id = "a", geometry = geometry)
  df <- aoi_from_sf(aoi(aster_summary_polygon(179.4, 180, 0.25, 0.75)), id_col = "id")
  summarize <- function(geometry) {
    fetched <- fetch_aster_ged(df, dir, src = src, aoi = aoi(geometry), progress = FALSE)
    tiles <- attr(fetched, "aster_ged_tiles")
    expect_setequal(tiles$tile_id, c("1.179", "1.-180"))
    expect_identical(tiles$status, rep("exists", 2))
    attach_aster_ged(fetched, bands = "elev", naming = "traktr", extras = "n_cells",
                     progress = FALSE, cache = FALSE)
  }

  ring <- rbind(c(179.4, 0.25), c(-179.7, 0.25), c(-179.7, 0.75),
                c(179.4, 0.75), c(179.4, 0.25))
  unsplit <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
  split <- sf::st_combine(c(aster_summary_polygon(179.4, 180, 0.25, 0.75),
                            aster_summary_polygon(-180, -179.7, 0.25, 0.75)))
  actual <- summarize(unsplit)
  expected <- summarize(split)
  expect_equal(actual$aster_elev_mean, 40 / 3, tolerance = 1e-8)
  expect_equal(actual$aster_elev_mean, expected$aster_elev_mean, tolerance = 1e-9)
  expect_identical(actual$aster_n_cells, expected$aster_n_cells)
  expect_equal(actual$aster_n_cells, 180)

  centre <- sf::st_sfc(sf::st_point(c(179.9995, 0.52)), crs = 4326)
  local_crs <- "+proj=aeqd +lat_0=0.52 +lon_0=179.9995 +datum=WGS84 +units=m"
  buffer <- sf::st_transform(sf::st_buffer(sf::st_transform(centre, local_crs), 200), 4326)
  wrapped <- sf::st_wrap_dateline(buffer, options = "WRAPDATELINE=YES", quiet = TRUE)
  actual <- summarize(buffer)
  expected <- summarize(wrapped)
  expect_identical(as.character(sf::st_geometry_type(wrapped)), "MULTIPOLYGON")
  expect_equal(actual$aster_elev_mean, expected$aster_elev_mean, tolerance = 1e-9)
  expect_identical(actual$aster_n_cells, expected$aster_n_cells)
  expect_equal(actual$aster_n_cells, 1)
})

test_that("the join preserves repeated IDs, row order, attributes and geometry", {
  x <- aster_summary_data(withr::local_tempdir())
  x <- x[c(2L, 1L, 1L), ]
  x$geometry <- sf::st_geometry(attr(x, "traktr_aoi"))[c(2L, 1L, 1L)]
  x$sequence <- 1:3
  attr(x, "metadata") <- "preserved"
  previous_s2 <- getOption("sf_use_s2")
  on.exit(options(sf_use_s2 = previous_s2))
  suppressMessages(sf::sf_use_s2(TRUE))
  expected_area <- as.numeric(sf::st_area(x$geometry)) / 1e6
  suppressMessages(sf::sf_use_s2(FALSE))

  result <- withVisible(attach_aster_ged(x, progress = FALSE))
  out <- result$value
  expect_false(result$visible)
  expect_false(sf::sf_use_s2())
  expect_identical(class(out), class(x))
  expect_identical(out$id, c("b", "a", "a"))
  expect_identical(out$sequence, 1:3)
  expect_equal(out$ASTER_county_elevation_mean, c(8, 10, 10))
  expect_identical(tail(names(out), 1), "geometry")
  expect_identical(out$geometry, x$geometry)
  expect_equal(out$county_area_km2_geodesic, expected_area)
  for (attr_name in c("traktr_aoi", "traktr_aoi_params", "aster_ged_tiles", "aster_ged_map", "metadata")) {
    expect_identical(attr(out, attr_name), attr(x, attr_name))
  }
  layer <- attr(out, "traktr_layers")$aster_ged_summary
  expect_identical(layer$params, list(
    bands = c("elev", "ndvi_mean", "ndvi_sd", "temp_mean", "temp_sd", "lwm", "nobs"),
    naming = "county", extras = c("median", "n_cells", "area"),
    convert = FALSE, ecolozip_na_quirk = FALSE
  ))
  expect_identical(attach_aster_ged(out, progress = FALSE), out)
})

test_that("cache is optional and produces identical results from stored moments", {
  dir <- withr::local_tempdir()
  x <- aster_summary_data(dir)
  cache_dir <- file.path(dir, "_moments")
  direct <- attach_aster_ged(x, cache = FALSE, progress = FALSE)
  expect_false(dir.exists(cache_dir))
  first <- attach_aster_ged(x, progress = FALSE)
  expect_length(list.files(cache_dir, pattern = "[.]rds$"), 1)
  local_mocked_bindings(.tile_moments = function(...) stop("cache was not read"))
  expect_identical(attach_aster_ged(x, progress = FALSE), first)
  expect_identical(first, direct)
})

test_that("the recorded fetch directory takes precedence for the cache", {
  dir <- withr::local_tempdir()
  x <- aster_summary_data(dir)
  out_dir <- withr::local_tempdir()
  x <- .set_traktr_layer(x, "aster_ged", attr(x, "aster_ged_tiles"), list(out_dir = out_dir))

  attach_aster_ged(x, progress = FALSE)
  expect_false(dir.exists(file.path(dir, "_moments")))
  expect_length(list.files(file.path(out_dir, "_moments"), pattern = "[.]rds$"), 1)
})

test_that("the area covers every AOI row of an ID, like the combined moments", {
  x <- aster_summary_data(withr::local_tempdir())
  aoi <- attr(x, "traktr_aoi")
  aoi[[attr(x, "traktr_aoi_params")$id_col]] <- "a"
  attr(x, "aster_ged_aoi") <- aoi
  attr(x, "aster_ged_map") <- data.frame(aoi_id = "a", tile_id = "one")
  out <- attach_aster_ged(x, bands = "elev", progress = FALSE)
  expect_equal(out$county_area_km2_geodesic[out$id == "a"],
               rep(sum(.aster_ged_area(aoi)), sum(out$id == "a")))
  expect_true(all(is.na(out$county_area_km2_geodesic[out$id == "b"])))
})
