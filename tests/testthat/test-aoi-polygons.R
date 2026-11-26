county_fixture <- function() testthat::test_path("fixtures", "counties5.gpkg")

test_that("FIPS codes are validated, padded and deduplicated", {
  expect_identical(.fips_normalize(c(1001, "01003", " 1005 ", 72127)),
                   c("01001", "01003", "01005", "72127"))
  expect_identical(.fips_normalize(1001.0), "01001")
  expect_error(.fips_normalize(1001.5), "1001.5", fixed = TRUE)
  expect_error(.fips_normalize(c("99001", "000001", "01abc")),
               "99001, 000001, 01abc", fixed = TRUE)
  expect_error(.fips_normalize(letters[1:11]), "and 1 more")
  expect_warning(out <- .fips_normalize(c(NA, " ", "1001")), "Removed 2")
  expect_identical(out, "01001")
  expect_warning(out <- .fips_normalize(NA), "Removed 1")
  expect_identical(out, character())
  expect_warning(out <- .fips_normalize(c("01001", "1001")), "1 duplicate")
  expect_identical(out, "01001")
})

test_that("read_fips_csv preserves codes and reports the available columns", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("x,county name", "8031,Denver", "02016,Aleutians West"), path)
  expect_identical(read_fips_csv(path), c("08031", "02016"))
  expect_error(read_fips_csv(path, col = "geoid"), "Available columns:.*x.*county name")
})

test_that("read_fips_csv warns about physically blank lines", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("x", "1001", "", "1003"), path)
  expect_warning(out <- read_fips_csv(path), "Removed 1 missing or blank FIPS code")
  expect_identical(out, c("01001", "01003"))
})

test_that("read_fips_csv skips blank lines before the header", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("", "x", "1001", "", "1003"), path)
  expect_warning(out <- read_fips_csv(path), "Removed 1 missing or blank FIPS code")
  expect_identical(out, c("01001", "01003"))

  empty <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("", " "), empty)
  expect_error(read_fips_csv(empty), "no header")
})

test_that("aoi_from_sf stores polygon AOIs and puts the IDs first", {
  counties <- sf::st_read(county_fixture(), quiet = TRUE)
  out <- aoi_from_sf(counties, "GEOID")
  expect_s3_class(out, "traktr_df")
  expect_s3_class(out, "tbl_df")
  expect_identical(names(out)[1], "GEOID")
  expect_length(names(out), ncol(counties) - 1L)
  aoi <- attr(out, "traktr_aoi")
  expect_identical(sf::st_crs(aoi)$epsg, 4326L)
  expect_identical(names(aoi), c("GEOID", "geom"))
  expect_identical(
    attr(out, "traktr_aoi_params"),
    list(method = "polygon", deduplicate = FALSE, crs_source = 4269L, id_col = "GEOID")
  )
  expect_named(aoi_from_sf(counties, "GEOID", keep = FALSE), "GEOID")
  expect_error(aoi_from_sf(counties[c(1, 1), ], "GEOID"), "unique")
  counties$GEOID[1] <- NA_character_
  expect_error(aoi_from_sf(counties, "GEOID"), "non-missing")
})

test_that("aoi_from_sf rejects wide outer rings and identifies their AOIs", {
  ring <- matrix(c(-100, 10, 100, 10, 100, 20, -100, 20, -100, 10),
                 ncol = 2L, byrow = TRUE)
  aoi <- sf::st_sf(GEOID = "wide-ring",
                   geom = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326))
  expect_error(aoi_from_sf(aoi, "GEOID"), "180 degrees.*wide-ring")
})

test_that("county AOIs preserve request order and Aleutian polygon parts", {
  out <- aoi_us_counties(c(8031, 2016, 15001), src = county_fixture())
  expect_s3_class(out, "traktr_df")
  expect_named(out, c("GEOID", "STATEFP", "COUNTYFP", "NAME", "NAMELSAD",
                      "STUSPS", "ALAND", "AWATER", "region"))
  expect_identical(out$GEOID, c("08031", "02016", "15001"))
  expect_identical(out$region, c("conus", "ak", "hi"))
  params <- attr(out, "traktr_aoi_params")
  expect_identical(params[c("source", "year", "resolution", "territories")],
                   list(source = "census_cb", year = 2020, resolution = "500k",
                        territories = "PR"))
  aleutians <- sf::st_geometry(attr(out, "traktr_aoi"))[[2L]]
  expect_s3_class(aleutians, "MULTIPOLYGON")
  expect_length(aleutians, 2L)
  x_ranges <- lapply(aleutians, function(part) range(part[[1L]][, 1L]))
  expect_true(all(vapply(x_ranges, diff, numeric(1)) < 10))
  expect_gt(x_ranges[[1L]][1L], 170)
  expect_lt(x_ranges[[2L]][2L], -170)
})

test_that("county AOIs filter territories and diagnose missing FIPS and vintage", {
  src <- county_fixture()
  out <- aoi_us_counties(territories = "none", src = src)
  expect_identical(out$GEOID, c("02016", "02185", "06071", "08031", "15001"))
  expect_error(aoi_us_counties(72127, territories = "none", src = src),
               "excluded.*72127")
  expect_warning(out <- aoi_us_counties(c(1001, 8031), src = src), "missing.*01001")
  expect_identical(attr(out, "fips_missing"), "01001")
  expect_identical(out$GEOID, "08031")
  expect_error(aoi_us_counties(1001, src = src), "No requested FIPS.*01001")
  expect_error(aoi_us_counties("09110", year = 2020, src = src), "2022")
  expect_error(aoi_us_counties("09015", year = 2022, src = src), "2020")
})

test_that("optional Census columns are filled when absent", {
  counties <- sf::st_read(county_fixture(), quiet = TRUE)
  counties$STUSPS <- NULL
  counties$NAMELSAD <- NULL
  src <- withr::local_tempfile(fileext = ".gpkg")
  sf::st_write(counties, src, layer = "counties", quiet = TRUE)

  out <- aoi_us_counties(c("08031", "72127"), src = src)
  expect_identical(out$STUSPS, c("CO", "PR"))
  expect_identical(out$NAMELSAD, rep(NA_character_, 2L))
  expect_identical(out$region, c("conus", "pr"))
})

test_that("the Census download helper is used only when src is NULL", {
  calls <- list()
  src <- county_fixture()
  local_mocked_bindings(.census_county_file = function(year, resolution, cache, quiet) {
    calls[[length(calls) + 1L]] <<- list(year, resolution, cache, quiet)
    src
  })
  remote <- aoi_us_counties(year = 2023, resolution = "5m", cache = FALSE, quiet = TRUE)
  expect_identical(calls, list(list(2023, "5m", FALSE, TRUE)))
  local <- aoi_us_counties(year = 2023, resolution = "5m", src = src)
  expect_identical(calls, list(list(2023, "5m", FALSE, TRUE)))
  expect_identical(local, remote)
})

test_that("an unreadable download is never added to the Census cache", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    writeLines("not a zip", destfile)
  }, .package = "curl")

  expect_error(aoi_us_counties(quiet = TRUE), "downloaded county archive is unreadable")
  directory <- file.path(tools::R_user_dir("traktr", "cache"), "census")
  expect_length(list.files(directory, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("a corrupt Census cache is refreshed and a good cache survives a bad refresh", {
  gdal <- sf::sf_extSoftVersion()[["GDAL"]]
  skip_if(utils::compareVersion(gdal, "3.1.0") < 0L,
          "Writing zipped shapefiles requires GDAL 3.1")
  archive <- withr::local_tempfile(fileext = ".shp.zip")
  counties <- sf::st_read(county_fixture(), quiet = TRUE)
  sf::st_write(counties, archive, layer = "counties",
               driver = "ESRI Shapefile", quiet = TRUE)

  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  directory <- file.path(tools::R_user_dir("traktr", "cache"), "census")
  dir.create(directory, recursive = TRUE)
  archive_name <- "cb_2020_us_county_500k.zip"
  cached <- file.path(directory, archive_name)
  writeLines("not a zip", cached)
  downloads <- 0L
  download_ok <- TRUE
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    downloads <<- downloads + 1L
    if (download_ok) {
      file.copy(archive, destfile)
    } else {
      writeLines("not a zip", destfile)
    }
  }, .package = "curl")

  expect_warning(out <- aoi_us_counties(8031, quiet = TRUE),
                 "cached county archive.*unreadable")
  expect_identical(out$GEOID, "08031")
  expect_identical(downloads, 1L)
  good_checksum <- unname(tools::md5sum(archive))
  expect_identical(unname(tools::md5sum(cached)), good_checksum)
  expect_identical(aoi_us_counties(8031, quiet = TRUE), out)
  expect_identical(downloads, 1L)

  download_ok <- FALSE
  expect_error(aoi_us_counties(cache = FALSE, quiet = TRUE),
               "downloaded county archive is unreadable")
  expect_identical(downloads, 2L)
  expect_identical(unname(tools::md5sum(cached)), good_checksum)
  expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), archive_name)
})

test_that("AOI deduplication uses the active geometry column", {
  counties <- sf::st_read(county_fixture(), quiet = TRUE)[c(1, 2, 1), ]
  counties$GEOID[3L] <- "test-copy"
  aoi <- attr(aoi_from_sf(counties, "GEOID"), "traktr_aoi")
  expect_identical(attr(aoi, "sf_column"), "geom")
  deduped <- .aoi_deduplicate(aoi, "GEOID")
  expect_identical(deduped$aoi_unique$GEOID, c("08031", "06071"))
  expect_identical(deduped$map_back$.aoi_row[match(aoi$GEOID, deduped$map_back$GEOID)],
                   c(1L, 2L, 1L))
  kept <- .aoi_deduplicate(aoi, "GEOID", deduplicate = FALSE)
  expect_identical(kept$aoi_unique$GEOID, aoi$GEOID)
  expect_identical(kept$map_back$GEOID, aoi$GEOID)
  expect_identical(kept$map_back$.aoi_row, 1:3)
})

test_that("geometry deduplication retains the previous S2 normalization of valid rings", {
  withr::local_options(sf_use_s2 = TRUE)
  ring <- matrix(c(1, 2, 3, 2, 3, 4, 1, 4, 1, 2), ncol = 2L, byrow = TRUE)
  rings <- list(ring, ring[5:1, ], ring[c(2:5, 2), ])
  polygons <- lapply(rings, function(ring) sf::st_polygon(list(ring)))
  aoi <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(polygons, crs = 4326))
  expect_true(all(sf::st_is_valid(aoi)))

  deduped <- .aoi_deduplicate(aoi, "id", deduplicate = TRUE)
  # Baseline from f86b17a with S2 enabled.
  expect_equal(nrow(deduped$aoi_unique), 2L)
})
