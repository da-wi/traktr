county_export_fixture <- function(dir, convert = FALSE, naming = "county") {
  polygon <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1, ymax = 1), crs = 4326))
  aoi <- sf::st_sf(GEOID = "01001", STATEFP = "01", COUNTYFP = "001", NAME = "Autauga",
                    NAMELSAD = "Autauga County", STUSPS = "AL", ALAND = 1000,
                    AWATER = 100, region = "conus", geometry = polygon)
  df <- aoi_from_sf(aoi, id_col = "GEOID")
  attr(df, "traktr_aoi_params") <- c(attr(df, "traktr_aoi_params"),
                                    list(source = "census_cb", year = 2020L, resolution = "500k"))
  attr(df, "fips_missing") <- "01999"
  df$custom <- 1.23456789012345
  df$shape <- polygon

  raster <- terra::rast(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                        nrows = 1, ncols = 2, nlyrs = 2, crs = "EPSG:4326")
  terra::values(raster) <- cbind(elev = c(0, 20), ndvi_mean = c(10, 30))
  names(raster) <- c("elev", "ndvi_mean")
  tile <- file.path(dir, "ag100.tif")
  terra::writeRaster(raster, tile, overwrite = TRUE)
  tiles <- data.frame(tile_id = "1.0", status = "converted", file = tile,
                      size_bytes = file.size(tile))
  attr(df, "aster_ged_tiles") <- tiles
  attr(df, "aster_ged_map") <- data.frame(aoi_id = df$GEOID, tile_id = tiles$tile_id)
  df <- .set_traktr_layer(df, "aster_ged", tiles, list(
    product = "AG100", version = "003", concept_id = "C2763266348-LPCLOUD",
    retrieved = as.POSIXct("2026-01-01 12:34:56", tz = "UTC"), out_dir = dir
  ))
  attach_aster_ged(df, bands = c("elev", "ndvi_mean"), convert = convert,
                  naming = naming, progress = FALSE, cache = FALSE)
}

test_that("CSV preserves column order, FIPS codes, precision and empty fields", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  stem <- file.path(dir, "new", "summary")
  result <- withVisible(write_county_summary(df, paste0(stem, ".csv"),
                                              dictionary = FALSE, provenance = FALSE))
  expect_false(result$visible)
  expect_identical(result$value, c(csv = paste0(stem, ".csv")))
  csv <- utils::read.csv(result$value[["csv"]], colClasses = "character")
  expect_named(csv, names(df)[names(df) != "shape"])
  expect_identical(csv$GEOID, "01001")
  expect_match(readLines(result$value[["csv"]])[2L], '^"01001",')
  expect_identical(csv$ASTER_county_elevation_sum, "")
  expect_equal(as.numeric(csv$custom), df$custom, tolerance = 1e-14)

  write_county_summary(df, paste0(stem, ".parquet"), na = "MISSING", overwrite = TRUE,
                       dictionary = FALSE, provenance = FALSE)
  csv <- utils::read.csv(result$value[["csv"]], colClasses = "character")
  expect_identical(csv$ASTER_county_elevation_sum, "MISSING")
  expect_false(file.exists(paste0(stem, ".parquet.csv")))
})

test_that("the dictionary covers every column and follows naming and units", {
  dir <- withr::local_tempdir()
  for (convert in c(FALSE, TRUE)) {
    df <- county_export_fixture(dir, convert = convert, naming = "traktr")
    paths <- write_county_summary(df, file.path(dir, paste0("summary-", convert)))
    dict <- utils::read.csv(paths[["dictionary"]], colClasses = "character")
    expect_named(dict, c("name", "source", "band", "stat", "unit", "description", "ecolozip_name"))
    expect_identical(dict$name, setdiff(names(df), "shape"))
    expect_identical(dict$source[dict$name == "custom"], "input")
    expect_identical(dict$source[dict$name == "ALAND"], "US Census cartographic boundary 2020")
    expect_identical(dict$unit[dict$name == "ALAND"], "m\u00b2")
    expect_identical(dict$ecolozip_name[dict$name == "aster_elev_mean"], "ASTER_zip_elevation_mean")
    elev_names <- paste0("aster_elev_", c("mean", "min", "max", "sd", "cin", "sum",
                                          "missing_perc", "median"))
    elev_rows <- dict$description[match(elev_names, dict$name)]
    expect_false(anyNA(elev_rows))
    expect_true(all(grepl("elevation band was entirely missing in the converted",
                          elev_rows, fixed = TRUE)))
    expect_true(all(grepl("no ASTER elevation for any Hawaii ZIP code",
                          elev_rows, fixed = TRUE)))
    expect_true(all(grepl("ASTER GDEM artefacts: isolated pits near water bodies",
                          dict$description[dict$name %in% c("aster_elev_min", "aster_elev_max")],
                          fixed = TRUE)))
    expect_true(all(grepl("includes county area in tiles without an AG100 granule",
                          dict$description[dict$stat %in% "missing_perc"], fixed = TRUE)))
    ndvi <- dict[dict$name == "aster_ndvi_mean_mean", ]
    expect_identical(ndvi$source, "ASTER GED AG100 v003")
    expect_identical(ndvi$band, "ndvi_mean")
    expect_identical(ndvi$stat, "mean")
    expect_identical(ndvi$unit, if (convert) "NDVI" else "NDVI\u00d7100")
  }
})

test_that("provenance describes the used tiles without exporting local paths", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  tile <- attr(df, "aster_ged_tiles")$file
  paths <- write_county_summary(df, file.path(dir, "summary"), formats = character(),
                                dictionary = FALSE)
  provenance <- readLines(paths[["provenance"]])
  expect_true("- CMR collection ID: C2763266348-LPCLOUD" %in% provenance)
  expect_true("- Retrieved (UTC): 2026-01-01T12:34:56Z" %in% provenance)
  expect_true("- Boundary source: US Census cartographic boundary" %in% provenance)
  expect_true("- Boundary year: 2020" %in% provenance)
  expect_true("- Boundary resolution: 500k" %in% provenance)
  expect_true("- Rows: 1" %in% provenance)
  expect_true("- Unique AOI IDs: 1" %in% provenance)
  expect_true("- Missing FIPS codes: 01999" %in% provenance)
  expect_true("- `bands`: elev, ndvi_mean" %in% provenance)
  expect_true(paste("| 1.0 | converted | ag100.tif |", file.size(tile), "|",
                    unname(tools::md5sum(tile)), "|") %in% provenance)
  expect_false(any(grepl(tempdir(), provenance, fixed = TRUE)))

  unlink(tile)
  write_county_summary(df, file.path(dir, "summary"), formats = character(),
                       dictionary = FALSE, overwrite = TRUE)
  expect_true(any(grepl("^\\| 1\\.0 \\| converted \\| ag100\\.tif \\| [0-9]+ \\|  \\|$",
                        readLines(paths[["provenance"]]))))
})

test_that("provenance lists tiles without data from the fetch index", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  fetch <- attr(df, "traktr_layers")$aster_ged
  unavailable <- data.frame(tile_id = c("1.1", "1.2", "1.3"),
                            status = c("missing_remote", "missing_local", "error"),
                            file = NA_character_, size_bytes = NA_real_)
  df <- .set_traktr_layer(df, "aster_ged", rbind(fetch$index, unavailable), fetch$params)
  attr(df, "aster_ged_map") <- rbind(attr(df, "aster_ged_map"), data.frame(
    aoi_id = c("01001", "01003", "01001"), tile_id = "1.1"
  ))
  stem <- file.path(dir, "summary")
  paths <- write_county_summary(df, stem, formats = character(), dictionary = FALSE)
  provenance <- readLines(paths[["provenance"]])
  expect_true("## Tiles without data" %in% provenance)
  expect_true("| 1.1 | missing_remote | 01001, 01003 |" %in% provenance)
  expect_true("| 1.2 | missing_local |  |" %in% provenance)
  expect_true("| 1.3 | error |  |" %in% provenance)

  attr(df, "aster_ged_map") <- NULL
  write_county_summary(df, stem, formats = character(), dictionary = FALSE, overwrite = TRUE)
  expect_true("| 1.1 | missing_remote |  |" %in% readLines(paths[["provenance"]]))
})

test_that("existing outputs abort before any files are published", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  stem <- file.path(dir, "result")
  dict_path <- paste0(stem, "_dictionary.csv")
  writeLines("keep me", dict_path)
  expect_error(write_county_summary(df, stem), "already exist")
  expect_identical(readLines(dict_path), "keep me")
  expect_false(file.exists(paste0(stem, ".csv")))
  expect_false(file.exists(paste0(stem, "_PROVENANCE.md")))

  paths <- write_county_summary(df, stem, overwrite = TRUE)
  expect_true(all(file.exists(paths)))
  expect_identical(utils::read.csv(dict_path)$name, setdiff(names(df), "shape"))
  expect_setequal(list.files(dir), c("ag100.tif", basename(paths)))
})

test_that("failure publishing a later file rolls back new files and restores originals", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  publish <- .publish_county_summary_file
  published <- character()
  local_mocked_bindings(.publish_county_summary_file = function(from, to, overwrite) {
    if (length(published)) return(FALSE)
    result <- publish(from, to, overwrite)
    if (result) published <<- to
    result
  })

  stem <- file.path(dir, "result")
  expect_error(write_county_summary(df, stem), "Could not publish.*result_dictionary")
  expect_identical(published, paste0(stem, ".csv"))
  expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), "ag100.tif")

  originals <- paste0(stem, c(".csv", "_PROVENANCE.md"))
  writeLines("old csv", originals[1L])
  writeLines("old provenance", originals[2L])
  before <- tools::md5sum(originals)
  published <- character()
  expect_error(write_county_summary(df, stem, overwrite = TRUE), "Could not publish.*result_dictionary")
  expect_identical(published, paste0(stem, ".csv"))
  expect_identical(tools::md5sum(originals), before)
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE),
                  c("ag100.tif", basename(originals)))
})

test_that("an output created during publication is preserved without overwrite", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  stem <- file.path(dir, "result")
  concurrent <- paste0(stem, "_dictionary.csv")
  publish <- .publish_county_summary_file
  local_mocked_bindings(.publish_county_summary_file = function(from, to, overwrite) {
    if (identical(to, concurrent)) writeLines("concurrent export", to)
    publish(from, to, overwrite)
  })

  expect_error(write_county_summary(df, stem), "Could not publish.*result_dictionary")
  expect_identical(readLines(concurrent), "concurrent export")
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE),
                  c("ag100.tif", basename(concurrent)))
})

test_that("a dangling symlink counts as an existing output", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  stem <- file.path(dir, "result")
  link <- paste0(stem, "_dictionary.csv")
  missing <- file.path(dir, "missing.csv")
  expect_true(file.symlink(missing, link))

  expect_error(write_county_summary(df, stem), "already exist.*result_dictionary")
  expect_identical(Sys.readlink(link), missing)
  expect_false(file.exists(missing))
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE), c("ag100.tif", basename(link)))
})

test_that("Parquet round-trips FIPS codes as character", {
  skip_if_not_installed("nanoparquet")
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  paths <- write_county_summary(df, file.path(dir, "result.parquet"),
                                formats = c("parquet", "csv"))
  expect_named(paths, c("parquet", "csv", "dictionary", "provenance"))
  result <- nanoparquet::read_parquet(paths[["parquet"]])
  expect_identical(result$GEOID, "01001")
  expect_named(result, setdiff(names(df), "shape"))
})

test_that("requesting unavailable Parquet support aborts without writing files", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  local_mocked_bindings(.has_nanoparquet = function() FALSE)
  stem <- file.path(dir, "new", "result")
  expect_error(write_county_summary(df, stem, formats = c("csv", "parquet")),
               "install.packages.*nanoparquet")
  expect_false(dir.exists(dirname(stem)))
})

test_that("export requires a summary layer and rejects list columns by name", {
  stem <- file.path(withr::local_tempdir(), "result")
  expect_error(write_county_summary(data.frame(GEOID = "01001"), stem), "run attach_aster_ged")
  df <- county_export_fixture(dirname(stem))
  df$nest <- list(c(1, 2))
  df$labels <- list(c("x", "y"))
  expect_error(write_county_summary(df, stem), "nest.*labels")
  expect_error(write_county_summary(df, stem, formats = "txt"), "formats.*csv.*parquet")
  expect_false(file.exists(paste0(stem, ".csv")))
})

test_that("rollback keeps an earlier output that a concurrent export replaced", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  stem <- file.path(dir, "result")
  csv <- paste0(stem, ".csv")
  dictionary <- paste0(stem, "_dictionary.csv")
  publish <- .publish_county_summary_file
  local_mocked_bindings(.publish_county_summary_file = function(from, to, overwrite) {
    if (identical(to, dictionary)) {
      writeLines("other export", csv)
      return(FALSE)
    }
    publish(from, to, overwrite)
  })

  expect_error(write_county_summary(df, stem), "Could not publish.*result_dictionary")
  expect_identical(readLines(csv), "other export")
})

test_that("provenance reports tiles that vanished between fetch and attach", {
  dir <- withr::local_tempdir()
  df <- county_export_fixture(dir)
  unlink(file.path(dir, "ag100.tif"))
  expect_warning(
    df <- attach_aster_ged(df, bands = c("elev", "ndvi_mean"), allow_missing = TRUE,
                           progress = FALSE, cache = FALSE),
    "Unavailable ASTER GED tiles"
  )
  expect_equal(df$ASTER_county_elevation_missing_perc, 100)
  paths <- write_county_summary(df, file.path(dir, "summary"), formats = character(),
                                dictionary = FALSE)
  provenance <- readLines(paths[["provenance"]])
  expect_true("| 1.0 | missing_local | 01001 |" %in% provenance)
})
