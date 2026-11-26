ag100_rectangle <- function(xmin, xmax, ymin, ymax) {
  ring <- matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
                 ncol = 2L, byrow = TRUE)
  sf::st_sf(id = "a", geom = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326))
}

cmr_response <- function(page, cursor = NULL) {
  path <- testthat::test_path("fixtures", paste0("cmr_ag100_page", page, ".json"))
  header <- if (length(cursor)) paste0("CMR-Search-After: ", cursor, "\r\n")
  list(status_code = 200L,
       content = readBin(path, "raw", n = file.size(path)),
       headers = charToRaw(paste0("HTTP/1.1 200 OK\r\n", header, "\r\n")))
}

test_that("Earthdata checks credentials locally without disclosing values", {
  netrc <- withr::local_tempfile()
  withr::local_envvar(c(EARTHDATA_TOKEN = "secret-xyz", EARTHDATA_NETRC = netrc))
  writeLines("machine urs.earthdata.nasa.gov login u password p", netrc)

  output <- capture.output(
    messages <- capture.output(token <- withVisible(earthdata_auth_check()), type = "message")
  )
  expect_identical(token, list(value = structure(TRUE, method = "token"), visible = FALSE))
  expect_match(paste(messages, collapse = "\n"), "token")
  expect_no_match(paste(c(output, messages), collapse = "\n"), "secret-xyz", fixed = TRUE)
  expect_no_message(earthdata_auth_check(verbose = FALSE))

  withr::local_envvar(c(EARTHDATA_TOKEN = ""))
  expect_message(netrc_result <- earthdata_auth_check(), "netrc")
  expect_identical(netrc_result, structure(TRUE, method = "netrc"))
  writeLines(c("# machine urs.earthdata.nasa.gov login u password p",
               "machine urs.earthdata.nasa.gov login u",
               "machine other.example.test password p"), netrc)
  expect_message(none <- earthdata_auth_check(), "EARTHDATA_TOKEN")
  expect_identical(none, structure(FALSE, method = "none"))
  unlink(netrc)
  expect_identical(earthdata_auth_check(FALSE), none)
  expect_error(.earthdata_handle(), "EARTHDATA_TOKEN")
})

test_that("Earthdata curl handles use a bearer header or netrc and a cookie jar", {
  netrc <- withr::local_tempfile()
  writeLines(c("machine urs.earthdata.nasa.gov", "  login u", "  password p"), netrc)
  cache <- withr::local_tempdir()
  withr::local_envvar(c(EARTHDATA_TOKEN = "secret-xyz", EARTHDATA_NETRC = netrc,
                       R_USER_CACHE_DIR = cache))

  headers <- options <- NULL
  new_handle <- curl::new_handle
  local_mocked_bindings(
    handle_setheaders = function(handle, ...) {
      headers <<- list(...)
      invisible(handle)
    },
    new_handle = function(...) {
      options <<- list(...)
      new_handle(...)
    },
    .package = "curl"
  )

  expect_s3_class(.earthdata_handle(), "curl_handle")
  expect_identical(headers, list(Authorization = "Bearer secret-xyz"))
  expect_identical(options, list(followlocation = TRUE))

  withr::local_envvar(c(EARTHDATA_TOKEN = ""))
  expect_s3_class(.earthdata_handle(), "curl_handle")
  jar <- file.path(tools::R_user_dir("traktr", "cache"), "earthdata_cookies.txt")
  expect_true(startsWith(jar, cache))
  expect_true(dir.exists(dirname(jar)))
  expect_identical(options, list(followlocation = TRUE, netrc = 1L, netrc_file = netrc,
                                 unrestricted_auth = 1L, cookiefile = jar, cookiejar = jar))
})

test_that("Earthdata mentions project Renviron precedence when the file exists", {
  withr::local_dir(withr::local_tempdir())
  withr::local_envvar(c(EARTHDATA_TOKEN = "", EARTHDATA_NETRC = "missing-netrc"))
  messages <- capture.output(earthdata_auth_check(), type = "message")
  expect_no_match(paste(messages, collapse = "\n"), "takes precedence")
  expect_error(.earthdata_handle(), "EARTHDATA_TOKEN")

  expect_true(file.create(".Renviron"))
  expect_message(earthdata_auth_check(), "takes precedence")
  expect_error(.earthdata_handle(), "takes precedence")
})

test_that("Earthdata matches netrc hostnames case-insensitively", {
  netrc <- withr::local_tempfile()
  withr::local_envvar(c(EARTHDATA_TOKEN = "", EARTHDATA_NETRC = netrc))
  writeLines("machine URS.EARTHDATA.NASA.GOV login u password p", netrc)

  expect_identical(earthdata_auth_check(FALSE), structure(TRUE, method = "netrc"))
})

test_that("AG100 tiles use north and west edges and maintain unique AOI pairs", {
  withr::local_options(sf_use_s2 = TRUE)
  aoi <- ag100_rectangle(-105.3, -104.2, 39.1, 39.9)
  south <- ag100_rectangle(10.2, 10.8, -0.8, -0.2)
  south$id <- "b"
  result <- .ag100_tiles_for_aoi(rbind(aoi, south, aoi), "id")
  expect_identical(result$tiles$tile_id, c("40.-106", "40.-105", "0.10"))
  expect_identical(result$tiles$lat_n, c(40L, 40L, 0L))
  expect_identical(result$tiles$lon_w, c(-106L, -105L, 10L))
  expect_identical(result$tiles$ymin, c(39, 39, -1))
  expect_identical(result$map$aoi_id, c("a", "a", "b"))
  expect_setequal(result$map$tile_id, result$tiles$tile_id)
  expect_true(sf::sf_use_s2())
  expect_identical(.ag100_tiles_for_aoi(aoi[0, ], "id")$tiles, result$tiles[0, ])
})

test_that("AG100 tile selection excludes boundary-only intersections", {
  aoi <- ag100_rectangle(-105, -104, 39, 40)
  expect_identical(.ag100_tiles_for_aoi(aoi, "id")$tiles$tile_id, "40.-105")

  outer <- sf::st_geometry(ag100_rectangle(0, 3, 0, 3))[[1L]][[1L]]
  hole <- sf::st_geometry(ag100_rectangle(1, 2, 1, 2))[[1L]][[1L]]
  aoi <- sf::st_sf(id = "hole", geom = sf::st_sfc(sf::st_polygon(list(outer, hole)), crs = 4326))
  tiles <- .ag100_tiles_for_aoi(aoi, "id")$tiles$tile_id
  expect_length(tiles, 8L)
  expect_false("2.1" %in% tiles)
})

test_that("AG100 county tiles respect parts on either side of the dateline", {
  counties <- aoi_us_counties(src = testthat::test_path("fixtures", "counties5.gpkg"))
  result <- .ag100_tiles_for_aoi(attr(counties, "traktr_aoi"), "GEOID")
  west <- result$map$tile_id[result$map$aoi_id == "02016"]
  longitude <- result$tiles$lon_w[match(west, result$tiles$tile_id)]
  expect_true(all(c(179L, -180L) %in% longitude))
  expect_true(all(longitude < -170 | longitude > 170))
  expect_setequal(result$map$aoi_id, counties$GEOID)
})

test_that("AG100 selects only adjoining longitude tiles for H3 and buffer AOIs at the dateline", {
  point <- tibble::tibble(lon = 179.9995, lat = 52)
  for (out in list(mutate_aoi_h3(point, res = 6), mutate_aoi_buffer(point, dist = 1000))) {
    aoi <- attr(out, "traktr_aoi")
    id_col <- attr(out, "traktr_aoi_params")$id_col
    result <- .ag100_tiles_for_aoi(aoi, id_col)
    expect_setequal(result$tiles$lon_w, c(179L, -180L))
  }
})

test_that("AG100 normalizes longitude bounds in both latitude bands of a crossing polygon", {
  aoi <- ag100_rectangle(179.5, -179.5, 51.5, 52.5)
  result <- .ag100_tiles_for_aoi(aoi, "id")
  expect_identical(result$tiles$tile_id, c("53.-180", "53.179", "52.-180", "52.179"))
  expect_identical(result$tiles$xmin, c(-180, 179, -180, 179))
  expect_identical(result$tiles$xmax, c(-179, 180, -179, 180))
  expect_setequal(result$map$tile_id, result$tiles$tile_id)
})

test_that("AG100 maps AOI IDs even when the source column is named tile_id", {
  aoi <- ag100_rectangle(1.2, 1.8, 1.2, 1.8)
  aoi$tile_id <- "area-1"
  out <- aoi_from_sf(aoi, id_col = "tile_id")
  result <- .ag100_tiles_for_aoi(attr(out, "traktr_aoi"), "tile_id")

  expect_identical(result$id_col, "tile_id")
  expect_named(result$map, c("aoi_id", "tile_id"))
  expect_identical(result$map$aoi_id, "area-1")
  expect_identical(result$map$tile_id, "2.1")
})

test_that("AG100 granule names and URLs are vectorized", {
  name <- "AG100.v003.40.-105.0001"
  expect_identical(.ag100_tile_name(40L, -105L), name)
  expect_identical(.ag100_tile_name(c(40L, 0L), c(-105L, 10L)),
                   c(name, "AG100.v003.0.10.0001"))
  expect_identical(
    .ag100_url(40L, -105L),
    paste0("https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/AG100.003/",
           name, "/", name, ".h5")
  )
})

test_that("CMR follows the search cursor and uses titles and HTTPS data links", {
  skip_if_not_installed("jsonlite")
  search_headers <- list()
  requests <- list()
  cursor <- "[1685681106164,1234]"
  withr::local_envvar(c(EARTHDATA_TOKEN = "secret-xyz"))
  local_mocked_bindings(
    handle_setheaders = function(handle, ...) {
      search_headers <<- list(...)
      invisible(handle)
    },
    curl_fetch_memory = function(url, handle) {
      page <- length(requests) + 1L
      requests[[page]] <<- list(url = url, headers = search_headers, handle = handle)
      cmr_response(page, if (page == 1L) cursor)
    },
    .package = "curl"
  )

  granules <- .cmr_granules("C2763266348-LPCLOUD", c(-106, 39, -104, 40), page_size = 2)
  expect_length(requests, 2L)
  expect_identical(requests[[1L]]$headers, list())
  expect_identical(requests[[2L]]$headers, list(`CMR-Search-After` = cursor))
  expect_identical(requests[[1L]]$handle, requests[[2L]]$handle)
  expect_identical(
    requests[[1L]]$url,
    paste0("https://cmr.earthdata.nasa.gov/search/granules.json?",
           "collection_concept_id=C2763266348-LPCLOUD&bounding_box=-106,39,-104,40&page_size=2")
  )
  expect_identical(granules$title[2L], "AG100.v003.05.-075.0001")
  expect_identical(granules$lat_n, c(40L, 5L, 52L))
  expect_identical(granules$lon_w, c(-105L, -75L, 179L))
  expect_identical(granules$url[1L], "https://example.test/AG100.v003.40.-105.0001.h5")
})

test_that("AG100 resolves each longitude group independently and marks absent tiles", {
  skip_if_not_installed("jsonlite")
  east <- ag100_rectangle(179.2, 179.8, 51.2, 51.8)
  west <- ag100_rectangle(-179.8, -179.2, 51.2, 51.8)
  tiles <- .ag100_tiles_for_aoi(rbind(east, west), "id")$tiles
  urls <- character()
  local_mocked_bindings(curl_fetch_memory = function(url, handle) {
    urls <<- c(urls, url)
    if (grepl("bounding_box=-180,", url, fixed = TRUE)) {
      list(status_code = 200L, content = charToRaw('{"feed":{"entry":[]}}'))
    } else {
      cmr_response(2)
    }
  }, .package = "curl")

  resolved <- .ag100_resolve(tiles)
  expect_length(urls, 2L)
  expect_true(any(grepl("bounding_box=179,51,180,52&", urls, fixed = TRUE)))
  expect_true(any(grepl("bounding_box=-180,51,-179,52&", urls, fixed = TRUE)))
  expect_identical(resolved$tile_id, c("52.-180", "52.179"))
  expect_identical(resolved$status, c("missing_remote", NA_character_))
  expect_identical(resolved$url, c(NA_character_, "https://example.test/AG100.v003.52.179.0001.h5"))
  expect_identical(resolved$url_source, rep("cmr", 2L))
})

test_that("AG100 prefers HTTPS links and constructs missing links from validated CMR titles", {
  skip_if_not_installed("jsonlite")
  aoi <- rbind(ag100_rectangle(-105, -104, 39, 40),
               ag100_rectangle(-75, -74, 4, 5),
               ag100_rectangle(-70, -69, 2, 3))
  tiles <- .ag100_tiles_for_aoi(aoi, "id")$tiles
  local_mocked_bindings(curl_fetch_memory = function(url, handle) cmr_response(3),
                        .package = "curl")

  resolved <- .ag100_resolve(tiles)
  padded <- "AG100.v003.05.-075.0001"
  expect_identical(resolved$granule, c("AG100.v003.40.-105.0001", padded, "invalid title"))
  expect_identical(resolved$status, c(NA_character_, NA_character_, "error"))
  expect_identical(resolved$url_source, c("cmr", "constructed", "cmr"))
  expect_identical(resolved$url, c(
    "https://example.test/direct.h5",
    paste0("https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/AG100.003/",
           padded, "/", padded, ".h5"),
    NA_character_
  ))
})

test_that("AG100 assigns precise CMR boxes before title coordinates", {
  skip_if_not_installed("jsonlite")
  aoi <- rbind(ag100_rectangle(-78, -77, 7, 8),
               ag100_rectangle(-77, -76, 6, 7),
               ag100_rectangle(-75, -74, 4, 5),
               ag100_rectangle(-70, -69, 2, 3))
  tiles <- .ag100_tiles_for_aoi(aoi, "id")$tiles
  local_mocked_bindings(curl_fetch_memory = function(url, handle) cmr_response(4),
                        .package = "curl")

  warnings <- character()
  resolved <- withCallingHandlers(.ag100_resolve(tiles), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(warnings, "Discarded 2 CMR granules without usable tile coordinates.")
  expect_identical(resolved$tile_id, c("8.-78", "7.-77", "5.-75", "3.-70"))
  expect_identical(resolved$granule, c("AG100.v003.07.-077.0001", NA_character_,
                                      NA_character_, "AG100.v003.05.-075.unknown"))
  expect_identical(resolved$status, c(NA_character_, "missing_remote", "missing_remote", "error"))
  expect_identical(resolved$url[4L], NA_character_)
})

test_that("AG100 picks the first HTTPS link across both title-based and box-based matches", {
  skip_if_not_installed("jsonlite")
  tiles <- .ag100_tiles_for_aoi(ag100_rectangle(-71, -70, 3, 4), "id")$tiles
  local_mocked_bindings(curl_fetch_memory = function(url, handle) cmr_response(3),
                        .package = "curl")

  resolved <- .ag100_resolve(tiles)
  expect_identical(resolved$granule, "noncanonical title with direct link")
  expect_identical(resolved$url, "https://example.test/box-match.h5")
  expect_identical(resolved$url_source, "cmr")
  expect_identical(resolved$status, NA_character_)
})

test_that("AG100 constructs every URL after an HTTP failure or when CMR is disabled", {
  skip_if_not_installed("jsonlite")
  east <- ag100_rectangle(179.2, 179.8, 51.2, 51.8)
  west <- ag100_rectangle(-179.8, -179.2, 51.2, 51.8)
  tiles <- .ag100_tiles_for_aoi(rbind(east, west), "id")$tiles
  calls <- 0L
  local_mocked_bindings(curl_fetch_memory = function(url, handle) {
    calls <<- calls + 1L
    if (grepl("bounding_box=-180,", url, fixed = TRUE)) {
      list(status_code = 500L)
    } else {
      cmr_response(2)
    }
  }, .package = "curl")
  expect_error(.cmr_granules("id", c(-180, 51, -179, 52)), "HTTP 500")
  expect_warning(fallback <- .ag100_resolve(tiles), "constructed AG100 URLs")
  expect_identical(calls, 3L)
  expect_identical(fallback$url_source, rep("constructed", 2L))
  expect_identical(fallback$status, rep(NA_character_, 2L))
  expect_identical(fallback$url, .ag100_url(tiles$lat_n, tiles$lon_w))
  expect_warning(disabled <- .ag100_resolve(tiles, use_cmr = FALSE), "constructed AG100 URLs")
  expect_identical(disabled, fallback)
  expect_identical(calls, 3L)
})

test_that("AG100 constructs URLs when jsonlite is unavailable", {
  tiles <- .ag100_tiles_for_aoi(ag100_rectangle(-105, -104, 39, 40), "id")$tiles
  local_mocked_bindings(.has_jsonlite = function() FALSE)
  local_mocked_bindings(curl_fetch_memory = function(...) {
    fail("CMR should be skipped without jsonlite")
  }, .package = "curl")

  expect_warning(resolved <- .ag100_resolve(tiles), "constructed AG100 URLs")
  expect_identical(resolved$url_source, "constructed")
  expect_identical(resolved$granule, "AG100.v003.40.-105.0001")
  expect_identical(resolved$status, NA_character_)
})

test_that("CMR boxes at the 1e-6 tolerance limit are accepted", {
  expect_identical(.cmr_ag100_box(list("2.000001 -70 3.000001 -69")), c(3L, -70L))
  expect_identical(.cmr_ag100_box(list("2 -70 3 -69")), c(3L, -70L))
  expect_identical(.cmr_ag100_box(list("2.00001 -70 3.00001 -69")), c(NA_integer_, NA_integer_))
})
