#' Check local Earthdata credentials
#'
#' Checks for an Earthdata Login token in `EARTHDATA_TOKEN` or a netrc
#' entry for `urs.earthdata.nasa.gov`. Set `EARTHDATA_NETRC` to use a
#' specific netrc file. The default is `~/.netrc`, with `~/_netrc`
#' also checked on Windows. The check works offline.
#'
#' @param verbose Whether to report the detected method or print setup instructions.
#' @return Invisibly, a logical scalar with a `method` attribute of
#'   `"token"`, `"netrc"`, or `"none"`.
#' @examples
#' earthdata_auth_check(verbose = FALSE)
#' @export
earthdata_auth_check <- function(verbose = TRUE) {
  method <- .earthdata_credentials()$method
  if (verbose) {
    if (method == "none") {
      cli::cli_inform(.earthdata_auth_help())
    } else {
      cli::cli_inform("Earthdata credentials found via {.val {method}}.")
    }
  }
  invisible(structure(method != "none", method = method))
}

#' @noRd
.earthdata_auth_help <- function() {
  c(
    "Earthdata credentials are missing.",
    "i" = "Create an Earthdata Login account and generate a token.",
    "i" = "Set {.envvar EARTHDATA_TOKEN} in {.path ~/.Renviron} and restart R."
  )
}

#' @noRd
.earthdata_credentials <- function() {
  if (nzchar(Sys.getenv("EARTHDATA_TOKEN"))) {
    return(list(method = "token"))
  }

  paths <- Sys.getenv("EARTHDATA_NETRC")
  if (!nzchar(paths)) {
    paths <- "~/.netrc"
    if (.Platform$OS.type == "windows") paths <- c(paths, "~/_netrc")
  }

  for (path in path.expand(paths)) {
    if (.earthdata_netrc_has_login(path)) {
      return(list(method = "netrc", path = path))
    }
  }
  list(method = "none")
}

#' @noRd
.earthdata_netrc_has_login <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)

  tokens <- tryCatch(
    suppressWarnings(scan(path, what = character(), quiet = TRUE,
                          comment.char = "#", na.strings = NULL)),
    error = function(e) character()
  )
  earthdata <- login <- password <- FALSE
  i <- 1L
  while (i <= length(tokens)) {
    key <- tokens[[i]]
    if (key %in% c("machine", "default")) {
      earthdata <- key == "machine" && i < length(tokens) &&
        tolower(tokens[[i + 1L]]) == "urs.earthdata.nasa.gov"
      login <- password <- FALSE
    }

    if (key %in% c("machine", "login", "password", "account") && i < length(tokens)) {
      value <- tokens[[i + 1L]]
      if (key == "login") login <- nzchar(value)
      if (key == "password") password <- nzchar(value)
      i <- i + 1L
    }
    if (earthdata && login && password) return(TRUE)
    i <- i + 1L
  }
  FALSE
}

#' @noRd
.earthdata_handle <- function(error_call = rlang::caller_env()) {
  credentials <- .earthdata_credentials()
  if (credentials$method == "token") {
    h <- curl::new_handle(followlocation = TRUE)
    curl::handle_setheaders(h, Authorization = paste("Bearer", Sys.getenv("EARTHDATA_TOKEN")))
    return(h)
  }

  if (credentials$method == "netrc") {
    cache_dir <- tools::R_user_dir("traktr", "cache")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    jar <- file.path(cache_dir, "earthdata_cookies.txt")
    return(curl::new_handle(followlocation = TRUE, netrc = 1L,
                            netrc_file = credentials$path, unrestricted_auth = 1L,
                            cookiefile = jar, cookiejar = jar))
  }

  cli::cli_abort(.earthdata_auth_help(), call = error_call)
}

#' @noRd
.ag100_tiles_for_aoi <- function(aoi, id_col) {
  previous_s2 <- getOption("sf_use_s2")
  on.exit(options(sf_use_s2 = previous_s2), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  parts <- sf::st_cast(aoi[, id_col, drop = FALSE], "POLYGON", warn = FALSE)
  ids <- parts[[id_col]]
  geometry <- sf::st_geometry(parts)
  tile_rows <- vector("list", length(geometry))
  map_rows <- vector("list", length(geometry))
  for (i in which(!sf::st_is_empty(geometry))) {
    polygon <- geometry[i]
    bbox <- sf::st_bbox(polygon)
    crosses_dateline <- bbox[["xmax"]] - bbox[["xmin"]] > 180
    if (crosses_dateline) {
      # Keep rings crossing the dateline contiguous during planar intersection.
      polygon <- sf::st_shift_longitude(polygon)
      bbox <- sf::st_bbox(polygon)
    }
    lon_w <- seq.int(floor(bbox[["xmin"]]), ceiling(bbox[["xmax"]]) - 1L)
    lat_n <- seq.int(floor(bbox[["ymin"]]) + 1L, ceiling(bbox[["ymax"]]))
    grid <- expand.grid(lon_w = lon_w, lat_n = lat_n, KEEP.OUT.ATTRS = FALSE)
    cells <- sf::st_make_grid(
      polygon, cellsize = 1,
      offset = c(lon_w[1L], lat_n[1L] - 1L),
      n = c(length(lon_w), length(lat_n))
    )
    hits <- suppressMessages(sf::st_intersects(polygon, cells))[[1L]]
    has_area <- lengths(suppressMessages(
      sf::st_relate(cells[hits], polygon, pattern = "2********")
    )) > 0L
    grid <- grid[hits[has_area], , drop = FALSE]
    if (crosses_dateline) {
      shifted <- grid$lon_w >= 180L
      grid$lon_w[shifted] <- grid$lon_w[shifted] - 360L
    }

    tile_id <- paste(grid$lat_n, grid$lon_w, sep = ".")
    tile_rows[[i]] <- tibble::tibble(
      tile_id = tile_id,
      lat_n = as.integer(grid$lat_n),
      lon_w = as.integer(grid$lon_w),
      xmin = as.double(grid$lon_w), ymin = grid$lat_n - 1,
      xmax = grid$lon_w + 1, ymax = as.double(grid$lat_n)
    )
    map_rows[[i]] <- tibble::tibble(aoi_id = rep(ids[i], length(tile_id)), tile_id = tile_id)
  }

  empty_tiles <- tibble::tibble(tile_id = character(), lat_n = integer(), lon_w = integer(),
                                xmin = numeric(), ymin = numeric(),
                                xmax = numeric(), ymax = numeric())
  empty_map <- tibble::tibble(aoi_id = ids[0], tile_id = character())
  tiles <- unique(do.call(rbind, c(list(empty_tiles), tile_rows)))
  map <- unique(do.call(rbind, c(list(empty_map), map_rows)))
  list(tiles = tiles[order(-tiles$lat_n, tiles$lon_w), ], map = map, id_col = id_col)
}

#' @noRd
.ag100_tile_name <- function(lat_n, lon_w) {
  paste0("AG100.v003.", lat_n, ".", lon_w, ".0001", recycle0 = TRUE)
}

#' @noRd
.ag100_url <- function(lat_n, lon_w) {
  .ag100_url_from_name(.ag100_tile_name(lat_n, lon_w))
}

#' @noRd
.ag100_url_from_name <- function(name) {
  paste0("https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/AG100.003/",
         name, "/", name, ".h5", recycle0 = TRUE)
}

#' @noRd
.has_jsonlite <- function() {
  requireNamespace("jsonlite", quietly = TRUE)
}

#' @noRd
.cmr_ag100_box <- function(boxes) {
  if (length(boxes) != 1L || !is.character(boxes[[1L]])) {
    return(c(NA_integer_, NA_integer_))
  }

  bounds <- suppressWarnings(as.numeric(strsplit(trimws(boxes[[1L]]), "[[:space:]]+")[[1L]]))
  edges <- round(bounds)
  if (length(bounds) != 4L || any(!is.finite(bounds)) ||
      # 1e-6 tolerance, plus slack for floating-point rounding at the limit
      any(abs(bounds - edges) > 1e-6 + 1e-9) || any(edges[3:4] - edges[1:2] != 1) ||
      edges[1L] < -90 || edges[2L] < -180 || edges[3L] > 90 || edges[4L] > 180) {
    return(c(NA_integer_, NA_integer_))
  }
  as.integer(edges[c(3L, 2L)])
}

#' @noRd
.cmr_granules <- function(concept_id, bbox, page_size = 2000, max_pages = 50) {
  if (!.has_jsonlite()) cli::cli_abort("CMR requires the {.pkg jsonlite} package.")

  bounds <- paste(format(bbox, trim = TRUE, scientific = FALSE, decimal.mark = "."), collapse = ",")
  url <- paste0("https://cmr.earthdata.nasa.gov/search/granules.json?",
                "collection_concept_id=", concept_id,
                "&bounding_box=", bounds, "&page_size=", page_size)
  h <- curl::new_handle()
  entries <- vector("list", max_pages)
  search_after <- NULL
  for (page in seq_len(max_pages)) {
    if (length(search_after)) {
      curl::handle_setheaders(h, `CMR-Search-After` = search_after)
    }
    response <- curl::curl_fetch_memory(url, handle = h)
    if (response$status_code != 200L) {
      cli::cli_abort("CMR returned HTTP {response$status_code}.")
    }
    payload <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
    page_entries <- payload$feed$entry
    if (!is.list(page_entries)) cli::cli_abort("CMR returned an invalid granule feed.")
    if (!length(page_entries)) break

    entries[[page]] <- page_entries
    search_after <- curl::parse_headers_list(response$headers)[["cmr-search-after"]]
    if (!length(search_after) || !nzchar(search_after)) break
    if (page == max_pages) cli::cli_abort("CMR exceeded the limit of {max_pages} pages.")
  }

  granules <- unlist(entries, recursive = FALSE)
  title <- vapply(granules, function(entry) {
    if (is.null(entry$title)) NA_character_ else entry$title
  }, character(1))
  href <- vapply(granules, function(entry) {
    links <- vapply(entry$links, function(link) {
      if (is.null(link$href)) NA_character_ else link$href
    }, character(1))
    links[grepl("^https://.*\\.h5$", links)][1L]
  }, character(1))

  title_coords <- utils::strcapture(
    "^AG100\\.v003\\.(-?\\d+)\\.(-?\\d+)\\.0001$", title,
    list(lat_n = numeric(), lon_w = numeric())
  )
  valid_title <- !is.na(title_coords$lat_n) & !is.na(title_coords$lon_w)
  boxes <- vapply(granules, function(entry) .cmr_ag100_box(entry$boxes), integer(2L))
  lat_n <- boxes[1L, ]
  lon_w <- boxes[2L, ]
  from_title <- is.na(lat_n) & valid_title &
    title_coords$lat_n > -90 & title_coords$lat_n <= 90 &
    title_coords$lon_w >= -180 & title_coords$lon_w < 180
  lat_n[from_title] <- as.integer(title_coords$lat_n[from_title])
  lon_w[from_title] <- as.integer(title_coords$lon_w[from_title])
  tibble::tibble(title = title, url = href, lat_n = lat_n, lon_w = lon_w,
                 valid_title = valid_title)
}

#' @noRd
.ag100_resolve <- function(tiles, concept_id = "C2763266348-LPCLOUD", use_cmr = TRUE) {
  tiles$granule <- tiles$url <- tiles$url_source <- tiles$status <- rep(NA_character_, nrow(tiles))
  if (!nrow(tiles)) return(tiles)

  remote <- NULL
  if (use_cmr && .has_jsonlite()) {
    groups <- split(seq_len(nrow(tiles)), tiles$lon_w < 0L)
    remote <- tryCatch(
      lapply(groups, function(i) {
        bbox <- c(min(tiles$xmin[i]), min(tiles$ymin[i]),
                  max(tiles$xmax[i]), max(tiles$ymax[i]))
        .cmr_granules(concept_id, bbox)
      }),
      error = function(e) NULL
    )
  }

  if (is.null(remote)) {
    cli::cli_warn("CMR lookup unavailable or disabled. Using constructed AG100 URLs.")
    tiles$granule <- .ag100_tile_name(tiles$lat_n, tiles$lon_w)
    tiles$url <- .ag100_url(tiles$lat_n, tiles$lon_w)
    tiles$url_source <- "constructed"
    return(tiles)
  }

  remote <- do.call(rbind, remote)
  assigned <- !is.na(remote$lat_n) & !is.na(remote$lon_w)
  discarded <- sum(!assigned)
  if (discarded) {
    cli::cli_warn("Discarded {discarded} CMR granule{?s} without usable tile coordinates.")
  }
  remote <- remote[assigned, ]
  priority <- rep(3L, nrow(remote))
  priority[remote$valid_title] <- 2L
  priority[!is.na(remote$url)] <- 1L
  candidates <- order(priority)
  remote_keys <- paste(remote$lat_n, remote$lon_w)[candidates]
  matches <- candidates[match(paste(tiles$lat_n, tiles$lon_w), remote_keys)]

  tiles$granule <- remote$title[matches]
  tiles$url <- remote$url[matches]
  tiles$url_source <- "cmr"
  tiles$status[is.na(matches)] <- "missing_remote"
  selected <- priority[matches]
  tiles$status[which(selected == 3L)] <- "error"
  constructed <- which(selected == 2L)
  tiles$url[constructed] <- .ag100_url_from_name(tiles$granule[constructed])
  tiles$url_source[constructed] <- "constructed"
  tiles
}
