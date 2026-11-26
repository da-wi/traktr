.state_fips <- c(
  "01" = "AL", "02" = "AK", "04" = "AZ", "05" = "AR", "06" = "CA",
  "08" = "CO", "09" = "CT", "10" = "DE", "11" = "DC", "12" = "FL",
  "13" = "GA", "15" = "HI", "16" = "ID", "17" = "IL", "18" = "IN",
  "19" = "IA", "20" = "KS", "21" = "KY", "22" = "LA", "23" = "ME",
  "24" = "MD", "25" = "MA", "26" = "MI", "27" = "MN", "28" = "MS",
  "29" = "MO", "30" = "MT", "31" = "NE", "32" = "NV", "33" = "NH",
  "34" = "NJ", "35" = "NM", "36" = "NY", "37" = "NC", "38" = "ND",
  "39" = "OH", "40" = "OK", "41" = "OR", "42" = "PA", "44" = "RI",
  "45" = "SC", "46" = "SD", "47" = "TN", "48" = "TX", "49" = "UT",
  "50" = "VT", "51" = "VA", "53" = "WA", "54" = "WV", "55" = "WI",
  "56" = "WY", "60" = "AS", "66" = "GU", "69" = "MP", "72" = "PR",
  "78" = "VI"
)

.fips_preview <- function(x) {
  preview <- paste(utils::head(x, 10L), collapse = ", ")
  if (length(x) > 10L) {
    preview <- paste(preview, "and", length(x) - 10L, "more")
  }
  preview
}

.fips_normalize <- function(x) {
  if (!is.character(x) && !is.numeric(x)) {
    if (!is.logical(x) || !all(is.na(x))) {
      cli::cli_abort("FIPS codes must be numeric or character.")
    }
  }

  missing <- is.na(x) | trimws(as.character(x)) == ""
  if (any(missing)) {
    cli::cli_warn("Removed {sum(missing)} missing or blank FIPS code{?s}.")
    x <- x[!missing]
  }

  if (is.numeric(x)) {
    valid <- is.finite(x) & x >= 0 & x <= 99999 & x == trunc(x)
    digits <- as.integer(x[valid])
  } else {
    x <- trimws(as.character(x))
    valid <- grepl("^[0-9]{1,5}$", x)
    digits <- as.integer(x[valid])
  }
  fips <- character(length(x))
  fips[valid] <- formatC(digits, width = 5L, flag = "0")
  valid <- valid & substr(fips, 1L, 2L) %in% names(.state_fips)
  if (!all(valid)) {
    preview <- .fips_preview(x[!valid])
    cli::cli_abort("Invalid FIPS codes: {.val {preview}}.")
  }

  duplicates <- duplicated(fips)
  if (any(duplicates)) {
    cli::cli_warn("Removed {sum(duplicates)} duplicate FIPS code{?s}.")
  }
  fips[!duplicates]
}

#' Read county FIPS codes from a CSV file
#'
#' @param path Path to a CSV file with a header.
#' @param col Name of the column containing county FIPS codes.
#'
#' @return A character vector of unique five-digit codes in input order.
#'   Missing, blank and duplicate entries are removed with warnings.
#' @seealso [aoi_us_counties()]
#' @export
read_fips_csv <- function(path, col = "x") {
  lines <- readLines(path, warn = FALSE)
  # Skip blank lines before the header; blank lines after it count as entries
  first <- which(nzchar(trimws(lines)))[1]
  if (is.na(first)) {
    cli::cli_abort("{.file {path}} contains no header.")
  }
  data <- utils::read.csv(text = lines[first:length(lines)],
                          colClasses = "character", check.names = FALSE,
                          blank.lines.skip = FALSE)
  if (!col %in% names(data)) {
    cli::cli_abort("Column {.val {col}} is missing. Available columns: {.val {names(data)}}.")
  }
  .fips_normalize(data[[col]])
}

#' Create polygon AOIs from an sf object
#'
#' @param x An `sf` object with polygon or multipolygon geometries and a CRS.
#' @param id_col Name of the column containing unique, non-missing AOI IDs.
#' @param keep Keep all non-geometry columns if `TRUE`, or only the ID column.
#' @param make_valid Repair invalid geometries before attaching them.
#'
#' @details Outer rings may span at most 180 degrees of longitude in WGS84.
#'   Represent polygons on both sides of the antimeridian as separate parts.
#'
#' @return A `traktr_df` tibble with the ID column first. The
#'   `traktr_aoi` attribute contains the ID and geometry columns in
#'   EPSG:4326. The `traktr_aoi_params` attribute records the method,
#'   deduplication setting, source CRS and ID column.
#' @seealso [aoi_us_counties()]
#' @export
aoi_from_sf <- function(x, id_col, keep = TRUE, make_valid = TRUE) {
  if (!inherits(x, "sf")) {
    cli::cli_abort("{.arg x} must be an sf object.")
  }
  source_crs <- sf::st_crs(x)
  if (is.na(source_crs)) {
    cli::cli_abort("{.arg x} must have a CRS.")
  }
  geometry_col <- attr(x, "sf_column")
  if (!is.character(id_col) || length(id_col) != 1L ||
      !id_col %in% setdiff(names(x), geometry_col)) {
    cli::cli_abort("{.arg id_col} must name a non-geometry column in {.arg x}.")
  }
  ids <- x[[id_col]]
  if (anyNA(ids) || anyDuplicated(ids)) {
    cli::cli_abort("AOI IDs in {.val {id_col}} must be unique and non-missing.")
  }
  polygon_types <- c("POLYGON", "MULTIPOLYGON")
  if (!all(sf::st_geometry_type(x) %in% polygon_types)) {
    cli::cli_abort("{.arg x} must contain only polygons and multipolygons.")
  }

  if (!isTRUE(source_crs$epsg == 4326L)) {
    x <- sf::st_transform(x, 4326)
  }
  if (make_valid) {
    valid <- sf::st_is_valid(x)
    invalid <- which(is.na(valid) | !valid)
    if (length(invalid)) {
      geometry <- sf::st_geometry(x)
      geometry[invalid] <- sf::st_make_valid(geometry[invalid])
      sf::st_geometry(x) <- geometry
      if (!all(sf::st_geometry_type(x) %in% polygon_types)) {
        cli::cli_abort("Geometry repair must retain polygon or multipolygon geometries.")
      }
    }
  }
  .assert_no_wide_rings(x[, c(id_col, geometry_col), drop = FALSE])

  data <- sf::st_drop_geometry(x)
  columns <- if (keep) c(id_col, setdiff(names(data), id_col)) else id_col
  .aoi_attach_attrs(
    df = tibble::as_tibble(data[, columns, drop = FALSE]),
    id_vec = ids,
    aoi = x,
    id_col = id_col,
    params = list(
      method = "polygon", deduplicate = FALSE,
      crs_source = if (is.na(source_crs$epsg)) source_crs$Name else source_crs$epsg
    )
  )
}

.assert_no_wide_rings <- function(aoi, max_span = 180) {
  wide <- vapply(sf::st_geometry(aoi), function(geometry) {
    parts <- if (inherits(geometry, "POLYGON")) list(geometry) else geometry
    any(vapply(parts, function(polygon) {
      if (!length(polygon)) return(FALSE)
      longitude <- polygon[[1L]][, 1L]
      length(longitude) > 0L && diff(range(longitude)) > max_span
    }, logical(1)))
  }, logical(1))
  if (any(wide)) {
    ids <- sf::st_drop_geometry(aoi)[[1L]]
    cli::cli_abort("Outer rings span more than {max_span} degrees of longitude for AOI IDs: {.val {ids[wide]}}.")
  }
}

#' Create AOIs from U.S. Census county boundaries
#'
#' @param fips Numeric or character county FIPS codes. Use `NULL` for
#'   all counties allowed by `territories`.
#' @param year Census boundary vintage, an integer from 2020 to 2024.
#'   Connecticut uses planning regions from 2022 onward.
#' @param resolution Boundary scale, one of `"500k"`, `"5m"` or `"20m"`.
#' @param territories Include the 50 states and DC plus Puerto Rico
#'   (`"PR"`), the 50 states and DC alone (`"none"`), or also include
#'   the Island Areas (`"all"`).
#' @param src Optional local Census-style shapefile, zipped shapefile or
#'   GeoPackage. A local source bypasses downloads.
#' @param cache Reuse a previously downloaded archive if available.
#' @param quiet Suppress download messages.
#'
#' @return A `traktr_df` with county identifiers, names, land and water
#'   areas, and region labels. Rows follow the requested FIPS order or,
#'   when `fips` is `NULL`, GEOID order. The `traktr_aoi` attribute stores
#'   county polygons in WGS84. Requested codes absent from a partial
#'   match are reported in the `fips_missing` attribute and a warning.
#' @seealso [read_fips_csv()], [aoi_from_sf()]
#' @examples
#' \dontrun{
#' counties <- aoi_us_counties(c("08031", "06071"))
#' }
#' @export
aoi_us_counties <- function(fips = NULL, year = 2020,
                            resolution = c("500k", "5m", "20m"),
                            territories = c("PR", "none", "all"),
                            src = NULL, cache = TRUE, quiet = FALSE) {
  resolution <- match.arg(resolution)
  territories <- match.arg(territories)
  if (!is.numeric(year) || length(year) != 1L ||
      !is.finite(year) || year < 2020 || year > 2024 || year != trunc(year)) {
    cli::cli_abort("{.arg year} must be an integer between 2020 and 2024.")
  }

  island_states <- c("60", "66", "69", "78")
  excluded <- switch(territories, none = c("72", island_states),
                     PR = island_states, all = character())
  if (!is.null(fips)) {
    fips <- .fips_normalize(fips)
    if (!length(fips)) {
      cli::cli_abort("No FIPS codes were provided after normalization.")
    }
    excluded_fips <- fips[substr(fips, 1L, 2L) %in% excluded]
    if (length(excluded_fips)) {
      preview <- .fips_preview(excluded_fips)
      cli::cli_abort("FIPS codes excluded by territories = {territories}: {.val {preview}}.")
    }
    if (year < 2022 && any(grepl("^091[1-9][0-9]$", fips))) {
      cli::cli_abort("Connecticut planning-region codes require {.code year = 2022} or later.")
    }
    if (year >= 2022 && any(grepl("^090(0[1-9]|1[0-5])$", fips))) {
      cli::cli_abort("Use an earlier vintage such as {.code year = 2020} for Connecticut county codes.")
    }
  }

  if (is.null(src)) {
    src <- .census_county_file(year, resolution, cache, quiet)
  }
  counties <- .read_county_boundaries(src)

  counties$GEOID <- as.character(counties$GEOID)
  counties <- counties[!counties$STATEFP %in% excluded, , drop = FALSE]
  missing <- character()
  if (is.null(fips)) {
    counties <- counties[order(counties$GEOID), , drop = FALSE]
  } else {
    rows <- match(fips, counties$GEOID)
    missing <- fips[is.na(rows)]
    if (length(missing)) {
      preview <- .fips_preview(missing)
      if (all(is.na(rows))) {
        cli::cli_abort("No requested FIPS codes were found: {.val {preview}}.")
      }
      cli::cli_warn("FIPS codes missing from the county file: {.val {preview}}.")
    }
    counties <- counties[rows[!is.na(rows)], , drop = FALSE]
  }

  state <- as.character(counties$STATEFP)
  if (!"NAMELSAD" %in% names(counties)) {
    counties$NAMELSAD <- rep(NA_character_, nrow(counties))
  }
  if (!"STUSPS" %in% names(counties)) {
    counties$STUSPS <- unname(.state_fips[state])
  }
  regions <- c("02" = "ak", "15" = "hi", "72" = "pr")
  regions[island_states] <- "island"
  counties$region <- unname(regions[state])
  counties$region[is.na(counties$region)] <- "conus"
  columns <- c("GEOID", "STATEFP", "COUNTYFP", "NAME", "NAMELSAD",
               "STUSPS", "ALAND", "AWATER", "region")
  out <- aoi_from_sf(counties[, columns, drop = FALSE], id_col = "GEOID")
  attr(out, "traktr_aoi_params") <- c(
    attr(out, "traktr_aoi_params"),
    list(source = "census_cb", year = year, resolution = resolution,
         territories = territories)
  )
  if (length(missing)) {
    attr(out, "fips_missing") <- missing
  }
  out
}

.read_county_boundaries <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("County boundary file does not exist: {.path {path}}.")
  }
  if (tolower(tools::file_ext(path)) == "zip") {
    unpacked <- tempfile("traktr-census-")
    on.exit(unlink(unpacked, recursive = TRUE), add = TRUE)
    tryCatch(utils::unzip(path, exdir = unpacked), warning = function(warning) {
      cli::cli_abort("Could not unpack the county archive at {.path {path}}.", parent = warning)
    })
    path <- list.files(unpacked, pattern = "\\.shp$", full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)
    if (length(path) != 1L) {
      cli::cli_abort("The county archive must contain one shapefile.")
    }
  }
  counties <- sf::st_read(path, quiet = TRUE)
  required <- c("GEOID", "STATEFP", "COUNTYFP", "NAME", "ALAND", "AWATER")
  absent <- setdiff(required, names(counties))
  if (length(absent)) {
    cli::cli_abort("County boundary data is missing required columns: {.val {absent}}.")
  }
  counties
}

.census_county_file <- function(year, resolution, cache, quiet) {
  directory <- file.path(tools::R_user_dir("traktr", "cache"), "census")
  archive <- sprintf("cb_%d_us_county_%s.zip", year, resolution)
  path <- file.path(directory, archive)
  if (cache && file.exists(path)) {
    cached <- tryCatch(.read_county_boundaries(path), error = identity)
    if (!inherits(cached, "error")) {
      return(path)
    }
    cli::cli_warn("The cached county archive at {.path {path}} is unreadable. Downloading it again.")
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  url <- sprintf("https://www2.census.gov/geo/tiger/GENZ%d/shp/%s", year, archive)
  if (!quiet) {
    cli::cli_inform("Downloading Census county boundaries ({year}, {resolution}).")
  }
  temporary <- tempfile(tmpdir = directory, fileext = ".zip")
  on.exit(unlink(temporary), add = TRUE)
  curl::curl_download(url, temporary, quiet = quiet, mode = "wb")
  tryCatch(.read_county_boundaries(temporary), error = function(error) {
    cli::cli_abort("The downloaded county archive is unreadable.", parent = error)
  })
  if (!file.rename(temporary, path)) {
    cli::cli_abort("Could not save the county archive to {.path {path}}.")
  }
  path
}
