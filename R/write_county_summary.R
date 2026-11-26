#' Export an ASTER GED county summary
#'
#' Writes the non-geometry columns in their original order, with an optional
#' data dictionary and provenance report.
#'
#' @param df Output of [attach_aster_ged()]. List columns are unsupported.
#' @param path Output file stem. A trailing `.csv` or `.parquet` is removed.
#'   Parent directories are created as needed.
#' @param formats Character vector selecting `"csv"`, `"parquet"` or both.
#'   Parquet output requires the nanoparquet package.
#' @param dictionary Write a `_dictionary.csv` file with one entry per column.
#' @param provenance Write a `_PROVENANCE.md` file with source and tile metadata.
#' @param na String used for missing values in the summary CSV.
#' @param overwrite Allow replacement of existing output files.
#'
#' @details All output files are staged in the destination directory before
#'   publication. If publication fails, this call removes the files it
#'   published and restores any originals from backups. Without overwrite,
#'   publication uses hard links and requires a filesystem that supports them.
#'   Files that another export replaced in the meantime are left in place;
#'   concurrent exports to the same `path` are not coordinated otherwise.
#'   CSV output is UTF-8 and quotes character values.
#'   When reading a CSV, import GEOID as character to retain leading zeros.
#'
#' @return Invisibly, a named character vector of the written paths.
#' @examples
#' \dontrun{
#' counties <- attach_aster_ged(counties_with_tiles)
#' write_county_summary(counties, "output/us_counties_aster_2020")
#' }
#' @export
write_county_summary <- function(df, path, formats = "csv", dictionary = TRUE, provenance = TRUE,
                                 na = "", overwrite = FALSE) {
  layers <- attr(df, "traktr_layers", exact = TRUE)
  summary <- layers$aster_ged_summary
  if (is.null(summary)) {
    cli::cli_abort("Missing ASTER GED summary: run attach_aster_ged() first.")
  }
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    cli::cli_abort("{.arg path} must be a non-empty output file stem.")
  }
  stem <- sub("\\.(csv|parquet)$", "", path, ignore.case = TRUE)
  if (!nzchar(stem) || grepl("[/\\\\]$", stem)) {
    cli::cli_abort("{.arg path} must be a non-empty output file stem.")
  }
  if (!is.character(formats) || !all(formats %in% c("csv", "parquet"))) {
    cli::cli_abort("{.arg formats} must contain only {.val csv} or {.val parquet}.")
  }
  if ("parquet" %in% formats && !.has_nanoparquet()) {
    cli::cli_abort("Parquet output requires {.pkg nanoparquet}. Install it with install.packages(\"nanoparquet\").")
  }

  geometry <- vapply(df, inherits, logical(1), what = "sfc") |
    names(df) %in% attr(df, "sf_column", exact = TRUE)
  columns <- as.list(df)[!geometry]
  list_cols <- names(columns)[vapply(columns, is.list, logical(1))]
  if (length(list_cols)) {
    cli::cli_abort("Cannot export list columns: {.val {list_cols}}.")
  }
  data <- structure(columns, class = "data.frame", row.names = seq_len(nrow(df)))

  suffixes <- c(csv = ".csv", parquet = ".parquet", dictionary = "_dictionary.csv",
                provenance = "_PROVENANCE.md")
  outputs <- c(unique(formats), if (isTRUE(dictionary)) "dictionary",
               if (isTRUE(provenance)) "provenance")
  if (!length(outputs)) return(invisible(character()))

  stem <- path.expand(stem)
  paths <- stats::setNames(paste0(stem, suffixes[outputs]), outputs)
  overwrite <- isTRUE(overwrite)
  .check_county_summary_targets(paths, overwrite)

  directory <- dirname(stem)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(directory)) cli::cli_abort("Could not create output directory {.path {directory}}.")

  staged <- stats::setNames(vapply(outputs, function(output) tempfile(tmpdir = directory),
                                   character(1)), outputs)
  on.exit(unlink(staged, expand = FALSE), add = TRUE)
  if ("csv" %in% outputs) {
    utils::write.csv(data, file = staged[["csv"]], row.names = FALSE, na = na, fileEncoding = "UTF-8")
  }
  if ("parquet" %in% outputs) nanoparquet::write_parquet(data, staged[["parquet"]])
  if ("dictionary" %in% outputs) {
    dictionary_data <- .county_summary_dictionary(names(data), summary$params,
                                                  attr(df, "traktr_aoi_params", exact = TRUE))
    utils::write.csv(dictionary_data, file = staged[["dictionary"]], row.names = FALSE,
                     na = "", fileEncoding = "UTF-8")
  }
  if ("provenance" %in% outputs) {
    writeLines(enc2utf8(.county_summary_provenance(df, layers)), staged[["provenance"]],
               useBytes = TRUE)
  }

  .publish_county_summary(staged, paths, overwrite)
  invisible(paths)
}

.check_county_summary_targets <- function(paths, overwrite) {
  links <- Sys.readlink(paths)
  symlink <- !is.na(links) & nzchar(links)
  present <- file.exists(paths) | symlink
  if (!overwrite && any(present)) {
    cli::cli_abort("Output files already exist: {.file {paths[present]}}. Set overwrite = TRUE to replace them.")
  }
  directories <- dir.exists(paths) & !symlink
  if (any(directories)) {
    cli::cli_abort("Output paths point to directories: {.path {paths[directories]}}.")
  }
  present
}

.publish_county_summary_file <- function(from, to, overwrite) {
  if (overwrite) file.rename(from, to) else file.link(from, to)
}

.publish_county_summary <- function(staged, paths, overwrite) {
  present <- .check_county_summary_targets(paths, overwrite)
  backups <- stats::setNames(character(), character())
  published <- character()
  published_md5 <- character()
  committed <- FALSE
  on.exit({
    if (committed) {
      unlink(backups, expand = FALSE)
    } else {
      # Only remove files that still hold this run's content; a concurrent
      # export may have replaced them in the meantime.
      ours <- published[(unname(tools::md5sum(published)) == published_md5) %in% TRUE]
      if (length(ours) && unlink(ours, expand = FALSE) != 0L) {
        cli::cli_inform("Could not remove all partially published files: {.file {published}}.")
      }
      for (output in names(backups)) {
        if (!suppressWarnings(file.rename(backups[[output]], paths[[output]]))) {
          cli::cli_inform("Could not restore {.file {paths[[output]]}}. Its backup is at {.file {backups[[output]]}}.")
        }
      }
    }
  }, add = TRUE)

  for (output in names(paths)[present]) {
    path <- paths[[output]]
    backup <- tempfile("county-backup-", tmpdir = dirname(path))
    message <- "Could not back up existing output file {.file {path}}."
    saved <- withCallingHandlers(file.rename(path, backup), warning = function(warning) {
      cli::cli_abort(c(message, i = conditionMessage(warning)))
    })
    if (!saved) cli::cli_abort(message)
    backups[[output]] <- backup
  }
  for (output in names(paths)) {
    path <- paths[[output]]
    message <- "Could not publish output file {.file {path}}."
    staged_md5 <- unname(tools::md5sum(staged[[output]]))
    saved <- withCallingHandlers(
      .publish_county_summary_file(staged[[output]], path, overwrite),
      warning = function(warning) cli::cli_abort(c(message, i = conditionMessage(warning)))
    )
    if (!saved) cli::cli_abort(message)
    published <- c(published, path)
    published_md5 <- c(published_md5, staged_md5)
  }
  committed <- TRUE
  invisible(NULL)
}

.has_nanoparquet <- function() requireNamespace("nanoparquet", quietly = TRUE)

.county_summary_dictionary <- function(columns, params, aoi_params) {
  blank <- rep("", length(columns))
  dict <- data.frame(
    name = columns, source = rep("input", length(columns)), band = blank, stat = blank,
    unit = blank, description = blank, ecolozip_name = blank, stringsAsFactors = FALSE
  )

  if (identical(aoi_params$source, "census_cb")) {
    descriptions <- c(
      GEOID = "Five-digit county FIPS code.", STATEFP = "Two-digit state FIPS code.",
      COUNTYFP = "Three-digit county FIPS code within the state.", NAME = "County name.",
      NAMELSAD = "County name with its legal or statistical area designation.",
      STUSPS = "Two-letter postal abbreviation for the state or territory.",
      ALAND = "Census land area.", AWATER = "Census water area.",
      region = "Region assigned by traktr."
    )
    identity <- columns %in% names(descriptions)
    dict$source[identity] <- paste("US Census cartographic boundary", aoi_params$year)
    dict$description[identity] <- unname(descriptions[columns[identity]])
    dict$unit[columns %in% c("ALAND", "AWATER")] <- "m\u00b2"
  }

  aster <- .aster_ged_dictionary(params$convert)
  selected <- (!aster$extra & aster$band %in% params$bands) |
    (aster$stat == "median" & aster$band %in% params$bands & "median" %in% params$extras) |
    (is.na(aster$band) & aster$stat %in% params$extras)
  aster <- aster[selected, , drop = FALSE]
  rows <- match(columns, aster[[paste0("name_", params$naming)]])
  matched <- !is.na(rows)
  dictionary_columns <- c("band", "stat", "unit", "description", "ecolozip_name")
  aster_columns <- c("band", "stat", "unit", "description", "name_ecolozip")
  dict[matched, dictionary_columns] <- aster[rows[matched], aster_columns, drop = FALSE]
  dict$source[matched] <- "ASTER GED AG100 v003"
  dict
}

.county_summary_provenance <- function(df, layers) {
  summary <- layers$aster_ged_summary
  fetch <- layers$aster_ged$params
  aoi <- attr(df, "traktr_aoi_params", exact = TRUE)
  fetch_aoi <- attr(df, "aster_ged_aoi", exact = TRUE)
  id_col <- if (is.null(fetch_aoi)) aoi$id_col else setdiff(names(fetch_aoi), attr(fetch_aoi, "sf_column"))
  missing <- attr(df, "fips_missing", exact = TRUE)
  iso_time <- function(x) format(x, tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")

  product <- if (is.null(fetch$product)) "AG100" else fetch$product
  version <- if (is.null(fetch$version)) "003" else fetch$version
  fields <- c("bands", "naming", "extras", "convert", "ecolozip_na_quirk")
  values <- vapply(summary$params[fields], paste, character(1), collapse = ", ")
  values[!nzchar(values)] <- "(none)"
  tiles <- summary$index
  exists <- file.exists(tiles$file) & !dir.exists(tiles$file)
  size <- if (is.null(tiles$size_bytes)) rep(NA_real_, nrow(tiles)) else tiles$size_bytes
  md5 <- rep("", nrow(tiles))
  size[exists] <- file.size(tiles$file[exists])
  md5[exists] <- unname(tools::md5sum(tiles$file[exists]))
  md5[is.na(md5)] <- ""
  tile_rows <- paste0("| ", tiles$tile_id, " | ", tiles$status, " | ",
                      ifelse(is.na(tiles$file), "", basename(tiles$file)), " | ",
                      ifelse(is.na(size), "", format(size, trim = TRUE, scientific = FALSE)),
                      " | ", md5, " |", recycle0 = TRUE)

  index <- layers$aster_ged$index
  unavailable <- as.data.frame(index)[index$status %in% c("missing_remote", "missing_local", "error"),
                                      c("tile_id", "status"), drop = FALSE]
  # Tiles found missing at extraction time take precedence over the fetch status
  at_attach <- layers$aster_ged_summary$params$unavailable
  if (!is.null(at_attach) && nrow(at_attach)) {
    unavailable <- unavailable[!unavailable$tile_id %in% at_attach$tile_id, , drop = FALSE]
    unavailable <- rbind(unavailable, at_attach[, c("tile_id", "status")])
  }
  map <- attr(df, "aster_ged_map", exact = TRUE)
  aoi_by_tile <- if (!is.null(map)) split(map$aoi_id, map$tile_id)
  affected <- vapply(as.character(unavailable$tile_id), function(tile_id) {
    paste(unique(aoi_by_tile[[tile_id]]), collapse = ", ")
  }, character(1))
  unavailable_rows <- paste0("| ", unavailable$tile_id, " | ", unavailable$status,
                             " | ", affected, " |", recycle0 = TRUE)

  c(
    "# ASTER GED county summary provenance", "",
    paste0("- traktr version: ", utils::packageVersion("traktr")),
    paste0("- R version: ", getRversion()),
    paste0("- Exported (UTC): ", iso_time(Sys.time())),
    paste0("- Product: ASTER GED ", product, " v", version),
    if (!is.null(fetch$concept_id)) paste0("- CMR collection ID: ", fetch$concept_id),
    if (!is.null(fetch$retrieved)) paste0("- Retrieved (UTC): ", iso_time(fetch$retrieved)),
    if (identical(aoi$source, "census_cb")) "- Boundary source: US Census cartographic boundary",
    if (!is.null(aoi$year)) paste0("- Boundary year: ", aoi$year),
    if (!is.null(aoi$resolution)) paste0("- Boundary resolution: ", aoi$resolution),
    paste0("- Rows: ", nrow(df)),
    if (length(id_col) == 1L) paste0("- Unique AOI IDs: ", length(unique(df[[id_col]]))),
    if (length(missing)) paste0("- Missing FIPS codes: ", paste(missing, collapse = ", ")),
    "", "## Summary parameters", "",
    paste0("- `", fields, "`: ", values),
    "", "## Tiles used", "",
    "| Tile ID | Status | File | Size (bytes) | MD5 |",
    "| --- | --- | --- | ---: | --- |",
    tile_rows,
    "", "## Tiles without data", "",
    "| Tile ID | Status | Affected AOI IDs |",
    "| --- | --- | --- |",
    unavailable_rows
  )
}
