#' Resolve AOI from input or object attribute
#'
#' Returns an AOI \code{sf} object from the explicit \code{aoi_obj} argument
#' or from \code{attr(df, "traktr_aoi")}. Validates that the AOI is an \code{sf}
#' polygon/multipolygon with a defined CRS. Aborts with a clear message if the
#' AOI is missing or has no CRS.
#'
#' @param df An object that may carry \code{attr(df, "traktr_aoi")}.
#' @param aoi_obj Optional \code{sf} polygon/multipolygon to use instead of the
#'   attribute.
#' @param error_call Passed to \code{cli::cli_abort()} to attribute the error
#'   to the caller; typically \code{rlang::caller_env()}.
#'
#' @return An \code{sf} object representing the AOI.
#'
#' @details
#' If \code{aoi_obj} is \code{NULL}, the function looks up
#' \code{attr(df, "traktr_aoi")}. The AOI must have a non-NA CRS.
#'
#' @seealso \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @keywords internal
#' @noRd
.resolve_aoi <- function(df, aoi_obj, error_call = rlang::caller_env()) {
  x <- if (!is.null(aoi_obj)) aoi_obj else attr(df, "traktr_aoi")
  if (is.null(x) || !inherits(x, "sf")) {
    cli::cli_abort(
      "No AOI found: pass `aoi` or call a mutate_aoi_*() function first.",
      call = error_call
    )
  }
  if (is.na(sf::st_crs(x))) {
    cli::cli_abort("`aoi` has no CRS; assign one first.", call = error_call)
  }
  x
}

#' Store layer outputs and parameters on a traktr object
#'
#' Inserts or updates a named entry in \code{attr(df, "traktr_layers")}, storing
#' the layer index (e.g., an \code{sf} tibble) and the parameters used to create
#' it. Also ensures class \code{"traktr_df"} is present on \code{df}.
#'
#' @param df A traktr-like object.
#' @param name Character scalar; the layer name (e.g., "ndvi", "worldcover").
#' @param index The per-layer index object to store (e.g., \code{sf} or tibble).
#' @param params A named list of parameters used to build the layer.
#'
#' @return \code{df} with updated \code{traktr_layers} attribute (and class).
#'
#' @details
#' The \code{traktr_layers} attribute is a list. Each element is a list with
#' fields \code{index} and \code{params}. Existing entries with the same
#' \code{name} are replaced.
#'
#' @seealso \code{\link{fetch_ndvi}}, \code{\link{fetch_landuse}}
#'
#' @keywords internal
#' @noRd
.set_traktr_layer <- function(df, name, index, params) {
  layers <- attr(df, "traktr_layers")
  if (is.null(layers)) layers <- list()
  layers[[name]] <- list(index = index, params = params)
  attr(df, "traktr_layers") <- layers
  class(df) <- union("traktr_df", class(df))
  df
}

#' Sanitize numeric values for stable filenames
#'
#' Formats numeric-like values without scientific notation and replaces hyphens
#' and dots with underscores, producing cache-friendly filename components.
#'
#' @param x A numeric or character vector to sanitize.
#'
#' @return A character vector with \code{"-"} and \code{"."} replaced by
#'   \code{"_"} after non-scientific formatting.
#'
#' @examples
#' sanitize_num(c(-7.5, 12.0, 1e-6))
#' # " _7_5", "12_0", "0_000001"
#'
#' @keywords internal
#' @noRd
.sanitize_num <- function(x) {
  gsub("[-.]", "_", format(x, scientific = FALSE))
}


#' AOI helpers (internal)
#'
#' Internal utilities used by AOI mutate/attach functions for validation,
#' raster file collection, deduplication, and attribute restoration.
#'
#' These functions are not exported and are intended for internal use only.
#'
#' @keywords internal
#' @rdname aoi-helpers
#' @noRd
.aoi_detect <- function(df) {
  aoi <- attr(df, "traktr_aoi")
  if (is.null(aoi) || !"sf" %in% class(aoi)) {
    cli::cli_abort("No AOI found: call a mutate_aoi_*() function first.")
  }
  params <- attr(df, "traktr_aoi_params")
  id_col <- if (!is.null(params) && !is.null(params$id_col)) {
    params$id_col
  } else {
    setdiff(names(aoi), attr(aoi, "sf_column"))[1]
  }
  if (is.null(id_col) || !id_col %in% names(aoi)) {
    cli::cli_abort("Cannot determine AOI id column; set attr(df, 'traktr_aoi_params')$id_col.")
  }
  if (!id_col %in% names(df)) {
    cli::cli_abort(paste0("Input df does not contain AOI id column '", id_col, "'."))
  }
  list(aoi = aoi, id_col = id_col)
}

#' @description
#' \code{.collect_raster_files()} finds raster files either from a recorded
#' attribute (e.g., \code{"ndvi_aoi"}, \code{"worldcover_tiles"}) or from a
#' directory on disk.
#'
#' @param df A data object carrying tile metadata as attributes.
#' @param attr_name Name of the attribute containing tile metadata.
#' @param out_dir Optional directory path where raster files are stored.
#'
#' @return A character vector of raster file paths.
#'
#' @rdname aoi-helpers
#' @noRd
.collect_raster_files <- function(df, attr_name, out_dir) {
  tiles <- attr(df, attr_name)
  if (is.null(tiles)) {
    if (is.null(out_dir)) {
      cli::cli_abort("No recorded tiles; please supply 'out_dir'.")
    }
    files <- list.files(out_dir, "\\.tif(f)?$", full.names = TRUE)
    if (!length(files)) cli::cli_abort("No raster files found in 'out_dir'.")
  } else {
    files <- tiles$file[!is.na(tiles$file)]
    if (!length(files)) cli::cli_abort(paste0("No files recorded in attr(df, '", attr_name, "')."))
  }
  files
}

#' @description
#' \code{.aoi_deduplicate()} reduces AOIs to unique geometries for efficient
#' extraction, with special handling for H3 cells (unique by ID). It also
#' creates a mapping table back to the full AOI set so results can be joined
#' correctly.
#'
#' @param aoi An \code{sf} object of AOI polygons.
#' @param id_col Column name containing AOI identifiers.
#' @param deduplicate Logical. If \code{TRUE}, geometrically identical AOIs
#'   are deduplicated before extraction.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{aoi_unique}: deduplicated AOIs (with \code{.aoi_row}).
#'   \item \code{map_back}: table mapping original AOIs to \code{.aoi_row}.
#' }
#'
#' @rdname aoi-helpers
#' @noRd
.aoi_deduplicate <- function(aoi, id_col, deduplicate = TRUE) {
  aoi <- sf::st_make_valid(aoi)
  aoi4326 <- sf::st_transform(aoi, 4326)
  
  if (identical(id_col, "h3") && "h3" %in% names(aoi4326)) {
    # H3 cells are intrinsically unique by id
    aoi_unique <- aoi4326[!duplicated(aoi4326$h3), , drop = FALSE]
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    map_back <- sf::st_drop_geometry(aoi4326)[, c(id_col), drop = FALSE]
    map_back$.aoi_row <- match(map_back[[id_col]], aoi_unique[[id_col]])
    return(list(aoi_unique = aoi_unique, map_back = map_back))
  }
  
  if (deduplicate) {
    aoi4326$.__wkb <- vapply(sf::st_as_binary(aoi4326$geometry), paste, "", collapse = "")
    aoi_unique <- aoi4326[!duplicated(aoi4326$.__wkb), , drop = FALSE]
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    
    df1 <- sf::st_drop_geometry(aoi4326)[, c(id_col, ".__wkb"), drop = FALSE]
    df2 <- sf::st_drop_geometry(aoi_unique)[, c(".__wkb", ".aoi_row"), drop = FALSE]
    map_back <- merge(df1, df2, by = ".__wkb", all.x = TRUE, sort = FALSE)
    map_back$.__wkb <- NULL
  } else {
    aoi_unique <- aoi4326
    aoi_unique$.aoi_row <- seq_len(nrow(aoi_unique))
    map_back <- sf::st_drop_geometry(aoi_unique)[, c(id_col, ".aoi_row"), drop = FALSE]
  }
  
  list(aoi_unique = aoi_unique, map_back = map_back)
}

#' @description
#' \code{.restore_attrs_class()} copies attributes and classes from the
#' original \code{df} to a modified version after joins or transformations.
#'
#' @param df The original input object.
#' @param out The modified output object.
#'
#' @return The output object with attributes and class restored.
#'
#' @rdname aoi-helpers
#' @noRd
.restore_attrs_class <- function(df, out) {
  attrs <- attributes(df)
  attrs <- attrs[!names(attrs) %in% c("names", "row.names", "class")]
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  class(out) <- class(df)
  out
}
