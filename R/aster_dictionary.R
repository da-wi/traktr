.aster_ged_dictionary <- function(convert = FALSE) {
  bases <- c(
    elev = "ASTER_county_elevation", ndvi_mean = "ASTER_ndvi_mean",
    ndvi_sd = "ASTER_ndvi_sd", temp_mean = "ASTER_temp_mean",
    temp_sd = "ASTER_temp_sd", nobs = "ASTER_observation_number",
    lwm = "ASTER_lwm_water"
  )
  continuous <- c("mean", "min", "max", "sd")
  common <- c(continuous, "missing_perc")
  stats <- list(
    elev = c(continuous, "cin", "sum", "missing_perc"),
    ndvi_mean = common, ndvi_sd = common, temp_mean = common, temp_sd = common,
    nobs = c(continuous, "lt_4"),
    lwm = c("perc", "cin", "missing_perc")
  )
  row <- function(band, stat, extra, county, ecolozip, traktr, unit, description) {
    data.frame(band = band, stat = stat, extra = extra, name_county = county,
               name_ecolozip = ecolozip, name_traktr = traktr, unit = unit,
               description = description, stringsAsFactors = FALSE)
  }

  units <- c(elev = "m", ndvi_mean = "NDVI\u00d7100", ndvi_sd = "NDVI\u00d7100",
             temp_mean = "cK", temp_sd = "cK", nobs = "count", lwm = "%")
  if (isTRUE(convert)) {
    units[c("ndvi_mean", "ndvi_sd")] <- "NDVI"
    units[c("temp_mean", "temp_sd")] <- c("\u00b0C", "K")
  }
  descriptions <- c(
    mean = "Area-weighted mean of valid cells.",
    min = "Minimum among cells with valid values and positive covered area.",
    max = "Maximum among cells with valid values and positive covered area.",
    sd = "Area-weighted population standard deviation of valid cells.",
    cin = "Share of valid area at exactly zero metres.",
    sum = "Always NA in polygon mode because the point sum depends on the sample size.",
    missing_perc = "Share of area with missing values; includes county area in tiles without an AG100 granule.",
    lt_4 = "One if the area-weighted mean observation count is below four, zero otherwise.",
    perc = "Water share of valid area.",
    median = "Area-weighted median of valid cells."
  )

  band_rows <- function(band, stat, extra = FALSE) {
    county <- paste(bases[[band]], stat, sep = "_")
    ecolozip <- sub("^ASTER_county_elevation_", "ASTER_zip_elevation_", county)
    traktr <- paste("aster", band, ifelse(stat == "perc", "water_perc", stat), sep = "_")
    unit <- rep(units[[band]], length(stat))
    unit[stat %in% c("cin", "missing_perc")] <- "%"
    unit[stat == "lt_4"] <- "flag 0/1"
    if (isTRUE(convert) && band == "temp_mean") unit[stat == "sd"] <- "K"

    description <- unname(descriptions[stat])
    if (band == "lwm") description[stat == "cin"] <- "Share of area with valid values (100 minus the missing percentage)."
    if (band == "elev") {
      affected <- stat %in% c("min", "max", "sd")
      description[affected] <- paste(
        description[affected],
        "Sensitive to county size (MAUP). Elevation range and spread tend to grow with",
        "county area. Use county_area_km2_geodesic as a covariate."
      )
      extremes <- stat %in% c("min", "max")
      description[extremes] <- paste(
        description[extremes],
        "Sensitive to ASTER GDEM artefacts (pits near water bodies); prefer mean/median."
      )
    }
    row(band, stat, extra, county, ecolozip, traktr, unit, description)
  }

  core <- Map(band_rows, names(stats), stats)
  medians <- lapply(setdiff(names(stats), "lwm"), band_rows, stat = "median", extra = TRUE)
  do.call(rbind, unname(c(core, medians, list(
    row(NA_character_, "n_cells", TRUE, "n_cells", "n_cells", "aster_n_cells",
        "count", "Number of valid elevation cells, or valid cells in the first requested band."),
    row(NA_character_, "area", TRUE, "county_area_km2_geodesic",
        "county_area_km2_geodesic", "aoi_area_km2", "km\u00b2", "Geodesic area of the AOI.")
  ))))
}
