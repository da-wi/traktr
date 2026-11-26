args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop("Usage: Rscript validation/validate_ecolozip_county.R <us_csv> [county_csv] [out_dir]")
}
# us_csv: ecolo-zip US.csv from OSF (https://osf.io/download/gysr9/)
paths <- c(
  NA_character_,
  "testoutput/county-run/us_counties_aster_2022.csv",
  "testoutput/validation"
)
paths[seq_along(args)] <- args
us_csv <- paths[1L]
county_csv <- paths[2L]
out_dir <- paths[3L]

state_codes <- c(
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
state_fips <- setNames(names(state_codes), state_codes)

zip_keys <- function(state, county) {
  prefix <- unname(state_fips[state])
  geoid <- paste0(prefix, sprintf("%03d", county))
  status <- rep("eligible", length(county))
  status[prefix %in% "09"] <- "connecticut"
  status[is.na(prefix)] <- "unknown_state"
  status[is.na(county)] <- "missing_county"
  geoid[status != "eligible"] <- NA_character_
  data.frame(GEOID = geoid, status = status)
}

icc_2_1 <- function(x) {
  n <- nrow(x)
  k <- ncol(x)
  if (n < 2L || k < 2L) {
    return(NA_real_)
  }
  row_mean <- rowMeans(x)
  col_mean <- colMeans(x)
  ms_row <- k * var(row_mean)
  ms_col <- n * var(col_mean)
  residual <- x - outer(row_mean, col_mean, "+") + mean(x)
  ms_error <- sum(residual^2) / ((n - 1L) * (k - 1L))
  denominator <- ms_row + (k - 1L) * ms_error + k * (ms_col - ms_error) / n
  if (denominator == 0) {
    return(NA_real_)
  }
  (ms_row - ms_error) / denominator
}

# Shrout and Fleiss (1979), Table 2, six subjects rated by four judges.
ratings <- rbind(
  c(9, 2, 5, 8), c(6, 1, 3, 2), c(8, 4, 6, 8),
  c(7, 1, 2, 6), c(10, 5, 6, 9), c(6, 2, 4, 7)
)
test_keys <- zip_keys(c("AK", "CT"), c(13L, 1L))
stopifnot(
  round(icc_2_1(ratings), 4L) == 0.2898,
  identical(test_keys$GEOID, c("02013", NA_character_)),
  identical(test_keys$status, c("eligible", "connecticut"))
)

primary <- c(
  "elevation_mean", "ndvi_mean_mean", "ndvi_sd_mean", "temp_mean_mean",
  "temp_sd_mean", "observation_number_mean", "lwm_water_perc"
)
variables <- c(primary, "elevation_sd", "elevation_max", "ndvi_mean_sd")
elevation <- startsWith(variables, "elevation_")
zip_columns <- paste0(ifelse(elevation, "ASTER_zip_", "ASTER_"), variables)
county_columns <- paste0(ifelse(elevation, "ASTER_county_", "ASTER_"), variables)
targets <- c(
  elevation_mean = 0.95, elevation_sd = 0.80, ndvi_mean_mean = 0.85,
  temp_mean_mean = 0.85, lwm_water_perc = 0.60
)

zip <- read.delim(
  us_csv, quote = "", colClasses = c(admin_code1 = "character", admin_code2 = "integer")
)
county <- read.csv(county_csv, colClasses = c(GEOID = "character"))
metadata <- c("GEOID", "NAME", "STUSPS", "county_area_km2_geodesic")
stopifnot(
  all(c("admin_code1", "admin_code2", zip_columns) %in% names(zip)),
  all(c(metadata, county_columns) %in% names(county)),
  !anyNA(county$GEOID), all(grepl("^[0-9]{5}$", county$GEOID)),
  !anyDuplicated(county$GEOID)
)
keys <- zip_keys(zip$admin_code1, zip$admin_code2)
keys$status[keys$status == "eligible" & !keys$GEOID %in% county$GEOID] <- "unmatched"
unmatched_zip <- table(keys$GEOID[keys$status == "unmatched"])
unknown_states <- table(zip$admin_code1[keys$status == "unknown_state"], useNA = "ifany")
included <- keys$status == "eligible"

zip_groups <- split(which(included), keys$GEOID[included])
zip_means <- vapply(zip_columns, function(column) {
  vapply(zip_groups, function(rows) mean(zip[[column]][rows], na.rm = TRUE), numeric(1L))
}, numeric(length(zip_groups)))
zip_means <- matrix(
  zip_means, nrow = length(zip_groups), ncol = length(variables),
  dimnames = list(NULL, variables)
)
matched <- county[match(names(zip_groups), county$GEOID), metadata]
names(matched) <- c("GEOID", "name", "state", "county_area_km2_geodesic")
matched$n_zip <- unname(lengths(zip_groups))
county_only <- sort(setdiff(county$GEOID[!startsWith(county$GEOID, "09")], matched$GEOID))

county_values <- county[match(matched$GEOID, county$GEOID), county_columns]
hawaii <- startsWith(matched$GEOID, "15")
pairs <- setNames(lapply(seq_along(variables), function(i) {
  ecolo <- zip_means[, i]
  traktr <- county_values[[i]]
  comparable <- if (elevation[i]) !hawaii else rep(TRUE, nrow(matched))
  keep <- comparable & is.finite(ecolo) & is.finite(traktr)
  data.frame(
    matched[keep, ], ecolo = ecolo[keep], traktr = traktr[keep],
    delta = traktr[keep] - ecolo[keep]
  )
}), variables)

correlation <- function(pair, method = "pearson") {
  if (nrow(pair) < 2L || sd(pair$ecolo) == 0 || sd(pair$traktr) == 0) {
    return(NA_real_)
  }
  cor(pair$ecolo, pair$traktr, method = method)
}

metrics <- do.call(rbind, lapply(variables, function(variable) {
  pair <- pairs[[variable]]
  target <- unname(targets[variable])
  pearson <- correlation(pair)
  bias <- if (nrow(pair)) mean(pair$delta) else NA_real_
  spread <- sd(pair$delta)
  data.frame(
    variable = variable, n = nrow(pair), pearson_r = pearson,
    spearman_rho = correlation(pair, "spearman"),
    icc_2_1 = icc_2_1(as.matrix(pair[c("ecolo", "traktr")])),
    bias = bias, loa_lower = bias - 1.96 * spread, loa_upper = bias + 1.96 * spread,
    median_abs_delta = median(abs(pair$delta)), target = target,
    meets_target = pearson >= target
  )
}))

area <- matched$county_area_km2_geodesic
quartile_breaks <- quantile(area, c(0.25, 0.5, 0.75), na.rm = TRUE)
area_quartile <- sprintf("Q%d", 1L + rowSums(outer(area, quartile_breaks, ">")))
area_quartile[!is.finite(area)] <- NA_character_
zip_band <- cut(matched$n_zip, c(0, 1, 4, Inf), labels = c("1", "2–4", "≥5"))
stratification <- do.call(rbind, lapply(c("elevation_mean", "ndvi_mean_mean"), function(variable) {
  pair <- pairs[[variable]]
  rows <- match(pair$GEOID, matched$GEOID)
  groups <- list(
    "Area quartile" = factor(area_quartile[rows], levels = paste0("Q", 1:4)),
    "ZIP count" = zip_band[rows]
  )
  do.call(rbind, lapply(names(groups), function(stratum) {
    group <- groups[[stratum]]
    do.call(rbind, lapply(levels(group), function(level) {
      subset <- pair[which(group == level), ]
      data.frame(
        variable = variable, stratum = stratum, group = level, n = nrow(subset),
        pearson_r = correlation(subset), median_abs_delta = median(abs(subset$delta))
      )
    }))
  }))
}))

outliers <- do.call(rbind, lapply(primary, function(variable) {
  pair <- pairs[[variable]]
  spread <- sd(pair$delta)
  score <- if (is.finite(spread) && spread > 0) abs(pair$delta) / spread else NA_real_
  # The denominator is constant within a variable, so ranking by |delta| is equivalent.
  pair$abs_delta_over_sd <- rep_len(score, nrow(pair))
  top <- head(order(-abs(pair$delta), pair$GEOID), 10L)
  pair$variable <- rep(variable, nrow(pair))
  pair[top, c("variable", names(matched), "ecolo", "traktr", "delta", "abs_delta_over_sd")]
}))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE, na = "")
write.csv(outliers, file.path(out_dir, "outliers.csv"), row.names = FALSE, na = "")

save_plot <- function(file, draw) {
  png(file.path(out_dir, file), width = 1600, height = 1200, res = 150)
  on.exit(dev.off())
  draw()
}

scatter_titles <- c(
  elevation_mean = "Mean elevation", ndvi_mean_mean = "Mean NDVI",
  temp_mean_mean = "Mean land surface temperature", lwm_water_perc = "Water percentage"
)
save_plot("scatter.png", function() {
  par(mfrow = c(2, 2))
  for (variable in names(scatter_titles)) {
    pair <- pairs[[variable]]
    limits <- if (nrow(pair)) range(pair$ecolo, pair$traktr) else c(0, 1)
    plot(
      pair$ecolo, pair$traktr, xlim = limits, ylim = limits, asp = 1,
      main = scatter_titles[[variable]], xlab = "ecolo-zip", ylab = "traktr",
      pch = 16, cex = 0.45, col = adjustcolor("#246985", alpha.f = 0.4)
    )
    abline(a = 0, b = 1, col = "#A8472C", lwd = 2)
  }
})

save_plot("bland_altman_elevation.png", function() {
  pair <- pairs$elevation_mean
  elevation_metrics <- metrics[metrics$variable == "elevation_mean", ]
  lines <- unlist(elevation_metrics[c("bias", "loa_lower", "loa_upper")], use.names = FALSE)
  pair_mean <- (pair$traktr + pair$ecolo) / 2
  x_limits <- if (nrow(pair)) range(pair_mean) else c(0, 1)
  y_limits <- if (nrow(pair)) range(c(pair$delta, lines), finite = TRUE) else c(-1, 1)
  plot(
    pair_mean, pair$delta, xlim = x_limits, ylim = y_limits,
    main = "Bland–Altman plot for mean elevation",
    xlab = "Mean of traktr and ecolo-zip (m)", ylab = "traktr minus ecolo-zip (m)",
    pch = 16, cex = 0.55, col = adjustcolor("#246985", alpha.f = 0.4)
  )
  abline(h = lines[1L], col = "#A8472C", lwd = 2)
  abline(h = lines[2:3], col = "#A8472C", lwd = 2, lty = 2)
  legend("topright", c("Bias", "95% limits of agreement"), col = "#A8472C", lty = 1:2, lwd = 2)
})

md_table <- function(data) {
  row <- function(values) paste0("| ", paste(values, collapse = " | "), " |")
  data[] <- lapply(data, function(column) {
    values <- if (is.numeric(column)) as.character(signif(column, 5L)) else as.character(column)
    values[is.na(column)] <- "—"
    gsub("|", "&#124;", values, fixed = TRUE)
  })
  c(row(names(data)), row(rep("---", ncol(data))), apply(data, 1L, row))
}

zip_exclusions <- data.frame(
  reason = c("Missing county code", "Unknown state code", "Connecticut", "No traktr GEOID"),
  zip_rows = as.integer(table(factor(
    keys$status, levels = c("missing_county", "unknown_state", "connecticut", "unmatched")
  )))
)
unknown_labels <- names(unknown_states)
unknown_labels[is.na(unknown_labels)] <- "NA"
unknown_details <- sprintf("%s (%d)", unknown_labels, as.integer(unknown_states))
no_match_details <- sprintf("%s (%d)", names(unmatched_zip), as.integer(unmatched_zip))
list_or_none <- function(values) {
  if (length(values)) paste(values, collapse = ", ") else "None"
}
missing_pairs <- data.frame(
  variable = variables,
  hawaii_excluded = ifelse(elevation, sum(hawaii), 0L),
  missing_pair_values = nrow(matched) - ifelse(elevation, sum(hawaii), 0L) - metrics$n
)
input_paths <- c(us_csv, county_csv)

writeLines(c(
  "# Gate 3 — County comparison with ecolo-zip", "",
  "## Inputs", "",
  md_table(data.frame(
    source = c("ecolo-zip", "traktr"), file = basename(input_paths),
    rows = c(nrow(zip), nrow(county)), MD5 = unname(tools::md5sum(input_paths))
  )), "",
  "County GEOIDs combine the state FIPS prefix and three-digit county code.",
  "Each ecolo-zip county value is the unweighted mean of its available ZIP values.",
  "ZIP counts include all eligible rows in a matched county.", "",
  "## Exclusions", "",
  "ZIP-row exclusions are counted once, in table order.", "",
  md_table(zip_exclusions), "",
  sprintf("- Unknown state codes among rows with a county code include %s.", list_or_none(unknown_details)),
  sprintf("- %d ecolo-zip GEOIDs lack a traktr match. The ZIP row count follows each GEOID in parentheses. %s.",
          length(unmatched_zip), list_or_none(no_match_details)),
  sprintf("- The comparison excludes %d traktr Connecticut GEOIDs.", sum(startsWith(county$GEOID, "09"))),
  sprintf("- Outside Connecticut, %d traktr GEOIDs lack an eligible ecolo-zip match. %s.",
          length(county_only), list_or_none(county_only)),
  sprintf("- The comparison matched %d counties covering %d ZIP rows.", nrow(matched), sum(matched$n_zip)),
  "", "The table shows county exclusions among matched GEOIDs for each variable.", "",
  md_table(missing_pairs), "",
  "## Metrics", "",
  "Differences are traktr minus ecolo-zip. Limits of agreement are bias ± 1.96 SD.",
  "ICC(2,1) measures single-measure absolute agreement in a two-way random-effects model.",
  "Targets apply to Pearson r and are descriptive.", "",
  md_table(metrics), "",
  "## Stratification", "",
  sprintf(
    "Area quartiles use all matched county areas (km²), with cut points %s. Ties stay in the lower quartile.",
    paste(signif(quartile_breaks, 5L), collapse = ", ")
  ),
  "Each row uses the available pairs for that variable.", "",
  md_table(stratification), "",
  "## Largest differences in mean elevation", "",
  "Ranking uses |difference| divided by the sample SD of paired differences.",
  "The CSV includes the ten largest differences for each primary variable.", "",
  md_table(outliers[outliers$variable == "elevation_mean", -1L]), "",
  "## Known differences", "",
  "- traktr weights ASTER cells by area across county polygons. ecolo-zip samples 140 spiderweb points within 20 miles of ZIP centroids, concentrating coverage near population centers.",
  "- The modifiable areal unit problem (MAUP) affects SD and maxima. County polygons cover more terrain and can capture GDEM sinks, so elevation minima are excluded.",
  "- Elevation sums are excluded because the traktr polygon output leaves them missing and totals depend on sampling support.",
  "- ecolo-zip reports elevation in Hawaii. The comparison excludes Hawaii from elevation metrics because traktr has missing elevation values there.",
  "- ecolo-zip uses historical Connecticut counties and traktr uses Census 2022 planning regions. Other boundary vintages can leave unmatched GEOIDs, including changes in Alaska and Virginia.",
  "- Missing percentages are excluded. ecolo-zip counts missing sampled points and can report NA instead of zero. traktr measures missing area, including area outside granule coverage.",
  "- Granule inventories differ between the ecolo-zip 2021 retrieval and the traktr 2026 retrieval of ASTER GED v003. This can affect observation counts and NDVI."
), file.path(out_dir, "report.md"))

print(metrics, row.names = FALSE)
