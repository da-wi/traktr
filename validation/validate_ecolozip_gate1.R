started <- proc.time()[["elapsed"]]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 4L) {
  stop("Usage: Rscript validation/validate_ecolozip_gate1.R <us_csv> [tile_dir] [out_dir] [workers]")
}
defaults <- c(NA_character_, "testoutput/county-run/aster_ged",
              "testoutput/validation/gate1", "1")
defaults[seq_along(args)] <- args
us_csv <- defaults[1L]
tile_dir <- defaults[2L]
out_dir <- defaults[3L]
workers <- suppressWarnings(as.integer(defaults[4L]))
if (is.na(workers) || workers < 1L || !grepl("^[0-9]+$", defaults[4L])) {
  stop("workers must be a positive integer")
}
if (workers > 1L && .Platform$OS.type == "windows") {
  stop("workers > 1 requires fork support")
}

devtools::load_all(quiet = TRUE)
spiderweb_ecolozip <- traktr:::.spiderweb_ecolozip
spiderweb_points <- traktr:::.spiderweb_points
spiderweb_tile_id <- traktr:::.spiderweb_tile_id
sample_points <- traktr:::.aster_ged_sample_points
aggregate_points <- traktr:::.spss_aggregate
web <- spiderweb_ecolozip()
points_per_zip <- 140L
seed <- 20231106L
sample_size <- 500L

covered_zips <- function(points, n_zips, tile_ids) {
  point_tiles <- spiderweb_tile_id(points$lat, points$lon)
  tabulate(points$id[point_tiles %in% tile_ids], nbins = n_zips) == points_per_zip
}

sample_zips <- function(n) {
  set.seed(seed, kind = "Mersenne-Twister", sample.kind = "Rejection")
  sample.int(n, sample_size)
}

elevation_pair_status <- function(state, published, reproduced) {
  finite <- is.finite(published) & is.finite(reproduced)
  documented <- state %in% "HI" &
    is.na(published) & !is.nan(published) & is.na(reproduced) & !is.nan(reproduced)
  ifelse(finite, "finite", ifelse(documented, "documented", "unexpected"))
}

passes_elevation_gate <- function(status, correlation, median_error) {
  expected <- sum(status != "documented")
  !any(status == "unexpected") &&
    100 * sum(status == "finite") >= 99 * expected &&
    is.finite(correlation) && correlation > 0.999 &&
    is.finite(median_error) && median_error < 2
}

test_points <- spiderweb_points(c(40.5, 40.5), c(-105.5, -105.05), web)
test_sample <- sample_zips(600L)
stopifnot(
  identical(tabulate(test_points$id), rep(points_per_zip, 2L)),
  identical(covered_zips(test_points, 2L, "41.-106"), c(TRUE, FALSE)),
  all(covered_zips(test_points, 2L, c("41.-106", "41.-105"))),
  length(unique(test_sample)) == sample_size,
  identical(test_sample, sample_zips(600L))
)
test_elevation <- c(100, 200, rep(NA_real_, 498L))
stopifnot(
  !passes_elevation_gate(elevation_pair_status("MA", test_elevation, test_elevation), 1, 0),
  passes_elevation_gate(elevation_pair_status(c("MA", "MA", "HI"),
                                               c(100, 200, NA), c(100, 200, NA)), 1, 0),
  !passes_elevation_gate(elevation_pair_status(c("MA", "MA", "HI"),
                                                c(100, 200, NA), c(100, 200, 300)), 1, 0),
  identical(elevation_pair_status("HI", c(NA, Inf), c(NaN, Inf)),
            rep("unexpected", 2L))
)

dictionary <- traktr:::.aster_ged_dictionary()
core <- dictionary[!dictionary$extra, ]
variables <- core$name_ecolozip
exact_variables <- variables[core$stat %in% c("min", "max", "cin", "sum", "lt_4")]
elevation <- "ASTER_zip_elevation_mean"
relief <- "ASTER_zip_elevation_sd"

# The published US.csv is tab-delimited.
reference <- read.delim(
  us_csv, quote = "",
  colClasses = c(postalcode = "character", latitude = "numeric", longitude = "numeric")
)
required <- c("postalcode", "place_name", "admin_code1", "latitude", "longitude", variables)
stopifnot(all(required %in% names(reference)))

in_states <- reference$admin_code1 %in% c(state.abb, "DC")
military <- !in_states &
  (grepl("^(APO|FPO|DPO)([[:space:]]|$)", reference$place_name, ignore.case = TRUE) |
   reference$admin_code1 %in% c("AA", "AE", "AP"))
finite_coords <- is.finite(reference$latitude) & is.finite(reference$longitude)
eligible <- in_states & finite_coords
zip <- reference[eligible, ]
stopifnot(!anyNA(zip$postalcode), all(grepl("^[0-9]{5}$", zip$postalcode)),
          !anyDuplicated(zip$postalcode))

if (!dir.exists(tile_dir)) stop("Tile directory does not exist")
tile_pattern <- "^ag100_v003_(-?[0-9]+)_(-?[0-9]+)[.]tif$"
tile_paths <- list.files(tile_dir, pattern = tile_pattern, full.names = TRUE)
tile_paths <- tile_paths[file.exists(tile_paths) & !dir.exists(tile_paths)]
if (!length(tile_paths)) stop("No AG100 v003 tiles found")
tile_ids <- sub(tile_pattern, "\\1.\\2", basename(tile_paths))

# Group ZIPs by centroid tile. V1 opens one raster at a time within each
# group, so the full point grid and sampled bands never accumulate in memory.
centroid_tiles <- spiderweb_tile_id(zip$latitude, zip$longitude)
centroid_tiles[is.na(centroid_tiles)] <- "outside"
groups <- split(seq_len(nrow(zip)), centroid_tiles)
process_group <- function(rows) {
  lat <- zip$latitude[rows]
  points <- spiderweb_points(lat, zip$longitude[rows], web, lat_ref = lat)
  covered <- covered_zips(points, length(rows), tile_ids)
  if (!any(covered)) return(data.frame(id = integer()))

  sampled <- sample_points(points[covered[points$id], ], tile_dir)
  stopifnot(!any(sampled$tile_missing))
  values <- aggregate_points(sampled, sampled$id)
  values$id <- rows[values$id]
  values
}

message(sprintf("Checking coverage and sampling %d ZIPs with %d worker(s).", nrow(zip), workers))
results <- if (workers > 1L) {
  parallel::mclapply(groups, process_group, mc.cores = workers)
} else {
  lapply(groups, process_group)
}
failed <- vapply(results, inherits, logical(1L), "try-error")
if (any(failed)) stop(results[[which(failed)[1L]]])
if (!all(vapply(results, is.data.frame, logical(1L)))) {
  stop("At least one sampling worker did not return a result")
}
values <- do.call(rbind, unname(results))
if (is.null(values) || nrow(values) < sample_size) {
  stop(sprintf("Gate 1 requires at least %d fully covered ZIPs.", sample_size))
}
values <- values[order(values$id), ]
zip <- zip[values$id, ]
values <- data.frame(ZIP = zip$postalcode, values[variables])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(values, file.path(out_dir, "zip_values.csv"), row.names = FALSE)

pearson <- function(x, y) {
  if (length(x) < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  cor(x, y)
}

compare_variables <- function(rows) {
  do.call(rbind, lapply(variables, function(variable) {
    published <- zip[[variable]][rows]
    reproduced <- values[[variable]][rows]
    paired <- is.finite(published) & is.finite(reproduced)
    published <- published[paired]
    reproduced <- reproduced[paired]
    delta <- abs(reproduced - published)
    data.frame(
      Variable = variable,
      n = length(delta),
      "Pearson r" = pearson(published, reproduced),
      "Median |delta|" = median(delta),
      "Share < 1e-9" = if (length(delta)) mean(delta < 1e-9) else NA_real_,
      "Max |delta|" = if (length(delta)) max(delta) else NA_real_,
      "Share exactly equal" = if (length(delta) && variable %in% exact_variables) {
        mean(reproduced == published)
      } else {
        NA_real_
      },
      check.names = FALSE
    )
  }))
}

sample_rows <- sample_zips(nrow(values))
all_rows <- seq_len(nrow(values))
sample_stats <- compare_variables(sample_rows)
all_stats <- compare_variables(all_rows)
gate_rows <- rbind(sample_stats[sample_stats$Variable == elevation, ],
                   all_stats[all_stats$Variable == elevation, ])
pair_status <- elevation_pair_status(zip$admin_code1, zip[[elevation]], values[[elevation]])
populations <- list(pair_status[sample_rows], pair_status)
passes <- mapply(passes_elevation_gate, populations, gate_rows[["Pearson r"]],
                 gate_rows[["Median |delta|"]])
documented_counts <- vapply(populations, function(status) sum(status == "documented"), integer(1L))
unexpected_counts <- vapply(populations, function(status) sum(status == "unexpected"), integer(1L))
verdict <- if (all(passes)) "PASS" else "FAIL"
gate <- data.frame(
  Population = c("500-ZIP sample", "All covered ZIPs"),
  Verdict = ifelse(passes, "PASS", "FAIL"),
  "Expected pairs" = lengths(populations) - documented_counts,
  "Unexpected missing or non-finite pairs" = unexpected_counts,
  gate_rows[c("n", "Pearson r", "Median |delta|", "Max |delta|")],
  check.names = FALSE
)

unexpected <- which(pair_status == "unexpected")
unexpected_pairs <- data.frame(
  ZIP = values$ZIP[unexpected], State = zip$admin_code1[unexpected],
  "Published mean (m)" = format(zip[[elevation]][unexpected], trim = TRUE, digits = 10L),
  "Reproduced mean (m)" = format(values[[elevation]][unexpected], trim = TRUE, digits = 10L),
  "In sample" = ifelse(unexpected %in% sample_rows, "Yes", "No"),
  check.names = FALSE
)
paired <- pair_status == "finite"
delta <- values[[elevation]] - zip[[elevation]]
outliers <- head(which(paired)[order(abs(delta[paired]), decreasing = TRUE)], 20L)
top <- data.frame(
  ZIP = values$ZIP[outliers], State = zip$admin_code1[outliers],
  Latitude = zip$latitude[outliers], Longitude = zip$longitude[outliers],
  "Published mean (m)" = zip[[elevation]][outliers],
  "Reproduced mean (m)" = values[[elevation]][outliers],
  "Delta (m)" = delta[outliers],
  check.names = FALSE
)
association <- function(predictor, label) {
  rows <- which(paired & is.finite(predictor))
  data.frame(Predictor = label, n = length(rows),
             "Pearson r" = pearson(abs(delta[rows]), predictor[rows]),
             check.names = FALSE)
}
associations <- rbind(
  association(zip$latitude, "ZIP latitude"),
  association(zip[[relief]], "Published elevation SD")
)

plot_elevation <- function() {
  png(file.path(out_dir, "gate1_elevation.png"), width = 1650, height = 800, res = 150)
  on.exit(dev.off())
  par(mfrow = c(1, 2))
  if (!any(paired)) {
    plot.new()
    text(0.5, 0.5, "No finite elevation pairs")
    return(invisible(NULL))
  }
  limits <- range(c(zip[[elevation]][paired], values[[elevation]][paired]))
  plot(zip[[elevation]][paired], values[[elevation]][paired],
       xlab = "Published elevation mean (m)", ylab = "Reproduced elevation mean (m)",
       main = "Mean elevation by ZIP", xlim = limits, ylim = limits, asp = 1,
       pch = 16, cex = 0.45, col = adjustcolor("#194D71", alpha.f = 0.25))
  abline(0, 1, col = "#B74D32", lwd = 2)
  hist(delta[paired], breaks = 60L, xlab = "Reproduced minus published (m)",
       main = "Mean elevation differences", col = "#689AB4", border = "white")
  abline(v = 0, col = "#B74D32", lwd = 2)
}
plot_elevation()

markdown_table <- function(data) {
  columns <- lapply(data, function(x) {
    cells <- if (is.numeric(x)) format(x, digits = 10L, trim = TRUE) else as.character(x)
    cells[is.na(x)] <- "—"
    gsub("|", "&#124;", cells, fixed = TRUE)
  })
  row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
  c(row(gsub("|", "&#124;", names(data), fixed = TRUE)), row(rep("---", ncol(data))),
    vapply(seq_len(nrow(data)), function(i) row(vapply(columns, function(x) x[i], "")), ""))
}

covered_n <- nrow(values)
eligible_n <- sum(eligible)
uncovered_n <- eligible_n - covered_n
coverage_summary <- data.frame(
  Population = c("Eligible ZIPs in the 50 states and DC", "Fully covered ZIPs",
                 "ZIPs excluded for incomplete tile coverage"),
  Observed = c(eligible_n, covered_n, uncovered_n),
  Expected = c(40970L, 40389L, 581L)
)
coverage_note <- if (all(coverage_summary$Observed == coverage_summary$Expected)) {
  "All three counts match the offline coverage estimate."
} else {
  paste("The table shows the differences from the offline estimate.",
        "Eligibility uses finite centroid coordinates in the 50 states and DC.",
        "Coverage is recomputed from the tile inventory below using all 140 point tile IDs.")
}
exclusions <- data.frame(
  Reason = c("Overseas military mail (APO, FPO, DPO or military state code)",
             "Other ZIPs outside the 50 states and DC", "Non-finite centroid coordinates in the 50 states or DC",
             "Incomplete tile coverage (including coastal and border ZIPs)"),
  ZIPs = c(sum(military), sum(!in_states & !military),
           sum(in_states & !finite_coords), uncovered_n)
)
checksums <- tools::md5sum(c(us_csv, tile_paths))
stopifnot(!anyNA(checksums))
input_files <- data.frame(File = basename(names(checksums)), MD5 = unname(checksums))
elapsed <- proc.time()[["elapsed"]] - started

report <- c(
  "# Gate 1 — AG100 georeferencing at ZIP centroids",
  "",
  sprintf("**Overall Gate 1 verdict: %s.** Both populations must have Pearson r > 0.999 and median absolute mean elevation difference < 2 m.", verdict),
  "Each population must have finite elevation pairs for at least 99% of ZIPs after documented exceptions. Any other missing or non-finite elevation pair fails the gate.",
  "",
  markdown_table(gate),
  "",
  sprintf("Documented exceptions: %d Hawaii ZIPs (%d in the sample), elevation missing on both sides.",
          documented_counts[2L], documented_counts[1L]),
  "",
  if (nrow(unexpected_pairs)) {
    c("### ZIPs with unexpected missing or non-finite elevation pairs", "",
      markdown_table(unexpected_pairs), "")
  },
  "## Inputs and coverage",
  "",
  sprintf("Reference file %s (MD5 %s) contains %d ZIP records.",
          input_files$File[1L], input_files$MD5[1L], nrow(reference)),
  sprintf("The local AG100 v003 inventory contains %d GeoTIFF files.", length(tile_paths)),
  "",
  markdown_table(coverage_summary),
  "",
  coverage_note,
  "",
  markdown_table(exclusions),
  "",
  "Sampling used four rings with 35 points each, V1 nearest-cell extraction and SPSS aggregation.",
  "Each ZIP used its own latitude to calculate longitude offsets.",
  sprintf("The sample was drawn without replacement from the covered ZIPs in input order with seed %d (Mersenne-Twister, Rejection).", seed),
  "",
  "<details>",
  "<summary>AG100 tile filenames and MD5 checksums</summary>",
  "",
  markdown_table(input_files[-1L, ]),
  "",
  "</details>",
  "",
  "## Agreement with published ASTER values",
  "",
  sprintf("The V1 dictionary yields %d ASTER variables across the published bands, including the mean and SD bands for both NDVI and temperature.", length(variables)),
  "Each row uses finite pairs. Missing percentages are compared directly and SPSS missing values remain NA.",
  "A dash indicates an undefined metric or an inapplicable exact-equality check. The last column counts strict equality for raster extrema, cin percentages, the elevation sum and the observation flag.",
  "",
  "### 500-ZIP sample",
  "",
  markdown_table(sample_stats),
  "",
  "### All covered ZIPs",
  "",
  markdown_table(all_stats),
  "",
  "## Mean elevation differences",
  "",
  "![Mean elevation scatter plot and difference histogram](gate1_elevation.png)",
  "",
  "Pearson correlations of absolute differences with latitude and local relief across covered ZIPs.",
  "",
  markdown_table(associations),
  "",
  "The 20 largest absolute differences. Signed differences are reproduced minus published.",
  "",
  markdown_table(top),
  "",
  "### Known differences",
  "",
  "- The original processing used the mean latitude of each batch for longitude offsets. Batch membership and mean latitudes are unavailable.",
  "- The published values used 2021 granules, while this run uses the local 2026 AG100 v003 store. Granule differences can affect pixel values and observation counts.",
  "- Nearest-cell choices at tile and cell seams can change the sampled value. V1 assigns latitude seams to the southern tile and longitude seams to the eastern tile.",
  "",
  sprintf("Elapsed time was %.1f seconds with %d worker(s).", elapsed, workers)
)
writeLines(report, file.path(out_dir, "report.md"), useBytes = TRUE)
print(gate, row.names = FALSE)
message(sprintf("Overall Gate 1 %s. %d of %d eligible ZIPs covered in %.1f s.",
                verdict, covered_n, eligible_n, elapsed))
if (!all(passes)) quit(status = 1L)
