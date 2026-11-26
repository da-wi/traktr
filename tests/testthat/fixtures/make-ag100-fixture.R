# Run from the package root to regenerate the offline AG100 fixtures.
write_ag100_fixture <- function(path, lat_n, lon_w, bad = FALSE, vectors = FALSE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path)) unlink(path)
  access <- rhdf5::H5Pcreate("H5P_FILE_ACCESS")
  on.exit(rhdf5::H5Pclose(access), add = TRUE)
  rhdf5::H5Pset_libver_bounds(access, "H5F_LIBVER_V18", "H5F_LIBVER_V18")
  rhdf5::H5Fclose(rhdf5::H5Fcreate(path, fapl = access))

  n <- 20L
  rows <- row(matrix(0, n, n))
  cols <- col(rows)
  lat <- lat_n - (seq_len(n) - 0.5) / n
  lon <- lon_w + (seq_len(n) - 0.5) / n
  if (bad) lat <- rev(lat)

  bands <- list(
    "ASTER GDEM/ASTGDEM" = 1500L + 10L * rows + cols,
    "NDVI/Mean" = 300L + rows + cols,
    "NDVI/SDev" = 20L + cols,
    "Temperature/Mean" = 33000L + rows,
    "Temperature/SDev" = 10L * rows,
    "Land Water Map/LWmap" = cols %% 2L,
    "Observations/NumObs" = 2L * rows + cols
  )
  bands <- lapply(bands, function(band) {
    band[5:6, 7:8] <- -9999L
    band
  })
  geolocation <- list(
    "Geolocation/Latitude" = if (vectors) lat else matrix(lat, n, n),
    "Geolocation/Longitude" = if (vectors) lon else matrix(lon, n, n, byrow = TRUE)
  )
  datasets <- c(bands, geolocation)
  for (group in unique(dirname(names(datasets)))) rhdf5::h5createGroup(path, group)
  for (name in names(datasets)) {
    data <- datasets[[name]]
    # HDF5's row-major order is the reverse of R's matrix storage order.
    if (is.matrix(data)) data <- t(data)
    rhdf5::h5write(data, path, name)
  }
}

root <- file.path("tests", "testthat", "fixtures", "ag100")
normal <- "AG100.v003.40.-105.0001.h5"
write_ag100_fixture(file.path(root, normal), 40, -105)
write_ag100_fixture(file.path(root, "bad", normal), 40, -105, bad = TRUE)
write_ag100_fixture(file.path(root, "AG100.v003.05.-075.0001.h5"), 5, -75,
                    vectors = TRUE)
