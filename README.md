# traktr <a href='https://da-wi.github.io/traktr/'><img src="man/figures/logo.png" align="right" height="139" /></a>

[![R-CMD-check](https://github.com/da-wi/traktr/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/da-wi/traktr/actions/workflows/R-CMD-check.yml)

## Overview

The `traktr` package provides a tidy workflow for enriching GPS tracking
data from Sentinel-2 imagery and related products. It supports both
H3-based and buffer-based areas of interest (AOI), handles tiling,
querying and rate-limited downloads from the Copernicus Data Space
Ecosystem (CDSE), and attaches contextual statistics back to your GPS
data.

## Installation

Simply install from github using the package `pak`.

``` r
# install.packages("pak")

# Install the latest version from github
pak::pak("da-wi/traktr")
```

## Quick start

This snippet can be run offline. It creates the expected variables and
generates buffers around GPS coordinates.

``` r
library(traktr)
library(dplyr)

# Demo points
gps <- tibble::tibble(
  subject = c("Z001","Z001","Z002"),
  timestamp = as.POSIXct(c("2025-06-01 08:10","2025-06-01 12:05","2025-06-02 09:30")),
  lat = c(47.38, 47.39, 47.385),
  lon = c(8.53, 8.54, 8.545)
)

# Define areas of interest
buf <- mutate_aoi_buffer(gps, lat, lon, dist = 100)    # 100 m circles
hex <- mutate_aoi_h3(gps, lat, lon, res = 10)          # H3 grid

# Quick look 
plot(sf::st_geometry(attr(buf, "traktr_aoi"))
```

## System requirements

- R ≥ 4.2
- For fetching Sentinel-2 data: Copernicus Data Space CDSE account.
- Optional: GPS_SMB_ROOT for your CSV data directory.

## Documentation

See the “Getting started” guide.

## Citation

Please cite `traktr` and data sources (Sentinel-2, ESA WorldCover) in
publications.

## Contributing

Issues and PRs welcome.

## License

MIT
