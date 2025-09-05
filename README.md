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

Install with `pak`.

``` r
# install.packages("pak")
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

hex <- gps |> 
  mutate_aoi_h3(res = 9) |> 
  fetch_ndvi(from = as.Date("2025-04-01"),
           to = as.Date("2025-07-31"),
           out_dir = "data/.ndvi",
           client = oauth, res_m = 10) |>
  attach_ndvi()


ggplot(hex, aes(fill = ndvi_mean))+
  geom_sf() + 
  scale_fill_gradient2(
    low = "brown",
    mid = "#E7D6A4",
    high = "darkgreen",
    midpoint = 0.5,
    limits = c(0,1)
  )+
  labs(fill = "NDVI")
```

![](man/figures/README-plot01.png)

## System requirements

- R ≥ 4.2
- For fetching Sentinel-2 data: [CDSE
  account](https://dataspace.copernicus.eu/).
- Optional: GPS_SMB_ROOT for your CSV data directory.

## Documentation

See the “Getting started” guide.

## Citation

Please cite `traktr` and data sources (Sentinel-2, ESA WorldCover) in
publications.

``` r
citation("traktr")
```

## Contributing

Contributions are welcome! Please open an issue for bug reports or
feature requests.

## License

MIT
