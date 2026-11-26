# traktr 0.0.1

## New features

- `read_fips_csv()` reads and normalizes county FIPS codes from CSV files.
- `aoi_from_sf()` creates polygon areas of interest from sf objects.
- `aoi_us_counties()` creates county areas of interest from US Census boundaries.
- `earthdata_auth_check()` detects available Earthdata Login credentials offline.
- `fetch_aster_ged()` downloads intersecting AG100 tiles and converts them to GeoTIFF.
- `attach_aster_ged()` computes area-weighted statistics for county polygons.
- `write_county_summary()` exports county summaries as CSV or Parquet with a data dictionary and provenance report.
- Known data gaps include missing elevation in Hawaii and ASTER GDEM pits near water bodies. County area in tiles without AG100 granules counts toward `missing_perc`.

## Internal

- Added a tiled extraction engine that combines area-weighted moments across tiles.
- Fixed `.aoi_deduplicate()` to use the active sf geometry column.
- Updated `.join_stats()` to use a temporary row-order column that avoids collisions with existing columns.
