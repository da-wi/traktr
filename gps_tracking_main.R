# 1) Install pak if needed
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# 2) Declare package *specs* (prefix with cran:: or github::)
specs <- c(
  "cran::rnaturalearth",
  "cran::sf",
  "cran::tidyterra",
  "cran::ggnewscale",
  "cran::cmocean",
  "cran::tidyverse",
  "cran::data.table",
  "github::crazycapivara/h3-r",
  "cran::cartographr",
  "cran::CDSE"
)

pak::pkg_install(specs, upgrade = FALSE, ask = FALSE)

# 3) Load what you actually use
to_load <- c("h3","rnaturalearth","CDSE","cartographr",
             "sf","tidyterra","ggnewscale","cmocean",
             "tidyverse","data.table")
invisible(lapply(to_load, library, character.only = TRUE))


# get cartography data
zurich_osm <- get_osmdata(lat = 47.376188350502574,
                          lon = 8.538240122970143,
                          x_distance = 1200,
                          y_distance = 1200)

sf::st_bbox(xmin=8.522321, ymin=47.365408, xmax=8.554159, ymax= 47.386968)
########################################################
# Copernicus OAuth client data
id <- "sh-09a6a314-1c14-4612-be75-6ea7354309eb"
secret <- "wKZwZZUtWzYCfwhgYYxz0aTBLeEvZPdL"
oauth <- GetOAuthClient(id = id, secret = secret)

aoi <- zurich_osm$bbox |> st_as_sfc() |> st_transform(4326)

images <- SearchCatalog(aoi = aoi,
                        from = "2022-01-01",
                        to = "2025-07-07", 
                        collection = "sentinel-2-l2a",
                        with_geometry = TRUE, 
                        client = OAuthClient) 

# Pick the least-cloudy day
day <- images |> 
  arrange(tileCloudCover) |> 
  head(1) |> 
  pull(acquisitionDate)


#script_file <- system.file("scripts", "TrueColor.js", package = "CDSE")
# day <- "2025-04-04"
script_file <- system.file("scripts", "NDVI_float32.js", package = "CDSE")

fp <- GetImage(
  aoi = aoi, time_range = day, script = script_file,
  collection = "sentinel-2-l2a", format = "image/tiff",
  mosaicking_order = "leastCC", resolution = 10, mask = TRUE,
  buffer = 0, client = oauth
)

# fp[fp < 0.5] <- NA            # single-band; for multi-band see below
plot(fp)


## h3 tiles 
zurich_h3 <- st_as_sfc(zurich_osm$bbox) %>%
  st_as_sf() |># sf::st_buffer(dist = 2000) %>%
  h3::polyfill(res = 10) 
##

# 1) H3 polygons as sf (one row per cell)
h3_sf <- h3::h3_to_geo_boundary_sf(zurich_h3) |>   # returns polygons
  st_as_sf() |>
  st_transform(terra::crs(ras))                    # match raster CRS

# 2) Exact area-weighted mean NDVI per polygon
h3_sf$ndvi_mean <- exactextractr::exact_extract(
  fp, h3_sf, function(v, w) weighted.mean(v, w, na.rm = TRUE)
)

# (optional) fraction of polygon covered by valid NDVI
h3_sf$valid_frac <- exactextractr::exact_extract(
  fp, h3_sf, function(v, w) sum(w[!is.na(v)]) / sum(w)
)

ggplot(h3_sf)+
  geom_sf(alpha = 1, aes(geometry= geometry), show.legend = F)+
  cartographr:::adjust_viewport(zurich_osm)

set_attribution(attribution = F)
set_output_size(c(210,210))

pp <- plot_map(zurich_osm, palette = "grey") +
  ggnewscale::new_scale_fill()+
  cmocean::scale_fill_cmocean(name = "speed", na.value = "transparent")+
  geom_spatraster(data=ras, alpha = 0.5, show.legend = F)+
  cartographr:::adjust_viewport(zurich_osm)+
  theme_infomap_barlow()

pp_sf <- plot_map(zurich_osm, palette = "grey") +
  ggnewscale::new_scale_fill()+
  cmocean::scale_fill_cmocean(name = "speed", na.value = "transparent")+
  geom_spatraster(data=ras, alpha = 0.5)+
  geom_sf(data=h3_sf,alpha = 0, aes(geometry= geometry), show.legend = F)+
  cartographr:::adjust_viewport(zurich_osm)+
  theme_infomap_barlow()

pp_sf_filled <- plot_map(zurich_osm, palette = "grey") +
  ggnewscale::new_scale_fill()+
  cmocean::scale_fill_cmocean(name = "speed", na.value = "transparent")+
  geom_sf(data=h3_sf,alpha = 0.5, color=NA,aes(geometry= geometry, fill = ndvi_mean), show.legend = F)+
  #geom_spatraster(data=ras, alpha = 0.5)+
  cartographr:::adjust_viewport(zurich_osm)+
  theme_infomap_barlow()+
  labs(title = "Vegetation in Zurich",
       fill = "NDVI")
pp_sf_filled
#save_map(pp, "pp.png", device = "png")
#save_map(pp_sf, "pp_sf.png", device = "png")
#save_map(pp_sf_filled, "pp_sf_filled.png", device = "png")

st_as_sfc(zurich_osm$bbox) |> #ne_states(country = "Switzerland", returnclass = "sf") |>
  #select(iso_a2, geometry) |>
  #filter(name_en == "Zürich") |>
  st_as_sf() |>
  h3::polyfill(res = 9)


images <- SearchCatalog(aoi = st_as_sfc(zurich_osm$bbox), from = "2023-06-01", to = "2023-06-30", 
                        collection = "sentinel-5p-l2", with_geometry = TRUE, client = OAuthClient)

day <- images %>% arrange(tileCloudCover) %>% head(1) %>% pull(acquisitionDate)
script_file <- "Z:/Reisen/Zurich25/S5P_NO2.js"
ras <- GetImage(aoi = st_as_sfc(zurich_osm$bbox),
                time_range = day,
                script = script_file, 
                collection = "sentinel-5p-l2", #"sentinel-2-l2a",
                format = "image/tiff",
                resolution = 10, mask = F, buffer = 100, client = OAuthClient)

#ras[ras < 0.5] <- NA
plot(ras)


save_map(plot = pp, filename = "zurich.pdf")
  

# Get country boundary for Switzerland
switzerland_h3 <- st_as_sfc(zurich_osm$bbox) %>% #ne_states(country = "Switzerland", returnclass = "sf") |>
  #select(iso_a2, geometry) |>
  #filter(name_en == "Zürich") |>
  st_as_sf() |>
  h3::polyfill(res = 9)

pm25 <- terra::crop(a, st_as_sfc(zurich_osm$bbox))

# Now extract values. terra will handle reprojection if needed.
#points_sv <- terra::vect(vienna_h3_latlon, geom=c("lon", "lat"), crs="EPSG:4326")
#extracted_values <- terra::extract(a, points_sv)

a <-  terra::rast("C:/Users/David.Willinger/ASTER/data/tiles/sdei-global-annual-gwr-pm2-5-modis-misr-seawifs-aod-v4-gl-03-2019.tif")
library(raster)
a_zurich_crop <- terra::crop(a, zurich_osm$bbox)
print(a_zurich_crop)

switzerland_h3_latlon <- data.table(h3_index = switzerland_h3) |>
  mutate(h3 = lapply(h3_index, coord_h3)) |>
  unnest(h3) |>
  mutate(GWRPM25 = a[terra::cellFromXY(a, xy = data.frame(x=lon, y=lat))]$`sdei-global-annual-gwr-pm2-5-modis-misr-seawifs-aod-v4-gl-03-2019`) |>
  group_by(h3_index) |>
  reframe(GWRPM25 = mean(GWRPM25, na.rm = T), lat = lat, lon = lon) |>
  mutate(geometry = h3_to_geo_boundary_sf(h3_index)$geometry)


ggplot(switzerland_h3_latlon) +
  geom_sf(aes(fill = GWRPM25, geometry = geometry), color = NA) +
  scale_fill_distiller(palette = "RdYlBu")+
  theme_void()

pp2 <- plot_map(zurich_osm, palette = "grey") +
  ggnewscale::new_scale_fill()+
  #cmocean::scale_fill_cmocean(name = "speed", na.value = "transparent")+
  #geom_spatraster(data=ras, alpha = 0.5)+
  scale_fill_distiller(palette = "RdYlBu")+
  geom_sf(data = switzerland_h3_latlon, 
          aes(fill = GWRPM25, geometry = geometry), color = NA, alpha = 0.2)+
  #geom_spatraster_rgb(data=ras)+
  cartographr:::adjust_viewport(zurich_osm)+
  theme_infomap_barlow()+
  labs(title = "Vegetation in Zurich",
       fill = "NDVI")
  
pp3 <- #plot_map(zurich_osm, palette = "imhof") +
  ggplot()+
  ggnewscale::new_scale_fill()+
  cmocean::scale_fill_cmocean(name = "speed", na.value = "transparent")+
  geom_spatraster(data=pm25, alpha = 0.5)+
  #geom_spatraster_rgb(data=ras)+
  cartographr:::adjust_viewport(zurich_osm)+
  theme_infomap_barlow()+
  labs(title = "PM in Zurich",
       fill = "PM2.5")
