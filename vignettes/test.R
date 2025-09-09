library(tidyverse)
devtools::load_all()
#devtools::install_github("da-wi/cartographr@develop")
library(cartographr)

#gps_data <- read_gps_csv(pattern = paste0("^Z31(2|3).*\\.csv$"))  |> 
#  select(subject, timestamp, lat, lon, direction, speed)
gps_data1 <- read_gps_csv(pattern = paste0("^Z.*\\.csv$"))  |> 
  select(subject, timestamp, lat, lon, direction, speed)
gps_data <- read_gps_csv(pattern = paste0("^Z.*\\.csv$"))  |> 
  select(subject, timestamp, lat, lon, direction, speed)

gps_data <- vroom::vroom("output/enhanced.csv") |> select(subject, timestamp,lat,lon)

aa <- gps_data |> mutate_aoi_h3()

ggplot() +
  ggspatial::annotation_map_tile(type = "cartolight", zoom = 7) +
  geom_sf(data = attr(aa, "traktr_aoi"), fill = NA, colour = "#D81B60", linewidth = 0.8) +
  coord_sf(expand = FALSE)


aab <- aa |> select(-c(timestamp,lat,lon)) |> data.table::as.data.table() |> unique()

attr(aa, "traktr_aoi")

aa <- gps_data #|> filter(subject=="Z313") |> select(subject, timestamp,lat,lon)
bb <- aa |> mutate_aoi_h3(lat = "lat", lon = "lon") |>
  fetch_ndvi(
    from = min(aa$timestamp),
    to   = max(aa$timestamp),
    out_dir = file.path("output/ndvi_10m_individual2", unique(aa$subject)),
    client = oauth,
    res_m = 10
  ) |>
  attach_ndvi() |>
  fetch_landuse(out_dir = "output/worldcover_individual2") |>
  attach_landuse(out_dir = "output/worldcover_individual2")

 oauth <- CDSE::GetOAuthClient(id = Sys.getenv("CDSE_CLIENT_ID", ""),
                              secret = Sys.getenv("CDSE_CLIENT_SECRET", ""))  

ahsg <- gps_data |>
  dplyr::group_split(subject) |>
  # subject specific 
  purrr::map(function(df_subj) {
    cli::cli_alert_info("→ Subject {.strong {first(df_subj$subject)}}")
    df_subj |>
      mutate_aoi_h3(lat = "lat", lon = "lon") |>
      fetch_ndvi(
        from = min(df_subj$timestamp),
        to   = max(df_subj$timestamp),
        out_dir = file.path("output/ndvi_10m_individual", unique(df_subj$subject)),
        client = oauth,
        res_m = 10
      ) |>
      attach_ndvi() |>
      fetch_landuse(out_dir = "output/worldcover_individual") |>
      attach_landuse(out_dir = "output/worldcover_individual")
  }) |>
  dplyr::bind_rows()

ahsg |> vroom::vroom_write("output/enhanced.csv")


sub <- "Z313"
osm <- cartographr::get_osmdata(bbox = ahsg |> 
                                  filter(subject == sub) |> 
                                  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
                                  sf::st_bbox())
set_attribution(F)

{
  p <- osm |> plot_map("swiss") 
  p +
    ggnewscale::new_scale_fill()+
    geom_sf(data= ahsg |> filter(subject == sub) |>
              select(geometry, ndvi_mean, frac_green_core, h3) |> 
              unique(),
            aes(fill = frac_green_core), 
            color = NA, alpha = 0.6)+
    scale_fill_gradient2(
      low = "brown",midpoint = 0.5,
      mid = "#E7D6A4",
      high = "darkgreen",
      name = "values",
      limits = c(0,1)
    )+
    cartographr:::adjust_viewport(p)+
    cartographr::theme_infomap()
}

sub <- "Z313"
ggplot() +
  geom_sf(data=  ahsg |> filter(subject == sub) |> 
            select(geometry, ndvi_mean, frac_green_core, h3) |> 
            pivot_longer(cols = c(ndvi_mean, frac_green_core)) |>
            unique(),
          #aes(fill = frac_green_core), 
          aes(fill = value), 
          color = NA, alpha = 0.6)+
  scale_fill_gradient2(
    low = "brown",midpoint = 0.5,
    mid = "#E7D6A4",
    high = "darkgreen",
    name = "values",
    limits = c(0,1)
  ) + 
  hrbrthemes::theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size =7),
        axis.text.y = element_text( size =7))+
  facet_wrap(~name)

ggplot(ahsg |> st_drop_geometry() |> select(subject,h3, ndvi_mean,frac_green_core) |> unique())+
  geom_point(aes(x = ndvi_mean, y = frac_green_core))+
  facet_wrap(~subject)

MuMIn::r.squaredGLMM(
  lmerTest::lmer(ndvi_mean ~ frac_green_core  + (1|subject), data = ahsg |> sf::st_drop_geometry() |> select(subject,h3, ndvi_mean,frac_green_core) |> unique())
)


# create p, need it afterwards for adjusting viewport
p <- osm |> plot_map("arctic")
p + 
  ggnewscale::new_scale_fill()+
  geom_sf(data= gps_data |>
            data.table::as.data.table() |>
            dplyr::group_by(h3, geometry, ndvi_mean) |>
            count() |> 
            ungroup() |>
            mutate(dwell_s   = n * 10,
                   dwell_prop = dwell_s / max(dwell_s, na.rm = TRUE)) |>
            filter(n < 8000) |>
            sf::st_as_sf(),
          aes(fill = ndvi_mean, alpha = dwell_prop), 
          color = "#565656")+
  scale_fill_gradient2(
    low = "brown",midpoint = 0.5,
    mid = "#E7D6A4",
    high = "darkgreen",
    name = "NDVI",
    limits = c(0,1)
  )+
  cartographr:::adjust_viewport(p)+
  cartographr::theme_infomap()

#######################
#
pak::pak('thomasp85/gganimate')
library(gganimate)

sub <- "Z313"
bb <- ahsg |> filter(subject == sub) |>  
  select(geometry, timestamp, ndvi_mean, frac_green_core, h3) |>
  mutate(time = as.integer((timestamp-min(timestamp))/3600),
         first = min(timestamp)) |> 
  arrange(time) |> 
  data.table::as.data.table() |>
  group_by(time, h3, geometry, first) |> 
  summarize(ndvi_mean = mean(ndvi_mean),
            frac_green_core = mean(frac_green_core)) |>
  mutate(actual_time = lubridate::floor_date(first + time * 3600, "hours")) |>
  sf::st_as_sf() 

p <- ggplot(bb) +
  geom_sf(data = bb, aes(fill = frac_green_core), color = NA, alpha = 0.6)+
  scale_fill_gradient2(
    low = "brown",
    midpoint = 0.5,
    mid = "#E7D6A4",
    high = "darkgreen",
    name = "values",
    limits = c(0,1)
  )+ 
  hrbrthemes::theme_ipsum()+
  labs(caption = 'Time: {as.integer(frame_time)}', fill = "green (%)") +
  transition_time(time) +
  ease_aes("linear")

anim <- animate(p, renderer = ffmpeg_renderer(
  format = "wmv",
  ffmpeg = NULL,
  options = list(vcodec = "wmv2", acodec = "wmav2")), 
  duration = 30, 
  height = 4, 
  width = 6,
  units = "in", res = 150, nframes = bb |> ungroup() |> distinct(time) |> nrow())
anim_save(animation = anim, filename = "animation.wmv", path = ".")



##################################
ahsgdt <- data.table::as.data.table(ahsg)

ahsgdt |> group_by(subject)


MuMIn::r.squaredGLMM(
  lmerTest::lmer(ndvi_mean ~ frac_green_core  + (1|subject), data = ahsgdt |> sf::st_drop_geometry() |> select(subject,h3, ndvi_mean,frac_green_core) |> unique())
)


ahsgdt |> select(subject, ndvi_mean, frac_green_core) |> distinct() |>
  ggplot()+
  geom_point(aes(x=ndvi_mean, y= frac_green_core))+
  facet_wrap(~subject,scale = "free")+
  hrbrthemes::theme_ipsum()


ahsgdt |> select(subject, ndvi_mean, frac_green_core) |> 
  distinct() |> filter(ndvi_mean> 0) |>
  ggplot()+
  geom_density2d_filled(aes(x=ndvi_mean, y= frac_green_core))+ 
  hrbrthemes::theme_ipsum()



bb <- aa |> mutate_aoi_h3(lat = "lat", lon = "lon") |>
  fetch_ndvi(
    from = min(aa$timestamp),
    to   = max(aa$timestamp),
    out_dir = file.path("output/ndvi_10m_individual2", unique(aa$subject)),
    client = oauth,
    res_m = 10
  ) |>
  attach_ndvi() |>
  attach_landuse(out_dir = "output/worldcover_individual2")  |> fetch_landuse(out_dir = "output/worldcover_individual2") 



gps_data_Z313 <- gps_data_Z313 |>
  mutate_aoi_h3(lat = "lat", lon = "lon") |>
  fetch_ndvi(
    from = min(gps_data_Z313$timestamp),
    to   = max(gps_data_Z313$timestamp),
    out_dir = "output/.ndvi",
    client = oauth,
    res_m = 10
  ) |>
  attach_ndvi() |>
  fetch_landuse(out_dir = "output/.landuse") |>
  attach_landuse()


library(traktr)
library(dplyr)

# Demo points
gps <- tibble::tibble(
  subject = c("Z001","Z001","Z002"),
  timestamp = as.POSIXct(c("2025-06-01 08:10","2025-06-01 12:05","2025-06-02 09:30")),
  lat = c(47.382, 47.384, 47.385),
  lon = c(8.545, 8.542, 8.545)
)



# Define areas of interest
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

ggsave(
  "man/figures/README-plot01.png",
  width = 5, height = 3.5, dpi = 150
)


gps_data_Z313 <- gps_data |> filter(subject == "Z313")
gps_data_Z313_rich <- gps_data |> filter(subject == "Z313") |>
  select(lat, lon, subject, timestamp) |>
  mutate_aoi_h3(lat = "lat", lon = "lon") |>
  fetch_ndvi(
    from = min(gps_data_Z313$timestamp),
    to   = max(gps_data_Z313$timestamp),
    out_dir = "output/.ndvi",
    client = oauth,
    res_m = 10
  ) |>
  attach_ndvi() |>
  fetch_landuse(out_dir = "output/.landuse") |>
  attach_landuse()
