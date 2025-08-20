
library(sf)
library(tidyverse)
library(lubridate)

set.seed(42)

# ---------------- Parameter ----------------
n_people   <- 20          # Anzahl Personen
days       <- 14          # Tage
step_mins  <- 10          # Intervall in Minuten
tz_local   <- "Europe/Zurich"

# Grobe BBOX für Zürich (WGS84). Für Mock-Daten reicht das.
zrh <- st_as_sfc(st_bbox(c(xmin = 8.45, ymin = 47.32,
                           xmax = 8.64, ymax = 47.45), crs = 4326))

crs_wgs84 <- 4326
crs_ch    <- 2056  # LV95 für metrische Distanzen

# ------------- Hilfsfunktionen (metrisch) -------------
to_ch  <- function(g) st_transform(g, crs_ch)
to_wgs <- function(g) st_transform(g, crs_wgs84)

sample_point <- function(poly) st_sample(poly, 1) %>% st_cast("POINT")

move_from <- function(pt, dist_m, angle_rad) {
  p  <- to_ch(pt)
  xy <- st_coordinates(p)
  new_xy <- cbind(xy[1] + dist_m * cos(angle_rad),
                  xy[2] + dist_m * sin(angle_rad))
  to_wgs(st_sfc(st_point(new_xy), crs = crs_ch))
}

# ------------- Orte pro Person (Home/Work/Parks) -------------
mk_places <- function(id, poly) {
  home      <- sample_point(poly)
  work      <- move_from(home,      runif(1, 2000, 7000), runif(1, 0, 2*pi))
  park_home <- move_from(home,      runif(1,  300, 1200), runif(1, 0, 2*pi))
  park_work <- move_from(work,      runif(1,  200,  800), runif(1, 0, 2*pi))
  tibble(id = id,
         place = c("home","work","park_home","park_work"),
         geometry = c(home, work, park_home, park_work)) |>
    st_as_sf(crs = 4326)
}

ids     <- sprintf("P%02d", 1:n_people)
places  <- map_dfr(ids, mk_places, poly = zrh)

# Für schnelle Berechnung: Koordinaten als metrische x/y pro Ort breit pivotieren
places_xy <- places |>
  st_transform(crs_ch) |>
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) |>
  st_drop_geometry() |>
  select(id, place, x, y) |>
  pivot_wider(names_from = place, values_from = c(x,y), names_sep = "_")

# ------------- Zeitraster (2 Wochen, 10-minütig) -------------
start_datetime <- as.POSIXct(paste(Sys.Date() - (days - 1), "00:00:00"), tz = tz_local)
end_datetime   <- start_datetime + days*24*3600 - step_mins*60
times          <- seq(start_datetime, end_datetime, by = paste(step_mins, "mins"))

grid <- tidyr::crossing(id = ids, time = times) |>
  mutate(date       = as.Date(time, tz = tz_local),
         tod_min    = hour(time) * 60 + minute(time),
         wday       = lubridate::wday(time, week_start = 1),  # 1=Mo
         is_weekend = wday >= 6)

# Tages-spezifische (und personen-spezifische) Zeitpläne mit etwas Zufall
sched_params <- tidyr::crossing(id = ids, date = unique(grid$date)) |>
  mutate(
    commute_start1 = 7*60 + 30 + rnorm(n(), 0, 10),
    commute_dur1   = pmax(40, rnorm(n(), 60, 8)),
    work_start     = commute_start1 + commute_dur1,
    lunch_start    = 12*60 + rnorm(n(), 0, 15),
    lunch_dur      = pmax(20, rnorm(n(), 45, 10)),
    commute_start2 = 17*60 + rnorm(n(), 0, 10),
    commute_dur2   = pmax(40, rnorm(n(), 60, 8))
  )

# ------------- State-Zuordnung pro Zeitstempel -------------
traj <- grid |>
  left_join(sched_params, by = c("id","date")) |>
  mutate(
    state = case_when(
      !is_weekend & tod_min <  commute_start1                          ~ "home",
      !is_weekend & tod_min >= commute_start1 &
        tod_min <  commute_start1 + commute_dur1                       ~ "commute_to_work",
      !is_weekend & tod_min >= work_start &
        tod_min <  lunch_start                                         ~ "work",
      !is_weekend & tod_min >= lunch_start &
        tod_min <  lunch_start + lunch_dur                             ~ "park_lunch",
      !is_weekend & tod_min >= lunch_start + lunch_dur &
        tod_min <  commute_start2                                      ~ "work",
      !is_weekend & tod_min >= commute_start2 &
        tod_min <  commute_start2 + commute_dur2                       ~ "commute_to_home",
      !is_weekend & tod_min >= commute_start2 + commute_dur2           ~ "home",
      is_weekend  & tod_min >= 10*60 & tod_min < 16*60                 ~ "park_weekend",
      TRUE                                                             ~ "home"
    ),
    frac = case_when(
      state == "commute_to_work" ~ (tod_min - commute_start1) / commute_dur1,
      state == "commute_to_home" ~ (tod_min - commute_start2) / commute_dur2,
      TRUE ~ NA_real_
    )
  )

# ------------- Positionen je State (metrisch) -------------
# leichte Positionsjitter, damit Punkte an Home/Work/Park nicht exakt gleich sind
sd_home  <- 20    # m
sd_work  <- 20
sd_park  <- 15

traj_xy <- traj |>
  left_join(places_xy, by = "id") |>
  mutate(
    x = case_when(
      state == "home"          ~ rnorm(n(), x_home, sd_home),
      state == "work"          ~ rnorm(n(), x_work, sd_work),
      state == "park_lunch"    ~ rnorm(n(), x_park_work, sd_park),
      state == "park_weekend"  ~ rnorm(n(), x_park_home, sd_park),
      state == "commute_to_work" ~ x_home + frac * (x_work - x_home),
      state == "commute_to_home" ~ x_work + frac * (x_home - x_work),
      TRUE ~ rnorm(n(), x_home, sd_home)
    ),
    y = case_when(
      state == "home"          ~ rnorm(n(), y_home, sd_home),
      state == "work"          ~ rnorm(n(), y_work, sd_work),
      state == "park_lunch"    ~ rnorm(n(), y_park_work, sd_park),
      state == "park_weekend"  ~ rnorm(n(), y_park_home, sd_park),
      state == "commute_to_work" ~ y_home + frac * (y_work - y_home),
      state == "commute_to_home" ~ y_work + frac * (y_home - y_work),
      TRUE ~ rnorm(n(), y_home, sd_home)
    )
  )

# ------------- sf-Objekt + Längen/Geschwindigkeit -------------
traj_sf <- traj_xy |>
  st_as_sf(coords = c("x","y"), crs = crs_ch, remove = FALSE) |>
  st_transform(crs_wgs84) |>
  mutate(
    lon = st_coordinates(geometry)[,1],
    lat = st_coordinates(geometry)[,2]
  ) |>
  arrange(id, time) |>
  group_by(id) |>
  mutate(
    # Distanz & Geschwindigkeit (m, m/s) über metrische x/y (LV95)
    dist_m   = sqrt((x - lag(x))^2 + (y - lag(y))^2),
    speed_ms = dist_m / (step_mins * 60)
  ) |>
  ungroup() |>
  select(id, time, date, state, lon, lat, dist_m, speed_ms, geometry)

# ----- Beispiel-Plot (optional) -----
library(dplyr)
library(ggplot2)

id_show <- "P01"

df <- traj_sf %>%
  filter(id == id_show) %>%
  arrange(time)

ggplot(df, aes(lon, lat)) +
  geom_path(aes(color = state), linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(color = state), size = 0.6, alpha = 0.6) +
  coord_equal() +
  labs(title = paste("Trajectory of", id_show),
       x = "Longitude", y = "Latitude", color = "State") +
  theme_minimal()


traj_sf
