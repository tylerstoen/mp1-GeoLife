library(tidyverse)
library(lubridate)
library(fs)
library(sf)
library(ggspatial)
library(markovchain)
library(prettymapr)

# helper functions live here (defines sample_geolife_same_locale_fast, etc.)
source("scripts/sample_geolife.R")

geolife_root <- "data/Geolife Trajectories 1.3/Data"

# ------------------------------------------------------------
# 0) User-tunable knobs
# ------------------------------------------------------------
n_traj   <- 100   # how many trajectories to pull
n_plot   <- 10    # how many to show on the map
cell_deg <- 0.02  # grid resolution (~2.2km in latitude); adjust as desired

# ------------------------------------------------------------
# 1) Sample trajectories (single locale, fast)
# ------------------------------------------------------------
pts_local <- sample_geolife_same_locale_fast(
  n = n_traj,
  radius_km = 60,
  n_candidates = 1500,
  head_pts = 50,
  seed = 42,
  root = geolife_root
)

# ------------------------------------------------------------
# Cluster Grid
# ------------------------------------------------------------

unique_pts <- pts_local |>
  count(lon, lat, name = "w")

pts_sf <- st_as_sf(unique_pts,
                   coords = c("lon","lat"),
                   crs = 4326) |>
  st_transform(3857)

coords <- st_coordinates(pts_sf)

km <- kmeans(coords, centers = 10, iter.max = 100)


unique_pts$clust <- km$cluster

pts_local$state <- km$cluster[
  match(
    paste(pts_local$lon, pts_local$lat),
    paste(unique_pts$lon, unique_pts$lat)
  )
]

# --------------------------------------------------
# Plot clusters
# --------------------------------------------------

centers <- as.data.frame(km$centers)
colnames(centers) <- c("x", "y")

centers_sf <- st_as_sf(
  centers,
  coords = c("x", "y"),
  crs = 3857
) |>
  mutate(clust = seq_len(nrow(centers)))

# union centroids to build shapes on map
bbox <- st_as_sfc(st_bbox(centers_sf)) # bounding box for clipping

# Voronoi polygons
voronoi <- st_voronoi(st_union(centers_sf))
voronoi_polys <- st_collection_extract(voronoi) # extract polygons

# turn into sf object
voronoi_sf <- st_sf(
  geometry = voronoi_polys
) %>%
  st_intersection(bbox) %>% # clip to bounding box
  mutate(clust = seq_len(nrow(.))) # assign cluster IDs


ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_rect(
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "white", alpha = 0.35
  ) +
  geom_sf(
    data = voronoi_sf,
    aes(fill = factor(clust)),
    alpha = 0.4,
    color = "black",
    linewidth = 0.3
  ) +
  geom_sf(
    data = centers_sf,
    color = "red",
    size = 2
  ) +
  coord_sf(crs = st_crs(3857), expand = FALSE) +
  theme_minimal() +
  labs(title = "K-means Clusters (Voronoi Polygons)", fill = "Cluster") +
  guides(fill = "none")

# -------------------------------------------------------
# Define States
# -------------------------------------------------------

pts_states <- pts_local |>
  group_by(traj_id) |>
  filter(state != lag(state) | is.na(lag(state))) |>
  ungroup()




