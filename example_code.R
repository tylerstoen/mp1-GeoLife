# ============================================================
# GeoLife -> discretized states -> visit heatmap -> lazy RW fit
# (NO thinning; uses path completion to enforce 4-neighbor steps)
# ============================================================

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
# 2) Plot a few example trajectories on an OSM basemap (raw points)
# ------------------------------------------------------------
set.seed(12)

traj_ids <- pts_local |>
  distinct(traj_id)

n_keep <- min(n_plot, nrow(traj_ids))

traj_keep <- traj_ids |>
  slice_sample(n = n_keep) |>
  pull(traj_id)

pts_plot <- pts_local |>
  filter(traj_id %in% traj_keep)

traj_lines <- st_as_sf(pts_plot, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  arrange(traj_id, t) |>
  group_by(traj_id) |>
  summarise(do_union = FALSE, .groups = "drop") |>
  st_cast("LINESTRING") |>
  st_transform(3857)

ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  # lighten the basemap
  geom_rect(
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "white", alpha = 0.35
  ) +
  geom_sf(
    data = traj_lines,
    aes(color = factor(traj_id)),
    linewidth = 1.0,
    alpha = 0.9
  ) +
  coord_sf(crs = st_crs(3857), expand = FALSE) +
  theme_minimal() +
  labs(title = paste0("GeoLife trajectories (", n_keep, " sampled)")) +
  guides(color = "none")

# ------------------------------------------------------------
# 3) Discretize to grid cells
# ------------------------------------------------------------
pts_states <- pts_local |>
  arrange(traj_id, t) |>
  mutate(
    gx = floor(lon / cell_deg),
    gy = floor(lat / cell_deg),
    state = paste(gx, gy, sep = "_")
  )

# Drop consecutive duplicates in cell space (removes repeated GPS points in same cell)
pts_states <- pts_states |>
  group_by(traj_id) |>
  filter(gx != lag(gx) | gy != lag(gy) | is.na(lag(gx))) |>
  ungroup()

# ------------------------------------------------------------
# 4) Path completion (enforce 4-neighbor moves on Z^2)
#    For each observed jump (gx,gy)->(gx',gy'), insert intermediate cells
#    along a Manhattan path (x-then-y). This guarantees |dx|+|dy| in {0,1}.
# ------------------------------------------------------------

# Complete a single trajectory given integer vectors gx, gy (same length)
complete_path_xy <- function(gx, gy) {
  stopifnot(length(gx) == length(gy))
  n <- length(gx)
  if (n == 0) return(list(gx = integer(), gy = integer()))
  if (n == 1) return(list(gx = gx, gy = gy))
  
  out_x <- integer(0)
  out_y <- integer(0)
  
  # start at the first observed cell
  x <- gx[1]
  y <- gy[1]
  out_x <- c(out_x, x)
  out_y <- c(out_y, y)
  
  for (i in 2:n) {
    x1 <- gx[i]
    y1 <- gy[i]
    
    dx <- x1 - x
    dy <- y1 - y
    
    # move in x direction one unit at a time
    if (dx != 0) {
      sx <- sign(dx)
      for (k in seq_len(abs(dx))) {
        x <- x + sx
        out_x <- c(out_x, x)
        out_y <- c(out_y, y)
      }
    }
    
    # then move in y direction one unit at a time
    if (dy != 0) {
      sy <- sign(dy)
      for (k in seq_len(abs(dy))) {
        y <- y + sy
        out_x <- c(out_x, x)
        out_y <- c(out_y, y)
      }
    }
  }
  
  list(gx = out_x, gy = out_y)
}

# Build completed state sequences (no stitching across trajectories)
seq_tbl <- pts_states |>
  group_by(user, traj_id, traj_file) |>
  summarise(
    gx = list(gx),
    gy = list(gy),
    .groups = "drop"
  ) |>
  mutate(
    # complete the path per trajectory
    completed = purrr::map2(gx, gy, complete_path_xy),
    gx_c = purrr::map(completed, "gx"),
    gy_c = purrr::map(completed, "gy"),
    state_seq = purrr::map2(gx_c, gy_c, ~ paste(.x, .y, sep = "_"))
  ) |>
  select(user, traj_id, traj_file, state_seq, gx_c, gy_c)

seqs <- seq_tbl$state_seq

# estimated drift (mean increment) for the discretized GeoLife random walk
# assumes you already have pts_states with gx, gy and completed paths / unit steps

# --- option 1: drift from completed unit steps (best if you used path completion) ---
drift_hat <- steps |>
  summarise(
    mu_dx = mean(dx),
    mu_dy = mean(dy)
  )

drift_hat

# ------------------------------------------------------------
# 5) Visit-count heatmap (based on COMPLETED state sequences)
# ------------------------------------------------------------
visit_counts <- seq_tbl |>
  tidyr::unnest(state_seq) |>
  count(state_seq, name = "visits") |>
  rename(state = state_seq)

heat_sf <- visit_counts |>
  tidyr::separate(state, into = c("gx", "gy"), sep = "_", convert = TRUE) |>
  mutate(
    xmin = gx * cell_deg,
    xmax = (gx + 1) * cell_deg,
    ymin = gy * cell_deg,
    ymax = (gy + 1) * cell_deg,
    geometry = purrr::pmap(
      list(xmin, xmax, ymin, ymax),
      function(xmin, xmax, ymin, ymax) {
        st_polygon(list(matrix(
          c(xmin, ymin,
            xmax, ymin,
            xmax, ymax,
            xmin, ymax,
            xmin, ymin),
          ncol = 2, byrow = TRUE
        )))
      }
    )
  ) |>
  st_as_sf(crs = 4326) |>
  st_transform(3857)

ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_rect(
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "white", alpha = 0.60
  ) +
  geom_sf(
    data = heat_sf,
    aes(fill = log1p(visits)),
    color = NA,
    alpha = 0.85
  ) +
  coord_sf(crs = st_crs(3857), expand = FALSE) +
  theme_minimal() +
  labs(
    title = paste0(
      "GeoLife state visit heatmap (completed paths; cell_deg = ", cell_deg, ")"
    ),
    fill = "log(1 + visits)"
  )

# ------------------------------------------------------------
# 6) Estimate a strict lazy 4-neighbor random walk from completed paths
#    Allowed steps: stay, N, S, E, W
# ------------------------------------------------------------

# Expand completed gx/gy into a long table of unit steps (no stitching)
steps <- seq_tbl |>
  transmute(traj_id, gx = gx_c, gy = gy_c) |>
  tidyr::unnest(c(gx, gy)) |>
  group_by(traj_id) |>
  mutate(
    dx = gx - lag(gx),
    dy = gy - lag(gy)
  ) |>
  filter(!is.na(dx), !is.na(dy)) |>
  ungroup()

step_counts <- steps |>
  mutate(step = case_when(
    dx ==  0 & dy ==  0 ~ "stay",
    dx ==  1 & dy ==  0 ~ "E",
    dx == -1 & dy ==  0 ~ "W",
    dx ==  0 & dy ==  1 ~ "N",
    dx ==  0 & dy == -1 ~ "S",
    TRUE                ~ "other"
  )) |>
  count(step, name = "n") |>
  arrange(desc(n))

step_counts

# Under path completion, "other" should be zero (or extremely close)
frac_other <- step_counts |>
  summarise(frac_other = sum(n[step == "other"], na.rm = TRUE) / sum(n)) |>
  pull(frac_other)

frac_other

# MLE step probabilities for the strict lazy RW
rw_hat <- step_counts |>
  filter(step %in% c("stay", "N", "S", "E", "W")) |>
  mutate(p = n / sum(n)) |>
  arrange(match(step, c("stay", "N", "S", "E", "W")))

rw_hat

# ------------------------------------------------------------
# 7) Markovchain object for the STEP process (useful for prediction/simulation)
#    Under the RW assumption, steps are i.i.d. with distribution rw_hat.
#    In Markovchain form: each row is the same distribution.
# ------------------------------------------------------------
# 4-neighbor step states
step_states <- c("N", "S", "E", "W")

# step probabilities (MLE) from your step_counts table
p_step <- step_counts |>
  filter(step %in% step_states) |>
  mutate(p = n / sum(n)) |>
  select(step, p) |>
  tibble::deframe()

# enforce ordering
p_step <- p_step[step_states]

# transition matrix: identical rows, each equal to p_step
P_step <- matrix(
  rep(as.numeric(p_step), times = length(step_states)),
  nrow = length(step_states),
  byrow = TRUE,
  dimnames = list(step_states, step_states)
)

P_step

mc_steps <- new('markovchain',
                states = step_states,
                byrow = TRUE,
                transitionMatrix = P_step)

mc_steps

# ------------------------------------------------------------
# Reconstruct a grid-cell path from simulated steps
# ------------------------------------------------------------

set.seed(1)
sim_steps <- markovchain::rmarkovchain(n = 200, object = mc_steps, t0 = "N")

# map step labels -> (dx, dy)
step_to_dxdy <- function(step) {
  switch(step,
         "N" = c(0L,  1L),
         "S" = c(0L, -1L),
         "E" = c(1L,  0L),
         "W" = c(-1L, 0L),
         stop("Unknown step: ", step))
}

reconstruct_path <- function(steps, gx0 = 0L, gy0 = 0L) {
  inc <- t(vapply(steps, step_to_dxdy, FUN.VALUE = c(0L, 0L)))
  dx <- inc[, 1]; dy <- inc[, 2]
  tibble(
    k  = 0:length(steps),
    gx = c(gx0, gx0 + cumsum(dx)),
    gy = c(gy0, gy0 + cumsum(dy))
  )
}

path_sim <- reconstruct_path(sim_steps, gx0 = 0L, gy0 = 0L)

ggplot(path_sim, aes(gx, gy)) +
  geom_path(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 0.8, alpha = 0.5) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Simulated trajectory",
       x = "gx", y = "gy")


# ------------------------------------------------------------
# Is the increment process Markov?
# ------------------------------------------------------------

step_states <- c("N", "S", "E", "W")

# steps is your long table of unit steps from completed paths:
# columns: traj_id, dx, dy (constructed earlier)
steps_labeled <- steps |>
  mutate(step = case_when(
    dx ==  1 & dy ==  0 ~ "E",
    dx == -1 & dy ==  0 ~ "W",
    dx ==  0 & dy ==  1 ~ "N",
    dx ==  0 & dy == -1 ~ "S",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(step)) |>
  arrange(traj_id)  # keep trajectory grouping stable

# build step-to-step transition counts (no stitching across trajectories)
N_step <- with(steps_labeled |>
                 group_by(traj_id) |>
                 mutate(step_next = lead(step)) |>
                 filter(!is.na(step_next)) |>
                 ungroup(),
               table(factor(step, levels = step_states),
                     factor(step_next, levels = step_states)))

# transition matrix
P_step_dep <- prop.table(N_step, margin = 1) |> unclass()

mc_steps_dep <- new("markovchain",
                    states = step_states,
                    byrow = TRUE,
                    transitionMatrix = P_step_dep,
                    name = "4-neighbor step process (step-to-step Markov)"
)

mc_steps_dep

# long-run step-direction frequencies
steadyStates(mc_steps_dep)           
