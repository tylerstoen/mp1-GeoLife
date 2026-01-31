# ============================================================
# GeoLife: helpers + sampling 100 trajectories from one locale
# ============================================================

# ---- packages ----
library(dplyr)
library(tibble)
library(fs)
library(readr)
library(lubridate)

# ---- set your dataset root (the "Data/" folder) ----
geolife_root <- "data/Geolife Trajectories 1.3/Data"  # change to your path


# ============================================================
# 1) File system helpers
# ============================================================

# List user folders (e.g., "000", "001", ...)
list_users <- function(root = geolife_root) {
  dir_ls(root, type = "directory", recurse = FALSE) |>
    path_file() |>
    sort()
}

# List .plt files for a user
list_user_plts <- function(user_id, root = geolife_root) {
  dir_ls(path(root, user_id, "Trajectory"), glob = "*.plt") |>
    sort()
}


# ============================================================
# 2) Reading helpers
# ============================================================

# Read a full trajectory (.plt)
# GeoLife .plt format: 6 header lines, then columns:
# lat, lon, 0, altitude(feet), days, date, time
read_plt <- function(plt_path) {
  read_csv(
    plt_path,
    skip = 6,
    col_names = c("lat","lon","zero","alt_ft","days","date","time"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    mutate(
      # parse timestamp (UTC is fine for most teaching purposes)
      t = ymd_hms(paste(date, time), tz = "UTC")
    ) |>
    select(t, lat, lon) |>
    arrange(t)
}

# Read only the first few points of a trajectory (FAST).
# This is used to approximate the "locale" of a trajectory without
# reading the entire file.
read_plt_head <- function(plt_path, n_pts = 50) {
  readr::read_csv(
    plt_path,
    skip = 6,
    n_max = n_pts,
    col_names = c("lat","lon","zero","alt_ft","days","date","time"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    dplyr::summarise(
      lat = stats::median(lat, na.rm = TRUE),
      lon = stats::median(lon, na.rm = TRUE)
    )
}


# ============================================================
# 3) Geometry helper: distance in km (haversine)
# ============================================================

# Haversine distance in kilometers (no extra packages)
haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371.0
  to_rad <- function(x) x * pi / 180
  
  lat1 <- to_rad(lat1); lon1 <- to_rad(lon1)
  lat2 <- to_rad(lat2); lon2 <- to_rad(lon2)
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}


# ============================================================
# 4) Main function: sample n trajectories from a single locale
# ============================================================

# Strategy:
#  - Build a table of all (user, file) pairs
#  - Sample a *candidate* set of files
#  - Approximate each candidate file's location using read_plt_head()
#  - Keep only candidates within radius_km of a chosen center
#  - Fully read the final sample (only n files are read completely)
sample_geolife_same_locale_fast <- function(
    n = 100,
    root = geolife_root,
    seed = 1,
    radius_km = 40,          # larger radius = easier to find enough trajectories
    center = NULL,           # c(lat=..., lon=...) or NULL to auto-pick
    n_candidates = 1200,     # how many files to scan cheaply for locale
    head_pts = 50            # how many points to read from each file in the scan
) {
  set.seed(seed)
  
  # ---- (a) build the file table (one row per trajectory file) ----
  users <- list_users(root)
  
  file_tbl <- bind_rows(lapply(users, function(u) {
    files <- list_user_plts(u, root)
    if (length(files) == 0) return(NULL)
    tibble(user = u, file = files)
  }))
  
  if (nrow(file_tbl) == 0) {
    warning("No .plt files found. Check geolife_root points at the GeoLife Data/ folder.")
    return(tibble())
  }
  
  # ---- (b) sample candidate files to find a consistent locale ----
  # We add an explicit index so we can merge candidate locations back
  # without relying on row order.
  cand <- file_tbl |>
    dplyr::slice_sample(n = min(n_candidates, nrow(file_tbl))) |>
    dplyr::mutate(idx = dplyr::row_number())
  
  # ---- (c) cheaply estimate lat/lon for each candidate ----
  # For each candidate file, read only the first `head_pts` points and
  # summarize them to a single (lat, lon) "center" (median).
  # If a read fails, record NA so we can drop it cleanly.
  heads <- dplyr::bind_rows(lapply(seq_len(nrow(cand)), function(k) {
    out <- try(read_plt_head(cand$file[k], n_pts = head_pts), silent = TRUE)
    
    if (inherits(out, "try-error") || !is.data.frame(out)) {
      tibble::tibble(idx = cand$idx[k], lat = NA_real_, lon = NA_real_)
    } else {
      tibble::tibble(idx = cand$idx[k], lat = out$lat[1], lon = out$lon[1])
    }
  }))
  
  # ---- (d) attach the approximated centers back onto the candidate table ----
  # Drop candidates whose "head read" failed.
  cand <- cand |>
    dplyr::left_join(heads, by = "idx") |>
    dplyr::filter(!is.na(lat), !is.na(lon))
  
  if (nrow(cand) == 0) {
    warning("All candidate reads failed. Check file permissions / file format / path.")
    return(tibble())
  }
  
  # ---- (d) pick a center trajectory if none is supplied ----
  if (is.null(center)) {
    k <- sample.int(nrow(cand), 1)
    center <- c(lat = cand$lat[k], lon = cand$lon[k])
  }
  
  # ---- (e) filter to the locale ----
  cand$dist_km <- haversine_km(cand$lat, cand$lon, center["lat"], center["lon"])
  local <- cand |> filter(dist_km <= radius_km)
  
  if (nrow(local) == 0) {
    warning("No candidates within radius. Try increasing radius_km or n_candidates.")
    return(tibble())
  }
  
  # ---- (f) choose the final sample from the local candidates ----
  n_take <- min(n, nrow(local))
  samp <- local |>
    slice_sample(n = n_take) |>
    mutate(
      traj_id = row_number(),        # numeric id 1..n_take within the sample
      traj_file = basename(file)     # file name only (not full path)
    )
  
  # ---- (g) read the FULL trajectories for the selected sample ----
  dat_list <- lapply(seq_len(nrow(samp)), function(k) {
    out <- try(read_plt(samp$file[k]), silent = TRUE)
    if (inherits(out, "try-error") || !is.data.frame(out)) return(NULL)
    
    out |>
      mutate(
        user = samp$user[k],
        traj_id = samp$traj_id[k],
        traj_file = samp$traj_file[k]
      )
  })
  
  pts <- bind_rows(dat_list)
  
  # If everything failed, return an empty tibble (but with expected columns)
  if (nrow(pts) == 0) {
    warning("Selected files could not be read (all failed).")
    return(tibble(user = character(), traj_id = integer(), traj_file = character(),
                  t = as.POSIXct(character()), lat = numeric(), lon = numeric()))
  }
  
  # ---- (h) return tidy output ----
  pts |>
    select(user, traj_id, traj_file, everything())
}