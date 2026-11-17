# Shared helper functions for the Zion multi-tracer scripts.
# Sourced by 00_setup.R, so every numbered script can use them.
#
# Contents
#   General        decay(), dec_year(), gc_km(), water_year()
#   Paths, checks  raw_path(), ext_path(), check_close()
#   Sites          SITE_RENAME, EMERGENCE_POINTS, apply_emergence_elevation(),
#                  read_tritium_raw(), site_locations()
#   Outlets        read_emergence_points(), contact_points(), contact_plane_fit()
#   Baseflow       recession_alpha(), bfimax_backward(), baseflow_eckhardt()


# ---- General ----------------------------------------------------------------

# Radioactive decay: concentration remaining after t years.
#   c0         starting concentration
#   t          elapsed time (years)
#   half_life  half-life (years)
decay <- function(c0, t, half_life) {
  c0 * 0.5^(t / half_life)
}

# Decimal year of a Date, taken at mid-day (e.g. 2 July 1988 -> 1988.50).
dec_year <- function(d) {
  lubridate::year(d) + (lubridate::yday(d) - 0.5) / 365.25
}

# Great-circle distance (km) between two points given in decimal degrees, on
# a sphere of radius 6,371 km (haversine formula).
gc_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180   # degrees to radians
  a <- sin((lat2 - lat1) * rad / 2)^2 +
    cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
  6371 * 2 * asin(sqrt(a))
}

# Water year (October-September), labeled by the calendar year in which it
# ends (e.g. 1 Oct 1987 - 30 Sep 1988 is WY1988).
water_year <- function(d) {
  as.integer(lubridate::year(d) + (lubridate::month(d) >= 10))
}


# ---- Paths and checks -------------------------------------------------------

# Path to a raw input copied by 01_copy_raw_data.R. Stops with a clear message
# if the copy is missing.
raw_path <- function(...) {
  p <- file.path(RAW, ...)
  if (!file.exists(p)) {
    stop("Missing raw input: ", p,
         " -- run paper/scripts/01_copy_raw_data.R first")
  }
  p
}

# Path to a cached download from 02_download_external.R. Stops with a clear
# message if the file is missing.
ext_path <- function(...) {
  p <- file.path(EXT, ...)
  if (!file.exists(p)) {
    stop("Missing cached download: ", p,
         " -- run paper/scripts/02_download_external.R first")
  }
  p
}

# Regression check: compare a new result with the value reported in the paper
# and print a line starting with [OK] or [DIFF]. run_all.R collects these
# lines.
#   new, ref  numeric vectors of the same length (new result, reported value)
#   tol       largest absolute difference that still counts as OK
#   label     what is being checked
check_close <- function(new, ref, tol, label) {
  max_diff <- max(abs(new - ref), na.rm = TRUE)
  ok <- isTRUE(all(abs(new - ref) <= tol, na.rm = TRUE))
  status <- dplyr::if_else(ok, "OK", "DIFF")
  message(sprintf("[%s] %s (max |diff| = %.4g, tol %.4g)",
                  status, label, max_diff, tol))
  invisible(ok)
}


# ---- Sites ------------------------------------------------------------------

# hard-coded: site names in Zion_Tritium.csv -> names used in the paper. The
# Birch Spring sample was collected at the upper sediment trap.
SITE_RENAME <- c("Birch Sediment Trap" = "Birch Spring")

# hard-coded: spring emergence points that replace the recorded site
# coordinates, picked from high-resolution Google Earth imagery (K. Cognac,
# November 2026; raw/Spring Emergence Points from google earth.kml), with the
# Google Earth elevation of the point. On cliff faces the 8-m DEM can differ
# from these elevations by more than 100 m within 15 m, so the Google Earth
# elevation is used for these sites (apply_emergence_elevation()).
#   Birch Spring: sampled at the upper sediment trap (37.241474, -112.974323),
#     which the 2025 field notes place about 1,000 ft from and 400 ft below
#     the discharge point.
#   Refrigerator Spring: the coordinates recorded for the 1988-89 samples
#     (37.282778, -112.953611; rounded to arc-seconds) plot 411 m from the
#     emergence point visible in the imagery.
EMERGENCE_POINTS <- tibble::tribble(
  ~Location,             ~latitude,  ~longitude,   ~elev_ge_m,
  "Birch Spring",        37.243310,  -112.974886,  1469,
  "Refrigerator Spring", 37.285953,  -112.951227,  1383
)

# Replace the elevation column `col` of a site table with the Google Earth
# elevation, for the sites listed in EMERGENCE_POINTS. Other sites keep their
# value. `d` needs a Location column.
apply_emergence_elevation <- function(d, col = "elev_m") {
  d %>%
    dplyr::left_join(
      EMERGENCE_POINTS %>% dplyr::select(Location, .elev_ge = elev_ge_m),
      by = "Location"
    ) %>%
    dplyr::mutate("{col}" := dplyr::coalesce(.elev_ge, .data[[col]])) %>%
    dplyr::select(-.elev_ge)
}

# Tritium and conductivity samples (one row per sample), with sites renamed
# as in the paper.
read_tritium_raw <- function() {
  readr::read_csv(raw_path("Tritium/Park_Samples/Zion_Tritium.csv"),
                  show_col_types = FALSE) %>%
    dplyr::mutate(Location = dplyr::recode(Location, !!!SITE_RENAME))
}

# Site locations: one row per site with latitude, longitude and the source of
# the coordinates. Sources, in order of priority:
#   1. emergence point picked from imagery (EMERGENCE_POINTS)
#   2. 2025 field GPS fix, for the sites sampled in 2025
#   3. coordinates in Zion_Tritium.csv. For the 1988-89 springs these are
#      approximate, as recorded by Kimball and Christensen (1996).
# Lower Grotto Sandtrap stays at its sand trap, which is downslope of the
# emergence point (a rockfall cave).
site_locations <- function() {
  # 3. coordinates recorded with the samples
  recorded <- read_tritium_raw() %>%
    dplyr::distinct(Location, latitude, longitude)

  # 2. 2025 field GPS fixes
  # hard-coded: tritium site name -> location name in the 2025 field summary
  gps_names <- tibble::tribble(
    ~Location,               ~location,
    "East Entrance Well",    "East Entrance Well",
    "Taylor Creek Well",     "Taylor Creek Well",
    "Lower Grotto Sandtrap", "Lower Grotto Spring",
    "Birch Spring",          "Birch Spring"
  )
  gps <- readr::read_csv(
    raw_path("site_visit/2025_July_Zion_CEBR_PIPE_Sampling_Summary.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(location = stringr::str_trim(location)) %>%
    dplyr::inner_join(gps_names, by = "location") %>%
    dplyr::select(Location, lat_gps = latitude, lon_gps = longitude)

  # 1. emergence points
  emergence <- EMERGENCE_POINTS %>%
    dplyr::select(Location, lat_em = latitude, lon_em = longitude)

  recorded %>%
    dplyr::left_join(gps, by = "Location") %>%
    dplyr::left_join(emergence, by = "Location") %>%
    dplyr::mutate(
      coord_source = dplyr::case_when(
        !is.na(lat_em)  ~ "emergence point (imagery)",
        !is.na(lat_gps) ~ "2025 field GPS",
        TRUE            ~ "Zion_Tritium.csv"
      ),
      latitude  = dplyr::coalesce(lat_em, lat_gps, latitude),
      longitude = dplyr::coalesce(lon_em, lon_gps, longitude)
    ) %>%
    dplyr::select(Location, latitude, longitude, coord_source)
}


# ---- Spring outlets relative to the Navajo-Kayenta contact (11, 16, 19) ----

# Spring emergence points picked on Google Earth imagery (K. Cognac, November
# 2026), read from the KML copied by 01_copy_raw_data.R. Returns an sf point
# layer in `crs` with id, label (placemark name), longitude, latitude,
# elev_ge_m (Google Earth elevation stored in the KML), generic (TRUE when the
# name is only "... emergence") and n_placemarks.
# Placemarks closer than 5 m are the same emergence labelled twice; they are
# merged into one point, with the specific name where one of the pair has it
# and the mean of the Google Earth elevations.
read_emergence_points <- function(crs) {
  kml <- sf::st_read(raw_path("Spring Emergence Points from google earth.kml"),
                     quiet = TRUE)
  xyz <- sf::st_coordinates(kml)
  em <- kml %>%
    sf::st_zm() %>%
    dplyr::transmute(id = dplyr::row_number(), label = Name,
                     longitude = xyz[, 1], latitude = xyz[, 2],
                     elev_ge_m = xyz[, 3]) %>%
    sf::st_transform(crs)
  # grp_id = first placemark within 5 m of each placemark
  dup <- sf::st_distance(em) %>% units::drop_units()
  grp_id <- purrr::map_int(seq_len(nrow(em)), ~ min(which(dup[.x, ] < 5)))
  out <- em %>%
    dplyr::mutate(generic = stringr::str_detect(
                    label, stringr::regex("emergence", ignore_case = TRUE)),
                  grp_id = grp_id) %>%
    dplyr::group_by(grp_id) %>%
    dplyr::mutate(elev_ge_m = mean(elev_ge_m), n_placemarks = dplyr::n()) %>%
    dplyr::arrange(generic, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-grp_id)
  message(nrow(kml) - nrow(out), " duplicate placemark(s) merged")
  out
}

# Points every `spacing` m along the mapped Navajo-Kayenta contact (an sf line,
# UTM), with the elevation `z` from `dem` (a SpatRaster with layer "elev")
contact_points <- function(contact, dem, spacing = 20) {
  pts <- sf::st_line_sample(sf::st_cast(contact, "LINESTRING"),
                            density = 1 / spacing) %>%
    sf::st_cast("POINT") %>%
    sf::st_sf()
  xy <- sf::st_coordinates(pts)
  pts %>%
    dplyr::mutate(x = xy[, "X"], y = xy[, "Y"],
                  z = terra::extract(dem, terra::vect(pts))$elev) %>%
    sf::st_drop_geometry()
}

# Local contact elevation at each point of `pts` (sf, UTM): a plane
# z = a + bx + cy fitted to the contact points (from contact_points()) within
# `radius` m, with coordinates centered on the point so the intercept is the
# estimate. The residual standard deviation is the uncertainty. Needs at least
# 10 contact points; otherwise NA.
# Returns one row per point: contact_elev_m, contact_sd_m, contact_n.
contact_plane_fit <- function(pts, cpts, radius = 1500) {
  xy <- sf::st_coordinates(pts)
  purrr::map2_dfr(xy[, "X"], xy[, "Y"], function(px, py) {
    near <- cpts %>%
      dplyr::filter(!is.na(z), sqrt((x - px)^2 + (y - py)^2) <= radius) %>%
      dplyr::mutate(x = x - px, y = y - py)
    if (nrow(near) < 10) {
      return(tibble::tibble(contact_elev_m = NA_real_, contact_sd_m = NA_real_,
                            contact_n = nrow(near)))
    }
    m <- lm(z ~ x + y, data = near)
    tibble::tibble(contact_elev_m = unname(coef(m)[1]),
                   contact_sd_m = sigma(m), contact_n = nrow(near))
  })
}


# ---- Baseflow (14_baseflow_climate.R) ---------------------------------------

# Daily recession constant alpha: the median ratio Q[t+1] / Q[t] over
# recession days.
#   q               daily flows, in date order
#   min_len         a recession is a run of at least this many declining days
#   skip            first days of each run to drop (quickflow tail)
#   q_max_quantile  keep only days with flow at or below this quantile of all
#                   flows
# Returns a list with alpha and the number of recession days used.
recession_alpha <- function(q, min_len = 7, skip = 3, q_max_quantile = 0.5) {
  q_max <- quantile(q, q_max_quantile, na.rm = TRUE)

  recession_days <- tibble::tibble(q = q, q_next = dplyr::lead(q)) %>%
    dplyr::mutate(
      # a day is "declining" if the next day's flow is lower
      declining = dplyr::coalesce(q_next < q, FALSE),
      # number the runs of consecutive declining (or not declining) days
      run = dplyr::consecutive_id(declining)
    ) %>%
    dplyr::group_by(run) %>%
    dplyr::mutate(run_len = dplyr::n(), day_in_run = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      declining, run_len >= min_len,
      day_in_run > skip,        # drop the first `skip` days of the run
      day_in_run < run_len,     # the last day's ratio leaves the run
      q > 0, q <= q_max
    )

  list(alpha = median(recession_days$q_next / recession_days$q),
       n_days = nrow(recession_days))
}

# Maximum baseflow index from the backward filter of Collischonn and Fan
# (2013). Working back from the last day, b[t] = min(b[t+1] / alpha, Q[t]);
# then BFImax = sum(b) / sum(Q).
bfimax_backward <- function(q, alpha) {
  n <- length(q)
  # accumulate() walks backward from the last day; `b_next` is b[t+1]. With
  # .dir = "backward" the running value is the second argument.
  b <- purrr::accumulate(q[-n], function(q_t, b_next) min(b_next / alpha, q_t),
                         .init = q[n], .dir = "backward")
  sum(b) / sum(q)
}

# Eckhardt (2005) two-parameter recursive digital filter for baseflow
# separation. Baseflow never exceeds streamflow (b[t] <= Q[t]), and the first
# day starts at b[1] = BFImax * Q[1].
#   q       daily flows, in date order
#   alpha   recession constant
#   bfimax  maximum baseflow index
# Returns the daily baseflow, the same length as q.
baseflow_eckhardt <- function(q, alpha, bfimax) {
  # accumulate() walks forward from day 2; `b_prev` is b[t-1]
  purrr::accumulate(
    q[-1],
    function(b_prev, q_t) {
      min(((1 - bfimax) * alpha * b_prev + (1 - alpha) * bfimax * q_t) /
            (1 - alpha * bfimax),
          q_t)
    },
    .init = bfimax * q[1]
  )
}
