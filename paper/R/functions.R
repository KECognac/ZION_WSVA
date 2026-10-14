# Shared helper functions for the Zion multi-tracer scripts.
# Sourced by 00_setup.R, so every numbered script can use them.


# ---- General ----------------------------------------------------------------

# Radioactive decay: concentration remaining after t years.
decay <- function(c0, t, half_life) {
  c0 * 0.5^(t / half_life)
}

# Decimal year of a Date (mid-day of the day of year).
dec_year <- function(d) {
  as.numeric(format(d, "%Y")) + (as.numeric(format(d, "%j")) - 0.5) / 365.25
}

# Great-circle distance (km) on a sphere of radius 6,371 km.
gc_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  6371 * 2 * asin(sqrt(a))
}

# Water year (October-September), labeled by the calendar year in which it
# ends (e.g. 1 Oct 1987 - 30 Sep 1988 is WY1988).
water_year <- function(d) {
  as.integer(format(d, "%Y")) + (as.integer(format(d, "%m")) >= 10)
}


# ---- Paths and checks -------------------------------------------------------

# Path to a raw input copied by 01_copy_raw_data.R; stops with a clear
# message if the copy is missing.
raw_path <- function(...) {
  p <- file.path(RAW, ...)
  if (!file.exists(p)) {
    stop("Missing raw input: ", p,
         " -- run paper/scripts/01_copy_raw_data.R first")
  }
  p
}

# Path to a cached download from 02_download_external.R.
ext_path <- function(...) {
  p <- file.path(EXT, ...)
  if (!file.exists(p)) {
    stop("Missing cached download: ", p,
         " -- run paper/scripts/02_download_external.R first")
  }
  p
}

# Regression check: compare a new result with the value reported in the paper
# and print [OK] or [DIFF]. run_all.R collects these lines.
check_close <- function(new, ref, tol, label) {
  ok <- isTRUE(all(abs(new - ref) <= tol, na.rm = TRUE))
  message(sprintf("[%s] %s (max |diff| = %.4g, tol %.4g)",
                  if (ok) "OK" else "DIFF", label,
                  max(abs(new - ref), na.rm = TRUE), tol))
  invisible(ok)
}


# ---- Sites ------------------------------------------------------------------

# Site locations: coordinates from Zion_Tritium.csv, replaced by the 2025
# field GPS fixes for the 2025 sites. The 1988-89 spring locations are
# approximate, as recorded by Kimball and Christensen (1996).
site_locations <- function() {
  tr <- readr::read_csv(raw_path("Tritium/Park_Samples/Zion_Tritium.csv"),
                        show_col_types = FALSE) %>%
    dplyr::distinct(Location, latitude, longitude)
  fs <- readr::read_csv(
    raw_path("site_visit/2025_July_Zion_CEBR_PIPE_Sampling_Summary.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(location = stringr::str_trim(location))

  # hard-coded: tritium site name -> location name in the 2025 field summary
  gps25 <- c(
    "East Entrance Well"    = "East Entrance Well",
    "Taylor Creek Well"     = "Taylor Creek Well",
    "Lower Grotto Sandtrap" = "Lower Grotto Spring",
    "Birch Sediment Trap"   = "Birch Spring"
  )

  tr$coord_source <- "Zion_Tritium.csv"
  for (k in names(gps25)) {
    r <- fs[fs$location == gps25[[k]], ]
    tr$latitude[tr$Location == k] <- r$latitude
    tr$longitude[tr$Location == k] <- r$longitude
    tr$coord_source[tr$Location == k] <- "2025 field GPS"
  }
  tr
}


# ---- Baseflow (14_baseflow_climate.R) ---------------------------------------

# Daily recession constant alpha: median Q[t+1] / Q[t] on days inside runs of
# at least `min_len` consecutive declining flows. The first `skip` days of each
# run (quickflow tail) are dropped, and only days with Q below the
# `q_max_quantile` of all flows are kept. Returns alpha and the number of
# recession days used.
recession_alpha <- function(q, min_len = 7, skip = 3, q_max_quantile = 0.5) {
  declining <- c(diff(q) < 0, FALSE)
  runs <- rle(declining)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1

  long_runs <- which(runs$values & runs$lengths >= min_len)
  idx <- unlist(lapply(long_runs, function(i) (starts[i] + skip):(ends[i] - 1)))
  idx <- idx[q[idx] > 0 & q[idx] <= quantile(q, q_max_quantile, na.rm = TRUE)]

  list(alpha = median(q[idx + 1] / q[idx]), n_days = length(idx))
}

# Maximum baseflow index from the backward filter of Collischonn and Fan
# (2013): starting from the last day, b[t] = min(b[t+1] / alpha, Q[t]);
# BFImax = sum(b) / sum(Q).
bfimax_backward <- function(q, alpha) {
  n <- length(q)
  b <- numeric(n)
  b[n] <- q[n]
  for (t in (n - 1):1) {
    b[t] <- min(b[t + 1] / alpha, q[t])
  }
  sum(b) / sum(q)
}

# Eckhardt (2005) two-parameter recursive digital filter for baseflow
# separation, with b[t] <= Q[t] and b[1] = BFImax * Q[1].
baseflow_eckhardt <- function(q, alpha, bfimax) {
  b <- numeric(length(q))
  b[1] <- bfimax * q[1]
  for (t in 2:length(q)) {
    b[t] <- min(
      ((1 - bfimax) * alpha * b[t - 1] + (1 - alpha) * bfimax * q[t]) /
        (1 - alpha * bfimax),
      q[t]
    )
  }
  b
}
