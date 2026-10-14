# 10_response_time.R -- hydraulic response time of the Navajo aquifer (paper
# sections 2.1.3, 4.4; ESM Fig. S4).
#
# Equations from Ferguson et al. (2026):
#   confined    t_NE = (12 / pi^2) S L^2 / T
#   unconfined  t_NE = a Sy L^2 / T, with a = 0.7-3.6
#   mixed       t_NE = 3 Sy Lu (Lc + Lu / 2) / T
#               (Lc, Lu = confined and unconfined lengths of the path)
# Parameters (hard-coded below in run_case()):
#   T  2,500-14,000 ft2/d, log-uniform (Cordova 1981)
#   S  1e-4 to 0.012, log-uniform (Cordova 1981; McGavock and Levings 1973)
#   Sy 0.04-0.17, uniform (Cordova 1978)
#   a  0.7-3.6, uniform
#   Sensitivity cases: "heilweil" (T 100-19,000 ft2/d, S 0.001-0.002;
#   Heilweil et al. 2000) and "blackmesa" (T 75-350 ft2/d; McGavock and
#   Levings 1973). 200,000 draws.
# Flow paths: canyon rim (Navajo-Temple Cap contact distance), nearest high
# ground >= 2,400 m east and north, regional paths (25th-75th percentile
# distance to that high ground) and the Markagunt Plateau (Cedar Breaks NM to
# Brian Head Peak).
# Confined share of each path: from mapped Temple Cap/Carmel/Cretaceous cover
# (low) to everything except exposed Navajo/Kayenta (high), sampled every
# 100 m along straight paths.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
})
sf_use_s2(FALSE)

site_tab <- read_csv(file.path(DER, "site_geology.csv"), show_col_types = FALSE)
springs <- site_tab %>%
  filter(Location %in% REPEAT_SPRINGS) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM)
dem <- rast(ext_path("dem_z9.tif"))
names(dem) <- "elev_m"

# Spring outlet elevations from the zoom-14 terrain tiles (same source as
# Figs. 5-6). The 127-m dem_z9 cells blend canyon floor, cliff and rim and
# read 71-222 m high at the canyon springs, so dem_z9 is used only for
# plateau-scale features.
sp_elev <- read_csv(ext_path("site_elevations_z14.csv"),
                    show_col_types = FALSE) %>%
  filter(Location %in% REPEAT_SPRINGS) %>%
  pull(elev_m)
origin <- st_centroid(st_union(springs))   # center of the five springs
o_xy <- st_coordinates(origin)
o_elev <- median(sp_elev)


# ---- Candidate recharge zones -----------------------------------------------
# hard-coded: DEM cells >= 2,400 m within 60 km of the springs, in the north
# (azimuth 315-45 deg) and east (45-135 deg) sectors

xy <- as.data.frame(dem, xy = TRUE, na.rm = TRUE) %>%
  mutate(
    dx = x - o_xy[1],
    dy = y - o_xy[2],
    dist_m = sqrt(dx^2 + dy^2),
    az = (atan2(dx, dy) * 180 / pi) %% 360,
    sector = case_when(
      az >= 315 | az < 45 ~ "North",
      az >= 45 & az < 135 ~ "East",
      TRUE ~ NA_character_
    )
  )
hi <- xy %>% filter(!is.na(sector), elev_m >= 2400, dist_m <= 60000)

zones <- hi %>%
  group_by(sector) %>%
  summarise(
    n_cells = n(),
    elev_median_m = median(elev_m),
    elev_max_m = max(elev_m),
    L_min_km = min(dist_m) / 1000,
    L_p25_km = quantile(dist_m, .25) / 1000,
    L_median_km = median(dist_m) / 1000,
    L_p75_km = quantile(dist_m, .75) / 1000,
    .groups = "drop"
  ) %>%
  # median spring outlet to median high ground
  mutate(topo_gradient_at_median = (elev_median_m - o_elev) /
           (L_median_km * 1000))

local <- site_tab %>%
  filter(Location %in% REPEAT_SPRINGS) %>%
  select(Location, dist_navajo_templecap_contact_m)


# ---- Confined and unconfined shares along straight paths --------------------

# hard-coded GRI map-unit classes: exposed Navajo/Kayenta (unconfined) vs
# Temple Cap/Carmel/Cretaceous cover (confining)
geo <- st_read(raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb"),
               "zionglg", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(UTM) %>%
  mutate(cls = case_when(
    GLG_SYM %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl",
                   "Jk", "Jkt", "Jku", "Jks") ~ "exposed",
    grepl("^Jt|^Jc|^K", GLG_SYM) ~ "cover",
    grepl("^Qb", GLG_SYM) ~ "basalt",
    grepl("^Q", GLG_SYM) ~ "quaternary",
    TRUE ~ "other"
  ))

# Fraction of a straight path (springs -> xy_end) on cover and on exposed
# Navajo/Kayenta, sampled every 100 m
split_to <- function(xy_end) {
  ln <- st_sfc(st_linestring(rbind(o_xy[1:2], xy_end)), crs = UTM)
  pts <- st_line_sample(ln, density = 1 / 100) %>%
    st_cast("POINT") %>%
    st_as_sf()
  j <- st_join(pts, geo %>% select(cls), join = st_within)
  j$cls[is.na(j$cls)] <- "unmapped"
  tibble(
    L_km = as.numeric(st_length(ln)) / 1000,
    cover = mean(j$cls == "cover"),
    exposed = mean(j$cls == "exposed")
  )
}

# Path to the high-ground cell at distance quantile q of a sector
path_split <- function(sec, q) {
  d <- hi %>% filter(sector == sec)
  t <- d[which.min(abs(d$dist_m - quantile(d$dist_m, q))), ]
  split_to(c(t$x, t$y)) %>% mutate(sector = sec, quantile = q)
}
splits <- bind_rows(lapply(c("North", "East"), function(s) {
  bind_rows(path_split(s, 0), path_split(s, 0.5))
}))


# ---- Markagunt Plateau ------------------------------------------------------
# Cedar Breaks NM (nearest boundary point) and Brian Head Peak (highest cell
# within 8 km of 37.68 N, 112.83 W; see 02_download_external.R)

nps <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  st_make_valid()
cebr <- nps %>%
  filter(UNIT_CODE == "CEBR") %>%
  st_transform(UTM) %>%
  summarise()
bh_dem <- rast(ext_path("dem_brian_head_z11.tif"))
bh <- as.data.frame(bh_dem, xy = TRUE, na.rm = TRUE)
names(bh)[3] <- "elev"
bh <- bh[which.max(bh$elev), ]
bh_pt <- st_sfc(st_point(c(bh$x, bh$y)), crs = UTM)

basin <- st_read(ext_path("gage_09405501_basin.gpkg"), quiet = TRUE) %>%
  st_transform(UTM)
huc8 <- st_read(ext_path("huc_units.gpkg"), quiet = TRUE) %>%
  filter(level == "huc08") %>%
  st_transform(UTM)
cebr_near <- st_nearest_points(origin, cebr) %>%
  st_cast("POINT") %>%
  .[2]
bh_lonlat <- st_coordinates(st_transform(bh_pt, 4326))

mk <- tibble(
  brian_head_lat = bh_lonlat[2],
  brian_head_lon = bh_lonlat[1],
  brian_head_elev_m = bh$elev,
  dist_brian_head_km = as.numeric(st_distance(origin, bh_pt)) / 1000,
  dist_cebr_nearest_km = as.numeric(st_distance(origin, cebr)) / 1000,
  bh_in_upper_virgin = as.vector(st_within(bh_pt, huc8, sparse = FALSE))[1],
  cebr_in_gage_basin_frac = sum(as.numeric(st_area(st_intersection(cebr, basin)))) /
    as.numeric(st_area(cebr))
)
msplit <- bind_rows(
  split_to(st_coordinates(cebr_near)[1, 1:2]) %>%
    mutate(target = "Cedar Breaks NM (nearest boundary)"),
  split_to(st_coordinates(bh_pt)[1, 1:2]) %>%
    mutate(target = "Brian Head Peak")
)

write_csv(zones, file.path(DER, "recharge_zone_stats.csv"))
write_csv(splits, file.path(DER, "path_confinement_split.csv"))
write_csv(mk, file.path(DER, "markagunt_targets.csv"))
write_csv(msplit, file.path(DER, "markagunt_confinement_split.csv"))


# ---- Monte Carlo ------------------------------------------------------------

FT <- 0.3048     # m per ft
N <- 200000      # draws per scenario

# Range of the confined share: mapped cover (low) to all but exposed (high)
lc <- function(s) c(min(s$cover), max(1 - s$exposed))

# Flow-path scenarios: L = range of path length (km), lc = confined share
scen <- list(
  list(name = "Local rim recharge (Path B)",
       L = range(local$dist_navajo_templecap_contact_m) / 1000,
       lc = c(0, 0)),
  list(name = "East plateau, nearest high ground",
       L = rep(zones$L_min_km[zones$sector == "East"], 2),
       lc = lc(filter(splits, sector == "East", quantile == 0))),
  list(name = "North (Kolob), nearest high ground",
       L = rep(zones$L_min_km[zones$sector == "North"], 2),
       lc = lc(filter(splits, sector == "North", quantile == 0))),
  list(name = "Regional, N and E high ground (median)",
       L = c(min(zones$L_p25_km), max(zones$L_p75_km)),
       lc = lc(filter(splits, quantile == 0.5))),
  list(name = "Markagunt Plateau (Cedar Breaks to Brian Head)",
       L = c(mk$dist_cebr_nearest_km, mk$dist_brian_head_km),
       lc = lc(msplit))
)

run_case <- function(case) {
  set.seed(SEED)
  # hard-coded parameter ranges (see header); T in ft2/d
  T <- switch(case,
              primary   = exp(runif(N, log(2500), log(14000))),
              heilweil  = exp(runif(N, log(100), log(19000))),
              blackmesa = exp(runif(N, log(75), log(350))))
  S <- if (case == "heilweil") {
    runif(N, 0.001, 0.002)
  } else {
    exp(runif(N, log(1e-4), log(0.012)))
  }
  Sy <- runif(N, 0.04, 0.17)
  a <- runif(N, 0.7, 3.6)

  map_dfr(scen, function(s) {
    L <- runif(N, s$L[1], s$L[2]) * 1000 / FT   # path length, ft
    fc <- runif(N, s$lc[1], s$lc[2])            # confined share
    Lc <- fc * L
    Lu <- (1 - fc) * L
    # response times in years
    t_conf <- (12 / pi^2) * S * L^2 / T / 365.25
    t_unc <- a * Sy * L^2 / T / 365.25
    t_mix <- 3 * Sy * Lu * (Lc + Lu / 2) / T / 365.25
    t_use <- if (s$lc[2] == 0) t_unc else t_mix

    map_dfr(
      list("mixed/used" = t_use, "all confined" = t_conf,
           "all unconfined" = t_unc),
      function(t) {
        q <- quantile(t, c(.05, .25, .5, .75, .95))
        tibble(p5 = q[1], p25 = q[2], median = q[3], p75 = q[4], p95 = q[5])
      },
      .id = "case_detail"
    ) %>%
      mutate(case = case, scenario = s$name, L_km_lo = s$L[1],
             L_km_hi = s$L[2], Lc_lo = s$lc[1], Lc_hi = s$lc[2], .before = 1)
  })
}

res <- map_dfr(c("primary", "heilweil", "blackmesa"), run_case)
write_csv(res, file.path(TAB, "response_time_results.csv"))
res %>%
  filter(case == "primary") %>%
  select(scenario, case_detail, p5, median, p95) %>%
  as.data.frame() %>%
  print(digits = 3)


# ---- Regression checks against the paper (section 4.4) ----------------------

# Result for one case / scenario / detail
md <- function(c, sc, d = "mixed/used", col = "median") {
  res[[col]][res$case == c & res$scenario == sc & res$case_detail == d]
}

check_close(c(zones$L_min_km[zones$sector == "East"],
              zones$L_min_km[zones$sector == "North"]),
            c(8.54, 13.88), 0.01, "nearest high ground east and north (km)")
check_close(c(mk$dist_cebr_nearest_km, mk$dist_brian_head_km,
              mk$brian_head_elev_m),
            c(37.72, 46.47, 3446.8), 0.05,
            "Markagunt distances and Brian Head elevation")

ref <- c(
  md("primary", scen[[1]]$name),
  md("primary", scen[[2]]$name),
  md("primary", scen[[3]]$name),
  md("primary", scen[[4]]$name),
  md("primary", scen[[4]]$name, "all confined"),
  md("primary", scen[[4]]$name, "all unconfined"),
  md("primary", scen[[5]]$name),
  md("primary", scen[[5]]$name, "all confined"),
  md("primary", scen[[5]]$name, "all unconfined")
)
# hard-coded: medians from the earlier Python Monte Carlo (reported values)
old <- c(0.1094, 34.58, 97.56, 461.1, 7.815, 1145.8, 761.1, 11.74, 1729.0)
message("relative change vs Python Monte Carlo: ",
        paste(sprintf("%+.1f%%", 100 * (ref / old - 1)), collapse = " "))
check_close(ref / old, rep(1, length(old)), 0.03,
            "response-time medians within 3% of the reported values")
message(sprintf(
  "Heilweil / primary median ratio (regional): %.1f; Black Mesa medians: rim %.1f yr, regional %.0f yr",
  md("heilweil", scen[[4]]$name) / md("primary", scen[[4]]$name),
  md("blackmesa", scen[[1]]$name),
  md("blackmesa", scen[[4]]$name)
))
