# 16_spring_catchments.R -- outlet position and local surface catchment of
# each canyon spring (paper sections 2.9, 4.1; ESM Section S5, Fig. S7).
#
# For springs discharging on the canyon walls, the topographic catchment
# above the outlet is used as a proxy for the local recharge area (an
# approximation: groundwater divides need not follow surface divides).
# Also computed:
#   - height of each spring above the North Fork Virgin River: lowest DEM
#     cell within 50 m of the nearest NHDPlus river point, so a river line
#     drawn slightly onto the canyon wall does not inflate the river
#     elevation;
#   - height above the Navajo-Kayenta contact: contact elevations sampled
#     every 20 m along the mapped contact (shared boundary of the GRI Navajo
#     and Kayenta polygons, as in 09) within 1.5 km of the spring, fitted with
#     a plane z = a + bx + cy and evaluated at the spring (residual SD
#     reported as uncertainty); see contact_plane_fit() in paper/R/functions.R.
# Spring elevations: the 8-m DEM, except at emergence points picked from
# imagery, which use their Google Earth elevation (apply_emergence_elevation()).
# Pour points: the 1988-89 spring coordinates are approximate, so the pour
# point is the highest-flow-accumulation cell within a search radius (50,
# 100, 200 m) whose elevation is at or above the spring's land-surface
# elevation. Water cannot drain uphill to the outlet, and this keeps pour
# points from snapping to the canyon floor or the river.
# DEM: AWS terrain tiles, zoom 14 (3DEP-sourced, resampled to 8 m, UTM 12N);
# see 02_download_external.R. Whitebox Tools intermediates are written to a
# temporary directory.
#
# Outputs
#   derived: spring_outlet_position.csv, spring_catchments.csv,
#            spring_catchment_stability.csv, spring_catchments.gpkg
#   figures: figS7_spring_catchments.png

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(whitebox)
})
sf_use_s2(FALSE)

CANYON_SPRINGS <- c(REPEAT_SPRINGS, "Lower Grotto Sandtrap",
                    "Birch Spring")
RADII <- c(50, 100, 200)       # hard-coded pour-point search radii (m)
LOC_BUFFER_M <- 15             # assumed spring location uncertainty (m)
CONTACT_RADIUS_M <- 1500       # contact points used for the local plane fit
CONTACT_SPACING_M <- 20        # sampling interval along the contact
RIVER_BUFFER_M <- 50           # search radius for the river elevation

# Whitebox works on files; keep intermediates in a temporary directory
wd <- file.path(tempdir(), "wbt16")
dir.create(wd, showWarnings = FALSE)
w <- function(f) file.path(wd, f)


# ---- Flow routing -----------------------------------------------------------

dem_f <- ext_path("dem_canyon_z14.tif")
dem <- rast(dem_f)
names(dem) <- "elev"
# hard-coded: breach depressions up to 200 cells, fill the rest
wbt_breach_depressions_least_cost(dem_f, w("breached.tif"), dist = 200,
                                  fill = TRUE, verbose_mode = FALSE)
wbt_d8_pointer(w("breached.tif"), w("d8.tif"), verbose_mode = FALSE)
wbt_d8_flow_accumulation(w("breached.tif"), w("acc.tif"), out_type = "cells",
                         verbose_mode = FALSE)
acc <- rast(w("acc.tif"))
names(acc) <- "acc"


# ---- Spring elevations ------------------------------------------------------

elev <- read_csv(ext_path("site_elevations_z14.csv"), show_col_types = FALSE)
sp <- site_locations() %>%
  filter(Location %in% CANYON_SPRINGS) %>%
  left_join(elev %>% select(Location, elev_tile_m = elev_m), by = "Location") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(UTM)

# Elevation range of the DEM within the location uncertainty (15 m)
rng <- terra::extract(dem, vect(st_buffer(sp, LOC_BUFFER_M)), fun = range,
                      na.rm = TRUE)

sp <- sp %>%
  # spring elevation on the analysis DEM ...
  mutate(elev_m = terra::extract(dem, vect(sp))$elev) %>%
  # ... replaced by the Google Earth elevation at emergence points picked
  # from imagery (the DEM is unreliable at the foot of cliffs)
  apply_emergence_elevation() %>%
  mutate(elev_min15_m = rng[, 2], elev_max15_m = rng[, 3])


# ---- Outlet position: height above the river --------------------------------

river <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  rename_with(tolower) %>%
  filter(gnis_name == "North Fork Virgin River") %>%
  st_transform(UTM) %>%
  st_union()
# Nearest river point to each spring: st_nearest_points() returns a two-point
# segment (spring, river) per spring; keep every second point
rp <- st_nearest_points(sp, river) %>%
  st_cast("POINT") %>%
  .[seq(2, 2 * nrow(sp), 2)]
river_min <- terra::extract(dem, vect(st_buffer(rp, RIVER_BUFFER_M)),
                            fun = min, na.rm = TRUE)
sp <- sp %>%
  mutate(river_elev_m = river_min[, 2],
         dist_river_m = as.numeric(st_distance(sp, river)))


# ---- Outlet position: height above the Navajo-Kayenta contact ---------------

geo <- st_read(raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb"),
               "zionglg", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(UTM)
# hard-coded GRI map-unit symbols (same groups as 09)
nav <- geo %>%
  filter(GLG_SYM %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl")) %>%
  st_union()
kay <- geo %>%
  filter(GLG_SYM %in% c("Jk", "Jkt", "Jku", "Jks")) %>%
  st_union()
nk <- st_intersection(st_boundary(nav), st_boundary(kay)) %>%
  st_collection_extract("LINESTRING") %>%
  st_union()

# Contact elevation every 20 m along the mapped contact, and a local plane
# fitted to it at each spring
cpts <- contact_points(nk, dem, CONTACT_SPACING_M)
sp <- sp %>%
  bind_cols(contact_plane_fit(sp, cpts, CONTACT_RADIUS_M)) %>%
  mutate(height_above_river_m = elev_m - river_elev_m,
         height_above_contact_m = elev_m - contact_elev_m)


# ---- Catchments for each search radius --------------------------------------

# Pour point for spring `s` (a one-row sf): the highest-accumulation cell
# within `radius` that is not below the spring; if no cell qualifies, the
# highest cell in the search area
pour_point <- function(s, radius) {
  b <- vect(st_buffer(s, radius))
  cells <- terra::extract(acc, b, cells = TRUE, xy = TRUE) %>%
    mutate(z = terra::extract(dem, b)$elev)
  cand <- cells %>% filter(!is.na(z), z >= s$elev_m)
  if (nrow(cand) == 0) cand <- cells %>% slice_max(z, n = 1, with_ties = FALSE)
  best <- cand %>% slice_max(acc, n = 1, with_ties = FALSE)
  st_sf(Location = s$Location, radius_m = radius, acc_cells = best$acc,
        geometry = st_sfc(st_point(c(best$x, best$y)), crs = UTM))
}

# Watershed draining to one pour point, as one polygon
watershed <- function(i, rad) {
  s <- sp %>% slice(i)
  f_pp <- w(sprintf("pp_%d_%d.shp", i, rad))
  f_ws <- w(sprintf("ws_%d_%d.tif", i, rad))
  st_write(pour_point(s, rad), f_pp, quiet = TRUE, delete_dsn = TRUE)
  wbt_watershed(w("d8.tif"), f_pp, f_ws, verbose_mode = FALSE)
  as.polygons(rast(f_ws)) %>%
    st_as_sf() %>%
    st_union() %>%
    st_sf() %>%
    mutate(Location = s$Location, radius_m = rad)
}

# One watershed per spring and radius (run separately so nested catchments
# are not split between pour points)
dem_edge <- st_as_sfc(st_bbox(dem)) %>% st_boundary()
catch <- expand_grid(rad = RADII, i = seq_len(nrow(sp))) %>%
  pmap_dfr(function(rad, i) watershed(i, rad)) %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6,
         touches_dem_edge = as.vector(st_intersects(., dem_edge,
                                                    sparse = FALSE)))

# Geology inside each catchment: cover (Temple Cap/Carmel/Cretaceous) and
# exposed Navajo/Kayenta, as in 10
is_cover <- function(sym) str_detect(sym, "^Jt|^Jc|^K")
is_exposed <- function(sym) {
  sym %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl", "Jk", "Jkt", "Jku", "Jks")
}
park <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  filter(UNIT_CODE == "ZION") %>%
  st_transform(UTM) %>%
  st_union()

# Statistics of one catchment (a one-row sf)
catchment_stats <- function(cj) {
  g <- suppressWarnings(st_intersection(geo %>% select(GLG_SYM), cj)) %>%
    mutate(a = as.numeric(st_area(.)))
  zs <- terra::extract(dem, vect(cj))$elev
  s <- sp %>% filter(Location == cj$Location)
  tibble(
    Location = cj$Location,
    radius_m = cj$radius_m,
    area_km2 = cj$area_km2,
    touches_dem_edge = cj$touches_dem_edge,
    relief_m = max(zs, na.rm = TRUE) - s$elev_m,
    max_elev_m = max(zs, na.rm = TRUE),
    frac_navajo_kayenta_exposed = sum(g$a[is_exposed(g$GLG_SYM)]) / sum(g$a),
    frac_cover = sum(g$a[is_cover(g$GLG_SYM)]) / sum(g$a),
    frac_in_park = as.numeric(st_area(st_intersection(cj, park))) /
      as.numeric(st_area(cj)),
    # farthest point of the catchment boundary from the spring
    max_dist_km = max(as.numeric(
      st_distance(s, st_cast(st_boundary(cj), "POINT"))
    )) / 1000
  )
}
catch_stats <- map_dfr(seq_len(nrow(catch)), ~ catchment_stats(slice(catch, .x)))


# ---- Catchment stability across search radii --------------------------------
# For springs on the canyon walls the outlet either sits on the wall
# (catchment of a few cells) or snaps into the hanging side canyon above it,
# so the area can change by orders of magnitude with the radius. Areas are
# treated as defined only where max/min across radii < 1.5 and the area
# exceeds 0.05 km2 (hard-coded; a stable area of a few cells on a cliff face
# is not a catchment).

stab <- catch_stats %>%
  select(Location, radius_m, area_km2) %>%
  pivot_wider(names_from = radius_m, values_from = area_km2,
              names_prefix = "area_") %>%
  mutate(ratio_max_min = pmax(area_50, area_100, area_200) /
           pmin(area_50, area_100, area_200),
         catchment_defined = ratio_max_min < 1.5 & area_100 > 0.05) %>%
  arrange(Location)
write_csv(stab, file.path(DER, "spring_catchment_stability.csv"))
print(as.data.frame(stab), digits = 3)

outlet <- sp %>%
  st_drop_geometry() %>%
  select(Location, elev_m, elev_tile_m, elev_min15_m, elev_max15_m,
         river_elev_m, dist_river_m, height_above_river_m, contact_elev_m,
         contact_sd_m, contact_n, height_above_contact_m)
write_csv(outlet, file.path(DER, "spring_outlet_position.csv"))
write_csv(catch_stats, file.path(DER, "spring_catchments.csv"))
st_write(catch, file.path(DER, "spring_catchments.gpkg"), quiet = TRUE,
         delete_dsn = TRUE)

outlet %>%
  transmute(Location,
            elev_m = round(elev_m),
            above_river = round(height_above_river_m),
            above_contact = round(height_above_contact_m),
            contact_sd = round(contact_sd_m),
            elev_range15 = str_c(round(elev_min15_m), "-", round(elev_max15_m))) %>%
  as.data.frame() %>%
  print()
catch_stats %>%
  select(Location, radius_m, area_km2, relief_m, max_dist_km,
         frac_navajo_kayenta_exposed, frac_cover, frac_in_park,
         touches_dem_edge) %>%
  as.data.frame() %>%
  print(digits = 3)


# ---- ESM map (Fig. S7) ------------------------------------------------------

# Hillshade (sun 40 deg above the horizon, from the northwest)
hs <- shade(terrain(dem, "slope", unit = "radians"),
            terrain(dem, "aspect", unit = "radians"), 40, 315)
c100 <- catch %>% filter(radius_m == 100)
bb <- st_bbox(st_buffer(st_union(c100), 1500))
hdf <- crop(hs, ext(bb[c(1, 3, 2, 4)])) %>%
  as.data.frame(xy = TRUE) %>%
  rename(hs = 3)

p <- ggplot() +
  geom_raster(data = hdf, aes(x, y, fill = hs), show.legend = FALSE) +
  scale_fill_gradient(low = "grey25", high = "white") +
  geom_sf(data = st_crop(st_sf(geometry = nk), bb), color = "#cc3030",
          linewidth = 0.4) +
  geom_sf(data = c100, aes(color = Location), fill = NA, linewidth = 0.6) +
  geom_sf(data = sp, shape = 21, fill = "white", color = "black", size = 1.8) +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], datum = UTM, expand = FALSE) +
  labs(x = NULL, y = NULL, color = NULL,
       caption = str_c("Surface catchments above each spring (100-m ",
                       "pour-point search); red: Navajo-Kayenta contact")) +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(legend.position = "right")
ggsave(file.path(FIG, "figS7_spring_catchments.png"), p, width = 7.2,
       height = 6, dpi = 300, bg = "white")


# ---- Checks -----------------------------------------------------------------

message(sprintf("[%s] no catchment touches the DEM edge",
                if_else(!any(catch$touches_dem_edge), "OK", "DIFF")))
check_close(outlet$height_above_river_m > 0, rep(TRUE, nrow(outlet)), 0,
            "every canyon spring discharges above the river")
wr <- outlet %>% filter(Location == "Weeping Rock")
check_close(wr$height_above_contact_m, 0, 2 * wr$contact_sd_m,
            "Weeping Rock discharges at the Navajo-Kayenta contact (within 2 SD of the fitted contact)")
