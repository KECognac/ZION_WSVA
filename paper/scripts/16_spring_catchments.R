# 16_spring_catchments.R -- outlet position and local surface catchment of
# each canyon spring.
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
#     reported as uncertainty).
# Pour points: the 1988-89 spring coordinates are approximate, so the pour
# point is the highest-flow-accumulation cell within a search radius (50,
# 100, 200 m) whose elevation is at or above the spring's land-surface
# elevation. Water cannot drain uphill to the outlet, and this keeps pour
# points from snapping to the canyon floor or the river.
# DEM: AWS terrain tiles, zoom 14 (3DEP-sourced, resampled to 8 m, UTM 12N);
# see 02_download_external.R. Whitebox Tools intermediates are written to a
# temporary directory.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(whitebox)
})
sf_use_s2(FALSE)

CANYON_SPRINGS <- c(REPEAT_SPRINGS, "Lower Grotto Sandtrap",
                    "Birch Sediment Trap")
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


# ---- Spring elevations ------------------------------------------------------

elev <- read_csv(ext_path("site_elevations_z14.csv"), show_col_types = FALSE)
sp <- site_locations() %>%
  filter(Location %in% CANYON_SPRINGS) %>%
  left_join(elev %>% select(Location, elev_tile_m = elev_m), by = "Location") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(UTM)

# Spring elevation on the analysis DEM, and the elevation range within the
# location uncertainty
sp$elev_m <- terra::extract(dem, vect(sp))$elev
rng <- terra::extract(dem, vect(st_buffer(sp, LOC_BUFFER_M)), fun = range,
                      na.rm = TRUE)
sp$elev_min15_m <- rng[, 2]
sp$elev_max15_m <- rng[, 3]


# ---- Outlet position: height above the river --------------------------------

fl <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE)
names(fl) <- tolower(names(fl))
river <- fl %>%
  filter(gnis_name == "North Fork Virgin River") %>%
  st_transform(UTM) %>%
  st_union()
# nearest river point to each spring (every second point of the connecting
# segments)
rp <- st_nearest_points(sp, river) %>%
  st_cast("POINT") %>%
  .[seq(2, 2 * nrow(sp), 2)]
sp$river_elev_m <- terra::extract(dem, vect(st_buffer(rp, RIVER_BUFFER_M)),
                                  fun = min, na.rm = TRUE)[, 2]
sp$dist_river_m <- as.numeric(st_distance(sp, river))


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

# Contact elevation every 20 m along the mapped contact
cpts <- st_line_sample(st_cast(nk, "LINESTRING"),
                       density = 1 / CONTACT_SPACING_M) %>%
  st_cast("POINT") %>%
  st_sf()
cxy <- st_coordinates(cpts)
cpts$z <- terra::extract(dem, vect(cpts))$elev

# Plane fit to the contact points near each spring, evaluated at the spring
# (coordinates centered on the spring, so the intercept is the estimate).
# Needs at least 10 points.
contact_fit <- map_dfr(seq_len(nrow(sp)), function(i) {
  p <- st_coordinates(sp[i, ])
  sel <- sqrt((cxy[, 1] - p[1])^2 + (cxy[, 2] - p[2])^2) <= CONTACT_RADIUS_M &
    !is.na(cpts$z)
  dfc <- tibble(z = cpts$z[sel], x = cxy[sel, 1] - p[1], y = cxy[sel, 2] - p[2])
  if (nrow(dfc) < 10) {
    return(tibble(contact_elev_m = NA_real_, contact_sd_m = NA_real_,
                  contact_n = nrow(dfc)))
  }
  m <- lm(z ~ x + y, data = dfc)
  tibble(contact_elev_m = unname(coef(m)[1]), contact_sd_m = sigma(m),
         contact_n = nrow(dfc))
})

sp <- bind_cols(sp, contact_fit) %>%
  mutate(height_above_river_m = elev_m - river_elev_m,
         height_above_contact_m = elev_m - contact_elev_m)


# ---- Catchments for each search radius --------------------------------------

# Pour point: highest-accumulation cell within `radius` that is not below the
# spring; if no cell qualifies, the highest cell in the search area
pour_point <- function(i, radius) {
  b <- vect(st_buffer(sp[i, ], radius))
  a <- terra::extract(acc, b, cells = TRUE, xy = TRUE)
  z <- terra::extract(dem, b)$elev
  cand <- a[!is.na(z) & z >= sp$elev_m[i], ]
  if (!nrow(cand)) cand <- a[which.max(z), ]
  cand <- cand[which.max(cand[[2]]), ]
  st_sf(Location = sp$Location[i], radius_m = radius, acc_cells = cand[[2]],
        geometry = st_sfc(st_point(c(cand$x, cand$y)), crs = UTM))
}

ext_poly <- st_as_sfc(st_bbox(dem)) %>% st_boundary()

# One watershed per spring and radius (run separately so nested catchments
# are not split between pour points)
catch <- map_dfr(RADII, function(rad) {
  map_dfr(seq_len(nrow(sp)), function(i) {
    pp <- pour_point(i, rad)
    f_pp <- w(sprintf("pp_%d_%d.shp", i, rad))
    f_ws <- w(sprintf("ws_%d_%d.tif", i, rad))
    st_write(pp, f_pp, quiet = TRUE, delete_dsn = TRUE)
    wbt_watershed(w("d8.tif"), f_pp, f_ws, verbose_mode = FALSE)
    as.polygons(rast(f_ws)) %>%
      st_as_sf() %>%
      st_union() %>%
      st_sf() %>%
      mutate(Location = sp$Location[i], radius_m = rad)
  })
})
catch <- catch %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6,
         touches_dem_edge = as.vector(st_intersects(., ext_poly,
                                                    sparse = FALSE)))

# Geology inside each catchment: cover (Temple Cap/Carmel/Cretaceous) and
# exposed Navajo/Kayenta, as in 10
cover_syms <- function(s) grepl("^Jt|^Jc|^K", s)
exposed_syms <- function(s) {
  s %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl", "Jk", "Jkt", "Jku", "Jks")
}
park <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  filter(UNIT_CODE == "ZION") %>%
  st_transform(UTM) %>%
  st_union()

catch_stats <- map_dfr(seq_len(nrow(catch)), function(j) {
  cj <- catch[j, ]
  g <- suppressWarnings(st_intersection(geo %>% select(GLG_SYM), cj))
  ga <- as.numeric(st_area(g))
  A <- sum(ga)
  zs <- terra::extract(dem, vect(cj))$elev
  s <- sp[sp$Location == cj$Location, ]
  tibble(
    Location = cj$Location,
    radius_m = cj$radius_m,
    area_km2 = cj$area_km2,
    touches_dem_edge = cj$touches_dem_edge,
    relief_m = max(zs, na.rm = TRUE) - s$elev_m,
    max_elev_m = max(zs, na.rm = TRUE),
    frac_navajo_kayenta_exposed = sum(ga[exposed_syms(g$GLG_SYM)]) / A,
    frac_cover = sum(ga[cover_syms(g$GLG_SYM)]) / A,
    frac_in_park = as.numeric(st_area(st_intersection(cj, park))) /
      as.numeric(st_area(cj)),
    max_dist_km = max(as.numeric(
      st_distance(s, st_cast(st_boundary(cj), "POINT"))
    )) / 1000
  )
})


# ---- Catchment stability across search radii --------------------------------
# For springs on the canyon walls the outlet either sits on the wall
# (catchment of a few cells) or snaps into the hanging side canyon above it,
# so the area can change by orders of magnitude with the radius. Areas are
# treated as defined only where max/min across radii < 1.5 and the area
# exceeds 0.05 km2 (hard-coded; a stable area of a few cells on a cliff face
# is not a catchment).

stab <- catch_stats %>%
  group_by(Location) %>%
  summarise(area_50 = area_km2[radius_m == 50],
            area_100 = area_km2[radius_m == 100],
            area_200 = area_km2[radius_m == 200],
            ratio_max_min = max(area_km2) / min(area_km2),
            .groups = "drop") %>%
  mutate(catchment_defined = ratio_max_min < 1.5 & area_100 > 0.05)
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
            elev_range15 = paste0(round(elev_min15_m), "-",
                                  round(elev_max15_m))) %>%
  as.data.frame() %>%
  print()
catch_stats %>%
  select(Location, radius_m, area_km2, relief_m, max_dist_km,
         frac_navajo_kayenta_exposed, frac_cover, frac_in_park,
         touches_dem_edge) %>%
  as.data.frame() %>%
  print(digits = 3)


# ---- ESM map ----------------------------------------------------------------

hs <- shade(terrain(dem, "slope", unit = "radians"),
            terrain(dem, "aspect", unit = "radians"), 40, 315)
c100 <- catch %>% filter(radius_m == 100)
bb <- st_bbox(st_buffer(st_union(c100), 1500))
hs_c <- crop(hs, ext(bb[c(1, 3, 2, 4)]))
hdf <- as.data.frame(hs_c, xy = TRUE)
names(hdf)[3] <- "hs"

p <- ggplot() +
  geom_raster(data = hdf, aes(x, y, fill = hs), show.legend = FALSE) +
  scale_fill_gradient(low = "grey25", high = "white") +
  geom_sf(data = st_crop(st_sf(geometry = nk), bb), color = "#cc3030",
          linewidth = 0.4) +
  geom_sf(data = c100, aes(color = Location), fill = NA, linewidth = 0.6) +
  geom_sf(data = sp, shape = 21, fill = "white", color = "black", size = 1.8) +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], datum = UTM, expand = FALSE) +
  labs(x = NULL, y = NULL, color = NULL,
       caption = paste0("Surface catchments above each spring (100-m ",
                        "pour-point search); red: Navajo-Kayenta contact")) +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(legend.position = "right")
ggsave(file.path(FIG, "figS7_spring_catchments.png"), p, width = 7.2,
       height = 6, dpi = 300, bg = "white")


# ---- Checks -----------------------------------------------------------------

message(sprintf("[%s] no catchment touches the DEM edge",
                if (!any(catch$touches_dem_edge)) "OK" else "DIFF"))
check_close(outlet$height_above_river_m > 0, rep(TRUE, nrow(outlet)), 0,
            "every canyon spring discharges above the river")
wr <- outlet %>% filter(Location == "Weeping Rock")
check_close(wr$height_above_contact_m, 0, 2 * wr$contact_sd_m,
            "Weeping Rock discharges at the Navajo-Kayenta contact (within 2 SD of the fitted contact)")
