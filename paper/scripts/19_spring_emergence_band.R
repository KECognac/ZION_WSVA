# 19_spring_emergence_band.R -- elevation band of spring emergence in Zion
# Canyon, from emergence points picked on Google Earth imagery (K. Cognac,
# November 2026; KML copied by 01_copy_raw_data.R).
#
# For each emergence point:
#   - land-surface elevation on the 8-m DEM (and the range within 15 m),
#     next to the Google Earth elevation stored in the KML;
#   - mapped geologic unit (GRI) and planar distance to the mapped
#     Navajo-Kayenta contact;
#   - height above the Navajo-Kayenta contact: plane fitted to contact
#     elevations within 1.5 km, as in 16_spring_catchments.R;
#   - height above the North Fork Virgin River (lowest DEM cell within 50 m
#     of the nearest river point), as in 16.
# Points are matched to the sampled sites of the paper by distance.
# Heights use the Google Earth elevation of each point (the 8-m DEM can be off
# by more than 100 m at the foot of cliffs).
#
# Outputs
#   derived: spring_emergence_points.csv
#   figures: fig7_flow_emergence.png (Fig. 7; panels a and b come from 11)
# Runs after 11.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
})
sf_use_s2(FALSE)

LOC_BUFFER_M <- 15             # assumed location uncertainty (m), as in 16
CONTACT_RADIUS_M <- 1500       # contact points used for the local plane fit
CONTACT_SPACING_M <- 20        # sampling interval along the contact
RIVER_BUFFER_M <- 50           # search radius for the river elevation
MATCH_M <- 150                 # an emergence point within this distance of a
                               # sampled site is treated as that site

dem <- rast(ext_path("dem_canyon_z14.tif"))
names(dem) <- "elev"


# ---- Emergence points -------------------------------------------------------

kml <- st_read(raw_path("Spring Emergence Points from google earth.kml"),
               quiet = TRUE)
xyz <- st_coordinates(kml)
em <- kml %>%
  st_zm() %>%
  transmute(id = row_number(), label = Name, longitude = xyz[, 1],
            latitude = xyz[, 2], elev_ge_m = xyz[, 3]) %>%
  st_transform(UTM)
# Placemarks closer than 5 m are the same emergence labelled twice; keep one
# point, with the specific name where one of the pair has it, and the mean of
# the two Google Earth elevations. grp_id = first placemark within 5 m.
dup <- st_distance(em) %>% units::drop_units()
grp_id <- map_int(seq_len(nrow(em)), ~ min(which(dup[.x, ] < 5)))
em <- em %>%
  mutate(generic = str_detect(label, regex("emergence", ignore_case = TRUE)),
         grp_id = grp_id) %>%
  group_by(grp_id) %>%
  mutate(elev_ge_m = mean(elev_ge_m), n_placemarks = n()) %>%
  arrange(generic, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp_id)
message(nrow(kml) - nrow(em), " duplicate placemark(s) merged")

# Nearest sampled site of the paper
sites <- site_locations() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM)
nearest <- sites[st_nearest_feature(em, sites), ]
rng <- terra::extract(dem, vect(st_buffer(em, LOC_BUFFER_M)), fun = range,
                      na.rm = TRUE)
em <- em %>%
  mutate(
    nearest_site = nearest$Location,
    dist_site_m = as.numeric(st_distance(em, nearest, by_element = TRUE)),
    site = if_else(dist_site_m <= MATCH_M, nearest_site, NA_character_),
    # DEM elevation at the point and its range within 15 m
    elev_m = terra::extract(dem, vect(em))$elev,
    elev_min15_m = rng[, 2],
    elev_max15_m = rng[, 3]
  )
stopifnot(!anyNA(em$elev_m))


# ---- Geology and the Navajo-Kayenta contact ---------------------------------

geo <- st_read(raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb"),
               "zionglg", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(UTM)
em <- st_join(em, geo %>% select(GLG_SYM), join = st_intersects,
              largest = TRUE)
# hard-coded GRI map-unit symbols (same groups as 09 and 16)
nav <- geo %>%
  filter(GLG_SYM %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl")) %>%
  st_union()
kay <- geo %>%
  filter(GLG_SYM %in% c("Jk", "Jkt", "Jku", "Jks")) %>%
  st_union()
nk <- st_intersection(st_boundary(nav), st_boundary(kay)) %>%
  st_collection_extract("LINESTRING") %>%
  st_union()
em <- em %>% mutate(dist_contact_m = as.numeric(st_distance(em, nk)))

# Local contact elevation at each point (plane fit, as in 16)
cpts <- contact_points(nk, dem, CONTACT_SPACING_M)
em <- em %>% bind_cols(contact_plane_fit(em, cpts, CONTACT_RADIUS_M))


# ---- Height above the river -------------------------------------------------

river <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  rename_with(tolower) %>%
  filter(gnis_name == "North Fork Virgin River") %>%
  st_transform(UTM) %>%
  st_union()
# Nearest river point to each emergence point (second point of each
# connecting segment) and the lowest DEM cell within 50 m of it
rp <- st_nearest_points(em, river) %>%
  st_cast("POINT") %>%
  .[seq(2, 2 * nrow(em), 2)]
river_min <- terra::extract(dem, vect(st_buffer(rp, RIVER_BUFFER_M)),
                            fun = min, na.rm = TRUE)
em <- em %>%
  mutate(river_elev_m = river_min[, 2],
         dist_river_m = as.numeric(st_distance(em, river)))

out <- em %>%
  mutate(height_above_contact_m = elev_ge_m - contact_elev_m,
         height_above_river_m = elev_ge_m - river_elev_m,
         utm_x = st_coordinates(.)[, 1], utm_y = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()
write_csv(out, file.path(DER, "spring_emergence_points.csv"))
print(as.data.frame(out %>%
  select(id, label, site, dist_site_m, elev_ge_m, elev_m, elev_min15_m,
         elev_max15_m, GLG_SYM, dist_contact_m, contact_elev_m, contact_sd_m,
         height_above_contact_m, height_above_river_m)), digits = 4)


# ---- Fig. 7: flow map, canyon map with emergence points, elevation profile --

suppressPackageStartupMessages(library(ggnewscale))
pan <- readRDS(file.path(DER, "fig_flow_panels.rds"))   # from 11

# Canyon-wall points only in the profile; the "Upper Spring Emergence" points
# lie about 600 m higher (hard-coded cut at 1,800 m)
out <- out %>%
  mutate(level = if_else(elev_ge_m >= 1800, "Upper emergence point",
                         "Canyon-wall emergence point"),
         sampled = !is.na(site))
check_close(sum(out$level == "Upper emergence point"), 3, 0,
            "three upper emergence points")
em_sf <- st_as_sf(out, coords = c("utm_x", "utm_y"), crs = UTM)

# (b) canyon map: add the emergence points that are not sampled sites, and
# widen the frame to hold them
# hard-coded: the map frame stops at 37.30 N; the two "likely" points farther
# north appear in panel c only
other <- em_sf %>% filter(!sampled, latitude < 37.30)
site_pts <- pan$pB$layers[[6]]$data
bb <- st_bbox(st_buffer(st_union(st_geometry(other), st_geometry(site_pts)), 900))
# site labels are placed again for the new frame
is_lab <- map_lgl(pan$pB$layers, ~ inherits(.x$geom, "GeomLabelRepel"))
stopifnot(sum(is_lab) == 1)
lab_df <- pan$pB$layers[[which(is_lab)]]$data
pB <- pan$pB
pB$layers <- pB$layers[!is_lab]
pB <- suppressMessages(
  pB +
    geom_sf(data = other, aes(color = level), shape = 18, size = 2.3) +
    ggrepel::geom_label_repel(
      data = lab_df, aes(x, y, label = lab), size = 2.2, family = "Arial",
      label.size = 0, fill = alpha("white", 0.85), label.padding = 0.08,
      lineheight = 0.9, min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.7, point.padding = 0.3, force = 8, max.overlaps = Inf,
      seed = 11
    ) +
    scale_color_manual(values = c("Canyon-wall emergence point" = "black",
                                  "Upper emergence point" = "#7a1fa2"),
                       name = "Imagery (not sampled)",
                       guide = guide_legend(order = 4)) +
    coord_sf(crs = UTM, xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE,
             datum = 4326) +
    labs(title = "b  Zion Canyon: young-water fraction, geology and emergence points")
)

# (c) elevation of each canyon-wall emergence point, south to north, with the
# fitted Navajo-Kayenta contact (+/- 1 SD) and the river
band <- out %>%
  filter(level == "Canyon-wall emergence point") %>%
  arrange(utm_y) %>%
  mutate(
    name = case_when(
      !is.na(site) ~ str_remove(site, " Sandtrap"),
      !generic ~ label,
      str_detect(label, regex("likely", ignore_case = TRUE)) ~ "Unnamed (likely)",
      TRUE ~ "Unnamed"
    ),
    x = row_number()
  )
k1 <- "Emergence point (Google Earth elevation)"
k2 <- "Navajo–Kayenta contact (fitted, ± 1 SD)"
k3 <- "North Fork Virgin River"
pC <- ggplot(band, aes(x = x)) +
  geom_linerange(aes(ymin = contact_elev_m - contact_sd_m,
                     ymax = contact_elev_m + contact_sd_m),
                 color = "#c9a27a", linewidth = 5, alpha = 0.6, na.rm = TRUE) +
  geom_point(aes(y = contact_elev_m, shape = k2), color = "#8a5a2b", size = 5,
             na.rm = TRUE) +
  geom_point(aes(y = river_elev_m, shape = k3), color = "#3a8fa0", size = 5) +
  geom_linerange(aes(ymin = elev_min15_m, ymax = elev_max15_m),
                 color = "grey45", linewidth = 0.35) +
  geom_point(aes(y = elev_ge_m, fill = sampled, shape = k1), size = 2.2,
             stroke = 0.5) +
  scale_shape_manual(values = setNames(c(21, 95, 95), c(k1, k2, k3)),
                     name = NULL) +
  scale_fill_manual(values = c(`TRUE` = "#2060cc", `FALSE` = "white"),
                    labels = c(`TRUE` = "Sampled site", `FALSE` = "Not sampled"),
                    name = NULL) +
  # legend order: sampled / not sampled first, then the symbols
  guides(fill = guide_legend(order = 1,
                             override.aes = list(shape = 21, size = 2.2)),
         shape = guide_legend(order = 2, ncol = 2, override.aes = list(
           size = c(2.2, 5, 5), color = c("black", "#8a5a2b", "#3a8fa0")))) +
  scale_x_continuous(breaks = band$x, labels = band$name) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Canyon-wall emergence points, south to north", y = "Elevation (m)",
       title = "c  Emergence elevation relative to the base of the Navajo Sandstone") +
  theme_paper(8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(), legend.position = "top",
        legend.box = "vertical", legend.box.just = "left",
        legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0, "cm"), plot.title = element_text(size = 10))

lr <- theme(legend.position = "right", legend.box = "vertical",
            legend.spacing.y = unit(0.1, "cm"))
fig <- (pan$pA + lr) / (pB + lr) / pC + plot_layout(heights = c(1, 1.05, 0.6))
ggsave(file.path(FIG, "fig7_flow_emergence.png"), fig, width = 7.2,
       height = 11, dpi = 400, bg = "white")
message("Fig. 7 written to ", FIG)

# Summary used in the text
b9 <- band %>% filter(!is.na(contact_elev_m))
message(sprintf(
  "canyon-wall points: %d, %0.f-%0.f m (GE), %0.f-%0.f m above the river; %d with a contact fit: %0.f to %+0.f m relative to the contact",
  nrow(band), min(band$elev_ge_m), max(band$elev_ge_m),
  min(band$elev_ge_m - band$river_elev_m), max(band$elev_ge_m - band$river_elev_m),
  nrow(b9), min(b9$elev_ge_m - b9$contact_elev_m),
  max(b9$elev_ge_m - b9$contact_elev_m)))
print(as.data.frame(band %>% transmute(name, ge = round(elev_ge_m), unit = GLG_SYM,
  h_contact = round(elev_ge_m - contact_elev_m), sd = round(contact_sd_m),
  h_river = round(elev_ge_m - river_elev_m), dem_lo = round(elev_min15_m),
  dem_hi = round(elev_max15_m))))

