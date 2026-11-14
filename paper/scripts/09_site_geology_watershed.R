# 09_site_geology_watershed.R -- watershed, geology and site setting (paper
# section 2.1; Table 1; Fig. 1).
#
# Geology: NPS Geologic Resources Inventory digital geologic map (layers
# zionglg, zionflt, zionjln, zionunit) grouped into formations.
# Watershed: NLDI basin of USGS streamgage 09405501 and WBD hydrologic units.
# Streams: NHDPlus.
# Also writes each site's side of the North Fork Virgin River and its
# land-surface elevation (used by Fig. 5).

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(nhdplusTools)
})
sf_use_s2(FALSE)


# ---- Geology (same GDB and formation groups as the Rmd, Section 4) ----------

gdb <- raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb")
geo <- st_read(gdb, "zionglg", quiet = TRUE) %>% st_make_valid()
flt <- st_read(gdb, "zionflt", quiet = TRUE) %>%
  mutate(FNAME = ifelse(is.na(FNAME) | FNAME %in% c("", "NA"),
                        "(unnamed)", FNAME))
jln <- st_read(gdb, "zionjln", quiet = TRUE)
units <- st_read(gdb, "zionunit", quiet = TRUE) %>% st_drop_geometry()

# hard-coded: GRI map-unit symbols (GLG_SYM) in each formation group
grp <- list(
  "Carmel Formation"     = c("Jc", "Jcw", "Jcp", "Jcx", "Jcc", "Jccu", "Jccl"),
  "Temple Cap Formation" = c("Jt", "Jtw", "Jts"),
  "Navajo Sandstone"     = c("Jn", "Jnw", "Jnp", "Jnb", "Jnl"),
  "Kayenta Formation"    = c("Jk", "Jkt", "Jku", "Jks")
)
geo <- geo %>%
  mutate(group = case_when(
    GLG_SYM %in% grp[[1]] ~ names(grp)[1],
    GLG_SYM %in% grp[[2]] ~ names(grp)[2],
    GLG_SYM %in% grp[[3]] ~ names(grp)[3],
    GLG_SYM %in% grp[[4]] ~ names(grp)[4],
    GLG_SYM %in% c("JTRkm", "JTRm", "Jmw", "JTRmd") |
      grepl("^TR|^P|^PZ", GLG_SYM) ~ "Moenave and older",
    grepl("^K", GLG_SYM) ~ "Cretaceous formations",
    grepl("^Qb", GLG_SYM) & !GLG_SYM %in% c("Qblcp") ~ "Quaternary basalt",
    grepl("^Q", GLG_SYM) ~ "Quaternary deposits",
    TRUE ~ "Other"
  ))
geo_g <- geo %>%
  group_by(group) %>%
  summarise(.groups = "drop") %>%
  st_make_valid()
map_extent <- st_union(geo) %>% st_make_valid()


# ---- Sites (tritium sample locations) ---------------------------------------

sites <- site_locations() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(UTM)


# ---- Watersheds, streams and park -------------------------------------------

gage <- read_csv(ext_path("nwis_gage_09405501.csv"), show_col_types = FALSE)
basin <- st_read(ext_path("gage_09405501_basin.gpkg"), quiet = TRUE) %>%
  st_transform(UTM)
gage_pt <- st_as_sf(gage, coords = c("dec_long_va", "dec_lat_va"),
                    crs = 4269) %>%
  st_transform(UTM)
huc <- st_read(ext_path("huc_units.gpkg"), quiet = TRUE) %>% st_transform(UTM)

site_huc <- st_join(sites %>% select(site = Location),
                    huc %>% transmute(type = level, id = huc, name, areasqkm),
                    join = st_within) %>%
  st_drop_geometry() %>%
  arrange(site, type)
write_csv(site_huc, file.path(DER, "site_huc.csv"))

flow_all <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  st_transform(UTM)
flow <- flow_all[st_intersects(flow_all, basin, sparse = FALSE)[, 1], ]

park <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  filter(UNIT_CODE == "ZION") %>%
  st_make_valid() %>%
  st_transform(UTM) %>%
  summarise()


# ---- 1. Watershed statistics ------------------------------------------------

a_km2 <- function(x) as.numeric(sum(st_area(x))) / 1e6

# Area (km2 and % of mapped area) of each formation group within `aoi`
group_area <- function(aoi) {
  st_intersection(geo_g, aoi) %>%
    mutate(km2 = as.numeric(st_area(.)) / 1e6) %>%
    st_drop_geometry() %>%
    group_by(group) %>%
    summarise(km2 = sum(km2), .groups = "drop") %>%
    mutate(pct_of_mapped = 100 * km2 / sum(km2))
}

basin_km2 <- a_km2(basin)
basin_mapped <- st_intersection(basin, map_extent)
geo_in_basin <- group_area(basin)
park_km2 <- a_km2(park)
geo_in_park <- group_area(park)

# Fault and joint lengths
flt_b <- st_intersection(flt, basin)
jln_b <- st_intersection(jln, basin)
fault_stats <- st_drop_geometry(flt_b) %>%
  mutate(len_km = as.numeric(st_length(flt_b)) / 1000) %>%
  group_by(FNAME) %>%
  summarise(len_km = sum(len_km), n = n(), .groups = "drop") %>%
  arrange(desc(len_km))
flt_park <- st_intersection(flt, park)
fault_park_named <- st_drop_geometry(flt_park) %>%
  mutate(len_km = as.numeric(st_length(flt_park)) / 1000) %>%
  group_by(FNAME) %>%
  summarise(len_km = sum(len_km), .groups = "drop") %>%
  arrange(desc(len_km))

joint_km_basin <- as.numeric(sum(st_length(jln_b))) / 1000
ws <- list(
  gage_name = gage$station_nm,
  gage_drain_area_mi2 = gage$drain_area_va,
  gage_huc = gage$huc_cd,
  basin_km2 = basin_km2,
  basin_mapped_km2 = a_km2(basin_mapped),
  basin_mapped_pct = 100 * a_km2(basin_mapped) / basin_km2,
  park_km2 = park_km2,
  park_in_basin_km2 = a_km2(st_intersection(park, basin)),
  fault_km_basin = sum(fault_stats$len_km),
  joint_km_basin = joint_km_basin,
  fault_density_km_per_km2_mapped = sum(fault_stats$len_km) /
    a_km2(basin_mapped),
  joint_density_km_per_km2_mapped = joint_km_basin / a_km2(basin_mapped)
)


# ---- 2. Site-level geology --------------------------------------------------

navajo <- geo_g %>% filter(group == "Navajo Sandstone")
kay <- geo_g %>% filter(group == "Kayenta Formation")
# Contacts = shared boundaries of the mapped formation polygons
nk_contact <- st_intersection(st_boundary(navajo), st_boundary(kay)) %>%
  st_collection_extract("LINESTRING")
nav_carm_tc <- st_intersection(
  st_boundary(navajo),
  st_boundary(geo_g %>% filter(group == "Temple Cap Formation"))
) %>%
  st_collection_extract("LINESTRING")

# Distance (m) from each point to the nearest line
nearest_d <- function(p, lines) as.numeric(apply(st_distance(p, lines), 1, min))

flt_named <- flt %>% filter(FNAME != "(unnamed)")
nf <- st_nearest_feature(sites, flt)

site_tab <- sites %>%
  st_join(geo %>% select(GLG_SYM, group), join = st_within) %>%
  left_join(units %>% select(GLG_SYM, GLG_NAME), by = "GLG_SYM") %>%
  mutate(
    in_gage_basin = as.vector(st_within(., basin, sparse = FALSE)),
    in_park = as.vector(st_within(., park, sparse = FALSE)),
    dist_fault_m = nearest_d(., flt),
    nearest_fault = flt$FNAME[nf],
    dist_named_fault_m = nearest_d(., flt_named),
    nearest_named_fault = flt_named %>%
      slice(st_nearest_feature(sites, .)) %>%
      pull(FNAME),
    dist_joint_m = nearest_d(., jln),
    dist_navajo_kayenta_contact_m = nearest_d(., nk_contact),
    dist_navajo_templecap_contact_m = nearest_d(., nav_carm_tc)
  ) %>%
  left_join(site_huc %>%
              filter(type == "huc12") %>%
              select(Location = site, huc12 = id, huc12_name = name),
            by = "Location") %>%
  left_join(site_huc %>%
              filter(type == "huc10") %>%
              select(Location = site, huc10_name = name),
            by = "Location") %>%
  st_drop_geometry()

write.csv(site_tab, file.path(DER, "site_geology.csv"), row.names = FALSE)
write.csv(site_tab %>%
            select(Location, GLG_NAME, dist_navajo_kayenta_contact_m,
                   dist_fault_m, nearest_fault, dist_joint_m),
          file.path(TAB, "table1_site_geology.csv"), row.names = FALSE)
write.csv(geo_in_basin, file.path(DER, "basin_geology_area.csv"),
          row.names = FALSE)
write.csv(geo_in_park, file.path(DER, "park_geology_area.csv"),
          row.names = FALSE)
write.csv(fault_stats, file.path(DER, "basin_faults.csv"), row.names = FALSE)
write.csv(fault_park_named, file.path(DER, "park_faults.csv"),
          row.names = FALSE)
jsonlite::write_json(ws, file.path(DER, "watershed_stats.json"),
                     auto_unbox = TRUE, digits = 6, pretty = TRUE)


# ---- 3. Site figure (Fig. 1) ------------------------------------------------

# hard-coded figure colors for the formation groups
cols <- c(
  "Cretaceous formations" = "#a9b38a", "Carmel Formation" = "#7f9fbf",
  "Temple Cap Formation"  = "#b9a3c9", "Navajo Sandstone" = "#efc9a4",
  "Kayenta Formation"     = "#b5654d", "Moenave and older" = "#8c6d5a",
  "Quaternary basalt"     = "#4d4d4d", "Quaternary deposits" = "#e8e3c9",
  "Other mapped units"    = "#cccccc"
)
geo_plot <- geo_g %>%
  mutate(group = recode(group, "Other" = "Other mapped units"),
         group = factor(group, levels = names(cols)))
sites_p <- site_tab %>%
  select(Location, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM) %>%
  mutate(kolob = Location == "Taylor Creek Well")
huc10 <- huc %>% filter(level == "huc10")
unmapped <- st_difference(basin, map_extent)

# Line features share one color/linetype legend
lcols <- c(
  "Gage basin (USGS 09405501)" = "#2060cc", "HUC-10 boundary" = "#2060cc",
  "Park boundary" = "#208040", "Fault" = "black", "Joint" = "grey45",
  "Navajo-Kayenta contact" = "#cc3030", "Stream (NHDPlus)" = "#3a8fa0"
)
lty <- c(
  "Gage basin (USGS 09405501)" = "solid", "HUC-10 boundary" = "dashed",
  "Park boundary" = "solid", "Fault" = "solid", "Joint" = "solid",
  "Navajo-Kayenta contact" = "solid", "Stream (NHDPlus)" = "solid"
)
line_scales <- list(
  scale_color_manual(values = lcols, name = NULL, breaks = names(lcols),
                     limits = names(lcols), drop = FALSE),
  scale_linetype_manual(values = lty, name = NULL, breaks = names(lcols),
                        limits = names(lcols), drop = FALSE)
)
# Off-map dummy lines that only produce the legend keys
legend_df <- data.frame(what = rep(names(lcols), each = 2),
                        x = rep(c(0, 1), 7), y = 0)

theme_map <- theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 10),
    legend.key.size = unit(0.38, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8.5, face = "bold")
  )
lab_pt <- function(x) {
  x %>% mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])
}
label_style <- list(family = "Arial", size = 2.6, hjust = 0, label.size = 0,
                    fill = alpha("white", 0.85),
                    label.padding = unit(0.06, "lines"))

# (a) North Fork Virgin River basin
bb_a <- st_bbox(st_union(st_buffer(basin, 1500), st_buffer(sites_p, 3000)))
gage_xy <- st_coordinates(gage_pt)
kolob_xy <- st_coordinates(sites_p %>% filter(kolob))

pA <- ggplot() +
  geom_blank() +
  geom_sf(data = geo_plot, aes(fill = group), color = NA) +
  geom_sf(data = unmapped, fill = "white", color = NA) +
  geom_sf(data = huc10,
          aes(color = "HUC-10 boundary", linetype = "HUC-10 boundary"),
          fill = NA, linewidth = 0.35, show.legend = FALSE) +
  geom_sf(data = basin,
          aes(color = "Gage basin (USGS 09405501)",
              linetype = "Gage basin (USGS 09405501)"),
          fill = NA, linewidth = 0.7, show.legend = FALSE) +
  geom_sf(data = flow %>% filter(streamorde >= 3),
          aes(color = "Stream (NHDPlus)", linetype = "Stream (NHDPlus)"),
          linewidth = 0.3, show.legend = FALSE) +
  geom_sf(data = flt, aes(color = "Fault", linetype = "Fault"),
          linewidth = 0.25, show.legend = FALSE) +
  geom_sf(data = park, aes(color = "Park boundary", linetype = "Park boundary"),
          fill = NA, linewidth = 0.5, show.legend = FALSE) +
  geom_path(data = legend_df,
            aes(x, y, color = what, linetype = what, group = what),
            linewidth = 0.6, key_glyph = "path") +
  geom_sf(data = gage_pt, shape = 24, fill = "white", color = "black",
          size = 1.8) +
  geom_sf(data = sites_p, aes(shape = kolob), fill = "#cc3030",
          color = "black", size = 1.6, stroke = 0.4) +
  # hard-coded label offsets (m) chosen by eye
  do.call(annotate, c(list("label", x = gage_xy[1] + 900, y = gage_xy[2] - 1800,
                           label = "USGS 09405501"), label_style)) +
  do.call(annotate, c(list("label", x = kolob_xy[1] + 2500,
                           y = kolob_xy[2] - 6500,
                           label = "Taylor Creek Well\n(Kolob Canyons)"),
                      label_style)) +
  annotate("text", x = mean(bb_a[c(1, 3)]) + 9000, y = bb_a[4] - 9000,
           label = "Not covered by\ngeologic map", family = "Arial",
           size = 2.6, color = "grey40", fontface = "italic") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  scale_fill_manual(values = cols, name = "Explanation", drop = FALSE) +
  line_scales +
  coord_sf(crs = UTM, xlim = bb_a[c(1, 3)], ylim = bb_a[c(2, 4)],
           expand = FALSE, datum = 4326) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25,
                              text_family = "Arial", text_cex = 0.6,
                              height = unit(0.12, "cm")) +
  ggspatial::annotation_north_arrow(
    location = "tr", height = unit(0.7, "cm"), width = unit(0.5, "cm"),
    style = ggspatial::north_arrow_orienteering(text_family = "Arial",
                                                text_size = 6)
  ) +
  theme_map +
  labs(title = "(a) North Fork Virgin River basin")

# (b) Zion Canyon
zc <- sites_p %>% filter(!kolob)
bb_b <- st_bbox(st_buffer(zc, 1200))

pB <- ggplot() +
  geom_sf(data = geo_plot, aes(fill = group), color = NA, show.legend = FALSE) +
  geom_sf(data = jln, aes(color = "Joint", linetype = "Joint"),
          linewidth = 0.15, show.legend = FALSE) +
  geom_sf(data = flt, aes(color = "Fault", linetype = "Fault"),
          linewidth = 0.4, show.legend = FALSE) +
  geom_sf(data = nk_contact,
          aes(color = "Navajo-Kayenta contact",
              linetype = "Navajo-Kayenta contact"),
          linewidth = 0.55, show.legend = FALSE) +
  geom_sf(data = flow,
          aes(color = "Stream (NHDPlus)", linetype = "Stream (NHDPlus)"),
          linewidth = 0.35, show.legend = FALSE) +
  geom_sf(data = zc, shape = 21, fill = "#cc3030", color = "black",
          size = 1.8, stroke = 0.4) +
  ggrepel::geom_label_repel(
    data = lab_pt(zc), aes(x, y, label = Location), size = 2.7,
    family = "Arial", label.size = 0, fill = alpha("white", 0.8),
    label.padding = 0.08, min.segment.length = 0, segment.size = 0.25,
    box.padding = 0.45, force = 3, seed = 4
  ) +
  scale_fill_manual(values = cols) +
  line_scales +
  coord_sf(crs = UTM, xlim = bb_b[c(1, 3)], ylim = bb_b[c(2, 4)],
           expand = FALSE, datum = 4326) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.2,
                              text_family = "Arial", text_cex = 0.6,
                              height = unit(0.12, "cm"),
                              bar_cols = c("black", "white")) +
  theme_map +
  labs(title = "(b) Zion Canyon")

# Location inset: western conterminous USA (Albers, EPSG:5070), Utah
# highlighted, with Mexico clipped to the inset frame
st_us <- st_read(ext_path("basemap_states.gpkg"), quiet = TRUE) %>%
  st_transform(5070)
zion_pt <- st_centroid(park) %>% st_transform(5070)
inset_bb <- st_bbox(st_us)
inset_bb["ymin"] <- inset_bb["ymin"] - 280000   # room for northern Mexico
mex <- st_read(ext_path("basemap_mexico.gpkg"), quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(5070) %>%
  st_crop(inset_bb)
lab_st <- st_point_on_surface(st_us) %>%
  mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])

hlm <- ggplot() +
  geom_sf(data = mex, fill = "grey88", color = "grey40", linewidth = 0.2) +
  geom_sf(data = st_us, aes(fill = STUSPS == "UT"), color = "grey40",
          linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#e8d4b8", `FALSE` = "white")) +
  geom_sf(data = zion_pt, shape = 22, size = 1.8, fill = "#cc3030",
          color = "black") +
  geom_text(data = lab_st, aes(x, y, label = STUSPS), size = 1.8,
            family = "Arial", color = "grey25") +
  # hard-coded country label positions (m, EPSG:5070)
  annotate("text", x = inset_bb["xmax"] - 60000, y = inset_bb["ymax"] - 120000,
           label = "USA", hjust = 1, size = 2, family = "Arial",
           fontface = "bold") +
  annotate("text",
           x = inset_bb["xmin"] + 0.62 * (inset_bb["xmax"] - inset_bb["xmin"]),
           y = inset_bb["ymin"] + 90000, label = "MEX", size = 2,
           family = "Arial", fontface = "bold") +
  coord_sf(xlim = inset_bb[c("xmin", "xmax")], ylim = inset_bb[c("ymin", "ymax")],
           datum = NA, expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = "grey30",
                                        linewidth = 0.3))
pA <- pA + inset_element(hlm, left = 0.0, bottom = 0.715, right = 0.255,
                         top = 1.0, align_to = "panel")

fig <- pA + pB + plot_layout(widths = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.margin = margin(0, 0, 0, 0))
fig <- fig & guides(fill = guide_legend(ncol = 3, order = 1),
                    color = guide_legend(ncol = 3, order = 2),
                    linetype = guide_legend(ncol = 3, order = 2))
ggsave(file.path(FIG, "fig1_site_map.png"), fig, width = 7.0, height = 5.0,
       dpi = 450, bg = "white")

cat("done\n")
print(ws)
print(geo_in_basin)
print(geo_in_park)
print(head(fault_stats, 10))
print(head(fault_park_named, 10))
print(site_tab %>%
        select(Location, GLG_SYM, GLG_NAME, in_gage_basin, in_park,
               dist_fault_m, nearest_named_fault, dist_named_fault_m,
               dist_navajo_kayenta_contact_m, dist_joint_m, huc12_name))


# ---- Side of the North Fork Virgin River and site elevation (Fig. 5) --------

nfv <- flow_all %>%
  filter(grepl("North Fork Virgin", gnis_name)) %>%
  st_union()
zs <- sites %>% filter(Location != KOLOB)
# nearest point on the river to each site (every second point of the
# connecting segments)
near_pt <- st_nearest_points(zs, nfv) %>% st_cast("POINT")
near_pt <- near_pt[seq(2, length(near_pt), 2)]
elev <- read_csv(ext_path("site_elevations_z14.csv"), show_col_types = FALSE)

side <- tibble(
  Location = zs$Location,
  side = ifelse(st_coordinates(zs)[, 1] > st_coordinates(near_pt)[, 1],
                "east", "west"),
  dist_river_m = round(as.numeric(st_distance(zs, nfv)))
) %>%
  left_join(elev %>% select(Location, elev_m), by = "Location")
write_csv(side, file.path(DER, "site_side_elevation.csv"))
print(side)


# ---- Regression checks against the paper (Table 1, section 2.1, Fig. 5) -----

stopifnot(identical(
  side$side[match(c("The Grotto", "Refrigerator Spring", "Birch Spring",
                    "Weeping Rock"), side$Location)],
  c("east", "west", "west", "east")
))
check_close(c(ws$basin_km2, ws$park_in_basin_km2, ws$basin_mapped_pct),
            c(897.5, 185.2, 60.1), 0.05,
            "basin area, park area in basin, mapped share")
check_close(c(ws$fault_km_basin, ws$joint_km_basin), c(44.3, 239.0), 0.05,
            "fault and joint length in basin")
check_close(site_tab$dist_navajo_kayenta_contact_m[
  match(c("Weeping Rock", "Refrigerator Spring", "Birch Spring"),
        site_tab$Location)
], c(3, 497, 55), 0.5, "distance to Navajo-Kayenta contact")
