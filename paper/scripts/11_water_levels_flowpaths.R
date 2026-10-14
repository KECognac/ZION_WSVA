# 11_water_levels_flowpaths.R -- Navajo aquifer water levels and interpreted
# flow (paper section 4.2; ESM Section S4, Fig. S5).
#
# Heads at NWIS wells completed in the Navajo Sandstone (land-surface
# altitude minus mean depth to water) and land-surface elevations of the
# canyon springs are interpolated by inverse-distance weighting (power 2) on
# a 250-m grid, as in the NWIS groundwater section of
# Zion_Isotope_Analysis.Rmd. Contours are drawn only within 8 km of a data
# point. Land above 2,400 m marks candidate recharge areas.
# Panel b shows the mean primary young-water fraction on the geology.

source("paper/scripts/00_setup.R")
suppressMessages({
  library(sf)
  library(terra)
})
sf_use_s2(FALSE)


# ---- Geology (same groups as 09_site_geology_watershed.R) -------------------

gdb <- raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb")
geo <- st_read(gdb, "zionglg", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(UTM)
flt <- st_read(gdb, "zionflt", quiet = TRUE) %>% st_transform(UTM)
jln <- st_read(gdb, "zionjln", quiet = TRUE) %>% st_transform(UTM)

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
    TRUE ~ "Other mapped units"
  ))
geo_g <- geo %>%
  group_by(group) %>%
  summarise(.groups = "drop") %>%
  st_make_valid()

# hard-coded figure colors (same as Fig. 1)
cols <- c(
  "Cretaceous formations" = "#a9b38a", "Carmel Formation" = "#7f9fbf",
  "Temple Cap Formation"  = "#b9a3c9", "Navajo Sandstone" = "#efc9a4",
  "Kayenta Formation"     = "#b5654d", "Moenave and older" = "#8c6d5a",
  "Quaternary basalt"     = "#4d4d4d", "Quaternary deposits" = "#e8e3c9",
  "Other mapped units"    = "#cccccc"
)
geo_g$group <- factor(geo_g$group, levels = names(cols))

navajo <- geo_g %>% filter(group == "Navajo Sandstone")
kay <- geo_g %>% filter(group == "Kayenta Formation")
nk_contact <- st_intersection(st_boundary(navajo), st_boundary(kay)) %>%
  st_collection_extract("LINESTRING")
hurricane <- flt %>% filter(grepl("Hurricane", FNAME))

park <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  filter(UNIT_CODE == "ZION") %>%
  st_make_valid() %>%
  st_transform(UTM) %>%
  summarise()


# ---- Sites and site-mean young fractions ------------------------------------

sg <- read.csv(file.path(DER, "site_geology.csv"))
zc_names <- c(REPEAT_SPRINGS, ZC_2025)   # Zion Canyon sites (no Taylor Creek)
s1 <- read.csv(file.path(DER, "tritium_samples.csv"))
yf <- s1 %>%
  group_by(Location) %>%
  summarise(yf = mean(yf, na.rm = TRUE), .groups = "drop")

# hard-coded check: rounded site means (%) must match Table 2
expect <- c("Echo Canyon Spring" = 100, "Weeping Rock" = 73, "The Grotto" = 48,
            "Refrigerator Spring" = 33, "Larry's Grotto" = 10,
            "Birch Sediment Trap" = 0, "East Entrance Well" = 1,
            "Lower Grotto Sandtrap" = 0)
stopifnot(all(round(100 * yf$yf[match(names(expect), yf$Location)]) == expect))

sites <- sg %>%
  filter(Location %in% zc_names) %>%
  left_join(yf, by = "Location") %>%
  mutate(path = ifelse(Location %in% PATH_A,
                       "Path A (high conductivity)",
                       "Path B (low conductivity)")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM)


# ---- (a) Water-level surface ------------------------------------------------

w <- read.csv(ext_path("nwis_navajo_wells.csv")) %>%
  filter(!is.na(gw_alt_ft)) %>%
  st_as_sf(coords = c("dec_long_va", "dec_lat_va"), crs = 4269) %>%
  st_transform(UTM) %>%
  mutate(head_m = gw_alt_ft * 0.3048, type = "Navajo well (NWIS)")

# Springs use their land-surface elevation as the head, as in the Rmd
springs <- sites %>% filter(!grepl("Well", Location))
elev <- read.csv(ext_path("site_elevations_z14.csv"))
springs$head_m <- elev$elev_m[match(springs$Location, elev$Location)]
pts <- bind_rows(
  w %>% select(head_m, type),
  springs %>% transmute(head_m, type = "Spring (land-surface elevation)")
)

# hard-coded: extent of the geologic map (UTM 12N, m)
bb_a <- st_bbox(c(xmin = 300700, ymin = 4110200, xmax = 345000, ymax = 4153700),
                crs = st_crs(UTM))

# Inverse-distance weighting (power 2) on a 250-m grid, 2 km beyond the map
g <- rast(ext(bb_a[1] - 2000, bb_a[3] + 2000, bb_a[2] - 2000, bb_a[4] + 2000),
          resolution = 250, crs = paste0("EPSG:", UTM))
xy <- crds(g, na.rm = FALSE)
pc <- st_coordinates(pts)
D <- sqrt(outer(xy[, 1], pc[, 1], "-")^2 + outer(xy[, 2], pc[, 2], "-")^2)
W <- 1 / pmax(D, 1)^2
values(g) <- as.vector((W %*% pts$head_m) / rowSums(W))

# Contours only within 8 km of a measurement (no extrapolation across gaps);
# 100-m interval, short pieces (< 2.5 km) dropped
near <- st_union(st_buffer(pts, 8000))
g_m <- mask(g, vect(near))
cnt <- st_as_sf(as.contour(g_m, levels = seq(1000, 2000, 100))) %>%
  st_set_crs(UTM) %>%
  st_intersection(st_as_sfc(bb_a)) %>%
  st_cast("MULTILINESTRING") %>%
  st_cast("LINESTRING", warn = FALSE) %>%
  mutate(len = as.numeric(st_length(.))) %>%
  filter(len > 2500)
cnt_lab <- cnt %>%
  group_by(level) %>%
  slice_max(len, n = 1) %>%
  ungroup() %>%
  mutate(geometry = st_line_sample(geometry, sample = 0.5)) %>%
  st_cast("POINT")

# Gradient check: mean down-gradient azimuth east and west of the North Fork
# near the canyon. Aspect = direction of steepest descent, degrees from north.
gx <- terrain(g, "aspect", unit = "degrees")
zone <- function(xmin, xmax, ymin, ymax) {
  v <- crop(gx, ext(xmin, xmax, ymin, ymax))
  a <- values(v, na.rm = TRUE) * pi / 180
  (atan2(mean(sin(a)), mean(cos(a))) * 180 / pi) %% 360
}
canyon_x <- mean(st_coordinates(springs)[, 1])
# hard-coded zones: 1.5-12 km east or west of the canyon; northings (m) span
# the springs and the wells on each side
az <- data.frame(
  zone = c("East of North Fork (to 12 km east)",
           "West of North Fork (to 12 km west)"),
  azimuth_deg = c(zone(canyon_x + 1500, canyon_x + 12000, 4119000, 4132000),
                  zone(canyon_x - 12000, canyon_x - 1500, 4122000, 4140000))
)
print(az)
write.csv(az, file.path(DER, "navajo_head_gradient_azimuth.csv"),
          row.names = FALSE)

used <- w %>%
  filter(st_within(., st_as_sfc(bb_a), sparse = FALSE)[, 1]) %>%
  st_drop_geometry() %>%
  select(site_no, station_nm, alt_va, dtw_ft_mean, n_meas, first, last,
         gw_alt_ft, head_m)
print(used)
write.csv(used, file.path(DER, "navajo_wells_in_fig6.csv"), row.names = FALSE)
write.csv(st_drop_geometry(springs) %>% select(Location, head_m),
          file.path(DER, "spring_heads.csv"), row.names = FALSE)

# Candidate recharge areas: land above 2,400 m (polygons > 1 km2)
dem <- rast(ext_path("dem_z9.tif"))
dem_a <- crop(dem, ext(bb_a[1] - 3000, bb_a[3] + 3000,
                       bb_a[2] - 3000, bb_a[4] + 3000))
hi <- as.polygons(dem_a >= 2400) %>%
  st_as_sf() %>%
  filter(.[[1]] == 1) %>%
  st_make_valid() %>%
  st_cast("POLYGON") %>%
  mutate(km2 = as.numeric(st_area(.)) / 1e6) %>%
  filter(km2 > 1)
hs <- shade(terrain(dem_a, "slope", unit = "radians"),
            terrain(dem_a, "aspect", unit = "radians"), 40, 315)
hs_df <- as.data.frame(hs, xy = TRUE)
names(hs_df)[3] <- "hs"

flow <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  st_transform(UTM)
nfv <- flow %>% filter(grepl("North Fork Virgin", gnis_name))


# ---- Schematic arrows (interpreted flow directions) -------------------------

# hard-coded arrow start/end offsets (m) from the center of the springs,
# placed by eye to show the interpreted directions in panel a
sx <- st_coordinates(springs)
cx <- mean(sx[, 1])
cy <- mean(sx[, 2])
arr <- data.frame(
  x    = c(cx + 16500, cx - 13000, cx + 1500),
  y    = c(cy + 9000,  cy + 18000, bb_a[4] - 600),
  xend = c(cx + 3500,  cx - 2800,  cx + 600),
  yend = c(cy + 1800,  cy + 5200,  cy + 9500),
  what = c("east", "west", "north")
)
mk <- read.csv(file.path(DER, "markagunt_targets.csv"))

w_in <- w[st_within(w, st_as_sfc(bb_a), sparse = FALSE)[, 1], ]


# ---- Plotting ---------------------------------------------------------------

theme_map <- theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 10),
    legend.key.size = unit(0.38, "cm"),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8, face = "bold")
  )
lab_pt <- function(x) {
  x %>% mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])
}
lcols <- c("Water-level contour (m)" = "#1f4e99",
           "Hurricane Fault Zone" = "#cc3030",
           "North Fork Virgin River" = "#3a8fa0",
           "Park boundary" = "#208040",
           "Interpreted flow direction" = "black")
# Off-map dummy lines that only produce the legend keys
legend_df <- data.frame(what = rep(names(lcols), each = 2),
                        x = rep(c(0, 1), length(lcols)), y = 0)
hurricane_xy <- st_coordinates(st_centroid(hurricane %>% summarise()))

pA <- ggplot() +
  geom_raster(data = hs_df, aes(x, y, alpha = hs), fill = "grey20",
              show.legend = FALSE) +
  scale_alpha(range = c(0.35, 0)) +
  geom_sf(data = hi,
          aes(fill = "Land above 2,400 m (candidate recharge area)"),
          color = "#7a5a00", linewidth = 0.2, alpha = 0.45) +
  geom_sf(data = park, fill = NA, color = lcols["Park boundary"],
          linewidth = 0.45) +
  geom_sf(data = nfv, color = lcols["North Fork Virgin River"],
          linewidth = 0.6) +
  geom_sf(data = hurricane, color = lcols["Hurricane Fault Zone"],
          linewidth = 0.6) +
  geom_sf(data = cnt, color = lcols["Water-level contour (m)"],
          linewidth = 0.35) +
  geom_sf_label(data = cnt_lab, aes(label = format(level, big.mark = ",")),
                size = 2.2, family = "Arial", color = "#1f4e99",
                label.size = 0, fill = alpha("white", 0.75),
                label.padding = unit(0.04, "lines")) +
  geom_segment(data = arr, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 1.1, color = "black",
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_path(data = legend_df, aes(x, y, color = what, group = what),
            linewidth = 0.6, key_glyph = "path") +
  geom_sf(data = pts, aes(shape = type), fill = "white", color = "black",
          size = 1.5, stroke = 0.4) +
  ggrepel::geom_text_repel(
    data = lab_pt(w_in), aes(x, y, label = format(round(head_m), big.mark = ",")),
    size = 2.3, family = "Arial", fontface = "bold", min.segment.length = 0.2,
    segment.size = 0.2, box.padding = 0.35, seed = 2
  ) +
  annotate("label", x = arr$x[3] - 900, y = bb_a[4] - 1200, hjust = 1,
           vjust = 1, size = 2.3, family = "Arial", label.size = 0,
           fill = alpha("white", 0.8),
           label = sprintf(paste0("To Markagunt Plateau (hypothesized):\n",
                                  "Cedar Breaks NM %.0f km, Brian Head %.0f km"),
                           mk$dist_cebr_nearest_km, mk$dist_brian_head_km)) +
  annotate("text", x = hurricane_xy[1] + 1200, y = hurricane_xy[2],
           label = "Hurricane\nFault Zone", size = 2.3, family = "Arial",
           color = "#cc3030", hjust = 0) +
  scale_fill_manual(
    values = c("Land above 2,400 m (candidate recharge area)" = "#e0b040"),
    name = NULL
  ) +
  scale_shape_manual(values = c("Navajo well (NWIS)" = 24,
                                "Spring (land-surface elevation)" = 21),
                     name = NULL) +
  scale_color_manual(values = lcols, name = NULL, breaks = names(lcols),
                     limits = names(lcols)) +
  coord_sf(crs = UTM, xlim = bb_a[c(1, 3)], ylim = bb_a[c(2, 4)],
           expand = FALSE, datum = 4326) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25,
                              text_family = "Arial", text_cex = 0.6,
                              height = unit(0.12, "cm")) +
  ggspatial::annotation_north_arrow(
    location = "br", height = unit(0.7, "cm"), width = unit(0.5, "cm"),
    style = ggspatial::north_arrow_orienteering(text_family = "Arial",
                                                text_size = 6)
  ) +
  theme_map +
  labs(title = "(a) Navajo aquifer water levels and interpreted flow")

# (b) Zion Canyon: young-water fraction and geology
zc <- sites
bb_b <- st_bbox(st_buffer(zc, 1200))
# percent label: one decimal for values between 0 and 1%
yf_label <- function(yf) {
  ifelse(yf > 0 & yf < 0.01, sprintf("%.1f%%", 100 * yf),
         sprintf("%.0f%%", 100 * yf))
}

pB <- ggplot() +
  geom_sf(data = geo_g, aes(fill = group), color = NA) +
  scale_fill_manual(values = cols, name = "Geologic unit", drop = TRUE,
                    guide = guide_legend(order = 3, ncol = 1)) +
  geom_sf(data = jln, color = "grey45", linewidth = 0.15) +
  geom_sf(data = flt, color = "black", linewidth = 0.4) +
  geom_sf(data = nk_contact, color = "#cc3030", linewidth = 0.55) +
  geom_sf(data = flow, color = "#3a8fa0", linewidth = 0.35) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = zc, aes(fill = 100 * yf, shape = path), size = 3.2,
          color = "black", stroke = 0.5) +
  ggrepel::geom_label_repel(
    data = lab_pt(zc), aes(x, y, label = sprintf("%s\n%s", Location, yf_label(yf))),
    size = 2.4, family = "Arial", label.size = 0, fill = alpha("white", 0.8),
    label.padding = 0.08, lineheight = 0.9, min.segment.length = 0,
    segment.size = 0.25, box.padding = 0.5, force = 3, seed = 4
  ) +
  scale_fill_gradient(
    low = "white", high = "#08306b", limits = c(0, 100),
    name = "Mean young-water\nfraction (%)",
    guide = guide_colorbar(barheight = unit(2.2, "cm"),
                           barwidth = unit(0.3, "cm"), order = 1)
  ) +
  scale_shape_manual(
    values = c("Path A (high conductivity)" = 21,
               "Path B (low conductivity)" = 22),
    name = NULL,
    guide = guide_legend(order = 2, override.aes = list(fill = "grey70"))
  ) +
  coord_sf(crs = UTM, xlim = bb_b[c(1, 3)], ylim = bb_b[c(2, 4)],
           expand = FALSE, datum = 4326) +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.2,
                              text_family = "Arial", text_cex = 0.6,
                              height = unit(0.12, "cm")) +
  theme_map +
  labs(title = "(b) Zion Canyon: young-water fraction and geology")

# Line key for panel b (faults, joints, contact, stream): invisible dummy
# paths whose legend keys are recolored
key_b <- data.frame(what = rep(c("Fault", "Joint", "Navajo-Kayenta contact",
                                 "Stream (NHDPlus)"), each = 2),
                    x = 0, y = 0)
pB <- pB +
  geom_path(data = key_b, aes(x, y, linetype = what, group = what), color = NA) +
  scale_linetype_manual(
    values = rep("solid", 4), name = NULL,
    guide = guide_legend(override.aes = list(
      color = c("black", "grey45", "#cc3030", "#3a8fa0"), linewidth = 0.6
    ))
  )

lr <- theme(legend.position = "right", legend.box = "vertical",
            legend.spacing.y = unit(0.1, "cm"))
fig <- (pA + lr) / (pB + lr) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(FIG, "figS5_water_levels.png"), fig, width = 7.2, height = 9.4,
       dpi = 450, bg = "white")
message("written ", file.path(FIG, "figS5_water_levels.png"))
