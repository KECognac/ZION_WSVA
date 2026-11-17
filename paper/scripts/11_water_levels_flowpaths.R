# 11_water_levels_flowpaths.R -- Navajo aquifer water levels and interpreted
# flow (paper section 4.2; Fig. 7a-b; ESM Section S4, Fig. S5).
#
# The water-level surface is a thin-plate smoothing spline (mgcv::gam,
# s(x, y, k = 20)) on a 250-m grid, fitted to these heads:
#   - NWIS wells completed in the Navajo Sandstone (land-surface altitude minus
#     mean depth to water);
#   - spring emergence points picked on Google Earth imagery, at their Google
#     Earth elevation (springs emerge where the water table meets the land
#     surface). The three points K. Cognac labelled "Upper Spring Emergence"
#     (about 1,980-2,050 m, near the Temple Cap) are left out: they discharge
#     perched water far above the canyon springs;
#   - sampled springs with no emergence point within 150 m, at their
#     land-surface elevation (8-m DEM);
#   - one point on the North Fork Virgin River where it meets the
#     Navajo-Kayenta contact, at the river elevation. This is the downstream
#     end of the first reach (500 m or longer) that the geologic map shows on
#     Navajo Sandstone at the river. Downstream of it the river runs on the
#     Kayenta and on alluvium, where the Navajo water table need not meet it,
#     so no other river points are used.
# Inverse-distance weighting, used in the NWIS groundwater section of
# Zion_Isotope_Analysis.Rmd and in earlier drafts, honors every point and so
# draws a closed ring of contours around each canyon spring; with these heads
# the 1,400 m contour crossed the North Fork 12 times. The spline is smooth
# and does not pass through every point. k = 20 is the most flexible of the
# basis sizes tried (10, 20, 30, 50) whose contours cross the river at most
# once; the fit at the springs is printed below.
# Contours are drawn only within 8 km of a data point, and none may cross the
# North Fork more than once (checked below). Land above 2,400 m marks
# candidate recharge areas.
# Panel b shows the mean primary young-water fraction on the geology.
#
# Outputs
#   derived: navajo_head_gradient_azimuth.csv, navajo_wells_in_fig6.csv,
#            spring_heads.csv (all non-well heads used),
#            fig_flow_panels.rds (panels a and b for Fig. 7,
#            assembled in 19_spring_emergence_band.R)
#   figures: figS5_water_levels.png,
#            figS6_contour_methods.png (inverse-distance weighting vs spline)
# Runs after 03 (tritium) and 09 (site geology).

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
CARMEL     <- c("Jc", "Jcw", "Jcp", "Jcx", "Jcc", "Jccu", "Jccl")
TEMPLE_CAP <- c("Jt", "Jtw", "Jts")
NAVAJO     <- c("Jn", "Jnw", "Jnp", "Jnb", "Jnl")
KAYENTA    <- c("Jk", "Jkt", "Jku", "Jks")
MOENAVE    <- c("JTRkm", "JTRm", "Jmw", "JTRmd")

geo <- geo %>%
  mutate(group = case_when(
    GLG_SYM %in% CARMEL ~ "Carmel Formation",
    GLG_SYM %in% TEMPLE_CAP ~ "Temple Cap Formation",
    GLG_SYM %in% NAVAJO ~ "Navajo Sandstone",
    GLG_SYM %in% KAYENTA ~ "Kayenta Formation",
    GLG_SYM %in% MOENAVE | str_detect(GLG_SYM, "^TR|^P|^PZ") ~
      "Moenave and older",
    str_detect(GLG_SYM, "^K") ~ "Cretaceous formations",
    str_detect(GLG_SYM, "^Qb") & GLG_SYM != "Qblcp" ~ "Quaternary basalt",
    str_detect(GLG_SYM, "^Q") ~ "Quaternary deposits",
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
geo_g <- geo_g %>% mutate(group = factor(group, levels = names(cols)))

navajo <- geo_g %>% filter(group == "Navajo Sandstone")
kay <- geo_g %>% filter(group == "Kayenta Formation")
nk_contact <- st_intersection(st_boundary(navajo), st_boundary(kay)) %>%
  st_collection_extract("LINESTRING")
hurricane <- flt %>% filter(str_detect(FNAME, "Hurricane"))

park <- st_read(ext_path("nps_boundary/nps_boundary.shp"), quiet = TRUE) %>%
  filter(UNIT_CODE == "ZION") %>%
  st_make_valid() %>%
  st_transform(UTM) %>%
  summarise()


# ---- Sites and site-mean young fractions ------------------------------------

sg <- read_csv(file.path(DER, "site_geology.csv"), show_col_types = FALSE)
zc_names <- c(REPEAT_SPRINGS, ZC_2025)   # Zion Canyon sites (no Taylor Creek)
yf <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  group_by(Location) %>%
  summarise(yf = mean(yf, na.rm = TRUE), .groups = "drop")

# hard-coded check: rounded site means (%) must match Table 2
expect <- tribble(
  ~Location,               ~pct,
  "Echo Canyon Spring",    100,
  "Weeping Rock",          73,
  "The Grotto",            48,
  "Refrigerator Spring",   33,
  "Larry's Grotto",        10,
  "Birch Spring",          0,
  "East Entrance Well",    1,
  "Lower Grotto Sandtrap", 0
) %>%
  left_join(yf, by = "Location")
stopifnot(all(round(100 * expect$yf) == expect$pct))

sites <- sg %>%
  filter(Location %in% zc_names) %>%
  left_join(yf, by = "Location") %>%
  mutate(path = if_else(Location %in% PATH_A,
                        "Path A (high conductivity)",
                        "Path B (low conductivity)")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM)


# ---- (a) Water-level surface ------------------------------------------------

# Navajo wells with a mean water level: head = land surface - depth to water
w <- read_csv(ext_path("nwis_navajo_wells.csv"), show_col_types = FALSE,
              col_types = cols(site_no = col_character())) %>%
  filter(!is.na(gw_alt_ft)) %>%
  st_as_sf(coords = c("dec_long_va", "dec_lat_va"), crs = 4269) %>%
  st_transform(UTM) %>%
  mutate(head_m = gw_alt_ft * 0.3048, type = "Navajo well (NWIS)")

# Sampled springs at their land-surface elevation (Birch and Refrigerator at
# their emergence points, with Google Earth elevations; see site_locations()).
# These are the spring symbols on the maps and the reference for the arrows.
elev <- read_csv(ext_path("site_elevations_z14.csv"), show_col_types = FALSE)
springs <- sites %>%
  filter(!str_detect(Location, "Well")) %>%
  left_join(elev %>% select(Location, head_m = elev_m), by = "Location")

# Emergence points from Google Earth, without the perched "Upper" points
EMERGENCE_MATCH_M <- 150   # a sampled spring this close is the same spring
em <- read_emergence_points(UTM) %>%
  filter(!str_detect(label, regex("^upper spring emergence", ignore_case = TRUE)))
# sampled springs with no emergence point nearby keep their own position
springs_unmatched <- springs %>%
  filter(map_lgl(st_is_within_distance(springs, em, EMERGENCE_MATCH_M),
                 ~ length(.x) == 0))

# River point at the Navajo-Kayenta contact. River elevation every 50 m along
# the North Fork: lowest cell of the canyon DEM within 50 m of the NHDPlus line
# (the mapped line sometimes sits on the canyon wall), as in 16 and 19.
dem_c <- rast(ext_path("dem_canyon_z14.tif"))
names(dem_c) <- "elev"
nfv_line <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  st_transform(UTM) %>%
  filter(str_detect(gnis_name, "North Fork Virgin")) %>%
  st_union() %>%
  st_line_merge()
stopifnot(st_geometry_type(nfv_line) == "LINESTRING")
riv <- st_line_sample(nfv_line, density = 1 / 50) %>%
  st_cast("POINT") %>%
  st_sf(geometry = .) %>%
  mutate(head_m = terra::extract(dem_c, vect(st_buffer(geometry, 50)),
                                 fun = min, na.rm = TRUE)[, 2],
         navajo = lengths(st_intersects(geometry, geo_g %>%
                                          filter(group == "Navajo Sandstone"))) > 0) %>%
  filter(is.finite(head_m))
# order the points upstream (rising elevation)
if (first(riv$head_m) > last(riv$head_m)) riv <- riv[rev(seq_len(nrow(riv))), ]
# first run of >= 10 consecutive points (500 m) on the Navajo: its first point
nav_runs <- riv %>%
  mutate(i = row_number(), run = consecutive_id(navajo)) %>%
  st_drop_geometry() %>%
  filter(navajo) %>%
  group_by(run) %>%
  summarise(start = min(i), n = n(), .groups = "drop") %>%
  filter(n >= 10)
river_nk <- riv[first(nav_runs$start), ]
message(sprintf("River at the Navajo-Kayenta contact: %.0f m (UTM %.0f E, %.0f N)",
                river_nk$head_m, st_coordinates(river_nk)[1],
                st_coordinates(river_nk)[2]))

spring_heads <- bind_rows(
  em %>% transmute(name = label, head_m = elev_ge_m,
                   source = "Emergence point (Google Earth elevation)"),
  springs_unmatched %>% transmute(name = Location, head_m,
                                  source = "Sampled spring (8-m DEM elevation)"),
  river_nk %>% transmute(name = "North Fork Virgin River at the Navajo-Kayenta contact",
                         head_m, source = "River (8-m DEM elevation)")
)
pts <- bind_rows(
  w %>% select(head_m, type),
  spring_heads %>% transmute(head_m, type = "Spring (land-surface elevation)")
)
# points drawn on the maps: wells and the sampled springs, as before
pts_map <- bind_rows(
  w %>% select(head_m, type),
  springs %>% transmute(head_m, type = "Spring (land-surface elevation)")
)

# hard-coded: extent of the geologic map (UTM 12N, m)
bb_a <- st_bbox(c(xmin = 300700, ymin = 4110200, xmax = 345000, ymax = 4153700),
                crs = st_crs(UTM))

# Thin-plate smoothing spline of head on map position, evaluated on a 250-m
# grid extending 2 km beyond the map (see the header for the choice of k)
SPLINE_K <- 20
pc <- st_coordinates(pts)
heads_xy <- tibble(x = pc[, 1], y = pc[, 2], head_m = pts$head_m,
                   type = pts$type)
fit <- mgcv::gam(head_m ~ s(x, y, k = SPLINE_K), data = heads_xy)
g <- rast(ext(bb_a[1] - 2000, bb_a[3] + 2000, bb_a[2] - 2000, bb_a[4] + 2000),
          resolution = 250, crs = paste0("EPSG:", UTM))
xy <- crds(g, na.rm = FALSE)
values(g) <- as.vector(predict(fit, newdata = tibble(x = xy[, 1], y = xy[, 2])))

# Fit of the surface: residual (observed - fitted head, m) by kind of point
heads_xy %>%
  mutate(resid = residuals(fit)) %>%
  group_by(type) %>%
  summarise(n = n(), rmse_m = sqrt(mean(resid^2)),
            max_abs_m = max(abs(resid)), .groups = "drop") %>%
  print()

# For comparison (ESM Fig. S6): inverse-distance weighting (power 2) of the
# same heads on the same grid, the method of earlier drafts. D is the distance
# from every grid cell (rows) to every head (columns); each cell gets the
# weighted mean head with weights 1 / D^2 (D at least 1 m).
g_idw <- g
D <- sqrt(outer(xy[, 1], pc[, 1], "-")^2 + outer(xy[, 2], pc[, 2], "-")^2)
W <- 1 / pmax(D, 1)^2
values(g_idw) <- as.vector((W %*% pts$head_m) / rowSums(W))

# Contours of surface r only within 8 km of a measurement (no extrapolation
# across gaps); 100-m interval, short pieces (< 2.5 km) dropped
near <- st_union(st_buffer(pts, 8000))
contour_lines <- function(r) {
  st_as_sf(as.contour(mask(r, vect(near)), levels = seq(1000, 2000, 100))) %>%
    st_set_crs(UTM) %>%
    st_intersection(st_as_sfc(bb_a)) %>%
    st_cast("MULTILINESTRING") %>%
    st_cast("LINESTRING", warn = FALSE) %>%
    mutate(len = as.numeric(st_length(.))) %>%
    filter(len > 2500)
}
# Points where the contours cross the North Fork, and their number by level
crossing_points <- function(cnt) {
  st_intersection(cnt, st_sf(geometry = nfv_line)) %>%
    st_cast("MULTIPOINT", warn = FALSE) %>%
    st_cast("POINT", warn = FALSE)
}
river_crossings <- function(cnt) {
  crossing_points(cnt) %>%
    st_drop_geometry() %>%
    count(level, name = "n_crossings")
}
cnt <- contour_lines(g)
cnt_idw <- contour_lines(g_idw)
xing <- river_crossings(cnt)
xing_idw <- river_crossings(cnt_idw)
message("River crossings by contour level, inverse-distance weighting:")
print(xing_idw)
message("River crossings by contour level, smoothing spline:")
print(xing)

# Check: a contour of the spline surface may cross the North Fork at most once
if (any(xing$n_crossings > 1)) {
  stop("contour(s) crossing the North Fork more than once: ",
       str_c(xing$level[xing$n_crossings > 1], collapse = ", "))
}
message("[OK] each contour crosses the North Fork at most once")

cnt_lab <- cnt %>%
  group_by(level) %>%
  slice_max(len, n = 1) %>%
  ungroup() %>%
  mutate(geometry = st_line_sample(geometry, sample = 0.5)) %>%
  st_cast("POINT")

# Gradient check: mean down-gradient azimuth east and west of the North Fork
# near the canyon. Aspect = direction of steepest descent, degrees from north.
gx <- terrain(g, "aspect", unit = "degrees")
# circular mean of the aspect inside a rectangle (UTM m)
zone <- function(xmin, xmax, ymin, ymax) {
  v <- crop(gx, ext(xmin, xmax, ymin, ymax))
  a <- values(v, na.rm = TRUE) * pi / 180
  (atan2(mean(sin(a)), mean(cos(a))) * 180 / pi) %% 360
}
canyon_x <- mean(st_coordinates(springs)[, 1])
# hard-coded zones: 1.5-12 km east or west of the canyon; northings (m) span
# the springs and the wells on each side
az <- tibble(
  zone = c("East of North Fork (to 12 km east)",
           "West of North Fork (to 12 km west)"),
  azimuth_deg = c(zone(canyon_x + 1500, canyon_x + 12000, 4119000, 4132000),
                  zone(canyon_x - 12000, canyon_x - 1500, 4122000, 4140000))
)
print(az)
write_csv(az, file.path(DER, "navajo_head_gradient_azimuth.csv"))

# Wells inside the map frame (the file name keeps the figure number of an
# earlier draft)
w_in <- w %>% st_filter(st_as_sfc(bb_a), .predicate = st_within)
used <- w_in %>%
  st_drop_geometry() %>%
  select(site_no, station_nm, alt_va, dtw_ft_mean, n_meas, first, last,
         gw_alt_ft, head_m)
print(used)
write_csv(used, file.path(DER, "navajo_wells_in_fig6.csv"))
spring_heads %>%
  mutate(utm_x = st_coordinates(.)[, 1], utm_y = st_coordinates(.)[, 2]) %>%
  st_drop_geometry() %>%
  write_csv(file.path(DER, "spring_heads.csv"))

# Candidate recharge areas: land above 2,400 m (polygons > 1 km2)
dem <- rast(ext_path("dem_z9.tif"))
dem_a <- crop(dem, ext(bb_a[1] - 3000, bb_a[3] + 3000,
                       bb_a[2] - 3000, bb_a[4] + 3000))
hi <- as.polygons(dem_a >= 2400) %>%
  st_as_sf() %>%
  rename(above = 1) %>%
  filter(above == 1) %>%
  st_make_valid() %>%
  st_cast("POLYGON") %>%
  mutate(km2 = as.numeric(st_area(.)) / 1e6) %>%
  filter(km2 > 1)
# Hillshade (sun 40 deg above the horizon, from the northwest)
hs <- shade(terrain(dem_a, "slope", unit = "radians"),
            terrain(dem_a, "aspect", unit = "radians"), 40, 315)
hs_df <- as.data.frame(hs, xy = TRUE) %>% rename(hs = 3)

flow <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  st_transform(UTM)
nfv <- flow %>% filter(str_detect(gnis_name, "North Fork Virgin"))


# ---- Schematic arrows (interpreted flow directions) -------------------------

# hard-coded arrow start/end offsets (m) from the center of the springs,
# placed by eye to show the interpreted directions in panel a
sx <- st_coordinates(springs)
cx <- mean(sx[, 1])
cy <- mean(sx[, 2])
arr <- tibble(
  x    = c(cx + 16500, cx - 13000, cx + 1500),
  y    = c(cy + 9000,  cy + 18000, bb_a[4] - 600),
  xend = c(cx + 3500,  cx - 2800,  cx + 600),
  yend = c(cy + 1800,  cy + 5200,  cy + 9500),
  what = c("east", "west", "north")
)
mk <- read_csv(file.path(DER, "markagunt_targets.csv"), show_col_types = FALSE)


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
# Add x and y columns (map coordinates) to an sf point layer, for labels
lab_pt <- function(x) {
  x %>% mutate(x = st_coordinates(.)[, 1], y = st_coordinates(.)[, 2])
}
lcols <- c("Water-level contour (m)" = "#1f4e99",
           "Hurricane Fault Zone" = "#cc3030",
           "North Fork Virgin River" = "#3a8fa0",
           "Park boundary" = "#208040",
           "Interpreted flow direction" = "black")
# Off-map dummy lines that only produce the legend keys
legend_df <- tibble(what = rep(names(lcols), each = 2),
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
  geom_sf(data = pts_map, aes(shape = type), fill = "white", color = "black",
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
  # fixed legend order: recharge area, then point symbols, then lines
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2),
         color = guide_legend(order = 3)) +
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
    data = lab_pt(zc) %>% mutate(lab = sprintf("%s\n%s", Location, yf_label(yf))),
    aes(x, y, label = lab),
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
key_b <- tibble(what = rep(c("Fault", "Joint", "Navajo-Kayenta contact",
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

# Panels for main-text Fig. 7 (assembled in 19_spring_emergence_band.R): panel
# a without the interpolated contours, which rest on five wells; both panels
# are saved as ggplot objects.
is_cnt <- map_lgl(pA$layers, ~ identical(.x$data, cnt) || identical(.x$data, cnt_lab))
stopifnot(sum(is_cnt) == 2)
pA_main <- pA
pA_main$layers <- pA$layers[!is_cnt]
pA_main <- suppressMessages(
  pA_main +
    scale_color_manual(values = lcols[-1], name = NULL, breaks = names(lcols)[-1],
                       limits = names(lcols)[-1], na.translate = FALSE)
) +
  # Fig. 7a legend order: lines, then the recharge area, then point symbols
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2),
         shape = guide_legend(order = 3)) +
  labs(title = "a  Navajo aquifer water levels (m) and interpreted flow")
saveRDS(list(pA = pA_main, pB = pB, theme_map = theme_map),
        file.path(DER, "fig_flow_panels.rds"))

lr <- theme(legend.position = "right", legend.box = "vertical",
            legend.spacing.y = unit(0.1, "cm"))
fig <- (pA + lr) / (pB + lr) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(FIG, "figS5_water_levels.png"), fig, width = 7.2, height = 9.4,
       dpi = 450, bg = "white")
message("written ", file.path(FIG, "figS5_water_levels.png"))


# ---- ESM Fig. S6: inverse-distance weighting vs smoothing spline ------------

# Panel of contours `cnt` over the river and the heads, within bounding box
# bb; red points mark where a contour crosses the river
method_panel <- function(cnt, title, subtitle, bb, lon_breaks = waiver()) {
  lab <- cnt %>%
    group_by(level) %>%
    slice_max(len, n = 1) %>%
    ungroup() %>%
    mutate(geometry = st_line_sample(geometry, sample = 0.5)) %>%
    st_cast("POINT")
  ggplot() +
    geom_raster(data = hs_df, aes(x, y, alpha = hs), fill = "grey20",
                show.legend = FALSE) +
    scale_alpha(range = c(0.35, 0)) +
    geom_sf(data = nfv, color = lcols["North Fork Virgin River"],
            linewidth = 0.6) +
    geom_sf(data = cnt, color = lcols["Water-level contour (m)"],
            linewidth = 0.35) +
    geom_sf_label(data = lab, aes(label = format(level, big.mark = ",")),
                  size = 2, family = "Arial", color = "#1f4e99",
                  label.size = 0, fill = alpha("white", 0.75),
                  label.padding = unit(0.04, "lines")) +
    geom_sf(data = pts %>% mutate(type = if_else(
              type == "Navajo well (NWIS)", "Navajo well (NWIS)",
              "Spring emergence or river at the Navajo-Kayenta contact")),
            aes(shape = type), fill = "white", color = "black", size = 1.3,
            stroke = 0.4) +
    geom_sf(data = crossing_points(cnt) %>%
              mutate(what = "Contour crossing the river"),
            aes(color = what), size = 1.6) +
    scale_shape_manual(values = c("Navajo well (NWIS)" = 24,
                                  "Spring emergence or river at the Navajo-Kayenta contact" = 21),
                       name = NULL) +
    scale_color_manual(values = c("Contour crossing the river" = "#d7301f"),
                       name = NULL) +
    scale_x_continuous(breaks = lon_breaks) +
    coord_sf(crs = UTM, xlim = bb[c(1, 3)], ylim = bb[c(2, 4)],
             expand = FALSE, datum = 4326) +
    theme_map +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 8, color = "grey30")) +
    labs(title = title, subtitle = subtitle)
}
# subtitle: the largest number of river crossings of any contour
crossing_note <- function(x) {
  sprintf("Up to %d river crossing%s per contour", max(x$n_crossings),
          if_else(max(x$n_crossings) == 1, "", "s"))
}
# hard-coded: canyon window (UTM 12N, m) for panels c and d
bb_c <- st_bbox(c(xmin = 321000, ymin = 4117000, xmax = 333000,
                  ymax = 4136000), crs = st_crs(UTM))
lon_c <- c(-112.98, -112.94, -112.90)   # longitude labels, panels c and d
fig_s6 <- (method_panel(cnt_idw, "a  Inverse-distance weighting",
                        crossing_note(xing_idw), bb_a) |
             method_panel(cnt, "b  Smoothing spline", crossing_note(xing),
                          bb_a)) /
  (method_panel(cnt_idw, "c  Inverse-distance weighting, Zion Canyon", NULL,
                bb_c, lon_c) |
     method_panel(cnt, "d  Smoothing spline, Zion Canyon", NULL, bb_c,
                  lon_c)) +
  plot_layout(guides = "collect", heights = c(1, 1.2)) &
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "figS6_contour_methods.png"), fig_s6, width = 7.2,
       height = 9.4, dpi = 450, bg = "white")
message("written ", file.path(FIG, "figS6_contour_methods.png"))
