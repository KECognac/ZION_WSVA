# 20_figure_3d_emergence.R -- 3D view of Zion Canyon showing
#   1. the spring emergence points picked from Google Earth imagery,
#   2. the elevation band in which the canyon-wall springs emerge, and
#   3. the reach where the North Fork Virgin River climbs into that band.
#
# How the band is defined
#   The band is the mean Google Earth elevation of the canyon-wall emergence
#   points (from 19_spring_emergence_band.R) plus or minus one standard
#   deviation. It is drawn wherever the land surface lies inside that range.
#   It is a level band, not a mapped water table: springs emerge near the base
#   of the Navajo Sandstone, and that surface is not level.
#
# How the river elevation is measured
#   Points are placed every 50 m along the NHDPlus river line. Each point takes
#   the lowest DEM cell within 50 m (the mapped line sometimes sits on the
#   canyon wall). Going upstream, the elevation is not allowed to decrease,
#   which removes DEM noise.
#
# Rendering
#   rayshader, through the helper functions in paper/R/map_3d.R (copied from
#   the CEBR_WSVA repository). A full render takes about two minutes, so this
#   script is not part of run_all.R.
#     Quick preview (default):   Rscript paper/scripts/20_figure_3d_emergence.R
#     Full-size figure:          set PREVIEW <- FALSE below, or run with the
#                                environment variable ZION_3D_FINAL=1
#
# Run after 09_site_geology_watershed.R and 19_spring_emergence_band.R.

options(rgl.useNULL = TRUE)   # render off-screen; must be set before rayshader
source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
})
source("paper/R/map_3d.R")
sf_use_s2(FALSE)


# ---- Settings ---------------------------------------------------------------

# TRUE = fast, half-size preview; FALSE = full-size path-traced figure
PREVIEW <- !identical(Sys.getenv("ZION_3D_FINAL"), "1")

# Terrain drape: "imagery" (Esri World Imagery tiles, downloaded by maptiles
# at zoom 15 for the model area) or "hillshade" (no download; uses the
# hypso_colors palette in the render call)
DRAPE <- "imagery"

RIVER_BUFFER_M <- 50     # search radius for the river elevation (m), as in 16
RIVER_STEP_M <- 50       # spacing of river elevation samples (m)
UPPER_CUT_M <- 1800      # hard-coded: emergence points above this elevation are
                         # the "upper" points, as in 19
CANYON_EAST_LON <- -112.92   # hard-coded: emergence points east of this lie on
                             # the east plateau, outside the model
AOI_BUFFER_M <- 1000     # margin around the springs for the model area (m)

# Colors
PINK   <- "gray90"#"#ff2e93"      # canyon-wall springs and their callouts
PURPLE <- "#8a2be2"      # upper emergence points
YELLOW <- "#ffd21f"      # river reach inside the band

# Thousands separator for labels, e.g. 1415 -> "1,415"
fmt_m <- function(x) scales::comma(round(x))


# ---- 1. Emergence points ----------------------------------------------------

# One row per emergence point, written by 19_spring_emergence_band.R
em <- read_csv(file.path(DER, "spring_emergence_points.csv"),
               show_col_types = FALSE) %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = UTM)

wall  <- em %>% filter(elev_ge_m < UPPER_CUT_M)    # springs on the canyon walls
upper <- em %>% filter(elev_ge_m >= UPPER_CUT_M)   # springs near the rim

# Sampled springs of the paper (used for the named callouts)
sites <- site_locations() %>%
  filter(Location %in% c(REPEAT_SPRINGS, "Birch Spring",
                         "Lower Grotto Sandtrap")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM)


# ---- 2. Emergence elevation band --------------------------------------------

# Mean +/- 1 SD of the Google Earth elevations of the canyon-wall points
band <- wall %>%
  st_drop_geometry() %>%
  summarise(n = n(), mid = mean(elev_ge_m), sd = sd(elev_ge_m),
            min = min(elev_ge_m), max = max(elev_ge_m)) %>%
  mutate(lo = mid - sd, hi = mid + sd)

message(str_glue(
  "emergence band: mean {round(band$mid)} m, SD {round(band$sd)} m ",
  "({round(band$lo)}-{round(band$hi)} m); ",
  "full range {round(band$min)}-{round(band$max)} m; n = {band$n}"
))


# ---- 3. Model area and elevation model --------------------------------------

# The canyon corridor: emergence points west of CANYON_EAST_LON and the
# sampled springs, plus a margin
aoi <- em %>%
  filter(longitude < CANYON_EAST_LON) %>%
  st_geometry() %>%
  c(st_geometry(sites)) %>%
  st_union() %>%
  st_buffer(AOI_BUFFER_M) %>%
  st_bbox() %>%
  st_as_sfc()

# Keep only the upper points that fall inside the model
upper <- upper %>% st_filter(aoi)

# 8-m elevation model (cached by 02_download_external.R), cropped to the model
dem <- rast(ext_path("dem_canyon_z14.tif")) %>%
  crop(vect(aoi))
names(dem) <- "elev"

# Land surface inside the band, as polygons. The comparison gives a raster of
# TRUE/FALSE; only the TRUE polygons are kept.
band_poly <- (dem >= band$lo & dem <= band$hi) %>%
  as.polygons() %>%
  st_as_sf() %>%
  rename(in_band = 1) %>%
  filter(in_band == 1) %>%
  st_make_valid()


# ---- 4. River: where it reaches the band ------------------------------------

# North Fork Virgin River inside the model, as one line (the longest piece)
river <- st_read(ext_path("nhdplus_flowlines.gpkg"), quiet = TRUE) %>%
  rename_with(tolower) %>%
  filter(gnis_name == "North Fork Virgin River") %>%
  st_transform(UTM) %>%
  st_union() %>%
  st_line_merge() %>%
  st_intersection(aoi) %>%
  st_cast("LINESTRING") %>%
  st_sf(geometry = .) %>%
  mutate(length_m = as.numeric(st_length(geometry))) %>%
  slice_max(length_m, n = 1)

# Sample points every RIVER_STEP_M along the river
river_pts <- river %>%
  st_line_sample(density = 1 / RIVER_STEP_M) %>%
  st_cast("POINT") %>%
  st_sf(geometry = .)

# River elevation at each point: lowest DEM cell within RIVER_BUFFER_M
river_min <- terra::extract(dem, vect(st_buffer(river_pts, RIVER_BUFFER_M)),
                            fun = min, na.rm = TRUE)

river_pts <- river_pts %>%
  mutate(z = river_min$elev,
         northing = st_coordinates(geometry)[, "Y"]) %>%
  # order upstream (south to north), then force the elevation to be
  # non-decreasing upstream to remove DEM noise
  arrange(northing) %>%
  mutate(z = cummax(z),
         km = (row_number() - 1) * RIVER_STEP_M / 1000,
         in_band = between(z, band$lo, band$hi))
stopifnot(any(river_pts$in_band), any(river_pts$z >= band$mid))

# The reach inside the band, as a line
reach <- river_pts %>%
  filter(in_band) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

# Three key points, going upstream: where the river enters the band, where it
# reaches the band mean, and the last point inside the band
key_pts <- bind_rows(
  enter = river_pts %>% filter(in_band) %>% slice_head(n = 1),
  mean  = river_pts %>% filter(z >= band$mid) %>% slice_head(n = 1),
  leave = river_pts %>% filter(in_band) %>% slice_tail(n = 1),
  .id = "where"
)
cross <- key_pts %>% filter(where == "mean")   # marker for the callout

# Latitude and longitude of the key points, for the log message and the table
key_ll <- key_pts %>%
  st_transform(4326) %>%
  mutate(lon = st_coordinates(geometry)[, "X"],
         lat = st_coordinates(geometry)[, "Y"]) %>%
  st_drop_geometry()
ll <- function(w, col) key_ll %>% filter(where == w) %>% pull({{ col }})

reach_km <- sum(river_pts$in_band) * RIVER_STEP_M / 1000
# does the band continue past the upstream (north) edge of the model?
leaves_at <- if_else(last(river_pts$in_band), "edge of the model",
                     "inside the model")

message(str_glue(
  "river enters the band at {round(ll('enter', lat), 4)} N, ",
  "{round(-ll('enter', lon), 4)} W; reaches the mean ({round(band$mid)} m) at ",
  "{round(ll('mean', lat), 4)} N, {round(-ll('mean', lon), 4)} W; leaves the ",
  "band at {round(ll('leave', lat), 4)} N ({leaves_at}); ",
  "reach in band {reach_km} km"
))

# Summary table for the text
tibble(
  band_mean_m = band$mid, band_sd_m = band$sd,
  band_lo_m = band$lo, band_hi_m = band$hi,
  river_enter_lat = ll("enter", lat), river_enter_lon = ll("enter", lon),
  river_mean_lat = ll("mean", lat), river_mean_lon = ll("mean", lon),
  river_reach_in_band_km = reach_km
) %>%
  write_csv(file.path(DER, "emergence_band_river.csv"))


# ---- 5. Navajo-Kayenta contact ----------------------------------------------

# Shared boundary of the mapped Navajo Sandstone and Kayenta Formation
# polygons (GRI geologic map), clipped to the model, as in 09 and 16
geo <- st_read(raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb"),
               "zionglg", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(UTM)

# hard-coded GRI map-unit symbols (same groups as 09 and 16)
navajo <- geo %>%
  filter(GLG_SYM %in% c("Jn", "Jnw", "Jnp", "Jnb", "Jnl")) %>%
  st_union()
kayenta <- geo %>%
  filter(GLG_SYM %in% c("Jk", "Jkt", "Jku", "Jks")) %>%
  st_union()

contact <- st_intersection(st_boundary(navajo), st_boundary(kayenta)) %>%
  st_collection_extract("LINESTRING") %>%
  st_intersection(aoi) %>%
  st_sf(geometry = .)


# ---- 6. Layers --------------------------------------------------------------

# Layers are drawn in the order listed (later layers on top). "fill", "line"
# and "square" layers are painted onto the terrain. "dot" layers are drawn on
# the finished image, so the markers stay round on the cliffs.

# One sampled spring, by name
one <- function(name) sites %>% filter(Location == name)

# A named spring: one marker plus a callout box. x and y give the center of
# the box as fractions of the image width and height, from the top left.
named <- function(pt, label, x, y) {
  map3d_layer(pt, label, "dot", color = PINK, outline = "white", legend = FALSE,
              callout = map3d_callout(label, x, y, color = PINK))
}

band_name <- str_glue(
  "Emergence elevation band ({fmt_m(band$mid)} ± {round(band$sd)} m)"
)

layers <- list(
  # the band, the contact and the river, painted on the terrain
  map3d_layer(band_poly, band_name, "fill", color = "#2f8fe0", width = 1,
              outline = NA, fill_alpha = 0.5),
  map3d_layer(contact, "Navajo–Kayenta contact (mapped)", "line",
              color = "#8c1c13", width = 3),
  map3d_layer(river, "North Fork Virgin River", "line", color = "#0b3c8c",
              width = 5, casing = "white"),
  map3d_layer(reach, "River within the emergence band", "line", color = YELLOW,
              width = 9, casing = "black"),

  # where the river reaches the mean emergence elevation
  map3d_layer(cross, "River reaches the band mean", "square", color = YELLOW,
              size_m = 110, outline = "black", legend = FALSE,
              callout = map3d_callout(
                str_glue("River reaches the mean\n",
                         "emergence elevation ({fmt_m(band$mid)} m)"),
                0.74, 0.20)),

  # all emergence points
  map3d_layer(wall, "Canyon-wall emergence point", "dot", color = PINK,
              outline = "white"),
  map3d_layer(upper, "Upper emergence point", "dot", color = PURPLE,
              outline = "white"),

  # named springs (callout box positions are tuned for the view below)
  named(one("Birch Spring"), "Birch Spring", 0.20, 0.70),
  named(one("The Grotto"), "The Grotto", 0.72, 0.62),
  named(one("Weeping Rock"), "Weeping Rock", 0.74, 0.44),
  # Refrigerator Spring's emergence point is beside the Temple of Sinawava
  named(one("Refrigerator Spring"),
        "Refrigerator Spring\n(near Temple of Sinawava)", 0.18, 0.40)
)


# ---- 7. Camera and render ---------------------------------------------------

# from: side the camera sits on; phi: angle above the horizon (90 = straight
# down); zoom: larger values move the camera out; fov = 0: no perspective.
# The image is wider than the model so the legend sits to its right.
views <- list(south = map3d_view(from = "S", phi = 52, zoom = 0.80, fov = 0))

# Notes printed under the image
notes <- c(
  str_glue(
    "Band: mean ± 1 standard deviation of the Google Earth elevations of ",
    "{band$n} canyon-wall emergence points (range {fmt_m(band$min)}–",
    "{fmt_m(band$max)} m), drawn where the land surface lies within it."
  ),
  str_c(
    "Elevation: USGS 3D Elevation Program, 8-m grid. ",
    if_else(DRAPE == "imagery", "Imagery: Esri World Imagery.",
            "Shaded relief from the elevation model.")
  )
)

render_map_3d(
  aoi = aoi, layers = layers, views = views,
  # terrain
  dem = dem, crop_to_aoi = FALSE,
  dem_agg = 2,          # average 2 x 2 cells (16-m grid) to speed up rendering
  z_exag = 1,           # vertical exaggeration
  drape = DRAPE, texture_zoom = 15,
  # hard-coded sandstone palette, low to high (used only with "hillshade")
  hypso_colors = c("#9c4a32", "#c97a52", "#e2a877", "#efd0a8", "#f8ecd8",
                   "seagreen"),
  # output
  preview = PREVIEW, out_dir = FIG, prefix = "fig_3d_zion_emergence",
  width = 3400, height = 2000, samples = 256,
  # map furniture; positions are fractions of the image from the top left
  elev_bar = NULL,      # the helper labels the elevation bar in feet
  scalebar = list(km = 2, edge = "S", start = 0.08, inset = 0.04,
                  shift = c(0, 0)),
  north_arrow = c(0.05, 0.86),
  legend_corner = c(0.995, 0.97),   # bottom-right corner of the legend box
  text_cex = 2.3,       # smaller text keeps the legend clear of the model
  notes = notes
)
message("3D figure written to ", FIG)
