# 02_download_external.R -- fetch and cache every input that comes from a web
# service, so later scripts never call the network.
#
# Each item is downloaded only if its cached file is missing; delete a file in
# data/park/ZION/paper/data/external/ to refresh it. Results can change if a
# service updates its data (for example new NWIS water levels), so the cache
# is the record of what the paper used. The initial cache was seeded from the
# October 2026 analysis run.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dataRetrieval)
  library(nhdplusTools)
  library(elevatr)
})
sf_use_s2(FALSE)

# Run `fun(path)` only if the cached file does not exist yet.
fetch <- function(file, fun) {
  p <- file.path(EXT, file)
  if (!file.exists(p)) {
    message("fetching ", file)
    fun(p)
  } else {
    message("cached   ", file)
  }
  invisible(p)
}


# ---- NPS administrative boundaries (Zion, Cedar Breaks) ---------------------

fetch("nps_boundary/nps_boundary.shp", function(p) {
  # hard-coded: NPS IRMA DataStore reference 733895 (NPS unit boundaries)
  z <- tempfile(fileext = ".zip")
  download.file("https://irma.nps.gov/DataStore/DownloadFile/733895", z,
                mode = "wb", quiet = TRUE)
  d <- tempfile()
  unzip(z, exdir = d)
  b <- list.files(d, "nps_boundary.shp$", recursive = TRUE,
                  full.names = TRUE)[1] %>%
    st_read(quiet = TRUE) %>%
    filter(UNIT_CODE %in% c("ZION", "CEBR"))
  dir.create(dirname(p), showWarnings = FALSE)
  st_write(b, p, quiet = TRUE, delete_dsn = TRUE)
})


# ---- USGS streamgage 09405501 and its NLDI drainage basin -------------------

fetch("nwis_gage_09405501.csv", function(p) {
  write_csv(readNWISsite("09405501"), p)
})
fetch("gage_09405501_basin.gpkg", function(p) {
  basin <- get_nldi_basin(list(featureSource = "nwissite",
                               featureID = "USGS-09405501"))
  st_write(basin, p, quiet = TRUE, delete_dsn = TRUE)
})


# ---- Watershed Boundary Dataset units containing the sampling sites ---------

fetch("huc_units.gpkg", function(p) {
  tr <- read_tritium_raw() %>%
    distinct(Location, latitude, longitude)
  pts <- st_as_sf(tr, coords = c("longitude", "latitude"), crs = 4326)

  # HUC-8, HUC-10 and HUC-12; the ID column name differs by level
  h <- bind_rows(lapply(c("huc08", "huc10", "huc12"), function(t) {
    u <- get_huc(AOI = pts, type = t)
    id_col <- intersect(c("huc8", "huc10", "huc12"), names(u))[1]
    u %>%
      transmute(huc = .data[[id_col]], name, areasqkm, level = t) %>%
      distinct(huc, .keep_all = TRUE)
  }))
  st_write(h, p, quiet = TRUE, delete_dsn = TRUE)
})


# ---- NHDPlus flowlines: gage basin plus an 8-km margin covering all sites ---

fetch("nhdplus_flowlines.gpkg", function(p) {
  basin <- st_read(file.path(EXT, "gage_09405501_basin.gpkg"), quiet = TRUE)
  aoi <- st_buffer(st_transform(basin, UTM), 8000) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_transform(4326)
  st_write(get_nhdplus(AOI = aoi, realization = "flowline"), p,
           quiet = TRUE, delete_dsn = TRUE)
})


# ---- NWIS wells completed in the Navajo Sandstone within 30 km of the park --

fetch("nwis_navajo_wells.csv", function(p) {
  park <- st_read(file.path(EXT, "nps_boundary/nps_boundary.shp"),
                  quiet = TRUE) %>%
    filter(UNIT_CODE == "ZION") %>%
    st_transform(UTM)
  bb <- st_bbox(st_transform(st_buffer(park, 30000), 4269))

  inv <- whatNWISdata(bBox = round(as.numeric(bb), 6))
  gw <- inv %>%
    filter(data_type_cd == "gw") %>%
    distinct(site_no, .keep_all = TRUE)
  meta <- readNWISsite(unique(gw$site_no)) %>%
    select(site_no, station_nm, dec_lat_va, dec_long_va, alt_va, alt_datum_cd,
           aqfr_cd, aqfr_type_cd, well_depth_va)
  # hard-coded: NWIS aquifer code 220NVJO = Navajo Sandstone
  nav <- meta %>% filter(aqfr_cd == "220NVJO")

  # depth to water below land surface (parameter 72019), summarized per well
  lev <- read_waterdata_field_measurements(
    monitoring_location_id = paste0("USGS-", nav$site_no),
    parameter_code = "72019"
  ) %>%
    st_drop_geometry() %>%
    mutate(site_no = sub("USGS-", "", monitoring_location_id)) %>%
    group_by(site_no) %>%
    summarise(
      n_meas = n(),
      dtw_ft_mean = mean(value, na.rm = TRUE),
      dtw_ft_min = min(value, na.rm = TRUE),
      dtw_ft_max = max(value, na.rm = TRUE),
      first = min(time),
      last = max(time),
      .groups = "drop"
    )

  nav %>%
    left_join(lev, by = "site_no") %>%
    mutate(gw_alt_ft = as.numeric(alt_va) - dtw_ft_mean) %>%
    write_csv(p)
})


# ---- Digital elevation models (AWS terrain tiles via elevatr) ---------------

# Response-time and recharge-zone DEM: geologic map extent plus 40 km around
# the springs, zoom 9 (127-m cells after reprojection)
fetch("dem_z9.tif", function(p) {
  gdb <- raw_path("geology/zion_geology_gdb_pro/zion_geology.gdb")
  geo <- st_read(gdb, "zionglg", quiet = TRUE) %>% st_transform(UTM)
  sites <- read_tritium_raw() %>%
    distinct(Location, latitude, longitude) %>%
    filter(Location %in% REPEAT_SPRINGS) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM)
  aoi <- st_union(st_union(geo), st_buffer(st_union(sites), 40000)) %>%
    st_as_sf()
  r <- get_elev_raster(st_transform(aoi, 4326), z = 9, src = "aws",
                       clip = "bbox") %>%
    rast() %>%
    project(paste0("EPSG:", UTM))
  names(r) <- "elev_m"
  writeRaster(r, p, overwrite = TRUE)
})

# Locates Brian Head Peak: highest cell within 8 km of 37.68 N, 112.83 W
# (hard-coded approximate peak location), zoom 11
fetch("dem_brian_head_z11.tif", function(p) {
  a <- st_sfc(st_point(c(-112.83, 37.68)), crs = 4326) %>% st_transform(UTM)
  r <- st_buffer(a, 8000) %>%
    st_transform(4326) %>%
    st_as_sf() %>%
    get_elev_raster(z = 11, src = "aws", clip = "locations") %>%
    rast() %>%
    project(paste0("EPSG:", UTM))
  writeRaster(r, p, overwrite = TRUE)
})

# Spring outlet position and catchments (script 16): zoom 14 (3DEP-sourced)
# over the canyon springs plus 8 km, reprojected to UTM 12N at 8 m
fetch("dem_canyon_z14.tif", function(p) {
  s <- site_locations() %>%
    filter(Location %in% c(REPEAT_SPRINGS, "Lower Grotto Sandtrap",
                           "Birch Spring")) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM)
  aoi <- st_as_sfc(st_bbox(st_buffer(st_union(s), 8000))) %>% st_sf()
  r <- get_elev_raster(st_transform(aoi, 4326), z = 14, src = "aws",
                       clip = "bbox") %>%
    rast()
  r <- project(r, paste0("EPSG:", UTM), res = 8, method = "bilinear") %>%
    crop(vect(aoi))
  names(r) <- "elev_m"
  writeRaster(r, p, overwrite = TRUE, gdal = "COMPRESS=LZW")
})

# Land-surface elevation at each sampling site, zoom 14
fetch("site_elevations_z14.csv", function(p) {
  s <- site_locations() %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  get_elev_point(s, src = "aws", z = 14) %>%
    st_drop_geometry() %>%
    select(Location, latitude, longitude, elev_m = elevation) %>%
    write_csv(p)
})
# The cache must hold the current site coordinates (site_locations() in
# paper/R/functions.R); delete the file to download it again if a site moves
local({
  cached <- read_csv(ext_path("site_elevations_z14.csv"), show_col_types = FALSE)
  now <- site_locations()
  m <- match(now$Location, cached$Location)
  moved <- is.na(m) | abs(cached$latitude[m] - now$latitude) > 1e-6 |
    abs(cached$longitude[m] - now$longitude) > 1e-6
  if (any(moved)) {
    stop("site_elevations_z14.csv is out of date for: ",
         paste(now$Location[moved], collapse = ", "))
  }
})


# ---- Base maps for the location inset ---------------------------------------

fetch("basemap_states.gpkg", function(p) {
  # hard-coded: western states shown in the Fig. 1 inset
  st <- tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    filter(STUSPS %in% c("UT", "NV", "AZ", "CO", "NM", "ID", "WY", "CA", "OR"))
  st_write(st, p, quiet = TRUE, delete_dsn = TRUE)
})

fetch("basemap_mexico.gpkg", function(p) {
  # Natural Earth 1:110m admin-0 countries
  z <- tempfile(fileext = ".zip")
  download.file(
    "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip",
    z, mode = "wb", quiet = TRUE
  )
  d <- tempfile()
  unzip(z, exdir = d)
  mex <- st_read(list.files(d, "\\.shp$", full.names = TRUE)[1],
                 quiet = TRUE) %>%
    filter(ADMIN == "Mexico")
  st_write(mex, p, quiet = TRUE, delete_dsn = TRUE)
})


# ---- USGS tritium in precipitation by quadrangle ----------------------------
# Jurgens (2018), ver. 2.0, doi:10.5066/P92CEFXN. The URL is the ScienceBase
# file link for Table 2.

fetch("usgs_tritium_quadrangles/Table_2_TritiumInQuadrangles.csv", function(p) {
  dir.create(dirname(p), showWarnings = FALSE)
  download.file(
    paste0("https://www.sciencebase.gov/catalog/file/get/",
           "641b2a9bd34e4c337f0c246a?f=__disk__85%2F86%2F0d%2F",
           "85860dbc2b0c4ff65a359aaf3f34e58ab11038f3"),
    p, mode = "wb", quiet = TRUE
  )
})


# ---- Long-term streamflow: USGS 09405500 (daily) ----------------------------
# North Fork Virgin River near Springdale. The record is cut at the end of
# water year 2025 (hard-coded end date) so results do not shift as new
# provisional data arrive.

fetch("nwis_dv_09405500.csv", function(p) {
  readNWISdv("09405500", "00060",            # 00060 = discharge, ft3/s
             startDate = "1923-01-01", endDate = "2025-09-30") %>%
    renameNWISColumns() %>%
    write_csv(p)
})
fetch("nwis_gage_09405500.csv", function(p) {
  write_csv(readNWISsite("09405500"), p)
})

# Springdale Canal (09405499) diverts at river mile 4.9, just above 09405500
# (mile 4.7). For WY1969-1989 the published 09405500 record equals river +
# canal; 09405501 is the river alone (1968-1991). Christensen et al. (2005)
# report the combined flow. Both are used in 14 to build a consistent series.
fetch("nwis_dv_09405501.csv", function(p) {
  readNWISdv("09405501", "00060") %>%
    renameNWISColumns() %>%
    write_csv(p)
})
fetch("nwis_dv_09405499.csv", function(p) {
  readNWISdv("09405499", "00060") %>%
    renameNWISColumns() %>%
    write_csv(p)
})
fetch("gage_09405500_basin.gpkg", function(p) {
  basin <- get_nldi_basin(list(featureSource = "nwissite",
                               featureID = "USGS-09405500"))
  st_write(basin, p, quiet = TRUE, delete_dsn = TRUE)
})


# ---- GridMET daily climate averaged over the 09405500 basin -----------------
# pr (mm), pet (grass reference ET, mm), tmmx/tmmn (K). Area-weighted basin
# mean of the 4-km cells (exactextractr), downloaded one year at a time
# through climateR, 1979 through the end of water year 2025.

fetch("gridmet_basin_09405500_daily.csv", function(p) {
  b <- st_read(file.path(EXT, "gage_09405500_basin.gpkg"), quiet = TRUE) %>%
    st_transform(4326)
  # 5-km margin so every edge cell is returned
  aoi <- st_buffer(st_transform(b, UTM), 5000) %>% st_transform(4326)

  # One variable at a time: GridMET names both temperatures "air_temperature",
  # so a multi-variable request could not be matched up by name
  one <- function(v, s, e) {
    r <- climateR::getGridMET(AOI = aoi, varname = v,
                              startDate = s, endDate = e)[[1]]
    tibble(
      date = as.Date(time(r)),
      value = as.numeric(exactextractr::exact_extract(r, b, "mean",
                                                      progress = FALSE))
    )
  }

  out <- map_dfr(1979:2025, function(y) {
    s <- sprintf("%d-01-01", y)
    e <- if (y == 2025) "2025-09-30" else sprintf("%d-12-31", y)
    v <- lapply(c(pr = "pr", pet = "pet", tmmx = "tmmx", tmmn = "tmmn"),
                one, s = s, e = e)
    stopifnot(all(sapply(v, function(d) identical(d$date, v$pr$date))))
    tibble(
      date = v$pr$date,
      pr_mm = v$pr$value,
      pet_mm = v$pet$value,
      tmax_c = v$tmmx$value - 273.15,
      tmin_c = v$tmmn$value - 273.15
    )
  })
  write_csv(out, p)
})


# ---- Cuthbert et al. (2019) water table ratio and response time rasters -----
# figshare doi:10.6084/m9.figshare.7393304 (CC BY 4.0). Global 1-km Eckert IV
# GeoTIFFs (349 MB to 2 GB each), internally tiled 128 x 128, so only the
# window covering the North Fork Virgin gage basin + 20 km and the Zion NP
# boundary is read over HTTP. Cropped on the native grid (no resampling), so
# cell values are exactly as published.

# hard-coded: figshare file IDs for each raster
CUTHBERT_FILES <- c(
  Log_GRT_a     = "13677749",
  LOG_WTR_NL_01 = "13677845",
  LOG_WTR_L_01  = "13677830",
  L01_m         = "28402347"
)
CUTHBERT_BUFFER_M <- 20000

for (nm in names(CUTHBERT_FILES)) {
  fetch(file.path("cuthbert2019", paste0(nm, "_zion_clip.tif")), function(p) {
    r <- rast(paste0("/vsicurl/https://ndownloader.figshare.com/files/",
                     CUTHBERT_FILES[[nm]]))
    basin <- st_read(file.path(EXT, "gage_09405501_basin.gpkg"),
                     quiet = TRUE) %>%
      st_transform(UTM)
    park <- st_read(file.path(EXT, "nps_boundary/nps_boundary.shp"),
                    quiet = TRUE) %>%
      filter(UNIT_CODE == "ZION") %>%
      st_transform(UTM)
    aoi <- st_union(st_union(st_buffer(basin, CUTHBERT_BUFFER_M)),
                    st_union(park)) %>%
      st_transform(crs(r))
    clip <- crop(r, ext(vect(aoi)), snap = "out")
    dir.create(dirname(p), showWarnings = FALSE)
    writeRaster(clip, p, overwrite = TRUE, gdal = "COMPRESS=LZW")
  })
}

message("external cache complete: ", EXT)
