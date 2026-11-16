# 03_tritium.R -- tritium thresholds, classification and young-water
# fractions (paper sections 2.3, 3.1, 4.3).
#
# Method (Lindsey et al. 2019; reproduces the "Moving modern & premodern
# thresholds" chunk of Zion_Isotope_Analysis.Rmd)
#   1. GNIP precipitation tritium (IAEA/WMO WISER) is averaged to half-years
#      at each of the four stations.
#   2. Modern threshold on a given date = mean across stations of each
#      station's minimum post-bomb-peak value, decayed to that date.
#   3. Premodern threshold on a given date = lowest station mean for
#      2004-2012, decayed from 1 January 1953 to that date.
#   4. Each sample is Modern (at or above the modern threshold), Premodern (at
#      or below the premodern threshold) or Mixed.
#   5. Young fraction (primary estimate) = (C - Cpre) / (Cmod - Cpre), limited
#      to 0-1.
# Robustness (section 4.3 only): the modern endmember is replaced by the mean
# across stations of each station's own bomb-peak value, decayed to the sample
# date.
#
# Outputs
#   derived: tritium_samples.csv, tritium_thresholds_by_period.csv,
#            tritium_thresholds_1988-11-24.csv, gnip_halfyear.csv,
#            tritium_threshold_curve.csv
#   tables:  table2_tritium_by_site.csv, tableS1_tritium_conductivity.csv

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages({
  library(readxl)
  library(lubridate)
})

# Years elapsed from date b to date a
yrs <- function(a, b) as.numeric(a - b) / 365.25

# hard-coded: start of bomb tritium in precipitation, the reference date for
# the premodern threshold (Lindsey et al. 2019)
BOMB_START <- as.Date("1953-01-01")


# ---- 1. GNIP precipitation tritium: half-year means per station -------------

gnip <- list.files(raw_path("Tritium/GNIP WISER"), pattern = "\\.xlsx$",
                   full.names = TRUE) %>%
  map_dfr(~ suppressWarnings(read_xlsx(.x))) %>%
  janitor::clean_names() %>%
  filter(measurand_symbol == "H3") %>%
  # half-year bin: 1 Jan for January-June, 1 Jul for July-December
  mutate(half = floor_date(as_date(sample_date), "halfyear")) %>%
  group_by(station = sample_site_name, half) %>%
  summarise(c = mean(measurand_amount, na.rm = TRUE), .groups = "drop")

# Bomb peak at each station: the highest half-year value and its date
peaks <- gnip %>%
  group_by(station) %>%
  slice_max(c, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(station, peak_date = half, peak_c = c)

# Flag the half-years after each station's bomb peak
gnip <- gnip %>%
  left_join(peaks, by = "station") %>%
  mutate(post_peak = half > peak_date)

# Mean 2004-2012 concentration at each station (hard-coded reference window
# for the premodern threshold, Lindsey et al. 2019). Only Albuquerque has data
# in this window.
ref_0412 <- gnip %>%
  filter(between(year(half), 2004, 2012)) %>%
  group_by(station) %>%
  summarise(c = mean(c), .groups = "drop")


# ---- 2-3. Thresholds as of a given date -------------------------------------

# Each station's minimum post-bomb-peak value, decayed to date d (one row per
# station, column m)
modern_by_station <- function(d) {
  gnip %>%
    filter(post_peak, half <= d) %>%
    group_by(station) %>%
    summarise(m = min(decay(c, yrs(d, half), H3_HALF_LIFE)), .groups = "drop")
}

# All three thresholds on date d, as a one-row tibble
thresholds <- function(d) {
  tibble(
    # mean across stations of the decayed post-peak minimum
    modern = mean(modern_by_station(d)$m),
    # lowest station 2004-2012 mean, decayed from 1953
    premodern = min(decay(ref_0412$c, yrs(d, BOMB_START), H3_HALF_LIFE)),
    # robustness endmember: mean across stations of the decayed bomb peak
    modern_bombpeak = peaks %>%
      mutate(v = decay(peak_c, yrs(d, peak_date), H3_HALF_LIFE)) %>%
      pull(v) %>%
      mean()
  )
}


# ---- 4-5. Samples: classification and young fraction -----------------------

tr <- read_tritium_raw() %>%
  mutate(date = mdy(Date), row = row_number())

# Thresholds on each sample's own date
tr <- tr %>%
  bind_cols(map_dfr(tr$date, thresholds)) %>%
  mutate(
    f_raw = (Tritium - premodern) / (modern - premodern),
    yf = pmin(pmax(f_raw, 0), 1),                     # limited to 0-1
    yf_bombpeak = pmin(
      pmax((Tritium - premodern) / (modern_bombpeak - premodern), 0), 1
    ),
    cls = case_when(
      is.na(Tritium) ~ NA_character_,
      Tritium >= modern ~ "Modern",
      Tritium <= premodern ~ "Premodern",
      TRUE ~ "Mixed"
    ),
    # hard-coded: decay-corrected to 12 Aug 2025, the latest sample date in
    # Zion_Tritium.csv
    tritium_dc = decay(Tritium, yrs(as.Date("2025-08-12"), date), H3_HALF_LIFE),
    period = if_else(year(date) < 2000, "1988-89", "2025")
  )
write_csv(tr, file.path(DER, "tritium_samples.csv"))


# ---- Table 2: tritium results by site ---------------------------------------

site_order <- c(REPEAT_SPRINGS, "Birch Spring", "East Entrance Well",
                "Lower Grotto Sandtrap", KOLOB)

t2 <- tr %>%
  filter(!is.na(Tritium)) %>%
  group_by(Location) %>%
  summarise(
    # first and last sample year, e.g. "1988-1989" (one year if the same)
    sampled = str_c(unique(range(format(date, "%Y"))), collapse = "-"),
    n = n(),
    tritium_min = min(Tritium),
    tritium_max = max(Tritium),
    n_modern = sum(cls == "Modern"),
    n_mixed = sum(cls == "Mixed"),
    n_premodern = sum(cls == "Premodern"),
    yf_mean_pct = 100 * mean(yf),
    yf_min_pct = 100 * min(yf),
    yf_max_pct = 100 * max(yf),
    .groups = "drop"
  ) %>%
  arrange(factor(Location, levels = site_order))
write_csv(t2, file.path(TAB, "table2_tritium_by_site.csv"))

# Range of the thresholds and of the bomb-peak endmember across the sample
# dates of each campaign (sections 2.3 and 4.3)
summ <- tr %>%
  filter(!is.na(Tritium)) %>%
  group_by(period) %>%
  summarise(
    modern_min = min(modern),
    modern_max = max(modern),
    premodern_min = min(premodern),
    premodern_max = max(premodern),
    bombpeak_min = min(modern_bombpeak),
    bombpeak_max = max(modern_bombpeak),
    .groups = "drop"
  )
write_csv(summ, file.path(DER, "tritium_thresholds_by_period.csv"))

print(as.data.frame(t2), digits = 3)
print(as.data.frame(summ), digits = 4)
class_counts <- tr %>% filter(!is.na(cls)) %>% count(cls)
message("class counts: ",
        str_c(class_counts$cls, class_counts$n, sep = " ", collapse = "; "))

# Lowest Weeping Rock young fraction under the bomb-peak endmember (4.3)
wr_min_bp <- tr %>%
  filter(Location == "Weeping Rock", !is.na(Tritium)) %>%
  summarise(v = min(yf_bombpeak)) %>%
  pull(v)
message(sprintf(
  "Weeping Rock minimum young fraction with bomb-peak endmember: %.1f%%",
  100 * wr_min_bp
))


# ---- Regression checks against the values in the paper ---------------------

# hard-coded reference values: site-mean young fractions (%) from Table 2
ref <- tribble(
  ~Location,               ~ref_pct,
  "Echo Canyon Spring",    100,
  "Weeping Rock",          73,
  "The Grotto",            48,
  "Refrigerator Spring",   33,
  "Larry's Grotto",        10,
  "Birch Spring",          0,
  "East Entrance Well",    0.8,
  "Lower Grotto Sandtrap", 0.4,
  "Taylor Creek Well",     100
) %>%
  left_join(t2 %>% select(Location, yf_mean_pct), by = "Location") %>%
  # values below 1% are reported to one decimal, the rest to whole percent
  mutate(new = round(yf_mean_pct, if_else(ref_pct < 1, 1, 0)))
check_close(ref$new, ref$ref_pct, 0, "site-mean young fractions (Table 2)")

# class counts: 4 Modern, 18 Mixed, 2 Premodern
stopifnot(identical(
  class_counts %>% arrange(factor(cls, c("Modern", "Mixed", "Premodern"))) %>%
    pull(n),
  c(4L, 18L, 2L)
))

s88 <- summ %>% filter(period == "1988-89")
s25 <- summ %>% filter(period == "2025")
check_close(c(s88$modern_min, s88$modern_max,
              s88$premodern_min, s88$premodern_max),
            c(9.83, 10.12, 0.99, 1.02), 0.005, "1988-89 thresholds")
check_close(c(s25$modern_min, s25$premodern_min),
            c(1.2626, 0.1270), 0.0001, "2025 thresholds")
check_close(c(s88$bombpeak_min, s88$bombpeak_max, s25$bombpeak_min),
            c(622, 639, 79.8), 0.5, "bomb-peak endmember")
check_close(100 * wr_min_bp, 0.9, 0.05, "Weeping Rock bomb-peak minimum")


# ---- Outputs for later scripts ----------------------------------------------

# Thresholds on 24 Nov 1988 (hard-coded: middle of the 1988-89 campaign),
# used by ESM Fig. S3
thresholds(as.Date("1988-11-24")) %>%
  write_csv(file.path(DER, "tritium_thresholds_1988-11-24.csv"))

# Half-year GNIP series and the moving thresholds through time, used by
# Fig. 2 (18_figure_tritium_precip.R). The modern threshold is the mean
# across stations, with the station-to-station range; quarterly from
# hard-coded 1 Jan 1966 (after the last station's bomb peak) to the latest
# sample date.
gnip %>%
  select(station, half, c, post_peak) %>%
  write_csv(file.path(DER, "gnip_halfyear.csv"))

curve <- tibble(date = seq(as.Date("1966-01-01"), max(tr$date),
                           by = "3 months")) %>%
  mutate(
    st = map(date, modern_by_station),     # one small table per date
    modern = map_dbl(st, ~ mean(.x$m)),
    modern_min = map_dbl(st, ~ min(.x$m)),
    modern_max = map_dbl(st, ~ max(.x$m)),
    n_station = map_int(st, nrow),
    premodern = map_dbl(date, ~ thresholds(.x)$premodern)
  ) %>%
  select(-st)
write_csv(curve, file.path(DER, "tritium_threshold_curve.csv"))

# ESM Table S1: every sample (tritium, decay-corrected tritium, class,
# primary young fraction, conductivity)
tr %>%
  select(Location, Date, period, Tritium, tritium_dc, cls, yf, Cond) %>%
  write_csv(file.path(TAB, "tableS1_tritium_conductivity.csv"))
