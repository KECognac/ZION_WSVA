# 03_tritium.R -- tritium thresholds, classification and young-water
# fractions (paper sections 2.3, 3.1, 4.3).
#
# Thresholds follow Lindsey et al. (2019) and reproduce the "Moving modern &
# premodern thresholds" chunk of Zion_Isotope_Analysis.Rmd:
#   - GNIP precipitation (IAEA/WMO WISER) averaged to half-years
#   - modern threshold = mean across stations of each station's minimum
#     post-bomb-peak value, decayed to the sample date
#   - premodern threshold = lowest station 2004-2012 mean, decayed from
#     1 January 1953 to the sample date
# Young fraction (primary estimate) = (C - Cpre) / (Cmod - Cpre), capped to
# 0-1.
# Robustness (section 4.3 only): modern endmember = mean across stations of
# each station's own bomb-peak value decayed to the sample date.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages(library(readxl))

# Years between two dates
yrs <- function(a, b) as.numeric(a - b) / 365.25


# ---- GNIP precipitation tritium: half-year means per station ----------------

gnip <- list.files(raw_path("Tritium/GNIP WISER"), pattern = "\\.xlsx$",
                   full.names = TRUE) %>%
  map_dfr(~ suppressWarnings(read_xlsx(.x))) %>%
  janitor::clean_names() %>%
  filter(measurand_symbol == "H3") %>%
  mutate(
    sample_date = as.Date(sample_date),
    # half-year bin: 1 Jan for Jan-Jun, 1 Jul for Jul-Dec
    half = as.Date(paste0(
      format(sample_date, "%Y"),
      ifelse(as.numeric(format(sample_date, "%m")) <= 6, "-01-01", "-07-01")
    ))
  ) %>%
  group_by(station = sample_site_name, half) %>%
  summarise(c = mean(measurand_amount, na.rm = TRUE), .groups = "drop")

# Bomb peak (highest half-year value) at each station
peaks <- gnip %>%
  group_by(station) %>%
  slice_max(c, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(station, peak_date = half, peak_c = c)

gnip <- gnip %>%
  left_join(peaks, by = "station") %>%
  mutate(post_peak = half > peak_date)

# hard-coded: 2004-2012 reference window for the premodern threshold
# (Lindsey et al. 2019)
ref_0412 <- gnip %>%
  filter(between(as.numeric(format(half, "%Y")), 2004, 2012)) %>%
  group_by(station) %>%
  summarise(c = mean(c), .groups = "drop")


# ---- Thresholds as of a given date ------------------------------------------

thresholds <- function(d) {
  mod <- gnip %>%
    filter(post_peak, half <= d) %>%
    group_by(station) %>%
    summarise(m = min(decay(c, yrs(d, half), H3_HALF_LIFE)), .groups = "drop") %>%
    pull(m) %>%
    mean()
  # hard-coded: premodern reference date 1 Jan 1953 (start of bomb tritium)
  pre <- min(decay(ref_0412$c, yrs(d, as.Date("1953-01-01")), H3_HALF_LIFE))
  peak <- peaks %>%
    mutate(v = decay(peak_c, yrs(d, peak_date), H3_HALF_LIFE)) %>%
    pull(v) %>%
    mean()
  c(modern = mod, premodern = pre, modern_bombpeak = peak)
}


# ---- Samples: classification and young fraction -----------------------------

tr <- read_csv(raw_path("Tritium/Park_Samples/Zion_Tritium.csv"),
               show_col_types = FALSE) %>%
  mutate(date = as.Date(Date, format = "%m/%d/%y"), row = row_number())

th <- t(sapply(tr$date, thresholds)) %>% as_tibble()

tr <- bind_cols(tr, th) %>%
  mutate(
    f_raw = (Tritium - premodern) / (modern - premodern),
    yf = pmin(pmax(f_raw, 0), 1),
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
    period = ifelse(as.numeric(format(date, "%Y")) < 2000, "1988-89", "2025")
  )
write_csv(tr, file.path(DER, "tritium_samples.csv"))


# ---- Table 2: tritium results by site ---------------------------------------

site_order <- c(REPEAT_SPRINGS, "Birch Sediment Trap", "East Entrance Well",
                "Lower Grotto Sandtrap", KOLOB)

t2 <- tr %>%
  filter(!is.na(Tritium)) %>%
  group_by(Location) %>%
  summarise(
    sampled = paste(unique(range(format(date, "%Y"))), collapse = "-"),
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
  arrange(match(Location, site_order))
write_csv(t2, file.path(TAB, "table2_tritium_by_site.csv"))

# Thresholds and bomb-peak endmember by campaign (sections 2.3 and 4.3)
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
message("class counts: ",
        paste(names(table(tr$cls)), table(tr$cls), collapse = "; "))

wr_min_bp <- tr %>%
  filter(Location == "Weeping Rock", !is.na(Tritium)) %>%
  summarise(min(yf_bombpeak)) %>%
  pull()
message(sprintf(
  "Weeping Rock minimum young fraction with bomb-peak endmember: %.1f%%",
  100 * wr_min_bp
))


# ---- Regression checks against the October 2026 values in the paper --------

# hard-coded reference values: site-mean young fractions (%) from Table 2
ref <- c(
  "Echo Canyon Spring" = 100, "Weeping Rock" = 73, "The Grotto" = 48,
  "Refrigerator Spring" = 33, "Larry's Grotto" = 10,
  "Birch Sediment Trap" = 0, "East Entrance Well" = 0.8,
  "Lower Grotto Sandtrap" = 0.4, "Taylor Creek Well" = 100
)
check_close(round(t2$yf_mean_pct[match(names(ref), t2$Location)],
                  ifelse(ref < 1, 1, 0)),
            ref, 0, "site-mean young fractions (Table 2)")
# class counts: 4 Modern, 18 Mixed, 2 Premodern
stopifnot(identical(
  as.integer(table(tr$cls)[c("Modern", "Mixed", "Premodern")]),
  c(4L, 18L, 2L)
))
check_close(c(summ$modern_min[1], summ$modern_max[1],
              summ$premodern_min[1], summ$premodern_max[1]),
            c(9.83, 10.12, 0.99, 1.02), 0.005, "1988-89 thresholds")
check_close(c(summ$modern_min[2], summ$premodern_min[2]),
            c(1.2626, 0.1270), 0.0001, "2025 thresholds")
check_close(c(summ$bombpeak_min[1], summ$bombpeak_max[1],
              summ$bombpeak_min[2]),
            c(622, 639, 79.8), 0.5, "bomb-peak endmember")
check_close(100 * wr_min_bp, 0.9, 0.05, "Weeping Rock bomb-peak minimum")


# ---- Outputs for later scripts ----------------------------------------------

# Thresholds on 24 Nov 1988 (hard-coded: middle of the 1988-89 campaign),
# used by ESM Fig. S3
write_csv(as_tibble(as.list(thresholds(as.Date("1988-11-24")))),
          file.path(DER, "tritium_thresholds_1988-11-24.csv"))

# ESM Table S1: every sample (tritium, decay-corrected tritium, class,
# primary young fraction, conductivity)
tr %>%
  select(Location, Date, period, Tritium, tritium_dc, cls, yf, Cond) %>%
  write_csv(file.path(TAB, "tableS1_tritium_conductivity.csv"))
