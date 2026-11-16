# 04_tritium_age_models.R -- tritium age-distribution (lumped-parameter)
# models (paper section 4.3; ESM Table S8, Fig. S3).
#
# Input: monthly tritium in precipitation for the 37-39 N, 110-115 W
# quadrangle containing Zion, Aug 1953 - Dec 2022 (Jurgens 2018, ver. 2.0;
# Michel et al. 2018).
#   - The one missing month (Feb 1954) is linearly interpolated.
#   - Months after Dec 2022 are held at the mean of the last 12 months.
#   - Ages older than the record use 5 TU (sensitivity 0-10 TU).
# Models
#   PFM  piston flow: all water has the same age.
#   EM   exponential: ages are exponentially distributed around a mean age.
#   BMM  binary mixture of tritium-free old water and a young EM component of
#        mean age Ty; the young fraction is observed / modeled young tritium.
# Input weighting: all months (annual) or November-April only (winter).
#
# Outputs
#   derived: tritium_age_models_by_sample.csv,
#            tritium_age_models_pre1953_sensitivity.csv,
#            tritium_input_quadrangle.csv, tritium_em_young_1988.csv
#   tables:  tableS8_tritium_age_models.csv
# Runs after 03_tritium.R.

source("paper/scripts/00_setup.R")

LAM <- log(2) / H3_HALF_LIFE   # tritium decay constant (1/yr)
PRE <- 5                       # hard-coded: tritium (TU) for pre-record ages
TY <- c(2, 5, 10, 20, 30)      # hard-coded: young-component mean ages (yr)


# ---- Precipitation input record ---------------------------------------------

# Column "39-37_115-110" is the quadrangle containing Zion. t is the decimal
# year at mid-month.
q <- read_csv(ext_path("usgs_tritium_quadrangles/Table_2_TritiumInQuadrangles.csv"),
              show_col_types = FALSE,
              col_types = cols(.default = col_character())) %>%
  transmute(
    Year = as.numeric(Year),
    Month = as.numeric(Month),
    c = suppressWarnings(as.numeric(`39-37_115-110`))   # 'na' becomes NA
  ) %>%
  mutate(t = Year + (Month - 0.5) / 12)

# Feb 1954 is recorded as 'na': fill it by linear interpolation in time
q <- q %>%
  mutate(c = approx(t[!is.na(c)], c[!is.na(c)], t)$y)

# Extend past the end of the record (Dec 2022) to Sep 2025, holding the mean
# of the last 12 months
tmax <- max(q$t)
held <- q %>% filter(t > tmax - 1) %>% pull(c) %>% mean()
extension <- tibble(t = seq(tmax + 1 / 12, 2025.75 - 1e-9, by = 1 / 12)) %>%
  mutate(Month = (floor((t * 12) %% 12 + 0.5) %% 12) + 1, c = held)
inp <- bind_rows(q %>% select(t, Month, c), extension)


# ---- Lumped-parameter models ------------------------------------------------

# Input series: all months, or winter months (November-April) only
series <- function(winter) {
  if (winter) inp %>% filter(Month %in% c(11, 12, 1, 2, 3, 4)) else inp
}

# Modeled tritium at sample time ts (decimal year) for mean age T (EM) or
# piston age T (PFM). `pre` is the tritium assumed for ages older than the
# record.
model <- function(ts, T, kind, winter = FALSE, pre = PRE) {
  s <- series(winter) %>%
    mutate(tau = ts - t) %>%          # age of each month's water at ts
    filter(tau >= 0)
  tau <- s$tau
  cc <- s$c
  dt <- if (winter) 2 / 12 else 1 / 12   # weight of each month

  if (kind == "PFM") {
    # older than the record: decayed pre-record concentration
    if (T > max(tau)) return(pre * exp(-LAM * T))
    i <- which.min(abs(tau - T))      # the month closest to age T
    return(cc[i] * exp(-LAM * tau[i]))
  }

  # Exponential model: weighted sum over the record, plus the tail of the age
  # distribution older than the record (at concentration `pre`)
  g <- exp(-tau / T) / T * dt
  tail <- exp(-max(tau) / T)
  tail_c <- pre * exp(-LAM * max(tau)) / (1 + LAM * T)
  sum(cc * g * exp(-LAM * tau)) + tail * tail_c
}

# hard-coded: age grid for solving, 0.5-99.5 yr by 0.5, then 100-3,000 yr by 10
GRID <- c(seq(0.5, 99.5, 0.5), seq(100, 3000, 10))

# All ages on GRID at which the model matches the observed tritium, found by
# linear interpolation between neighboring grid points where the modeled
# value crosses the observation (there can be several)
solve_ages <- function(ts, obs, kind, winter = FALSE) {
  v <- map_dbl(GRID, ~ model(ts, .x, kind, winter))
  n <- length(GRID)
  i <- which((v[-n] - obs) * (v[-1] - obs) <= 0)   # crossing between i, i+1
  GRID[i] + (obs - v[i]) * (GRID[i + 1] - GRID[i]) / (v[i + 1] - v[i])
}

# Ages as text, e.g. "304; 363", or "none" if there is no solution
fmt <- function(x) {
  if (length(x)) str_c(sprintf("%.0f", x), collapse = "; ") else "none"
}

# Binary-mixture young fraction for young-component mean age Ty, at most 1
bmm_fraction <- function(ts, obs, Ty) min(obs / model(ts, Ty, "EM"), 1)


# ---- Apply to every sample --------------------------------------------------

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Tritium)) %>%
  mutate(ts = dec_year(date))

# BMM young fraction for each Ty: one column per Ty, BMM_f_Ty2 ... BMM_f_Ty30
bmm <- TY %>%
  set_names(str_c("BMM_f_Ty", TY)) %>%
  map(function(Ty) map2_dbl(tr$ts, tr$Tritium, bmm_fraction, Ty = Ty)) %>%
  as_tibble()

res <- tr %>%
  mutate(
    EM_annual = map2_chr(ts, Tritium, ~ fmt(solve_ages(.x, .y, "EM"))),
    EM_winter = map2_chr(ts, Tritium, ~ fmt(solve_ages(.x, .y, "EM", TRUE))),
    PFM_annual = map2_chr(ts, Tritium, ~ fmt(solve_ages(.x, .y, "PFM")))
  ) %>%
  bind_cols(bmm) %>%
  select(Location, date, tritium = Tritium, EM_annual, EM_winter, PFM_annual,
         starts_with("BMM"))
write_csv(res, file.path(DER, "tritium_age_models_by_sample.csv"))

# ESM Table S8: site means of the BMM fractions and every EM age
site <- res %>%
  group_by(Location) %>%
  summarise(
    across(starts_with("BMM"), mean),
    EM_ages = str_c(EM_annual, collapse = "; "),
    .groups = "drop"
  )
write_csv(site, file.path(TAB, "tableS8_tritium_age_models.csv"))

# Sensitivity to the pre-record tritium (0, 5, 10 TU) for two sample times
# (late 1988 and mid 2025) and two mean ages (50 and 200 yr)
sens <- expand_grid(ts = c(1988.75, 2025.6), T = c(50, 200),
                    pre = c(0, 5, 10)) %>%
  mutate(tritium = pmap_dbl(list(ts, T, pre),
                            ~ model(..1, ..2, "EM", pre = ..3)))
write_csv(sens, file.path(DER, "tritium_age_models_pre1953_sensitivity.csv"))

# The input record as used, with the held months flagged
inp %>%
  transmute(t, c, held = t > tmax) %>%
  write_csv(file.path(DER, "tritium_input_quadrangle.csv"))
print(as.data.frame(site %>% select(-EM_ages)), digits = 3)


# ---- Regression checks (section 4.3 / Table S8; Python version, Oct 2026) ---

ref <- tribble(
  ~Location,             ~Ty2,  ~Ty30,
  "Weeping Rock",        0.509, 0.183,
  "The Grotto",          0.353, 0.127,
  "Refrigerator Spring", 0.270, 0.095,
  "Larry's Grotto",      0.118, 0.043,
  "Echo Canyon Spring",  1.000, 0.868,
  "Taylor Creek Well",   0.365, 0.549
) %>%
  left_join(site, by = "Location")
check_close(c(ref$BMM_f_Ty2, ref$BMM_f_Ty30), c(ref$Ty2, ref$Ty30), 0.0015,
            "binary-mixture young fractions")
check_close(sens %>% filter(ts == 1988.75, T == 50) %>% pull(tritium),
            c(33.373, 33.463, 33.553), 0.002, "pre-1953 sensitivity")
check_close(held, 7.25, 0.005, "held value after 2022")

# Young-component (exponential) tritium in late 1988 for Ty = 2, 10, 30 yr,
# used by ESM Fig. S3
tibble(Ty = c(2, 10, 30)) %>%
  mutate(tritium = map_dbl(Ty, ~ model(1988.9, .x, "EM"))) %>%
  write_csv(file.path(DER, "tritium_em_young_1988.csv"))
