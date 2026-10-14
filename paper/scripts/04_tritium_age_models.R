# 04_tritium_age_models.R -- tritium age-distribution (lumped-parameter)
# models (paper section 4.3; ESM Table S8, Fig. S3).
#
# Input: monthly tritium in precipitation for the 37-39 N, 110-115 W
# quadrangle containing Zion, Aug 1953 - Dec 2022 (Jurgens 2018, ver. 2.0;
# Michel et al. 2018).
#   - The one missing month (Feb 1954) is linearly interpolated.
#   - Months after Dec 2022 are held at the mean of the last 12 months.
#   - Ages older than the record use 5 TU (sensitivity 0-10 TU).
# Models: piston flow (PFM); exponential (EM); binary mixture (BMM) of
# tritium-free old water and a young EM component of mean age Ty.
# Input weighting: all months (annual) or November-April only (winter).

source("paper/scripts/00_setup.R")

LAM <- log(2) / H3_HALF_LIFE   # tritium decay constant (1/yr)
PRE <- 5                       # hard-coded: tritium (TU) for pre-record ages


# ---- Precipitation input record ---------------------------------------------

# Column "39-37_115-110" is the quadrangle containing Zion
q <- read_csv(ext_path("usgs_tritium_quadrangles/Table_2_TritiumInQuadrangles.csv"),
              show_col_types = FALSE,
              col_types = cols(.default = col_character())) %>%
  transmute(
    Year = as.numeric(Year),
    Month = as.numeric(Month),
    c = suppressWarnings(as.numeric(`39-37_115-110`))
  ) %>%
  mutate(t = Year + (Month - 0.5) / 12)

# Feb 1954 is recorded as 'na': linear interpolation
q$c <- approx(q$t[!is.na(q$c)], q$c[!is.na(q$c)], q$t)$y

# Extend past the end of the record (Dec 2022) to Sep 2025, holding the mean
# of the last 12 months
tmax <- max(q$t)
held <- mean(q$c[q$t > tmax - 1])
ext <- tibble(t = seq(tmax + 1 / 12, 2025.75 - 1e-9, by = 1 / 12)) %>%
  mutate(Month = (floor((t * 12) %% 12 + 0.5) %% 12) + 1, c = held)
inp <- bind_rows(q %>% select(t, Month, c), ext)


# ---- Lumped-parameter models ------------------------------------------------

# Input series: all months, or winter months (Nov-Apr) only
series <- function(winter) {
  s <- if (winter) filter(inp, Month %in% c(11, 12, 1, 2, 3, 4)) else inp
  list(t = s$t, c = s$c)
}

# Modeled tritium at sample time ts for mean age (or piston age) T
model <- function(ts, T, kind, winter = FALSE, pre = PRE) {
  s <- series(winter)
  tau <- ts - s$t
  ok <- tau >= 0
  tau <- tau[ok]
  cc <- s$c[ok]
  dt <- if (winter) 2 / 12 else 1 / 12

  if (kind == "PFM") {
    if (T > max(tau)) return(pre * exp(-LAM * T))
    i <- which.min(abs(tau - T))
    return(cc[i] * exp(-LAM * tau[i]))
  }

  # exponential model: weighted sum over the record, plus the tail of the
  # age distribution older than the record (at concentration `pre`)
  g <- exp(-tau / T) / T * dt
  tail <- exp(-max(tau) / T)
  tail_c <- pre * exp(-LAM * max(tau)) / (1 + LAM * T)
  sum(cc * g * exp(-LAM * tau)) + tail * tail_c
}

# hard-coded: age grid for solving, 0.5-99.5 yr by 0.5, then 100-3,000 yr by 10
GRID <- c(seq(0.5, 99.5, 0.5), seq(100, 3000, 10))

# All ages on GRID at which the model matches the observed tritium (linear
# interpolation between grid points; there can be several)
solve_ages <- function(ts, obs, kind, winter = FALSE) {
  v <- sapply(GRID, function(T) model(ts, T, kind, winter))
  out <- c()
  for (i in seq_len(length(GRID) - 1)) {
    if ((v[i] - obs) * (v[i + 1] - obs) <= 0) {
      out <- c(out, GRID[i] + (obs - v[i]) * (GRID[i + 1] - GRID[i]) /
                 (v[i + 1] - v[i]))
    }
  }
  out
}

fmt <- function(x) {
  if (length(x)) paste(sprintf("%.0f", x), collapse = "; ") else "none"
}


# ---- Apply to every sample --------------------------------------------------

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Tritium)) %>%
  mutate(ts = dec_year(date))

# BMM young fraction for young-component mean ages Ty = 2, 5, 10, 20, 30 yr
res <- tr %>%
  rowwise() %>%
  mutate(
    EM_annual = fmt(solve_ages(ts, Tritium, "EM")),
    EM_winter = fmt(solve_ages(ts, Tritium, "EM", TRUE)),
    PFM_annual = fmt(solve_ages(ts, Tritium, "PFM")),
    BMM_f_Ty2 = min(Tritium / model(ts, 2, "EM"), 1),
    BMM_f_Ty5 = min(Tritium / model(ts, 5, "EM"), 1),
    BMM_f_Ty10 = min(Tritium / model(ts, 10, "EM"), 1),
    BMM_f_Ty20 = min(Tritium / model(ts, 20, "EM"), 1),
    BMM_f_Ty30 = min(Tritium / model(ts, 30, "EM"), 1)
  ) %>%
  ungroup() %>%
  select(Location, date, tritium = Tritium, EM_annual, EM_winter, PFM_annual,
         starts_with("BMM"))
write_csv(res, file.path(DER, "tritium_age_models_by_sample.csv"))

site <- res %>%
  group_by(Location) %>%
  summarise(
    across(starts_with("BMM"), mean),
    EM_ages = paste(EM_annual, collapse = "; "),
    .groups = "drop"
  )
write_csv(site, file.path(TAB, "tableS8_tritium_age_models.csv"))

# Sensitivity to the pre-record tritium (0, 5, 10 TU) for two sample times
# (late 1988 and mid 2025) and two mean ages (50 and 200 yr)
sens <- expand_grid(ts = c(1988.75, 2025.6), T = c(50, 200),
                    pre = c(0, 5, 10)) %>%
  rowwise() %>%
  mutate(tritium = model(ts, T, "EM", pre = pre)) %>%
  ungroup()
write_csv(sens, file.path(DER, "tritium_age_models_pre1953_sensitivity.csv"))
write_csv(tibble(t = inp$t, c = inp$c, held = inp$t > tmax),
          file.path(DER, "tritium_input_quadrangle.csv"))
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
)
m <- ref %>% left_join(site, by = "Location")
check_close(c(m$BMM_f_Ty2, m$BMM_f_Ty30), c(m$Ty2, m$Ty30), 0.0015,
            "binary-mixture young fractions")
check_close(sens$tritium[sens$ts == 1988.75 & sens$T == 50],
            c(33.373, 33.463, 33.553), 0.002, "pre-1953 sensitivity")
check_close(held, 7.25, 0.005, "held value after 2022")

# Young-component (exponential) tritium in late 1988 for Ty = 2, 10, 30 yr,
# used by ESM Fig. S3
write_csv(tibble(Ty = c(2, 10, 30),
                 tritium = sapply(c(2, 10, 30),
                                  function(T) model(1988.9, T, "EM"))),
          file.path(DER, "tritium_em_young_1988.csv"))
