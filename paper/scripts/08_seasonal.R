# 08_seasonal.R -- seasonal means and within-campaign trends at the five
# springs sampled in 1988-89 (paper section 3.1; ESM Table S5).
#
# Seasons: Dec-Feb winter, Mar-May spring, Jun-Aug summer, Sep-Nov fall.
# Trend = least-squares slope per year over the campaign (three or more
# samples with at least two distinct values).
#
# Outputs
#   derived: later_samples_repeat_sites.csv
#   tables:  tableS5_seasonal.csv
# Runs after 03 (tritium) and 06 (isotopes).

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages(library(lubridate))

# Season of each date
season <- function(d) {
  case_when(
    month(d) %in% c(12, 1, 2) ~ "Winter",
    month(d) %in% 3:5 ~ "Spring",
    month(d) %in% 6:8 ~ "Summer",
    month(d) %in% 9:11 ~ "Fall"
  )
}

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  mutate(season = season(date))
isd <- read_csv(file.path(TAB, "tableS2_isotopes.csv"), show_col_types = FALSE) %>%
  mutate(season = season(date))

# Linear trend (per year) of column `col` in data d: a one-row tibble of the
# number of samples, the slope and its p-value (NA with too few samples)
trend <- function(d, col) {
  d <- d %>% filter(!is.na(.data[[col]]))
  if (nrow(d) < 3 || n_distinct(d[[col]]) < 2) {
    return(tibble(n = nrow(d), slope = NA_real_, p = NA_real_))
  }
  d <- d %>% mutate(x = as.numeric(date - min(date)) / 365.25,   # years
                    y = .data[[col]])
  f <- summary(lm(y ~ x, data = d))$coefficients
  tibble(n = nrow(d), slope = f["x", "Estimate"], p = f["x", "Pr(>|t|)"])
}

# Variables: label in the table and column name
vars <- tribble(
  ~label,                            ~col,
  "Tritium, decay-corrected (TU)",   "tritium_dc",
  "Primary young fraction",          "yf",
  "Specific conductivity (uS/cm)",   "Cond",
  "Snow fraction",                   "fsnow"
)

# One row of Table S5: seasonal means and trend of one variable at one spring
seasonal_row <- function(s, label, col) {
  # 1988-89 samples only (isotope samples before 2000)
  d <- if (col == "fsnow") {
    isd %>% filter(group == s, date < as.Date("2000-01-01"))
  } else {
    tr %>% filter(Location == s, period == "1988-89")
  }
  m <- d %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(season) %>%
    summarise(v = mean(.data[[col]]), .groups = "drop")
  mean_in <- function(x) m %>% filter(season == x) %>% pull(v) %>% first()
  t <- trend(d, col)
  tibble(
    site = s,
    variable = label,
    Summer = mean_in("Summer"),
    Fall = mean_in("Fall"),
    Winter = mean_in("Winter"),
    Spring = mean_in("Spring"),
    n = t$n,
    slope_per_yr = t$slope,
    p = t$p,
    # season with the lowest and highest mean (only with two or more seasons)
    lowest = if (nrow(m) > 1) m$season[which.min(m$v)] else NA,
    highest = if (nrow(m) > 1) m$season[which.max(m$v)] else NA
  )
}

# every spring x every variable, springs in the order of REPEAT_SPRINGS
s5 <- expand_grid(site = REPEAT_SPRINGS, vars) %>%
  pmap_dfr(function(site, label, col) seasonal_row(site, label, col))
write_csv(s5, file.path(TAB, "tableS5_seasonal.csv"))

# Later (2005-2025) samples at the repeat-sampled springs
later <- bind_rows(
  tr %>%
    filter(Location %in% REPEAT_SPRINGS, period == "2025") %>%
    transmute(Location, date, value = Cond, variable = "Conductivity"),
  isd %>%
    filter(group %in% REPEAT_SPRINGS, date >= as.Date("2000-01-01")) %>%
    transmute(Location = group, date, value = fsnow, variable = "Snow fraction")
)
write_csv(later, file.path(DER, "later_samples_repeat_sites.csv"))

s5 %>%
  filter(variable %in% c("Tritium, decay-corrected (TU)",
                         "Specific conductivity (uS/cm)")) %>%
  as.data.frame() %>%
  print(digits = 3)


# ---- Regression checks ------------------------------------------------------

tdc <- s5 %>% filter(variable == "Tritium, decay-corrected (TU)")
stopifnot(all(tdc$lowest == "Fall"))   # lowest tritium in fall at every spring
rs <- s5 %>%
  filter(site == "Refrigerator Spring",
         variable == "Specific conductivity (uS/cm)")
check_close(c(rs$slope_per_yr, rs$p), c(22.9, 0.019), 0.05,
            "Refrigerator Spring conductivity trend")
check_close(range(tdc$p), c(0.16, 0.88), 0.005, "tritium trend p-values")
