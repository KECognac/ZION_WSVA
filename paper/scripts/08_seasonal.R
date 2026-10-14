# 08_seasonal.R -- seasonal means and within-campaign trends at the five
# springs sampled in 1988-89 (paper section 3.1; ESM Table S5).
#
# Seasons: Dec-Feb winter, Mar-May spring, Jun-Aug summer, Sep-Nov fall.
# Trend = least-squares slope per year over the campaign (three or more
# samples with at least two distinct values).

source("paper/scripts/00_setup.R")

# Season of each date (index = month)
season <- function(d) {
  c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer",
    "Summer", "Summer", "Fall", "Fall", "Fall", "Winter")[
      as.numeric(format(d, "%m"))
    ]
}

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  mutate(season = season(date))
isd <- read_csv(file.path(TAB, "tableS2_isotopes.csv"), show_col_types = FALSE) %>%
  mutate(season = season(date))

# Linear trend (per year) of column `col`; NA when there are too few samples
trend <- function(d, col) {
  d <- d[!is.na(d[[col]]), ]
  if (nrow(d) < 3 || n_distinct(d[[col]]) < 2) {
    return(c(n = nrow(d), slope = NA, p = NA))
  }
  x <- as.numeric(d$date - min(d$date)) / 365.25
  f <- summary(lm(d[[col]] ~ x))$coefficients
  c(n = nrow(d), slope = f[2, 1], p = f[2, 4])
}

# Variables: c(label, column name)
vars <- list(
  c("Tritium, decay-corrected (TU)", "tritium_dc"),
  c("Primary young fraction", "yf"),
  c("Specific conductivity (uS/cm)", "Cond"),
  c("Snow fraction", "fsnow")
)

s5 <- map_dfr(REPEAT_SPRINGS, function(s) {
  map_dfr(vars, function(v) {
    # 1988-89 samples only (isotope samples before 2000)
    d <- if (v[2] == "fsnow") {
      filter(isd, group == s, date < as.Date("2000-01-01"))
    } else {
      filter(tr, Location == s, period == "1988-89")
    }
    m <- d %>%
      filter(!is.na(.data[[v[2]]])) %>%
      group_by(season) %>%
      summarise(v = mean(.data[[v[2]]]), .groups = "drop")
    t <- trend(d, v[2])
    tibble(
      site = s,
      variable = v[1],
      Summer = m$v[m$season == "Summer"][1],
      Fall = m$v[m$season == "Fall"][1],
      Winter = m$v[m$season == "Winter"][1],
      Spring = m$v[m$season == "Spring"][1],
      n = t[["n"]],
      slope_per_yr = t[["slope"]],
      p = t[["p"]],
      lowest = if (nrow(m) > 1) m$season[which.min(m$v)] else NA,
      highest = if (nrow(m) > 1) m$season[which.max(m$v)] else NA
    )
  })
})
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
