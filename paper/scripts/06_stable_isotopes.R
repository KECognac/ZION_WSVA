# 06_stable_isotopes.R -- deuterium excess, snow fraction and isotopic
# variability (paper sections 2.5, 3.3; ESM Tables S2, S6).
#
# Primary snow fraction: f_snow = (d + 10.3) / 27.1, where d is deuterium
# excess (d = d2H - 8 d18O), from two hard-coded endmembers in isotope.csv:
#   - snow near Weeping Rock, 27 Jan 1987, d-excess = 16.8
#   - rain at Birch Creek, 12 Sep 2023, d-excess = -10.3
#   (27.1 = 16.8 - (-10.3))
# Sensitivity: winter and summer mean d-excess of Cedar City event
# precipitation, 1991-93 (Friedman et al. 2002); outliers with
# |d - median| > 4 scaled MAD removed; Monte Carlo on the endmember means
# (normal, standard error), 20,000 draws.
#
# Outputs
#   derived: snow_fraction_by_group.csv, d18O_variability_1988_89.csv
#   tables:  tableS2_isotopes.csv, tableS6_snow_fraction_sensitivity.csv

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages(library(lubridate))

iso <- read_csv(raw_path("isotope.csv"), show_col_types = FALSE) %>%
  mutate(
    date = mdy(Collection_Date),
    # site name = text before the first comma of the Location field
    loc = str_trim(str_remove(replace_na(Location, ""), ",.*"))
  )


# ---- Group isotope samples by site ------------------------------------------

# hard-coded: analysis group -> Location labels in isotope.csv (the source
# spells Refrigerator "Refridgerator")
groups <- tribble(
  ~group,                       ~label,
  "Echo Canyon Spring",         "Echo Canyon Spring",
  "Weeping Rock",               "Weeping Rock (spring)",
  "Weeping Rock",               "Weeping Rock",
  "The Grotto",                 "The Grotto",
  "Refrigerator Spring",        "Refridgerator Spring",
  "Larry's Grotto",             "Larry's Grotto",
  "Birch Spring",               "Birch Spring Box",
  "Birch Spring",               "Birch Spring upper sediment trap (spring box)",
  "East Entrance Well",         "East Entrance Well",
  "Lower/Upper Grotto Complex", "Lower Grotto Spring Box",
  "Lower/Upper Grotto Complex", "Upper Grotto Spring Box",
  "Lower/Upper Grotto Complex", "Upper Grotto Spring; ZION_S0002",
  "Lower/Upper Grotto Complex", "Lower Grotto sediment trap (spring box)",
  "Lower/Upper Grotto Complex", "Upper Grotto sediment trap (spring box)",
  "Lower/Upper Grotto Complex", "Upper Grotto hanging garden",
  "Taylor Creek Well",          "Taylor Creek Well"
)
group_order <- unique(groups$group)

isd <- iso %>%
  inner_join(groups, by = c("loc" = "label"), keep = FALSE) %>%
  mutate(
    label = loc,
    dex = delta_2H - 8 * delta_18O,   # deuterium excess
    fsnow = (dex + 10.3) / 27.1       # primary endmembers (see header)
  ) %>%
  arrange(factor(group, levels = group_order), date)

# ESM Table S2: every sample with its group
isd %>%
  select(group, label, Collection_Date, date, delta_18O, delta_2H, dex, fsnow) %>%
  write_csv(file.path(TAB, "tableS2_isotopes.csv"))

# Mean d-excess and snow fraction by group
site <- isd %>%
  group_by(group) %>%
  summarise(
    n = n(),
    dex = mean(dex),
    f_primary = mean(fsnow),
    f_min = min(fsnow),
    f_max = max(fsnow),
    .groups = "drop"
  )

# d18O variability within the 1988-89 campaign (samples before 2000)
sd88 <- isd %>%
  filter(group %in% REPEAT_SPRINGS, date < as.Date("2000-01-01")) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    sd_d18O = sd(delta_18O),
    min_d18O = min(delta_18O),
    max_d18O = max(delta_18O),
    .groups = "drop"
  )
write_csv(site, file.path(DER, "snow_fraction_by_group.csv"))
write_csv(sd88, file.path(DER, "d18O_variability_1988_89.csv"))


# ---- Cedar City endmember sensitivity ---------------------------------------

# Cedar City event precipitation, with d-excess and month
cc_all <- iso %>%
  filter(Location == "Cedar City(Utah)", Sample_Type == "P") %>%
  mutate(d = delta_2H - 8 * delta_18O, m = month(date))
# outlier screen: keep samples within 4 scaled MAD of the median
cc <- cc_all %>% filter(abs(d - median(d)) <= 4 * mad(d))

# Mean, standard error and n of d-excess for the given months
em <- function(months) {
  cc %>%
    filter(m %in% months) %>%
    summarise(mean = mean(d), se = sd(d) / sqrt(n()), n = n())
}

# hard-coded season definitions: list(winter months, summer months)
defs <- list(
  "Nov-Mar vs Jun-Sep" = list(c(11, 12, 1, 2, 3), 6:9),
  "Oct-Apr vs May-Sep" = list(c(10:12, 1:4), 5:9)
)
N_DRAWS <- 20000

# Snow fraction for d-excess d between endmembers summer s and winter w,
# limited to 0-1 (works on vectors of draws)
f_snow <- function(d, w, s) pmin(pmax((d - s) / (w - s), 0), 1)

set.seed(SEED)
sens <- imap_dfr(defs, function(def, name) {
  w <- em(def[[1]])   # winter endmember
  s <- em(def[[2]])   # summer endmember
  dw <- rnorm(N_DRAWS, w$mean, w$se)   # winter endmember draws
  ds <- rnorm(N_DRAWS, s$mean, s$se)   # summer endmember draws
  site %>%
    mutate(
      definition = name,
      winter = w$mean, winter_se = w$se, n_winter = w$n,
      summer = s$mean, summer_se = s$se, n_summer = s$n,
      f_point = f_snow(dex, w$mean, s$mean),
      # 95% interval and median of the snow fraction over the draws
      f_lo = map_dbl(dex, ~ quantile(f_snow(.x, dw, ds), 0.025)),
      f_med = map_dbl(dex, ~ median(f_snow(.x, dw, ds))),
      f_hi = map_dbl(dex, ~ quantile(f_snow(.x, dw, ds), 0.975))
    )
})
write_csv(sens, file.path(TAB, "tableS6_snow_fraction_sensitivity.csv"))

print(as.data.frame(site), digits = 3)
print(as.data.frame(sd88), digits = 3)
sens %>%
  filter(definition == "Nov-Mar vs Jun-Sep") %>%
  select(group, f_primary, f_med, f_lo, f_hi) %>%
  as.data.frame() %>%
  print(digits = 2)


# ---- Regression checks against the paper ------------------------------------

message("Cedar City samples used: ", nrow(cc),
        " (outliers removed: ", nrow(cc_all) - nrow(cc), ")")

# value of column `col` for the given groups, in the order given
pick_groups <- function(d, groups, col) {
  tibble(group = groups) %>% left_join(d, by = "group") %>% pull({{ col }})
}

check_close(pick_groups(site, c("Weeping Rock", "Echo Canyon Spring",
                                "Taylor Creek Well"), f_primary),
            c(0.815, 0.576, 0.697), 0.0005, "primary snow fractions")
check_close(pick_groups(sd88, c("Larry's Grotto", "Refrigerator Spring",
                                "The Grotto", "Weeping Rock",
                                "Echo Canyon Spring"), sd_d18O),
            c(0.04, 0.05, 0.08, 0.10, 2.93), 0.005,
            "1988-89 d18O standard deviations")

a <- sens %>% filter(definition == "Nov-Mar vs Jun-Sep")
check_close(c(first(a$winter), first(a$summer)), c(11.22, 4.34), 0.005,
            "Cedar City endmembers")
check_close(pick_groups(a, c("Echo Canyon Spring", "Refrigerator Spring",
                             "Weeping Rock"), f_med),
            c(0.14, 0.92, 1.00), 0.01, "Cedar City snow fractions (medians)")
check_close(a %>%
              filter(!group %in% c("Echo Canyon Spring", "Taylor Creek Well")) %>%
              pull(f_lo) %>%
              min(),
            0.51, 0.01, "lowest lower bound, Zion Canyon sites")
