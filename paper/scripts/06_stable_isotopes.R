# 06_stable_isotopes.R -- deuterium excess, snow fraction and isotopic
# variability (paper sections 2.5, 3.3; ESM Tables S2, S6).
#
# Primary snow fraction: f_snow = (d + 10.3) / 27.1, from two hard-coded
# endmembers in isotope.csv:
#   - snow near Weeping Rock, 27 Jan 1987, d-excess = 16.8
#   - rain at Birch Creek, 12 Sep 2023, d-excess = -10.3
#   (27.1 = 16.8 - (-10.3))
# Sensitivity: winter and summer mean d-excess of Cedar City event
# precipitation, 1991-93 (Friedman et al. 2002); outliers with
# |d - median| > 4 scaled MAD removed; Monte Carlo on the endmember means
# (normal, standard error), 20,000 draws.

source("paper/scripts/00_setup.R")

iso <- read_csv(raw_path("isotope.csv"), show_col_types = FALSE) %>%
  mutate(date = as.Date(Collection_Date, format = "%m/%d/%y"))

# Site name = text before the first comma of the Location field
loc <- str_trim(str_split_fixed(replace_na(iso$Location, ""), ",", 2)[, 1])


# ---- Group isotope samples by site ------------------------------------------

# hard-coded: analysis group -> Location labels in isotope.csv (the source
# spells Refrigerator "Refridgerator")
groups <- list(
  "Echo Canyon Spring" = "Echo Canyon Spring",
  "Weeping Rock" = c("Weeping Rock (spring)", "Weeping Rock"),
  "The Grotto" = "The Grotto",
  "Refrigerator Spring" = "Refridgerator Spring",
  "Larry's Grotto" = "Larry's Grotto",
  "Birch Spring" = c("Birch Spring Box",
                     "Birch Spring upper sediment trap (spring box)"),
  "East Entrance Well" = "East Entrance Well",
  "Lower/Upper Grotto Complex" = c(
    "Lower Grotto Spring Box", "Upper Grotto Spring Box",
    "Upper Grotto Spring; ZION_S0002",
    "Lower Grotto sediment trap (spring box)",
    "Upper Grotto sediment trap (spring box)",
    "Upper Grotto hanging garden"
  ),
  "Taylor Creek Well" = "Taylor Creek Well"
)

isd <- imap_dfr(groups, ~ iso[loc %in% .x, ] %>%
                  mutate(group = .y, label = loc[loc %in% .x])) %>%
  mutate(
    dex = delta_2H - 8 * delta_18O,
    fsnow = (dex + 10.3) / 27.1   # primary endmembers (see header)
  ) %>%
  arrange(match(group, names(groups)), date)

isd %>%
  select(group, label, Collection_Date, date, delta_18O, delta_2H, dex, fsnow) %>%
  write_csv(file.path(TAB, "tableS2_isotopes.csv"))

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

cc <- iso %>%
  filter(Location == "Cedar City(Utah)", Sample_Type == "P") %>%
  mutate(d = delta_2H - 8 * delta_18O, m = as.numeric(format(date, "%m")))
keep <- abs(cc$d - median(cc$d)) <= 4 * mad(cc$d)   # outlier screen
cc <- cc[keep, ]

# Mean, standard error and n of d-excess for the given months
em <- function(months) {
  s <- cc$d[cc$m %in% months]
  c(mean = mean(s), se = sd(s) / sqrt(length(s)), n = length(s))
}

# hard-coded season definitions: list(winter months, summer months)
defs <- list(
  "Nov-Mar vs Jun-Sep" = list(c(11, 12, 1, 2, 3), 6:9),
  "Oct-Apr vs May-Sep" = list(c(10:12, 1:4), 5:9)
)
N_DRAWS <- 20000

set.seed(SEED)
sens <- imap_dfr(defs, function(d, name) {
  w <- em(d[[1]])
  s <- em(d[[2]])
  dw <- rnorm(N_DRAWS, w[["mean"]], w[["se"]])   # winter endmember draws
  ds <- rnorm(N_DRAWS, s[["mean"]], s[["se"]])   # summer endmember draws
  site %>%
    rowwise() %>%
    mutate(
      definition = name,
      winter = w[["mean"]],
      winter_se = w[["se"]],
      n_winter = w[["n"]],
      summer = s[["mean"]],
      summer_se = s[["se"]],
      n_summer = s[["n"]],
      f_point = min(max((dex - s[["mean"]]) / (w[["mean"]] - s[["mean"]]), 0), 1),
      f_lo = quantile(pmin(pmax((dex - ds) / (dw - ds), 0), 1), 0.025),
      f_med = median(pmin(pmax((dex - ds) / (dw - ds), 0), 1)),
      f_hi = quantile(pmin(pmax((dex - ds) / (dw - ds), 0), 1), 0.975)
    ) %>%
    ungroup()
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
        " (outliers removed: ", sum(!keep), ")")

ck <- function(g, col) site[[col]][site$group == g]
check_close(c(ck("Weeping Rock", "f_primary"),
              ck("Echo Canyon Spring", "f_primary"),
              ck("Taylor Creek Well", "f_primary")),
            c(0.815, 0.576, 0.697), 0.0005, "primary snow fractions")
check_close(sd88$sd_d18O[match(c("Larry's Grotto", "Refrigerator Spring",
                                 "The Grotto", "Weeping Rock",
                                 "Echo Canyon Spring"), sd88$group)],
            c(0.04, 0.05, 0.08, 0.10, 2.93), 0.005,
            "1988-89 d18O standard deviations")

a <- sens %>% filter(definition == "Nov-Mar vs Jun-Sep")
check_close(c(a$winter[1], a$summer[1]), c(11.22, 4.34), 0.005,
            "Cedar City endmembers")
check_close(a$f_med[match(c("Echo Canyon Spring", "Refrigerator Spring",
                            "Weeping Rock"), a$group)],
            c(0.14, 0.92, 1.00), 0.01, "Cedar City snow fractions (medians)")
check_close(min(a$f_lo[!a$group %in% c("Echo Canyon Spring",
                                       "Taylor Creek Well")]),
            0.51, 0.01, "lowest lower bound, Zion Canyon sites")
