# 14_baseflow_climate.R -- baseflow of the North Fork Virgin River and its
# relation to climate (paper sections 2.8, 3.6, 4.2; Fig. 5b; ESM Fig. S6,
# Tables S10-S11).
#
# Streamflow: USGS 09405500, North Fork Virgin River near Springdale (daily,
# approved; complete water years 1926 and 1928-2025, with a gap from Nov 1923
# to Sep 1925 and a partial WY1927; drainage area 344 mi2).
# Diversion: the Springdale Canal (09405499) diverts at river mile 4.9, just
# above the gage (mile 4.7). For WY1969-1989 the published 09405500 record is
# river + canal (it equals 09405501 + 09405499 on ~90% of days), whereas
# 09405501 (1968-1991) is the river alone; after 1989 the canal is not
# gaged. Two series are therefore analysed:
#   "published"  09405500 as published (canal included only in WY1969-1989)
#   "river-only" 09405501 for WY1969-1991, 09405500 otherwise (canal excluded
#                wherever it was measured)
# and WY1990-2025, when both series are the same, as a third check.
# About 24-55 ft3/s of North Fork flow (up to 100% of the minimum and 85% of
# the median) is groundwater discharge from the Navajo Sandstone
# (Christensen et al. 2005), so baseflow here is mostly Navajo discharge.
#
# Climate: GridMET basin means, WY1980-WY2025: precipitation P, grass
# reference ET (PET), and P - PET (climatic water balance). "Cold-day
# precipitation" = P on days with mean air temperature <= 0 C, an index of
# snowfall.
# Baseflow: Eckhardt (2005) filter with alpha from recession analysis and
# BFImax from the backward filter of Collischonn and Fan (2013); BFImax
# sensitivity cases 0.5 and 0.8. The annual 7-day minimum flow is a
# filter-free check.
# Memory: Spearman correlation of water-year baseflow with climate averaged
# over the preceding k water years (k = 1-10).
#
# Baseflow integrates groundwater discharge from the whole basin (including
# the Narrows), not individual springs. Overlapping windows and serial
# correlation make p-values approximate, so correlations and trends are
# reported descriptively.
#
# Outputs
#   derived: baseflow_summary.csv, baseflow_annual.csv,
#            sampling_year_hydroclimate.csv
#   tables:  tableS10_baseflow_climate_memory.csv, tableS11_baseflow_trends.csv
#   figures: figS6_baseflow_climate.png

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages(library(trend))

# hard-coded: BFImax sensitivity bounds and sampling years highlighted below
BFIMAX_LO <- 0.5
BFIMAX_HI <- 0.8
SAMPLING_WY <- c(1988, 1989, 2025)


# ---- Streamflow -------------------------------------------------------------

site <- read_csv(ext_path("nwis_gage_09405500.csv"), show_col_types = FALSE)
area_km2 <- as.numeric(site$drain_area_va) * 2.589988   # mi2 -> km2

# Daily mean ft3/s -> mm/day over the basin (0.0283168 m3 per ft3)
CFS_TO_MM_DAY <- 0.0283168 * 86400 / (area_km2 * 1e6) * 1000

# Daily discharge (ft3/s) from a cached NWIS file
read_dv <- function(f) {
  read_csv(ext_path(f), show_col_types = FALSE) %>%
    transmute(date = Date, Q = Flow)
}
q_pub <- read_dv("nwis_dv_09405500.csv")
q_501 <- read_dv("nwis_dv_09405501.csv")
q_canal <- read_dv("nwis_dv_09405499.csv")

# Check the gage relation: published 09405500 = 09405501 + canal (WY1969-89)
chk <- q_pub %>%
  inner_join(q_501, by = "date", suffix = c("_pub", "_501")) %>%
  inner_join(q_canal %>% rename(Q_canal = Q), by = "date")
share_sum <- mean(abs(chk$Q_pub - chk$Q_501 - chk$Q_canal) < 0.5)

# River-only series: 09405501 where it exists, published 09405500 otherwise
q_riv <- q_pub %>%
  left_join(q_501 %>% rename(Q501 = Q), by = "date") %>%
  mutate(Q = coalesce(Q501, Q)) %>%
  select(date, Q)

series <- list(published = q_pub, `river-only` = q_riv) %>%
  map(~ .x %>% mutate(wy = water_year(date)) %>% arrange(date))

# Complete water years only (WY1926, WY1928-2025)
keep_wy <- series$published %>%
  count(wy) %>%
  filter(n >= 365) %>%
  pull(wy)


# ---- Baseflow separation ----------------------------------------------------
# One continuous pass over each record; gaps are bridged by the filter's
# memory. alpha and BFImax are estimated from the published record and used
# for both series.

rec <- recession_alpha(series$published$Q)
bfi_max <- bfimax_backward(series$published$Q, rec$alpha)
k_days <- -1 / log(rec$alpha)   # linear-reservoir time constant (days)
message(sprintf(paste0("recession alpha = %.4f (%d recession days); ",
                       "linear-reservoir time constant = %.0f days; ",
                       "BFImax (backward filter) = %.2f"),
                rec$alpha, rec$n_days, k_days, bfi_max))

# Annual (water-year) flow, baseflow and 7-day minimum flow
annual <- function(q) {
  q %>%
    mutate(
      bf = baseflow_eckhardt(Q, rec$alpha, bfi_max),
      bf_lo = baseflow_eckhardt(Q, rec$alpha, BFIMAX_LO),
      bf_hi = baseflow_eckhardt(Q, rec$alpha, BFIMAX_HI)
    ) %>%
    filter(wy %in% keep_wy) %>%
    group_by(wy) %>%
    summarise(
      Q_mm = sum(Q) * CFS_TO_MM_DAY,
      bf_mm = sum(bf) * CFS_TO_MM_DAY,
      bf_lo_mm = sum(bf_lo) * CFS_TO_MM_DAY,
      bf_hi_mm = sum(bf_hi) * CFS_TO_MM_DAY,
      # lowest centered 7-day moving average of daily flow
      q7min_cfs = min(stats::filter(Q, rep(1 / 7, 7), sides = 2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(bfi = bf_mm / Q_mm)
}
ann_all <- map(series, annual)
ann <- ann_all$`river-only`   # primary series (canal excluded where gaged)


# ---- Climate (GridMET basin means) ------------------------------------------

clim <- read_csv(ext_path("gridmet_basin_09405500_daily.csv"),
                 show_col_types = FALSE) %>%
  mutate(wy = water_year(date), tmean = (tmax_c + tmin_c) / 2) %>%
  filter(wy >= 1980, wy <= 2025) %>%
  group_by(wy) %>%
  summarise(
    n = n(),
    P = sum(pr_mm),
    PET = sum(pet_mm),
    P_cold = sum(pr_mm[tmean <= 0]),   # hard-coded 0 C snowfall threshold
    .groups = "drop"
  ) %>%
  filter(n >= 365) %>%
  mutate(PmPET = P - PET) %>%
  select(-n)


# ---- Memory: correlation with climate over the preceding k water years ------

# Trailing mean of the current and k - 1 previous values (NA for the first
# k - 1 values)
roll <- function(x, k) as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))

clim_vars <- c(P = "P", `P - PET` = "PmPET", `Cold-day P` = "P_cold",
               PET = "PET")
resp_vars <- c(`Baseflow (Eckhardt)` = "bf_mm",
               `7-day minimum flow` = "q7min_cfs",
               `Total streamflow` = "Q_mm")

# Spearman correlations of each response with each climate variable averaged
# over k = 1-10 water years, for one annual series over water years `wys`
memory <- function(a, wys) {
  d <- a %>%
    inner_join(clim, by = "wy") %>%
    filter(wy %in% wys)
  expand_grid(climate = names(clim_vars), response = names(resp_vars),
              k = 1:10) %>%
    mutate(stats = pmap(list(climate, response, k), function(cv, rv, k) {
      # mean of WY t-k+1 ... t (climate from WY1980 on), aligned to WY t
      clim_k <- tibble(wy = clim$wy, x = roll(clim[[clim_vars[[cv]]]], k))
      xy <- d %>%
        select(wy, y = all_of(resp_vars[[rv]])) %>%
        left_join(clim_k, by = "wy") %>%
        filter(!is.na(x))
      tibble(rho = cor(xy$x, xy$y, method = "spearman"), n = nrow(xy))
    })) %>%
    unnest(stats)
}

mem <- bind_rows(
  memory(ann_all$`river-only`, 1980:2025) %>%
    mutate(series = "river-only, WY1980-2025"),
  memory(ann_all$published, 1980:2025) %>%
    mutate(series = "published, WY1980-2025"),
  memory(ann_all$`river-only`, 1990:2025) %>%
    mutate(series = "WY1990-2025 (both series)")
) %>%
  select(series, everything())

best <- mem %>%
  group_by(series, climate, response) %>%
  slice_max(rho, n = 1) %>%
  ungroup() %>%
  filter(response == "Baseflow (Eckhardt)") %>%
  arrange(series, desc(rho))
print(as.data.frame(best), digits = 3)
write_csv(mem, file.path(TAB, "tableS10_baseflow_climate_memory.csv"))


# ---- Trends (Sen slope, Mann-Kendall; descriptive) --------------------------

# Sen slope (per decade) and Mann-Kendall p-value of a series
trend_stats <- function(y) {
  s <- sens.slope(y)
  m <- mk.test(y)
  tibble(n = length(y), sen_per_decade = 10 * unname(s$estimates),
         p_mk = m$p.value)
}
d <- ann %>% inner_join(clim, by = "wy")
d90 <- d %>% filter(wy >= 1990)

trends <- bind_rows(
  trend_stats(d$bf_mm) %>%
    mutate(series = "Baseflow (mm), river-only, WY1980-2025"),
  trend_stats(ann_all$published %>% filter(wy >= 1980) %>% pull(bf_mm)) %>%
    mutate(series = "Baseflow (mm), published, WY1980-2025"),
  trend_stats(d90$bf_mm) %>%
    mutate(series = "Baseflow (mm), WY1990-2025"),
  trend_stats(d$q7min_cfs) %>%
    mutate(series = "7-day minimum flow (ft3/s), river-only, WY1980-2025"),
  trend_stats(d90$q7min_cfs) %>%
    mutate(series = "7-day minimum flow (ft3/s), WY1990-2025"),
  trend_stats(d$PmPET) %>%
    mutate(series = "P - PET (mm), WY1980-2025"),
  trend_stats(d90$PmPET) %>%
    mutate(series = "P - PET (mm), WY1990-2025"),
  trend_stats(d$P_cold) %>%
    mutate(series = "Cold-day P (mm), WY1980-2025")
) %>%
  select(series, everything())
print(as.data.frame(trends), digits = 3)


# ---- Hydroclimate context of the sampling years -----------------------------

# Percentile of water year y within WY1980-2025 for column v (share of years
# with a value at or below that year's)
pct <- function(v, y) {
  this_year <- d %>% filter(wy == y) %>% pull(all_of(v))
  round(100 * mean(d[[v]] <= this_year))
}

ctx <- tibble(wy = SAMPLING_WY) %>%
  filter(wy %in% d$wy) %>%
  mutate(
    P_pct = map_dbl(wy, ~ pct("P", .x)),
    PmPET_pct = map_dbl(wy, ~ pct("PmPET", .x)),
    Pcold_pct = map_dbl(wy, ~ pct("P_cold", .x)),
    bf_pct = map_dbl(wy, ~ pct("bf_mm", .x))
  )
print(as.data.frame(ctx))

summ <- tibble(
  area_km2 = area_km2,
  alpha = rec$alpha,
  recession_days = rec$n_days,
  k_days = k_days,
  bfi_max = bfi_max,
  bfi_mean = mean(ann$bfi),
  bf_mm_median = median(ann$bf_mm),
  wy_first = min(ann$wy),
  wy_last = max(ann$wy),
  n_wy = nrow(ann),
  n_wy_climate = nrow(d),
  canal_check_share = share_sum,
  canal_median_cfs = median(q_canal$Q, na.rm = TRUE)
)
write_csv(summ, file.path(DER, "baseflow_summary.csv"))
write_csv(bind_rows(ann_all, .id = "series"),
          file.path(DER, "baseflow_annual.csv"))
write_csv(trends, file.path(TAB, "tableS11_baseflow_trends.csv"))
write_csv(ctx, file.path(DER, "sampling_year_hydroclimate.csv"))


# ---- ESM figure: annual baseflow, P - PET, and memory curves ----------------

# Standardized anomalies for panel b
dd <- d %>%
  mutate(z_bf = (bf_mm - mean(bf_mm)) / sd(bf_mm),
         z_c = (PmPET - mean(PmPET)) / sd(PmPET))

pa <- ggplot(ann, aes(wy, bf_mm)) +
  geom_ribbon(aes(ymin = bf_lo_mm, ymax = bf_hi_mm), fill = "grey85") +
  geom_line(data = ann_all$published %>% filter(wy >= 1969, wy <= 1991),
            color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_line(color = "#3070ad", linewidth = 0.5) +
  geom_point(data = filter(ann, wy %in% SAMPLING_WY), color = "#c66426",
             size = 1.6) +
  labs(x = NULL, y = "Baseflow (mm/yr)",
       title = "a  Annual baseflow, North Fork Virgin River near Springdale",
       caption = paste0("Blue: river only; dashed: published record including ",
                        "the Springdale Canal (WY1969-1989); band: BFImax ",
                        "0.5-0.8")) +
  theme_paper()

pb <- ggplot(dd, aes(wy)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(aes(y = z_c, color = "P - PET (GridMET)"), linewidth = 0.5) +
  geom_line(aes(y = z_bf, color = "Baseflow"), linewidth = 0.5) +
  scale_color_manual(values = c("P - PET (GridMET)" = "#dca237",
                                "Baseflow" = "#3070ad"),
                     name = NULL) +
  labs(x = "Water year", y = "Standardized anomaly",
       title = "b  Baseflow and climatic water balance, WY1980-2025") +
  theme_paper() +
  theme(legend.position = "bottom")

pc <- ggplot(mem %>% filter(response == "Baseflow (Eckhardt)",
                            series == "river-only, WY1980-2025"),
             aes(k, rho, color = climate)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c(P = "#6eb2e4", `P - PET` = "#dca237",
                                `Cold-day P` = "#459a74", PET = "#c07da5"),
                     name = NULL) +
  labs(x = "Climate averaged over the preceding k water years",
       y = "Spearman correlation with annual baseflow",
       title = "c  Memory of baseflow") +
  theme_paper() +
  theme(legend.position = "bottom")

ggsave(file.path(FIG, "figS6_baseflow_climate.png"), pa / (pb | pc),
       width = 7.2, height = 6.8, dpi = 300, bg = "white")


# ---- Checks -----------------------------------------------------------------

check_close(area_km2, 891, 1, "09405500 drainage area (km2)")
check_close(nrow(ann), 99, 0, "complete water years (WY1926, WY1928-2025)")
message(sprintf(
  "[%s] alpha within the 0.90-0.995 range typical of daily baseflow recessions (%.4f)",
  if_else(between(rec$alpha, 0.90, 0.995), "OK", "DIFF"), rec$alpha
))
message(sprintf(
  "[%s] published 09405500 = 09405501 + canal on most overlapping days (%.0f%%)",
  if_else(share_sum > 0.85, "OK", "DIFF"), 100 * share_sum
))
