# 15_spring_sensitivity_index.R -- first-order discharge sensitivity of each
# spring to a decline in recharge.
#
# On timescales shorter than the response time of the premodern component, a
# fractional change in recharge (dR/R) acts mainly through the modern
# component, so the initial fractional change in discharge is
# dQ/Q ~ f_young * dR/R, and dQ/Q approaches dR/R once the whole flow system
# has re-equilibrated (section 4.4 response times).
# This is a scenario index, not a prediction: it assumes the modern
# component's discharge scales with recharge.
#
# Recharge scenarios: 10% and 20% declines, the range of average projected
# declines in total recharge for aquifers in the southern part of the western
# USA (Meixner et al. 2016).
#
# Main table (Table 4): PRIMARY young fraction only (an upper bound on the
# modern fraction, section 4.3). ESM table (Table S9): the same index from
# the binary-mixing age models (young-component mean age 2-30 years).
# Taylor Creek Well (Kolob Canyons) is excluded, as from the four-tracer
# synthesis.

source("paper/scripts/00_setup.R")

DR <- c(0.10, 0.20)   # hard-coded recharge declines (Meixner et al. 2016)
sites <- c(REPEAT_SPRINGS, ZC_2025)


# ---- Table 4: primary young fraction ----------------------------------------

prim <- read_csv(file.path(TAB, "table2_tritium_by_site.csv"),
                 show_col_types = FALSE) %>%
  filter(Location %in% sites) %>%
  transmute(Location,
            n,
            f_mean = yf_mean_pct / 100,
            f_min = yf_min_pct / 100,
            f_max = yf_max_pct / 100)

# Initial discharge decline (%) for each recharge scenario, using the site
# mean and the range across sample dates
tab4 <- prim %>%
  mutate(dQ_10_mean = 100 * f_mean * DR[1],
         dQ_20_mean = 100 * f_mean * DR[2],
         dQ_10_min = 100 * f_min * DR[1],
         dQ_10_max = 100 * f_max * DR[1],
         dQ_20_min = 100 * f_min * DR[2],
         dQ_20_max = 100 * f_max * DR[2],
         dQ_long_term = paste0(100 * DR[1], "-", 100 * DR[2])) %>%
  arrange(desc(f_mean))

write_csv(tab4, file.path(TAB, "table4_sensitivity_index.csv"))

tab4 %>%
  transmute(
    Location,
    n,
    young_pct = round(100 * f_mean),
    `dQ (10% dR)` = sprintf("%.1f (%.1f-%.1f)", dQ_10_mean, dQ_10_min, dQ_10_max),
    `dQ (20% dR)` = sprintf("%.1f (%.1f-%.1f)", dQ_20_mean, dQ_20_min, dQ_20_max)
  ) %>%
  as.data.frame() %>%
  print()


# ---- ESM Table S9: age-model young fractions --------------------------------
# Binary mixing model, young-component mean age 2-30 years (04)

am <- read_csv(file.path(TAB, "tableS8_tritium_age_models.csv"),
               show_col_types = FALSE) %>%
  filter(Location %in% sites) %>%
  select(Location, starts_with("BMM_f_Ty")) %>%
  pivot_longer(-Location, names_to = "young_age", values_to = "f") %>%
  group_by(Location) %>%
  summarise(f_lo = min(f), f_hi = max(f), .groups = "drop") %>%
  mutate(dQ_10_lo = 100 * f_lo * DR[1],
         dQ_10_hi = 100 * f_hi * DR[1],
         dQ_20_lo = 100 * f_lo * DR[2],
         dQ_20_hi = 100 * f_hi * DR[2]) %>%
  arrange(desc(f_hi))
write_csv(am, file.path(TAB, "tableS9_sensitivity_index_age_models.csv"))


# ---- Checks: the index is a direct product of the primary fractions ---------

wr <- tab4 %>% filter(Location == "Weeping Rock")
check_close(c(wr$dQ_10_mean, wr$dQ_20_mean), c(7.35, 14.69), 0.01,
            "Weeping Rock initial discharge decline for 10% and 20% recharge decline (%)")
check_close(tab4$dQ_20_mean[tab4$Location == "Echo Canyon Spring"], 20, 1e-9,
            "Echo Canyon Spring (entirely modern) declines with recharge")
