# 05_carbon14.R -- carbon-14 (paper sections 2.4, 3.2, 4.1, 4.3; ESM Tables
# S3-S4).
#
# 1. Path A ages from the published Kimball and Christensen (1996) scenarios:
#    central age = mean of Plateau 3 and 4; range = all published scenarios.
# 2. Rebuilt mass-balance correction (from Carbon_14/
#    kimball_christensen_c14_correction.R and Zion_Carbon14_Methods.docx):
#    initial water in equilibrium with calcite at each scenario's PCO2; mass
#    transfer of calcite, dolomite, albite, halite, gypsum and CO2; modern
#    carbon pool at 50 pmc; half-life 5,730 yr. Validated against the
#    published A_nd for the four 1988-89 springs, then applied to the 2025
#    samples with final-water chemistry bracketed by the Path A and Path B
#    group means.
# 3. Binary mixing check: A_obs = f A_modern + (1 - f) A_nd exp(-lambda t_old),
#    with A_modern = 100 pmc and f = site-mean young fraction (primary
#    estimate, and the age-model fractions of section 4.3).

source("paper/scripts/00_setup.R")


# ---- Mass-balance model -----------------------------------------------------

# hard-coded equilibrium constants (standard 25 C values), copied from
# raw/Carbon_14/kimball_christensen_c14_correction.R: CO2 solubility (K0),
# carbonic acid (K1, K2), calcite solubility product (Ksp)
K0 <- 10^-1.47
K1 <- 10^-6.35
K2 <- 10^-10.33
Ksp <- 10^-8.48

# Initial water (mmol/L) in equilibrium with calcite at log PCO2 = l
initial_water <- function(l) {
  p <- 10^l
  co2 <- K0 * p
  h <- (2 * Ksp * K1 * co2 / K2)^(1 / 3)
  list(Ca = h / 2 * 1000, C = (h + co2) * 1000)
}

# Stoichiometry: rows = Ca, Mg, Na, Cl, SO4, C; columns = mass transfer of
# calcite, dolomite, albite, halite, gypsum, CO2
S <- matrix(c(1, 1, 0, 0, 1, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 1, 1, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              1, 2, 0, 0, 0, 1), 6, byrow = TRUE)

# Corrected initial activity (And), modeled d13C, and age for final-water
# chemistry fc, log PCO2 l, d13C of minerals and gas, and measured activity
mass_balance <- function(fc, l, c13_min, c13_gas, a14, pmc_modern = 50) {
  iw <- initial_water(l)
  dm <- c(fc[["Ca_mM"]] - iw$Ca, fc[["Mg_mM"]], fc[["Na_mM"]],
          fc[["Cl_mM"]], fc[["SO4_mM"]], fc[["C_mM"]] - iw$C)
  a <- solve(S, dm)
  mineral <- a[1] + 2 * a[2]   # carbon from calcite and dolomite (dead)
  modern <- iw$C + a[6]        # carbon from soil CO2 (modern)
  and <- modern * pmc_modern / (mineral + modern)
  c(
    And = and,
    C13 = (modern * c13_gas + mineral * c13_min) / (mineral + modern),
    age = if (and > 0) C14_HALF_LIFE / log(2) * log(and / a14) else NA_real_
  )
}

scen <- read_csv(raw_path("Carbon_14/scenario_definitions.csv"),
                 show_col_types = FALSE)
pub <- read_csv(raw_path("Carbon_14/kimball_christensen_published_table.csv"),
                show_col_types = FALSE)
fchem <- read_csv(raw_path("Carbon_14/final_water_chemistry_mmolL.csv"),
                  show_col_types = FALSE)
proxy <- read_csv(raw_path("Carbon_14/group_average_chemistry_mmolL.csv"),
                  show_col_types = FALSE)

# Final-water chemistry vector for row i of a chemistry table
chem <- function(df, i) {
  unlist(df[i, c("Ca_mM", "Mg_mM", "Na_mM", "Cl_mM", "SO4_mM", "C_mM")])
}


# ---- 1. Path A ages (published scenarios) -----------------------------------

pathA <- pub %>%
  filter(Site %in% PATH_A) %>%
  group_by(Site) %>%
  summarise(
    central = mean(C14_age_corrected[Scenario %in% c("Plateau 3", "Plateau 4")]),
    min = min(C14_age_corrected),
    max = max(C14_age_corrected),
    .groups = "drop"
  )


# ---- 2a. Validate the rebuilt model against published A_nd (1988-89) -------

val <- expand_grid(i = seq_len(nrow(fchem)), j = seq_len(nrow(scen))) %>%
  rowwise() %>%
  mutate(
    site = fchem$site[i],
    scenario = scen$scenario[j],
    a14 = pub$C14_pmc[pub$Site == site][1],
    And_model = mass_balance(chem(fchem, i), scen$logPCO2atm[j],
                             scen$c13_dissolving_minerals[j],
                             scen$c13_gas[j], a14)[["And"]],
    And_published = pub$And_pmc[pub$Site == site &
                                  pub$Scenario == scenario][1]
  ) %>%
  ungroup() %>%
  select(-i, -j)
write_csv(val, file.path(DER, "c14_model_validation.csv"))

dv <- val %>%
  filter(!is.na(And_published)) %>%
  mutate(diff = round(And_model, 1) - And_published)
v14 <- dv %>%
  filter(scenario %in% c("Plateau 1", "Outcrop")) %>%
  pull(diff) %>%
  abs()
v24 <- dv %>%
  filter(scenario %in% c("Plateau 2", "Plateau 3", "Plateau 4")) %>%
  pull(diff)


# ---- 2b. 2025 samples -------------------------------------------------------

# hard-coded from the Beta Analytic reports: Beta-764442 Lower Grotto
# Sandtrap (22 Jul 2025) and Beta-764272 East Entrance Well (24 Jul 2025)
new <- tribble(
  ~site,                ~lab_no,       ~pmc,  ~pmc_sd, ~d13C,  ~apparent_age,
  "Lower Grotto",       "Beta-764442", 57.82, 0.21,    -9.15,  4400,
  "East Entrance Well", "Beta-764272", 83.29, 0.30,    -10.24, 1470
)
# the transcribed table must agree with the reports
stopifnot(
  all.equal(pub$C14_pmc[match(new$site, pub$Site)], new$pmc),
  all.equal(pub$C13_permil[match(new$site, pub$Site)], new$d13C)
)

nb <- expand_grid(k = seq_len(nrow(new)), j = seq_len(nrow(scen)),
                  g = seq_len(nrow(proxy))) %>%
  rowwise() %>%
  mutate(
    site = new$site[k],
    scenario = scen$scenario[j],
    proxy_chem = proxy$proxy_group[g],
    r = list(mass_balance(chem(proxy, g), scen$logPCO2atm[j],
                          scen$c13_dissolving_minerals[j], scen$c13_gas[j],
                          new$pmc[k])),
    And_pmc = round(r[["And"]], 1),
    C14_age = round(r[["age"]])
  ) %>%
  ungroup() %>%
  select(site, scenario, proxy_chem, And_pmc, C14_age)
write_csv(nb, file.path(TAB, "tableS4_c14_2025_scenarios.csv"))

ref <- read_csv(raw_path("Carbon_14/new_samples_c14_correction.csv"),
                show_col_types = FALSE)
mm <- nb %>% inner_join(ref, by = c("site", "scenario", "proxy_chem"))
check_close(c(mm$And_pmc, mm$C14_age),
            c(mm$And_pmc_massbalance, mm$C14_age_massbalance), 3,
            "2025 corrections vs new_samples_c14_correction.csv")


# ---- 3. Mixing check --------------------------------------------------------

yf <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Tritium)) %>%
  group_by(Location) %>%
  summarise(f = mean(yf), .groups = "drop")
bmm <- read_csv(file.path(TAB, "tableS8_tritium_age_models.csv"),
                show_col_types = FALSE)

# Old-water age for modern fraction f; NA when the modern contribution alone
# exceeds the measured activity. hard-coded: A_modern = 100 pmc
t_old <- function(and, aobs, f, am = 100) {
  d <- aobs - f * am
  if (d <= 0) NA_real_ else C14_HALF_LIFE / log(2) * log(and * (1 - f) / d)
}

mix <- pub %>%
  filter(Site %in% PATH_A, Scenario %in% c("Plateau 3", "Plateau 4")) %>%
  group_by(Site) %>%
  summarise(And = mean(And_pmc), Aobs = first(C14_pmc), .groups = "drop") %>%
  left_join(yf, by = c("Site" = "Location")) %>%
  left_join(bmm, by = c("Site" = "Location")) %>%
  rowwise() %>%
  mutate(
    age_no_mixing = t_old(And, Aobs, 0),
    age_primary_f = t_old(And, Aobs, f),
    age_Ty2 = t_old(And, Aobs, BMM_f_Ty2),
    age_Ty10 = t_old(And, Aobs, BMM_f_Ty10),
    age_Ty20 = t_old(And, Aobs, BMM_f_Ty20),
    age_Ty30 = t_old(And, Aobs, BMM_f_Ty30)
  ) %>%
  ungroup() %>%
  select(Site, And, Aobs, f, starts_with("age"))
write_csv(mix, file.path(DER, "c14_mixing_check.csv"))


# ---- ESM Table S3: site-level carbon-14 -------------------------------------

s3 <- pub %>%
  group_by(Site) %>%
  summarise(
    pmc = first(C14_pmc),
    d13C = first(C13_permil),
    apparent_age = first(Uncorrected_Age),
    .groups = "drop"
  ) %>%
  left_join(pathA, by = "Site")
write_csv(s3, file.path(TAB, "tableS3_c14_by_site.csv"))
print(as.data.frame(pathA))
print(as.data.frame(mix), digits = 5)


# ---- Regression checks against the paper ------------------------------------

check_close(c(pathA$central, pathA$min, pathA$max),
            c(4200, 1280, 3390, 300, 6290, 3690), 0.5,
            "Path A central ages and ranges")
check_close(c(min(v14), max(v14)), c(1.0, 1.6), 0.05,
            "rebuilt-model match, Plateau 1 and Outcrop (pmc)")
check_close(max(v24), 27.8, 0.05,
            "rebuilt-model maximum overestimate, Plateau 2-4 (pmc)")

lg <- mix %>% filter(Site == "Larry's Grotto")
wr <- mix %>% filter(Site == "Weeping Rock")
check_close(c(lg$age_no_mixing, lg$age_primary_f), c(4239, 9277), 3,
            "Larry's Grotto mixing check")
stopifnot(is.na(wr$age_primary_f))   # undefined at Weeping Rock
check_close(c(wr$age_Ty10, wr$age_Ty20, wr$age_Ty30), c(10677, 7281, 8010), 1,
            "Weeping Rock old-water age with age-model fractions")
check_close(c(lg$age_Ty20, lg$age_Ty10, lg$age_Ty2), c(5756, 6207, 10548), 1,
            "Larry's Grotto old-water age with age-model fractions")
