# 05_carbon14.R -- carbon-14 (paper sections 2.4, 3.2, 4.1, 4.3; ESM Tables
# S3-S4).
#
# 1. Path A ages from the published Kimball and Christensen (1996) scenarios:
#    central age = mean of Plateau 3 and 4; range = all published scenarios.
# 2. Rebuilt mass-balance correction (from Carbon_14/
#    kimball_christensen_c14_correction.R and Zion_Carbon14_Methods.docx):
#    initial water in equilibrium with calcite at each scenario's PCO2; mass
#    transfer of calcite, dolomite, albite, halite, gypsum and CO2; modern
#    carbon pool at 50 pmc; half-life 5,730 yr.
#    a. Validated against the published A_nd for the four 1988-89 springs.
#    b. Applied to the 2025 samples, with final-water chemistry bracketed by
#       the Path A and Path B group means (no major ions were measured).
# 3. Binary mixing check: A_obs = f A_modern + (1 - f) A_nd exp(-lambda t_old),
#    with A_modern = 100 pmc and f = site-mean young fraction (primary
#    estimate, and the age-model fractions of section 4.3).
#
# Outputs
#   derived: c14_model_validation.csv, c14_mixing_check.csv
#   tables:  tableS3_c14_by_site.csv, tableS4_c14_2025_scenarios.csv
# Runs after 03 and 04 (young fractions).

source("paper/scripts/00_setup.R")


# ---- Mass-balance model -----------------------------------------------------

# hard-coded equilibrium constants (standard 25 C values), copied from
# raw/Carbon_14/kimball_christensen_c14_correction.R: CO2 solubility (K0),
# carbonic acid (K1, K2), calcite solubility product (Ksp)
K0 <- 10^-1.47
K1 <- 10^-6.35
K2 <- 10^-10.33
Ksp <- 10^-8.48

# Initial water (mmol/L of Ca and of total carbon) in equilibrium with calcite
# at log PCO2 = l
initial_water <- function(l) {
  co2 <- K0 * 10^l
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

# Chemistry columns (mmol/L) of the final water, in the order of the rows of S
CHEM <- c("Ca_mM", "Mg_mM", "Na_mM", "Cl_mM", "SO4_mM", "C_mM")

# Corrected initial activity (And, pmc), modeled d13C and age (years).
#   fc        named final-water chemistry (CHEM, mmol/L)
#   l         log PCO2 of the scenario
#   c13_min   d13C of the dissolving minerals; c13_gas: d13C of soil CO2
#   a14       measured carbon-14 activity (pmc)
mass_balance <- function(fc, l, c13_min, c13_gas, a14, pmc_modern = 50) {
  iw <- initial_water(l)
  # change in each solute from the initial to the final water
  dm <- c(fc[["Ca_mM"]] - iw$Ca, fc[["Mg_mM"]], fc[["Na_mM"]],
          fc[["Cl_mM"]], fc[["SO4_mM"]], fc[["C_mM"]] - iw$C)
  a <- solve(S, dm)            # mass transfer of each phase
  mineral <- a[1] + 2 * a[2]   # carbon from calcite and dolomite (dead)
  modern <- iw$C + a[6]        # carbon from soil CO2 (modern)
  and <- modern * pmc_modern / (mineral + modern)
  c(
    And = and,
    C13 = (modern * c13_gas + mineral * c13_min) / (mineral + modern),
    age = if (and > 0) C14_HALF_LIFE / log(2) * log(and / a14) else NA_real_
  )
}

# Run the mass balance for every row of a table that holds the chemistry
# columns, the scenario columns and a14; returns one column of the result
# ("And" or "age")
mass_balance_col <- function(d, what) {
  pmap_dbl(d, function(...) {
    r <- list(...)
    mass_balance(unlist(r[CHEM]), r$logPCO2atm, r$c13_dissolving_minerals,
                 r$c13_gas, r$a14)[[what]]
  })
}

scen <- read_csv(raw_path("Carbon_14/scenario_definitions.csv"),
                 show_col_types = FALSE)
pub <- read_csv(raw_path("Carbon_14/kimball_christensen_published_table.csv"),
                show_col_types = FALSE)
fchem <- read_csv(raw_path("Carbon_14/final_water_chemistry_mmolL.csv"),
                  show_col_types = FALSE)
proxy <- read_csv(raw_path("Carbon_14/group_average_chemistry_mmolL.csv"),
                  show_col_types = FALSE)

# Measured activity and d13C of each site (one row per site)
pub_site <- pub %>%
  group_by(Site) %>%
  summarise(a14 = first(C14_pmc), d13C = first(C13_permil),
            apparent_age = first(Uncorrected_Age), .groups = "drop")


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

# Every site with final-water chemistry x every scenario
val <- fchem %>%
  select(site, all_of(CHEM)) %>%
  cross_join(scen) %>%
  left_join(pub_site %>% select(site = Site, a14), by = "site") %>%
  mutate(And_model = mass_balance_col(pick(everything()), "And")) %>%
  left_join(pub %>%
              distinct(Site, Scenario, .keep_all = TRUE) %>%
              select(site = Site, scenario = Scenario, And_published = And_pmc),
            by = c("site", "scenario")) %>%
  select(site, scenario, a14, And_model, And_published)
write_csv(val, file.path(DER, "c14_model_validation.csv"))

# Rebuilt minus published A_nd (rebuilt rounded as published, to 0.1 pmc)
dv <- val %>%
  filter(!is.na(And_published)) %>%
  mutate(diff = round(And_model, 1) - And_published)
v14 <- dv %>% filter(scenario %in% c("Plateau 1", "Outcrop")) %>%
  pull(diff) %>% abs()
v24 <- dv %>% filter(scenario %in% c("Plateau 2", "Plateau 3", "Plateau 4")) %>%
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
check_new <- new %>% left_join(pub_site, by = c("site" = "Site"),
                               suffix = c("", "_pub"))
stopifnot(all.equal(check_new$a14, check_new$pmc),
          all.equal(check_new$d13C_pub, check_new$d13C))

# Every 2025 sample x every scenario x both bracketing chemistries
nb <- new %>%
  select(site, a14 = pmc) %>%
  cross_join(scen) %>%
  cross_join(proxy %>% select(proxy_chem = proxy_group, all_of(CHEM))) %>%
  mutate(And_pmc = round(mass_balance_col(pick(everything()), "And"), 1),
         C14_age = round(mass_balance_col(pick(everything()), "age"))) %>%
  select(site, scenario, proxy_chem, And_pmc, C14_age)
write_csv(nb, file.path(TAB, "tableS4_c14_2025_scenarios.csv"))

ref <- read_csv(raw_path("Carbon_14/new_samples_c14_correction.csv"),
                show_col_types = FALSE)
mm <- nb %>% inner_join(ref, by = c("site", "scenario", "proxy_chem"))
check_close(c(mm$And_pmc, mm$C14_age),
            c(mm$And_pmc_massbalance, mm$C14_age_massbalance), 3,
            "2025 corrections vs new_samples_c14_correction.csv")


# ---- 3. Mixing check --------------------------------------------------------

# Site-mean primary young fraction (03) and age-model fractions (04)
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
  mutate(
    age_no_mixing = map2_dbl(And, Aobs, t_old, f = 0),
    age_primary_f = pmap_dbl(list(And, Aobs, f), t_old),
    age_Ty2 = pmap_dbl(list(And, Aobs, BMM_f_Ty2), t_old),
    age_Ty10 = pmap_dbl(list(And, Aobs, BMM_f_Ty10), t_old),
    age_Ty20 = pmap_dbl(list(And, Aobs, BMM_f_Ty20), t_old),
    age_Ty30 = pmap_dbl(list(And, Aobs, BMM_f_Ty30), t_old)
  ) %>%
  select(Site, And, Aobs, f, starts_with("age"))
write_csv(mix, file.path(DER, "c14_mixing_check.csv"))


# ---- ESM Table S3: site-level carbon-14 -------------------------------------

s3 <- pub_site %>%
  select(Site, pmc = a14, d13C, apparent_age) %>%
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
