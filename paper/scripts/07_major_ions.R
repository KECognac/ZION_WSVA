# 07_major_ions.R -- major-ion composition of the 1988-89 springs by flow path
# (paper sections 2.5, 3.4, 4.1; ESM Table S7).
#
# Data: Kimball and Christensen (1996) Table 2, as transcribed. Alkalinity
# (mg/L as CaCO3) is converted to HCO3 by multiplying by 1.2193
# (= 2 x 61.017 / 100.09).

source("paper/scripts/00_setup.R")

# hard-coded: molar mass (g/mol) and charge of each ion
M <- c(Ca = 40.078, Mg = 24.305, Na = 22.990, K = 39.098, Cl = 35.453,
       SO4 = 96.06, HCO3 = 61.017)
Z <- c(Ca = 2, Mg = 2, Na = 1, K = 1, Cl = 1, SO4 = 2, HCO3 = 1)

c0 <- read_csv(raw_path("Carbon_14/kimball_christensen_table2_chemistry.csv"),
               show_col_types = FALSE) %>%
  mutate(HCO3_mgL = Alk_mgL_CaCO3 * 1.2193)

# mg/L -> mmol/L and meq/L for every ion
for (k in names(M)) {
  c0[[paste0(k, "_mmol")]] <- c0[[paste0(k, "_mgL")]] / M[[k]]
  c0[[paste0(k, "_meq")]] <- c0[[paste0(k, "_mmol")]] * Z[[k]]
}

c0 <- c0 %>%
  mutate(
    cat = Ca_meq + Mg_meq + Na_meq + K_meq,
    an = Cl_meq + SO4_meq + HCO3_meq,
    charge_balance_pct = 100 * (cat - an) / (cat + an),
    NaClSO4_share = (Na_meq + Cl_meq + SO4_meq) / (cat + an),
    Na_Cl = Na_mmol / Cl_mmol,
    SO4_Cl = SO4_mmol / Cl_mmol,
    CaMg_HCO3 = (Ca_meq + Mg_meq) / HCO3_meq,
    path = ifelse(site %in% PATH_A, "A", "B")
  )

vars <- c("cond_uScm", "Na_mgL", "Cl_mgL", "SO4_mgL", "HCO3_mgL", "Ca_mgL",
          "Mg_mgL", "NaClSO4_share", "Na_Cl", "SO4_Cl", "CaMg_HCO3",
          "charge_balance_pct")
s7 <- c0 %>%
  group_by(path, site) %>%
  summarise(across(all_of(vars), list(mean = mean, min = min, max = max)),
            .groups = "drop")
write_csv(s7, file.path(TAB, "tableS7_major_ions.csv"))

# Ranges by flow path, for the text
bp <- c0 %>%
  group_by(path) %>%
  summarise(across(c(NaClSO4_share, Cl_mgL, CaMg_HCO3, HCO3_mgL,
                     charge_balance_pct),
                   list(min = min, max = max)),
            .groups = "drop")
print(as.data.frame(bp), digits = 3)


# ---- Regression checks (values in order: Path A, Path B) --------------------

check_close(c(bp$NaClSO4_share_min, bp$NaClSO4_share_max),
            c(0.62, 0.05, 0.69, 0.08), 0.006,
            "Na+Cl+SO4 share of ionic charge (A, B)")
check_close(c(bp$Cl_mgL_min, bp$Cl_mgL_max), c(180, 1.4, 293, 3.2), 0.01,
            "chloride ranges (A, B)")
check_close(c(bp$CaMg_HCO3_min, bp$CaMg_HCO3_max),
            c(1.46, 0.96, 2.16, 1.08), 0.006, "(Ca+Mg)/HCO3 (A, B)")
