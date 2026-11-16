# 13_esm_figures.R -- ESM Figs. S1-S3, ESM Table S6 rows, and Table 3.
# (The response-time panel formerly drawn here as Fig. S4 is now main-text
# Fig. 5a, drawn by 17_figure_timescales.R.)
#   S1  seasonal plots by month of collection: five springs sampled in
#       1988-89 (rows) x four tracers (columns)
#   S2  snow fraction by month for the isotope groups sampled in more than
#       one season across years
#   S3  tritium input record and dependence of the young fraction on the
#       assumed young-water tritium
#
# Outputs
#   figures: figS1_seasonal.png, figS2_seasonal_isotope_groups.png,
#            figS3_tritium_models.png
#   tables:  tableS6_rows.tex, table3_synthesis.csv
# Runs after 03-06.

source("paper/scripts/00_setup.R")
suppressPackageStartupMessages(library(lubridate))

# Position of a date on a 1-13 month axis (e.g. 15 Mar -> about 3.45)
mpos <- function(d) {
  month(d) + (day(d) - 1) / as.numeric(days_in_month(d))
}

# hard-coded site colors (same as main-text figures)
SITE_COL <- c(
  "Echo Canyon Spring" = "#2060cc", "Weeping Rock" = "#cc3030",
  "The Grotto" = "#208040", "Refrigerator Spring" = "#cc7020",
  "Larry's Grotto" = "#8040cc", "Lower/Upper Grotto Complex" = "#7a7a7a",
  "Birch Spring" = "#b08020", "East Entrance Well" = "#3a8fa0"
)
MLAB <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE)
isd <- read_csv(file.path(TAB, "tableS2_isotopes.csv"), show_col_types = FALSE)


# ---- S1: seasonal plots, five springs x four tracers ------------------------

vars <- c(tritium_dc = "tritium, decay-corr. (TU)",
          yf = "young-water fraction",
          Cond = "conductivity (µS/cm)",
          fsnow = "snow fraction")

long <- bind_rows(
  tr %>%
    filter(Location %in% REPEAT_SPRINGS) %>%
    select(site = Location, date, tritium_dc, yf, Cond) %>%
    pivot_longer(-c(site, date), names_to = "var"),
  isd %>%
    filter(group %in% REPEAT_SPRINGS) %>%
    transmute(site = group, date, var = "fsnow", value = fsnow)
) %>%
  filter(!is.na(value)) %>%
  mutate(x = mpos(date),
         late = year(date) >= 2000,             # 2005-2025 samples
         site = factor(site, levels = REPEAT_SPRINGS),
         var = factor(var, levels = names(vars), labels = vars))

# monthly means, joined by a line in each panel
means <- long %>%
  group_by(site, var, m = month(date)) %>%
  summarise(x = mean(x), y = mean(value), .groups = "drop")

# hard-coded fixed y-limits for the bounded tracers (blank points)
lims <- expand_grid(site = factor(REPEAT_SPRINGS, levels = REPEAT_SPRINGS),
                    var = factor(vars, levels = vars)) %>%
  mutate(
    lo = case_when(var == vars["yf"] ~ -0.05,
                   var == vars["fsnow"] ~ 0.3,
                   var == vars["tritium_dc"] ~ 0,
                   TRUE ~ NA_real_),
    hi = case_when(var == vars["yf"] ~ 1.05,
                   var == vars["fsnow"] ~ 0.97,
                   TRUE ~ NA_real_)
  ) %>%
  pivot_longer(c(lo, hi), values_to = "y") %>%
  filter(!is.na(y))

s1 <- ggplot(long, aes(x, value)) +
  # shaded November, the month with the lowest tritium at every spring
  annotate("rect", xmin = 11, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "grey90") +
  geom_blank(data = lims, aes(x = 6, y = y)) +
  geom_line(data = means, aes(x, y, color = site), alpha = 0.6,
            linewidth = 0.4) +
  geom_point(aes(color = site, fill = site, shape = late), size = 1.3,
             stroke = 0.5) +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 1),
    labels = c("1988-89 sample", "2005-2025 sample"), name = NULL,
    guide = guide_legend(override.aes = list(fill = c("grey40", NA),
                                             color = "grey40"))
  ) +
  scale_fill_manual(values = SITE_COL, guide = "none") +
  scale_color_manual(values = SITE_COL, guide = "none") +
  scale_x_continuous(limits = c(1, 13), breaks = 1:12 + 0.5, labels = MLAB,
                     minor_breaks = 1:13) +
  facet_wrap(~ site + var, ncol = 4, scales = "free_y",
             labeller = function(l) list(paste0(l$site, ":\n", l$var))) +
  labs(x = "Month of sample collection (all years)", y = NULL) +
  theme_paper(7) +
  theme(strip.text = element_text(size = 6), legend.position = "top",
        panel.grid.major.x = element_blank())
ggsave(file.path(FIG, "figS1_seasonal.png"), s1, width = 6.8, height = 8.2,
       dpi = 450, bg = "white")


# ---- S2: snow fraction by month, other isotope groups -----------------------

extra <- c("Lower/Upper Grotto Complex", "Birch Spring", "East Entrance Well")
g2 <- isd %>%
  filter(group %in% extra) %>%
  mutate(x = mpos(date), group = factor(group, levels = extra))
m2 <- g2 %>%
  group_by(group, m = month(date)) %>%
  summarise(x = mean(x), y = mean(fsnow), .groups = "drop")

s2 <- ggplot(g2, aes(x, fsnow, color = group)) +
  annotate("rect", xmin = 11, xmax = 12, ymin = -Inf, ymax = Inf,
           fill = "grey90") +
  geom_line(data = m2, aes(x, y), alpha = 0.6, linewidth = 0.4) +
  geom_point(shape = 1, size = 1.6, stroke = 0.7) +
  scale_color_manual(values = SITE_COL, guide = "none") +
  facet_wrap(~group, nrow = 1) +
  scale_x_continuous(limits = c(1, 13), breaks = 1:12 + 0.5, labels = MLAB) +
  coord_cartesian(ylim = c(0.3, 0.97)) +
  labs(x = "Month of sample collection (all years)", y = "Snow fraction") +
  theme_paper(7) +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(FIG, "figS2_seasonal_isotope_groups.png"), s2, width = 6.8,
       height = 2.0, dpi = 450, bg = "white")


# ---- S3: tritium input and young-fraction models ----------------------------

inp <- read_csv(file.path(DER, "tritium_input_quadrangle.csv"),
                show_col_types = FALSE)
th <- read_csv(file.path(DER, "tritium_thresholds_1988-11-24.csv"),
               show_col_types = FALSE)
em <- read_csv(file.path(DER, "tritium_em_young_1988.csv"),
               show_col_types = FALSE)

# (a) precipitation record; shaded bands = sampling campaigns (hard-coded
# decimal years)
pa <- ggplot(inp, aes(t, c)) +
  annotate("rect", xmin = c(1988.6, 2025.5), xmax = c(1989.75, 2026.05),
           ymin = 0, ymax = Inf, fill = "#cc3030", alpha = 0.25) +
  geom_line(aes(linetype = held), linewidth = 0.3, color = "#333333") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dotted"),
                        guide = "none") +
  scale_y_log10(labels = scales::comma) +
  annotate("text", x = c(1987.5, 2024.5), y = 3000,
           label = c("1988-89\nsamples", "2025\nsamples"), hjust = 1,
           size = 2.2, color = "#cc3030") +
  labs(x = NULL, y = "Tritium in precipitation (TU)",
       title = "(a) Precipitation, quadrangle containing Zion") +
  theme_paper(7.5)

# (b) tritium expected for each young fraction: primary two-component line
# and the binary-mixture lines for young-water ages 2, 10 and 30 yr
f <- seq(0, 1, 0.01)
PRIMARY <- "Primary estimate (minimum post-peak)"
lines_b <- bind_rows(
  tibble(f = f, C = th$premodern + f * (th$modern - th$premodern),
         model = PRIMARY),
  # one line per young-water age Ty: tritium = f x young-water tritium
  pmap_dfr(em, function(Ty, tritium) {
    tibble(f = f, C = f * tritium,
           model = sprintf("Age model, young water %d yr", Ty))
  })
)
# site means of the 1988-89 samples, excluding Echo Canyon (> 20 TU)
sm <- tr %>%
  filter(!is.na(Tritium), period == "1988-89") %>%
  group_by(Location) %>%
  summarise(t = mean(Tritium)) %>%
  filter(t < 20)
model_cols <- c(
  "Primary estimate (minimum post-peak)" = "black",
  "Age model, young water 2 yr" = "#2060cc",
  "Age model, young water 10 yr" = "#208040",
  "Age model, young water 30 yr" = "#cc7020"
)

pb <- ggplot(lines_b, aes(100 * f, C, color = model)) +
  geom_hline(data = sm, aes(yintercept = t), linetype = "dashed",
             color = "grey60", linewidth = 0.3) +
  geom_text(data = sm, aes(x = 101, y = t, label = Location),
            inherit.aes = FALSE, hjust = 0, size = 2, color = "grey35") +
  geom_line(aes(linewidth = model == PRIMARY)) +
  scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.5),
                         guide = "none") +
  scale_color_manual(values = model_cols, breaks = names(model_cols),
                     name = NULL) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 20), clip = "off") +
  labs(x = "Young-water fraction (%)", y = "Tritium in spring water (TU)",
       title = "(b) 1988-89 samples") +
  theme_paper(7.5) +
  theme(legend.position = "bottom", plot.margin = margin(5, 60, 5, 5),
        legend.key.height = unit(0.3, "cm")) +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(FIG, "figS3_tritium_models.png"),
       pa + pb + plot_layout(widths = c(1.1, 1)),
       width = 6.8, height = 3.4, dpi = 450, bg = "white")


message("ESM figures written to ", FIG)


# ---- ESM Table S6 rows (LaTeX), formatted from the R results ----------------

s6 <- read_csv(file.path(TAB, "tableS6_snow_fraction_sensitivity.csv"),
               show_col_types = FALSE)
ord6 <- c("Weeping Rock", "Refrigerator Spring", "Larry's Grotto",
          "Birch Spring", "The Grotto", "East Entrance Well",
          "Lower/Upper Grotto Complex", "Taylor Creek Well",
          "Echo Canyon Spring")
# "median (lower--upper)" of the snow fraction
fmt <- function(d) sprintf("%.2f (%.2f--%.2f)", d$f_med, d$f_lo, d$f_hi)
# One LaTeX table row per group: both season definitions side by side
rows <- map_chr(ord6, function(g) {
  a <- s6 %>% filter(group == g, definition == "Nov-Mar vs Jun-Sep")
  b <- s6 %>% filter(group == g, definition == "Oct-Apr vs May-Sep")
  sprintf("%s%s & %d & %.1f & %.2f & %s & %s \\\\",
          g, if_else(g == KOLOB, "$^{a}$", ""), a$n, a$dex, a$f_primary,
          fmt(a), fmt(b))
})
writeLines(rows, file.path(TAB, "tableS6_rows.tex"))


# ---- Table 3: four-tracer synthesis (Taylor Creek Well excluded) ------------

s3 <- read_csv(file.path(TAB, "tableS3_c14_by_site.csv"), show_col_types = FALSE)
fs <- read_csv(file.path(DER, "snow_fraction_by_group.csv"),
               show_col_types = FALSE)
t3 <- tr %>%
  filter(Location != KOLOB) %>%
  group_by(Location) %>%
  summarise(cond_mean = mean(Cond, na.rm = TRUE),
            yf_pct = 100 * mean(yf, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(s3 %>% select(Location = Site, c14_central_age = central),
            by = "Location") %>%
  mutate(path = if_else(Location %in% PATH_A, "A", "B"),
         # isotope group used for each site's snow fraction (the Lower
         # Grotto sample uses the whole Grotto complex)
         fsnow_group = recode(Location,
                              "Lower Grotto Sandtrap" =
                                "Lower/Upper Grotto Complex")) %>%
  relocate(path, .after = yf_pct) %>%
  left_join(fs %>% select(fsnow_group = group, f_snow = f_primary,
                          n_isotope = n),
            by = "fsnow_group") %>%
  arrange(desc(cond_mean))
write_csv(t3, file.path(TAB, "table3_synthesis.csv"))
print(as.data.frame(t3), digits = 4)
