# 12_figures_main.R -- main-text Figs. 3, 4 and 6, and ESM Fig. S4.
# Manuscript figure numbers changed when the tritium-in-precipitation figure
# became Fig. 2; the output file names keep the earlier numbers.
#   Fig. 3   (fig2_tritium_young_fraction_by_site.png) primary young-water
#            fraction for each tritium sample, by site
#   Fig. 4   (fig3_conductivity_vs_young_fraction.png) specific conductivity
#            versus young-water fraction (flow-path discriminant)
#   Fig. 6   (fig5_conceptual_model.png) conceptual model (schematic
#            west-east section; spring sides and elevations from 09)
#   Fig. S4  (figS4_dual_isotope.png) dual-isotope plot with the snow and rain
#            endmembers (ESM)
# Other figures: Fig. 1 by 09, Fig. 2 by 18, Fig. 5 by 17, Fig. 7 by 19,
# ESM Fig. S5 by 11.
# Runs after 03, 06 and 09.

source("paper/scripts/00_setup.R")

# hard-coded colors for the tritium classes
CLS <- c(Modern = "#2060cc", Mixed = "#b08020", Premodern = "#cc3030")

tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Tritium))
site_order <- c(REPEAT_SPRINGS, "Birch Spring", "East Entrance Well",
                "Lower Grotto Sandtrap", KOLOB)

# Percent label: one decimal for values between 0 and 1%
pct <- function(m) {
  if_else(m > 0 & m < 1, sprintf("%.1f%%", m), sprintf("%.0f%%", m))
}


# ---- Fig. 3: young fraction by sample and site ------------------------------

# One row per site (Taylor Creek Well set 0.6 lower, below a dashed line);
# repeat samples spread vertically within the row (`jit`)
d2 <- tr %>%
  mutate(Location = factor(Location, levels = rev(site_order)),
         y = as.numeric(Location) - if_else(Location == KOLOB, 0.6, 0)) %>%
  group_by(Location) %>%
  mutate(jit = if (n() > 1) seq(-0.18, 0.18, length.out = n()) else 0) %>%
  ungroup()
m2 <- d2 %>%
  group_by(Location, y) %>%
  summarise(m = 100 * mean(yf), n = n(), .groups = "drop")
labs2 <- m2 %>%
  arrange(desc(y)) %>%
  mutate(lab = if_else(
    Location == KOLOB,
    "Taylor Creek Well (n=1)\nKolob Canyons, outside Zion Canyon",
    sprintf("%s (n=%d)", Location, n)
  ))

f2 <- ggplot(d2) +
  geom_hline(yintercept = min(d2$y) + 0.8, color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(100 * yf, y + jit, color = cls, shape = Location == KOLOB),
             size = 1.9, stroke = 0.7) +
  geom_segment(data = m2, aes(x = m, xend = m, y = y - 0.32, yend = y + 0.32),
               color = "#222222", linewidth = 0.7) +
  geom_text(data = m2, aes(x = 103, y = y, label = pct(m)), hjust = 0,
            size = 2.6, family = "Arial") +
  scale_color_manual(values = CLS, name = "Explanation", breaks = names(CLS)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_y_continuous(breaks = labs2$y, labels = labs2$lab) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(-3, 112)) +
  labs(x = "Primary young-water fraction (%)", y = NULL) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")
ggsave(file.path(FIG, "fig2_tritium_young_fraction_by_site.png"), f2,
       width = 6.8, height = 2.6, dpi = 450, bg = "white")


# ---- Fig. 4: conductivity vs young fraction ---------------------------------

# hard-coded point shapes and colors by site
sp <- c("Echo Canyon Spring" = 21, "Weeping Rock" = 22, "The Grotto" = 24,
        "Refrigerator Spring" = 25, "Larry's Grotto" = 23,
        "Birch Spring" = 3, "East Entrance Well" = 4,
        "Lower Grotto Sandtrap" = 8, "Taylor Creek Well" = 1)
pal <- c("#2060cc", "#cc3030", "#208040", "#cc7020", "#8040cc", "#b08020",
         "#3a8fa0", "#a0507a", "#666666")
names(pal) <- names(sp)
sp_labels <- c(head(names(sp), -1), "Taylor Creek Well (Kolob,\nnot assigned)")

f3 <- ggplot(tr %>% mutate(Location = factor(Location, levels = names(sp))),
             aes(100 * yf, Cond, shape = Location, color = Location,
                 fill = Location)) +
  # 1,000 uS/cm: Path A / Path B boundary (Kimball and Christensen 1996)
  geom_hline(yintercept = 1000, linetype = "dashed", color = "grey35",
             linewidth = 0.35) +
  annotate("text", x = 38, y = 1080, label = "Path A", size = 2.5,
           color = "grey30", family = "Arial") +
  annotate("text", x = 38, y = 920, label = "Path B", size = 2.5,
           color = "grey30", family = "Arial") +
  geom_point(size = 1.9, stroke = 0.6) +
  scale_y_log10(limits = c(100, 2500), breaks = c(100, 200, 500, 1000, 2000),
                labels = scales::comma) +
  scale_shape_manual(values = sp, labels = sp_labels, name = "Explanation") +
  scale_color_manual(values = pal, labels = sp_labels, name = "Explanation") +
  scale_fill_manual(values = pal, labels = sp_labels, name = "Explanation") +
  scale_x_continuous(limits = c(-5, 105)) +
  labs(x = "Primary young-water fraction (%)",
       y = "Specific conductivity (µS/cm)") +
  theme_paper(8) +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(0.32, "cm"))
ggsave(file.path(FIG, "fig3_conductivity_vs_young_fraction.png"), f3,
       width = 5.2, height = 2.8, dpi = 450, bg = "white")


# ---- ESM Fig. S4: dual-isotope plot ----------------------------------------------

iso <- read_csv(raw_path("isotope.csv"), show_col_types = FALSE) %>%
  # site name = text before the first comma of the Location field
  mutate(loc = str_trim(str_remove(replace_na(Location, ""), ",.*")))
s2 <- read_csv(file.path(TAB, "tableS2_isotopes.csv"), show_col_types = FALSE)

# Endmember samples (see 06_stable_isotopes.R)
snow <- iso %>% filter(loc == "SNOW PILE NEAR WEEPING ROCK")
rain <- iso %>% filter(str_starts(loc, "Birch Creek Precipitation"))
# hard-coded: the one rain-leaning Echo Canyon sample, labeled on the plot
echo3 <- s2 %>% filter(group == "Echo Canyon Spring", Collection_Date == "11/3/88")

gp <- c("Weeping Rock" = "#cc3030", "Refrigerator Spring" = "#cc7020",
        "Larry's Grotto" = "#8040cc", "Birch Spring" = "#b08020",
        "The Grotto" = "#208040", "East Entrance Well" = "#3a8fa0",
        "Lower/Upper Grotto Complex" = "#a0507a",
        "Taylor Creek Well" = "#666666", "Echo Canyon Spring" = "#2060cc")
# Small text note at (x, y); hard-coded offsets (per mil) chosen by eye
note <- function(x, y, label, hjust) {
  annotate("text", x = x, y = y, label = label, hjust = hjust, size = 2.2,
           family = "Arial")
}

f4 <- ggplot(s2 %>% mutate(group = factor(group, levels = names(gp))),
             aes(delta_18O, delta_2H)) +
  # Global Meteoric Water Line, d2H = 8.02 d18O + 9.47
  geom_abline(intercept = 9.47, slope = 8.02, color = "grey35",
              linewidth = 0.4) +
  geom_point(aes(color = group, shape = group == KOLOB), size = 1.4,
             stroke = 0.6) +
  geom_point(data = bind_rows(snow, rain), shape = 8, size = 2.2,
             color = "#222222") +
  note(snow$delta_18O + 0.3, snow$delta_2H - 3,
       "Snow endmember (d-excess 16.8)", hjust = 0) +
  note(rain$delta_18O - 0.4, rain$delta_2H + 7,
       "Rain endmember\n(d-excess −10.3)", hjust = 1) +
  note(echo3$delta_18O + 0.4, echo3$delta_2H - 8,
       "Echo Canyon, 3 Nov 1988\n(snow fraction 0.336)", hjust = 0) +
  scale_color_manual(
    values = gp, name = "Explanation",
    guide = guide_legend(override.aes = list(
      shape = if_else(names(gp) == KOLOB, 1, 16)
    ))
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  coord_cartesian(xlim = c(-17, -0.5), ylim = c(-118, 2)) +
  labs(x = "δ¹⁸O (‰ VSMOW)", y = "δ²H (‰ VSMOW)") +
  theme_paper(8) +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))
ggsave(file.path(FIG, "figS4_dual_isotope.png"), f4, width = 5.2, height = 3.0,
       dpi = 450, bg = "white")


# ---- Fig. 6: conceptual model -----------------------------------------------
# Schematic west-east section. x runs 0-100 (arbitrary units); y is elevation
# (m). Spring elevations and sides of the river come from 09; everything else
# below (layer contacts, plateau surfaces, flow arrows, label positions) is
# hard-coded for illustration only.

se <- read_csv(file.path(DER, "site_side_elevation.csv"), show_col_types = FALSE)
yfs <- tr %>%
  group_by(Location) %>%
  summarise(m = 100 * mean(yf), .groups = "drop")

# hard-coded schematic elevations (m) of the layer tops
NK <- 1390      # Navajo-Kayenta contact
NT <- 1990      # Navajo-Temple Cap contact
TC <- 2030      # top of Temple Cap
CK <- 2120      # top of Carmel / plateau surface
FLOOR <- 1320   # canyon floor
XW <- 33        # west canyon wall (x)
XE <- 67        # east canyon wall (x)

xs <- seq(0, 100, length.out = 801)

# Land surface: plateau (slightly higher at both ends) cut by the canyon.
# `cany` rises from 0 on the plateau to 1 inside the canyon over 1.2 x-units
# at each wall.
surf <- function(x) {
  y <- case_when(x < 10 ~ CK + 45, x > 88 ~ CK + 80, TRUE ~ CK)
  cany <- pmin(pmax((x - XW) / 1.2, 0, na.rm = TRUE), pmax((XE - x) / 1.2, 0))
  cany <- pmin(pmax(cany, 0), 1)
  y - cany * (y - FLOOR)
}
top <- surf(xs)

# Polygon for a layer between elevations lo and hi, clipped at the surface
lay <- function(lo, hi, unit) {
  yb <- pmin(lo, top)
  yt <- pmin(hi, top)
  tibble(x = c(xs, rev(xs)), y = c(yb, rev(yt)), unit = unit)
}
LC <- c(K = "#a9b38a", Jc = "#7f9fbf", Jt = "#b9a3c9", Jn = "#efc9a4",
        Jk = "#b5654d")
layers <- bind_rows(lay(1240, NK, "Jk"), lay(NK, NT, "Jn"), lay(NT, TC, "Jt"),
                    lay(TC, CK, "Jc"), lay(CK, CK + 90, "K"))

# Spring positions: on the wall of their side of the river, at their
# elevation; ly = hard-coded label heights (m), staggered to avoid overlap
label_y <- tribble(
  ~Location,               ~ly,
  "Echo Canyon Spring",    1600,
  "Larry's Grotto",        1530,
  "The Grotto",            1465,
  "Weeping Rock",          1405,
  "Lower Grotto Sandtrap", 1345,
  "Refrigerator Spring",   1780,
  "Birch Spring",          1430
)
pos <- se %>%
  filter(Location != "East Entrance Well") %>%
  mutate(x = case_when(Location == "Lower Grotto Sandtrap" ~ XE - 0.4,
                       side == "east" ~ XE + 0.2,
                       TRUE ~ XW - 0.2),
         y = elev_m) %>%
  left_join(yfs, by = "Location") %>%
  left_join(label_y, by = "Location") %>%
  mutate(name = str_remove(Location, " Sediment Trap| Sandtrap"),
         lx = if_else(side == "east", XE - 1.4, XW + 1.4),
         lab = sprintf("%s %s\n(%s m)", name, pct(m), scales::comma(round(y))))

# Arrow along a path of points
P <- function(x, y, col, lw = 0.75, lty = "solid") {
  geom_path(data = tibble(x = x, y = y), aes(x, y), color = col,
            linewidth = lw, linetype = lty,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed"))
}
# x and y of one spring (the end point of its arrow)
pp <- function(l) pos %>% filter(Location == l) %>% select(x, y)
PA <- "#8040cc"   # Path A color
PB <- "#208040"   # Path B color
PE <- "#2060cc"   # Echo Canyon / side-canyon color
w <- pp("Weeping Rock")
l <- pp("Larry's Grotto")
g <- pp("The Grotto")
lg <- pp("Lower Grotto Sandtrap")
e <- pp("Echo Canyon Spring")
r <- pp("Refrigerator Spring")
b <- pp("Birch Spring")

sec <- ggplot() +
  geom_polygon(data = layers, aes(x, y, group = unit, fill = unit), color = NA) +
  scale_fill_manual(values = LC, guide = "none") +
  geom_line(data = tibble(x = xs, y = top), aes(x, y), color = "#333333",
            linewidth = 0.35) +
  annotate("segment", x = 0, xend = XW + 0.2, y = NK, yend = NK,
           color = "#cc3030", linewidth = 0.6) +
  annotate("segment", x = XE - 0.2, xend = 100, y = NK, yend = NK,
           color = "#cc3030", linewidth = 0.6) +
  annotate("rect", xmin = 47, xmax = 53, ymin = FLOOR - 4, ymax = FLOOR + 4,
           fill = "#3a8fa0") +
  # flow arrows (hard-coded schematic paths ending at each spring)
  P(c(97, 93, 87, 79, w$x), c(2195, 1960, 1720, 1490, w$y), PA) +
  P(c(87, 78, l$x), c(1720, 1600, l$y), PA, 0.6) +
  P(c(103.5, 95), c(1830, 1805), PA, 0.5, "dashed") +
  P(c(73.5, 72.5, 69.5, g$x), c(2110, 1640, 1480, g$y), PB) +
  P(c(71, 70.5, 68.6, lg$x), c(2110, 1560, 1400, lg$y), PB, 0.6) +
  P(c(25, 29.5, r$x), c(2110, 1920, r$y), PB) +
  P(c(19, 27, b$x), c(2110, 1700, b$y), PB, 0.6) +
  P(c(80, 76, 70, e$x), c(2110, 1800, 1640, e$y), PE, 0.6) +
  P(c(76, 75.5, w$x + 0.4), c(1800, 1520, w$y + 8), PE, 0.5, "dashed") +
  geom_segment(data = pos, aes(x = x, y = y, xend = lx, yend = ly),
               color = "grey40", linewidth = 0.2) +
  geom_point(data = pos, aes(x, y), shape = 21, fill = "white", size = 2,
             stroke = 0.6) +
  geom_label(data = pos,
             aes(lx, ly, label = lab, hjust = if_else(side == "east", 1, 0)),
             size = 1.85, family = "Arial", label.size = 0,
             fill = alpha("white", 0.85), label.padding = unit(0.08, "lines"),
             lineheight = 0.9) +
  # text labels (hard-coded positions)
  annotate("text", x = 50, y = 2262,
           label = "Snowmelt-dominated recharge on the plateaus", color = PE,
           size = 2.3, family = "Arial") +
  annotate("text", x = c(1, 98.5), y = 2262, label = c("W", "E"),
           fontface = "bold", size = 3, hjust = 0) +
  annotate("text", x = c(1.5, 1.5, 11, 89.5),
           y = c(1700, NT + 20, CK - 45, CK + 60),
           label = c("Navajo Sandstone\n(principal aquifer)", "Temple Cap Fm.",
                     "Carmel Formation", "Cretaceous"),
           hjust = 0, size = c(2.5, 2, 2, 2), family = "Arial") +
  annotate("text", x = 15, y = 1262,
           label = "Kayenta Formation (less permeable)", color = "white",
           hjust = 0, size = 2.5, family = "Arial") +
  annotate("text", x = 99, y = NK - 14,
           label = "Navajo-Kayenta contact (schematic)", color = "#cc3030",
           hjust = 1, vjust = 1, size = 2, family = "Arial") +
  annotate("text", x = 50, y = FLOOR - 12, label = "North Fork\nVirgin River",
           color = "white", vjust = 1, size = 2, family = "Arial") +
  annotate("text", x = 77.2, y = 1830, label = "side-canyon\nrecharge",
           color = PE, hjust = 0, vjust = 0, size = 2, family = "Arial") +
  annotate("text", x = 99.5, y = 1785,
           label = "regional flow from\nhigher plateaus\n(hypothesized)",
           color = PA, hjust = 1, vjust = 1, size = 2, family = "Arial") +
  scale_y_continuous(breaks = seq(1300, 2200, 100), labels = scales::comma,
                     name = "Elevation (m)", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(1240, 2290), clip = "off") +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Text card below the section
card <- function(title, body, col) {
  ggplot() +
    xlim(0, 10) +
    ylim(0, 10) +
    theme_void() +
    annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.3, ymax = 9.7,
             fill = "white", color = col, linewidth = 0.5) +
    annotate("text", x = 5, y = 8.8, label = title, fontface = "bold",
             color = col, size = 2.4, family = "Arial") +
    annotate("text", x = 5, y = 4.3, label = body, size = 2.05,
             family = "Arial", lineheight = 1.05)
}
cB <- card(
  "Path B (all other sites)",
  paste0("Local recharge on Navajo outcrop near\n",
         "the rims, on both sides of the canyon;\n",
         "short flow paths\n",
         "Conductivity 139-394 µS/cm\n",
         "Ca-Mg-HCO₃ water (carbonate cement)\n",
         "Substantial modern fractions at the\n",
         "springs sampled in 1988-89"),
  PB
)
cA <- card(
  "Path A (Weeping Rock, Larry's Grotto)",
  paste0("Recharge on the plateau, through or along\n",
         "the Carmel Formation; long flow paths\n",
         "Conductivity >1,000 µS/cm\n",
         "Na-Cl-SO₄ water\n",
         "Carbon-14 ages about 1,300-4,200 yr,\n",
         "mixed with modern water; Weeping Rock\n",
         "may also receive side-canyon recharge\n",
         "(hypothesized)"),
  PA
)

# Spring-box schematic (park staff observation, July 2025 site visit)
box <- ggplot() +
  xlim(0, 10) +
  ylim(0, 10) +
  theme_void() +
  annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.3, ymax = 9.7,
           fill = "white", color = "grey35", linewidth = 0.5) +
  annotate("text", x = 5, y = 8.9, label = "Spring box (schematic)",
           fontface = "bold", size = 2.4, family = "Arial") +
  annotate("rect", xmin = 1.2, xmax = 4.2, ymin = 3.6, ymax = 6.8,
           fill = "#cfe3f5", color = "#333333", linewidth = 0.4) +
  annotate("text", x = 2.7, y = 5.2, label = "full\nyear-round", size = 2,
           family = "Arial") +
  annotate("segment", x = 0.5, y = 7.4, xend = 1.2, yend = 5.6,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = PE) +
  annotate("text", x = 0.6, y = 7.6, label = "spring flow", hjust = 0,
           vjust = 0, size = 1.9, color = PE, family = "Arial") +
  annotate("segment", x = 4.2, y = 5.2, xend = 6.0, yend = 5.2,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = 6.2, y = 5.2, label = "captured:\npark water\nsupply",
           hjust = 0, size = 1.95, family = "Arial") +
  annotate("segment", x = 2.7, y = 3.6, xend = 2.7, yend = 1.3,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = PB) +
  annotate("text", x = 3.1, y = 2.1,
           label = "uncaptured: varies\nseasonally; available to\nspring ecosystems",
           hjust = 0, size = 1.95, color = "#145a28", family = "Arial")

caption <- str_c(
  "Schematic west-east section; spring elevations and sides of the river ",
  "are from site coordinates; layer contacts are schematic.\n",
  "Percentages are site-mean primary young-water fractions; Lower Grotto ",
  "Sandtrap and Birch Spring are single 2025 samples."
)
f5 <- sec / (cB | box | cA) +
  plot_layout(heights = c(3.3, 1.15)) +
  plot_annotation(
    caption = caption,
    theme = theme(plot.caption = element_text(size = 6, color = "grey35",
                                              hjust = 0.5, family = "Arial"))
  )
ggsave(file.path(FIG, "fig5_conceptual_model.png"), f5, width = 7.2,
       height = 6.0, dpi = 450, bg = "white")
message("Figs. 3, 4, 6 and S4 written to ", FIG)
