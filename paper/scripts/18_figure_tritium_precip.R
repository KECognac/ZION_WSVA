# 18_figure_tritium_precip.R -- tritium in precipitation, the moving
# thresholds and the samples (remake of the "tritium_plot" figure of
# Zion_Isotope_Analysis.Rmd).
#   a  half-year mean tritium in precipitation at the four GNIP stations, the
#      modern and premodern thresholds through time (Lindsey et al. 2019) and
#      all 24 samples
#   b  the 1988-89 samples by site against the thresholds on their dates
#   c  the 2025 samples by site against the 2025 thresholds
# Runs after 03_tritium.R, whose outputs it reads.

source("paper/scripts/00_setup.R")

gnip <- read_csv(file.path(DER, "gnip_halfyear.csv"), show_col_types = FALSE)
curve <- read_csv(file.path(DER, "tritium_threshold_curve.csv"),
                  show_col_types = FALSE)
tr <- read_csv(file.path(DER, "tritium_samples.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Tritium))

# hard-coded colors: threshold classes as in Fig. 2; one color per station
CLS <- c(Modern = "#2060cc", Mixed = "#b08020", Premodern = "#cc3030")
ST <- c("Albuquerque, New Mexico" = "#3a8fa0", "Flagstaff, Arizona" = "#8a6fc0",
        "Salt Lake City, Utah" = "#5a9a5a", "Santa Maria, California" = "#999999")

# "ALBUQUERQUE (NEW MEXICO)" -> "Albuquerque, New Mexico"
gnip <- gnip %>%
  mutate(station = str_to_title(station) %>%
           str_replace(" \\((.*)\\)$", ", \\1") %>%
           factor(levels = names(ST)))
stopifnot(!anyNA(gnip$station))

tu_labels <- function(x) format(x, big.mark = ",", drop0trailing = TRUE,
                                scientific = FALSE, trim = TRUE)


# ---- (a) Precipitation record, thresholds and samples -----------------------

# threshold labels sit on the curves at a hard-coded date (1 Jan 2010)
lab <- curve %>% filter(date == as.Date("2010-01-01"))

pa <- ggplot() +
  geom_ribbon(data = curve, aes(date, ymin = modern_min, ymax = modern_max),
              fill = "grey55", alpha = 0.3) +
  geom_line(data = gnip, aes(half, c, color = station), linewidth = 0.45) +
  geom_line(data = curve, aes(date, modern), linetype = "dashed",
            linewidth = 0.45) +
  geom_line(data = curve, aes(date, premodern), linetype = "dotted",
            linewidth = 0.55) +
  geom_point(data = tr, aes(date, Tritium), shape = 21, fill = "white",
             size = 1.3, stroke = 0.45) +
  annotate("text", x = lab$date, y = c(lab$modern * 2.6, lab$premodern / 1.9),
           label = c("Modern threshold", "Premodern threshold"), size = 2.3,
           family = "Arial") +
  annotate("text", x = as.Date(c("1988-11-15", "2024-06-01")),
           y = c(0.09, 12), label = c("1988–1989\nsamples (b)",
                                      "2025\nsamples (c)"),
           size = 2.3, lineheight = 0.9, family = "Arial") +
  scale_y_log10(limits = c(0.008, 6000), breaks = 10^(-2:3),
                labels = tu_labels) +
  scale_x_date(limits = as.Date(c("1961-01-01", "2027-01-01")),
               date_breaks = "10 years", date_labels = "%Y", expand = c(0, 0)) +
  scale_color_manual(values = ST, name = "Precipitation") +
  labs(x = NULL, y = "Tritium (TU)",
       title = "a  Tritium in precipitation, thresholds and samples") +
  theme_paper(7.5) +
  theme(legend.position = "top", legend.key.width = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0))


# ---- (b, c) Samples by site against the thresholds on their dates -----------

# Sampling round of the 1988-89 samples
tr <- tr %>%
  mutate(
    round = case_when(
      period == "2025" ~ "2025",
      date < as.Date("1988-09-15") ~ "August 1988",
      date < as.Date("1989-01-01") ~ "October–November 1988",
      TRUE ~ "February 1989"
    ),
    cls = factor(cls, levels = names(CLS))
  )
stopifnot(!anyNA(tr$round))

# Thresholds span a narrow range across the sample dates of each campaign
th <- tr %>%
  group_by(period) %>%
  summarise(mod_lo = min(modern), mod_hi = max(modern),
            pre_lo = min(premodern), pre_hi = max(premodern))

site_panel <- function(per, sites, labels, ylim, title, show_y = TRUE) {
  d <- tr %>%
    filter(period == per) %>%
    mutate(x = match(Location, sites)) %>%
    group_by(x) %>%
    mutate(dx = if (n() > 1) seq(-0.24, 0.24, length.out = n())[rank(date,
      ties.method = "first")] else 0) %>%
    ungroup()
  t1 <- th %>% filter(period == per)
  ggplot(d) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = t1$mod_lo,
             ymax = t1$mod_hi, fill = "grey30") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = t1$pre_lo,
             ymax = t1$pre_hi, fill = "grey30") +
    geom_hline(yintercept = t1$mod_lo, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = t1$pre_lo, linetype = "dotted", linewidth = 0.5) +
    geom_point(aes(x + dx, Tritium, color = cls, shape = round), size = 1.9,
               stroke = 0.7) +
    scale_x_continuous(breaks = seq_along(sites), labels = labels,
                       limits = c(0.5, length(sites) + 0.5)) +
    scale_y_log10(limits = ylim, breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10,
                                            30),
                  labels = tu_labels) +
    scale_color_manual(values = CLS, name = "Classification", drop = FALSE) +
    scale_shape_manual(
      values = c("August 1988" = 16, "October–November 1988" = 17,
                 "February 1989" = 15, "2025" = 18),
      name = "Sampled",
      breaks = c("August 1988", "October–November 1988", "February 1989"),
      guide = if (per == "2025") "none" else "legend"
    ) +
    labs(x = NULL, y = if (show_y) "Tritium (TU)" else NULL, title = title) +
    theme_paper(7.5) +
    theme(panel.grid.major.x = element_blank())
}

pb <- site_panel(
  "1988-89", REPEAT_SPRINGS,
  c("Echo Canyon\nSpring", "Weeping\nRock", "The\nGrotto",
    "Refrigerator\nSpring", "Larry's\nGrotto"),
  c(0.2, 60), "b  1988–1989 samples"
)
pc <- site_panel(
  "2025", c(ZC_2025, KOLOB),
  c("Birch\nSpring", "East Entrance\nWell", "Lower Grotto\nSandtrap",
    "Taylor Creek\nWell (Kolob)"),
  c(0.008, 3), "c  2025 samples"
)

fig <- pa / (pb + pc + plot_layout(widths = c(1.25, 1), guides = "collect") &
               theme(legend.position = "bottom", legend.box = "vertical",
                     legend.margin = margin(0, 0, 0, 0),
                     legend.spacing.y = unit(0.05, "cm"))) +
  plot_layout(heights = c(1.15, 1))
ggsave(file.path(FIG, "fig_tritium_precip_thresholds.png"), fig, width = 7.2,
       height = 5.8, dpi = 450, bg = "white")
message("Tritium precipitation figure written to ", FIG)
