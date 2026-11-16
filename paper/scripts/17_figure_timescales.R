# 17_figure_timescales.R -- main-text Fig. 5: response time versus water age
# (output file fig4_timescales.png keeps the earlier figure number).
#   a  hydraulic response time by flow-path scenario (10_response_time.R)
#      against the age of the water (modern band; carbon-14 ages, 05)
#   b  memory of North Fork Virgin River baseflow: Spearman correlation of
#      annual baseflow with climate averaged over the preceding k water years
#      (14_baseflow_climate.R, river-only series, WY1980-2025)
# Runs after 10 and 14, whose outputs it reads.
# Output: figures/fig4_timescales.png

source("paper/scripts/00_setup.R")


# ---- (a) Response time vs water age -----------------------------------------

rt <- read_csv(file.path(TAB, "response_time_results.csv"),
               show_col_types = FALSE) %>%
  filter(case == "primary")
# scenario name (from 10_response_time.R) -> axis label
sc <- c(
  "Local rim recharge (Path B)" = "Local rim recharge\n(0.16-0.73 km)",
  "East plateau, nearest high ground" = "East high ground\n(8.5 km)",
  "North (Kolob), nearest high ground" = "North (Kolob) high ground\n(13.9 km)",
  "Regional, N and E high ground (median)" = "Regional high ground\n(29-40 km)",
  "Markagunt Plateau (Cedar Breaks to Brian Head)" =
    "Markagunt Plateau\n(38-46 km; hypothesized)"
)
# y position: first scenario at the top
rt <- rt %>%
  mutate(y = length(sc) + 1 - as.integer(factor(scenario, levels = names(sc))))
mx <- rt %>% filter(case_detail == "mixed/used")

# hard-coded water-age bands: modern water (< 72 yr, recharged since about
# 1953) and the carbon-14 ages and ranges of Weeping Rock (1,280 yr;
# 300-3,690) and Larry's Grotto (4,200 yr; 3,390-6,290) from 05_carbon14.R
pa <- ggplot() +
  annotate("rect", xmin = 1e-2, xmax = 72, ymin = -Inf, ymax = Inf,
           fill = "#d6e6f7") +
  annotate("rect", xmin = c(300, 3390), xmax = c(3690, 6290), ymin = -Inf,
           ymax = Inf, fill = c("#b08020", "#8040cc"), alpha = 0.13) +
  geom_vline(xintercept = c(1280, 4200), color = c("#b08020", "#8040cc"),
             linetype = "dashed", linewidth = 0.35) +
  annotate("text", x = c(0.4, 1280 / 1.12, 4200 * 1.12), y = length(sc) + 0.75,
           label = c("Modern water (recharged since about 1953)",
                     "Weeping Rock\ncarbon-14 age",
                     "Larry's Grotto\ncarbon-14 age"),
           hjust = c(0.5, 1, 0), size = 2.1,
           color = c("#2060cc", "#b08020", "#8040cc")) +
  geom_segment(data = mx, aes(x = p5, xend = p95, y = y, yend = y),
               color = "#7fa0c4", linewidth = 1.2) +
  geom_segment(data = mx, aes(x = p25, xend = p75, y = y, yend = y),
               color = "#2060cc", linewidth = 2.8) +
  geom_point(data = mx, aes(median, y), shape = 21, fill = "white", size = 2) +
  # "<" and ">" symbols: medians if fully confined / fully unconfined
  geom_point(data = rt %>% filter(case_detail == "all confined", y < length(sc)),
             aes(median, y), shape = 60, size = 3, color = "#555555") +
  geom_point(data = rt %>% filter(case_detail == "all unconfined",
                                  y < length(sc)),
             aes(median, y), shape = 62, size = 3, color = "#555555") +
  scale_x_log10(limits = c(1e-2, 1e5), breaks = 10^(-2:4),
                labels = scales::label_math(10^.x, format = log10)) +
  scale_y_continuous(breaks = rev(seq_along(sc)), labels = sc,
                     limits = c(0.4, length(sc) + 1.1)) +
  labs(x = "Time (years)", y = NULL,
       title = "a  Hydraulic response time and water age") +
  theme_paper(7.5) +
  theme(panel.grid.major.y = element_blank())


# ---- (b) Baseflow memory ----------------------------------------------------

mem <- read_csv(file.path(TAB, "tableS10_baseflow_climate_memory.csv"),
                show_col_types = FALSE) %>%
  filter(response == "Baseflow (Eckhardt)",
         series == "river-only, WY1980-2025")

pb <- ggplot(mem, aes(k, rho, color = climate)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c(P = "#6eb2e4", `P - PET` = "#dca237",
                                `Cold-day P` = "#459a74", PET = "#c07da5"),
                     breaks = c("P", "P - PET", "Cold-day P", "PET"),
                     name = NULL) +
  labs(x = "Climate averaged over the current and preceding\nwater years (k years)",
       y = "Spearman correlation with\nannual baseflow",
       title = "b  Memory of North Fork Virgin River baseflow") +
  theme_paper(7.5) +
  theme(legend.position = "bottom")

fig <- pa + pb + plot_layout(widths = c(1.45, 1))
ggsave(file.path(FIG, "fig4_timescales.png"), fig, width = 7.2, height = 3.6,
       dpi = 450, bg = "white")
message("Fig. 5 written to ", FIG)
