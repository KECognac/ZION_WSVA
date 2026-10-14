# 00_setup.R -- packages, paths and shared settings for the Zion multi-tracer
# paper. Every numbered script sources this file. Run all scripts from the
# repository root (ZION_WSVA/), for example:
#   Rscript paper/scripts/03_tritium.R
#
# Data layout (all under the symlinked data directory):
#   data/park/ZION/paper/data/raw/       copies of the original inputs
#                                        (read-only)
#   data/park/ZION/paper/data/external/  cached downloads from web services
#   data/park/ZION/paper/data/derived/   intermediate results written by the
#                                        scripts
# Outputs (figures and tables for the manuscript) go to paper/outputs/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})


# ---- Paths ------------------------------------------------------------------

DATA <- "data/park/ZION/paper/data"
RAW  <- file.path(DATA, "raw")
EXT  <- file.path(DATA, "external")
DER  <- file.path(DATA, "derived")
OUT  <- "paper/outputs"
FIG  <- file.path(OUT, "figures")
TAB  <- file.path(OUT, "tables")
for (d in c(RAW, EXT, DER, FIG, TAB)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ---- Constants --------------------------------------------------------------

UTM <- 26912            # NAD83 / UTM zone 12N, used for all areas and distances
H3_HALF_LIFE <- 12.32   # tritium half-life, years (Lucas and Unterweger 2000)
C14_HALF_LIFE <- 5730   # carbon-14 half-life used by Kimball and Christensen
                        # (1996)
SEED <- 20261006        # fixed seed for every Monte Carlo step


# ---- Site groups used throughout --------------------------------------------

# hard-coded site names; they must match the Location column of
# Zion_Tritium.csv
REPEAT_SPRINGS <- c("Echo Canyon Spring", "Weeping Rock", "The Grotto",
                    "Refrigerator Spring", "Larry's Grotto")   # sampled 1988-89
ZC_2025 <- c("Birch Sediment Trap", "East Entrance Well",
             "Lower Grotto Sandtrap")                          # sampled 2025
# Outside Zion Canyon: reported, never assigned a flow path
KOLOB <- "Taylor Creek Well"
# Conductivity > 1,000 uS/cm (Kimball and Christensen 1996)
PATH_A <- c("Larry's Grotto", "Weeping Rock")


# ---- Figure theme shared by all non-map figures -----------------------------

theme_paper <- function(base = 8.5) {
  theme_minimal(base_family = "Arial", base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      axis.line = element_line(color = "grey30", linewidth = 0.3),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(size = base + 0.5)
    )
}

source("paper/R/functions.R")
