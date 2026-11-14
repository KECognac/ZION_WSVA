# run_all.R -- run the full Zion multi-tracer analysis in order, from the
# repository root:
#   Rscript paper/scripts/run_all.R
#
# Steps 01 (copy raw data) and 02 (cache downloads) only need to run once;
# 01 is skipped when the raw manifest exists, and 02 skips every cached file.
# Each script prints [OK]/[DIFF] regression checks against the values
# reported in the paper; they are collected here.

scripts <- c(
  "01_copy_raw_data.R",
  "02_download_external.R",
  "03_tritium.R",
  "04_tritium_age_models.R",
  "05_carbon14.R",
  "06_stable_isotopes.R",
  "07_major_ions.R",
  "08_seasonal.R",
  "09_site_geology_watershed.R",
  "10_response_time.R",
  "11_water_levels_flowpaths.R",
  "12_figures_main.R",
  "13_esm_figures.R",
  "14_baseflow_climate.R",
  "15_spring_sensitivity_index.R",
  "16_spring_catchments.R",
  "17_figure_timescales.R",
  "18_figure_tritium_precip.R"
)

skip_setup <- file.exists("data/park/ZION/paper/data/raw_manifest.csv")

for (s in scripts) {
  if (skip_setup && s == "01_copy_raw_data.R") {
    message("== skip ", s, " (raw data already copied)")
    next
  }
  message("== ", s)
  t0 <- Sys.time()

  # each script runs in its own R session; stdout and stderr are captured
  out <- system2("Rscript", file.path("paper/scripts", s),
                 stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status")
  checks <- grep("^\\[(OK|DIFF)\\]", out, value = TRUE)
  if (length(checks)) message(paste(checks, collapse = "\n"))
  if (!is.null(st) && st != 0) {
    message(paste(tail(out, 20), collapse = "\n"))
    stop(s, " failed")
  }
  message(sprintf("   done in %.0f s",
                  as.numeric(Sys.time() - t0, units = "secs")))
}

writeLines(capture.output(sessionInfo()), "paper/sessionInfo.txt")
message("All steps complete. Figures and tables: paper/outputs/; ",
        "intermediate data: data/park/ZION/paper/data/derived/")
