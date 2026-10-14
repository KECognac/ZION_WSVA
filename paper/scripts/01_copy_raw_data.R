# 01_copy_raw_data.R -- copy the original inputs from the shared ZION data
# folder into data/park/ZION/paper/data/raw/, so every later script reads only
# from the paper data directory. The originals are never modified.
# Re-running overwrites the copies (set OVERWRITE <- FALSE to keep them).

source("paper/scripts/00_setup.R")

SRC <- "data/park/ZION"   # shared folder holding the original inputs
OVERWRITE <- TRUE


# ---- Files to copy (paths relative to SRC) ----------------------------------

files <- c(
  # tritium: sample results, laboratory summary, and GNIP precipitation
  # records (IAEA/WMO WISER)
  "Tritium/Park_Samples/Zion_Tritium.csv",
  "Tritium/Park_Samples/0235 - Data Summary - Tritium.xlsm",
  file.path("Tritium/GNIP WISER",
            list.files(file.path(SRC, "Tritium/GNIP WISER"),
                       pattern = "\\.xlsx$")),
  # carbon-14: Kimball and Christensen (1996) transcriptions, 2025 laboratory
  # reports, and the correction method
  "Carbon_14/kimball_christensen_published_table.csv",
  "Carbon_14/carbon_14_kimball_christenson.csv",
  "Carbon_14/kimball_christensen_table2_chemistry.csv",
  "Carbon_14/final_water_chemistry_mmolL.csv",
  "Carbon_14/group_average_chemistry_mmolL.csv",
  "Carbon_14/scenario_definitions.csv",
  "Carbon_14/new_samples_c14_correction.csv",
  "Carbon_14/ZION East Entrance Well-Radiocarbon Dating Report.pdf",
  "Carbon_14/ZION Lower Grotto-Radiocarbon Dating Report.pdf",
  "Carbon_14/Zion_Carbon14_Methods.docx",
  "Carbon_14/kimball_christensen_c14_correction.R",
  # stable isotopes: springs, wells and precipitation, including Cedar City
  # (Friedman et al. 2002)
  "isotope.csv",
  # 2025 field visit summary: GPS fixes, field parameters, sampling notes
  "site_visit/2025_July_Zion_CEBR_PIPE_Sampling_Summary.csv"
)

for (f in files) {
  dest <- file.path(RAW, f)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (OVERWRITE || !file.exists(dest)) {
    stopifnot(file.copy(file.path(SRC, f), dest, overwrite = TRUE))
  }
}

# NPS Geologic Resources Inventory digital geologic map. A file geodatabase
# is a directory, so it is copied recursively.
gdb_src <- file.path(SRC, "geology/zion_geology_gdb_pro")
dir.create(file.path(RAW, "geology"), recursive = TRUE, showWarnings = FALSE)
if (OVERWRITE || !dir.exists(file.path(RAW, "geology/zion_geology_gdb_pro"))) {
  stopifnot(file.copy(gdb_src, file.path(RAW, "geology"),
                      recursive = TRUE, overwrite = TRUE))
}


# ---- Manifest ---------------------------------------------------------------

# Sizes and checksums, so the copies can be verified against the originals
man <- tibble(file = list.files(RAW, recursive = TRUE)) %>%
  mutate(
    bytes = file.size(file.path(RAW, file)),
    md5 = unname(tools::md5sum(file.path(RAW, file)))
  )
write_csv(man, file.path(DATA, "raw_manifest.csv"))
message(nrow(man), " raw files copied; manifest written to ",
        file.path(DATA, "raw_manifest.csv"))
