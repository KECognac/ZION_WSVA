# 01_copy_raw_data.R -- copy the original inputs from the shared ZION data
# folder into data/park/ZION/paper/data/raw/, so every later script reads only
# from the paper data directory. The originals are never modified.
#
# A file is copied only if its copy is missing or differs from the original
# (MD5 checksum). Rewriting unchanged files in a OneDrive-synced folder can
# leave empty or conflicting copies, so unchanged copies are left alone.
# Office files (.xlsx, .xlsm, .docx) are the exception: OneDrive adds
# sensitivity-label metadata to copies written into the folder, so their
# checksums never match the originals even though the data are identical.
# They are copied only if missing or if the original is newer than the copy.
# Each copy is written to a temporary file next to its destination, checked
# against the original, and only then renamed into place.
#
# Outputs: the copies under RAW, and DATA/raw_manifest.csv (size and checksum
# of every raw file).

source("paper/scripts/00_setup.R")

SRC <- "data/park/ZION"   # shared folder holding the original inputs


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
  "site_visit/2025_July_Zion_CEBR_PIPE_Sampling_Summary.csv",
  # spring emergence points picked from Google Earth imagery (K. Cognac,
  # November 2026)
  "Spring Emergence Points from google earth.kml"
)

# NPS Geologic Resources Inventory digital geologic map. A file geodatabase is
# a folder of many files; every file in it is listed, except ArcGIS lock files
# (*.lock), which only mark an open session.
gdb <- "geology/zion_geology_gdb_pro"
gdb_files <- list.files(file.path(SRC, gdb), recursive = TRUE) %>%
  discard(~ str_detect(.x, "\\.lock$")) %>%
  file.path(gdb, .)

to_copy <- tibble(file = c(files, gdb_files)) %>%
  mutate(src = file.path(SRC, file), dest = file.path(RAW, file))
stopifnot(all(file.exists(to_copy$src)))


# ---- Copy only missing or changed files -------------------------------------

md5 <- function(p) unname(tools::md5sum(p))

to_copy <- to_copy %>%
  mutate(
    md5_src = md5(src),
    md5_dest = if_else(file.exists(dest), md5(dest), NA_character_),
    office = str_detect(file, "\\.(xlsx|xlsm|docx)$"),
    needs_copy = case_when(
      is.na(md5_dest) ~ TRUE,                         # missing
      office ~ file.mtime(src) > file.mtime(dest),    # original is newer
      TRUE ~ md5_dest != md5_src                      # contents differ
    )
  )

# Copy one file through a temporary file, then check it before renaming
safe_copy <- function(src, dest, md5_src) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".tmp_copy")   # same folder, so the rename is atomic
  ok <- file.copy(src, tmp, overwrite = TRUE) && md5(tmp) == md5_src
  if (!ok) {
    unlink(tmp)
    stop("Copy failed or is incomplete: ", src)
  }
  stopifnot(file.rename(tmp, dest))
}

changed <- to_copy %>% filter(needs_copy)
pwalk(list(changed$src, changed$dest, changed$md5_src), safe_copy)
message(nrow(changed), " of ", nrow(to_copy), " files copied (",
        nrow(to_copy) - nrow(changed), " already up to date)")
if (nrow(changed) > 0) message("  ", paste(changed$file, collapse = "\n  "))


# ---- Manifest ---------------------------------------------------------------

# Size and checksum of every raw file, so the copies can be verified against
# the originals
manifest <- tibble(file = list.files(RAW, recursive = TRUE)) %>%
  mutate(path = file.path(RAW, file),
         bytes = file.size(path),
         md5 = md5(path)) %>%
  select(-path)
write_csv(manifest, file.path(DATA, "raw_manifest.csv"))
message("manifest of ", nrow(manifest), " raw files written to ",
        file.path(DATA, "raw_manifest.csv"))
