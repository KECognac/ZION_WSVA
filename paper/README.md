# Modern-water inputs and recharge sensitivity of springs in Zion Canyon, Utah, USA

Reproducible analysis for the manuscript by K.E. Cognac and K. Willi (target journal: *Hydrogeology Journal*).
All code is R. Run every script from the repository root (`ZION_WSVA/`).

## Quick start

```r
# install packages once (versions used are listed in paper/sessionInfo.txt)
install.packages(c("dplyr", "tidyr", "purrr", "readr", "stringr", "tibble", "lubridate", "ggplot2", "patchwork",
                   "readxl", "janitor", "sf", "terra", "units", "elevatr", "nhdplusTools", "dataRetrieval", "tigris",
                   "ggspatial", "ggrepel", "ggnewscale", "scales", "jsonlite", "exactextractr", "trend", "whitebox",
                   "climateR", "magick"))
# only for the 3D figure (script 20)
install.packages(c("rayshader", "rgl", "maptiles"))
```

```bash
Rscript paper/scripts/run_all.R
```

`run_all.R` runs scripts 01-19 in order (script 20, the 3D figure, is run separately; see below). Each script prints `[OK]`/`[DIFF]` checks against the values reported in the paper.

## Run order

| Script | What it does | Paper |
|---|---|---|
| `00_setup.R` | Packages, paths, constants, shared plot theme (sourced by every script) | |
| `01_copy_raw_data.R` | Copies the original inputs into `data/park/ZION/paper/data/raw/` and writes a checksum manifest | |
| `02_download_external.R` | Downloads and caches all web inputs to `.../external/` (skipped when cached) | |
| `03_tritium.R` | GNIP thresholds, classification, primary young fractions, bomb-peak endmember | 2.3, 3.1, 4.3; Table 2, S1 |
| `04_tritium_age_models.R` | Exponential, piston-flow and binary-mixture age models | 4.3; Table S8 |
| `05_carbon14.R` | Kimball and Christensen (1996) ages, rebuilt mass-balance correction, 2025 samples, mixing check | 2.4, 3.2, 4.1, 4.3; Tables S3, S4 |
| `06_stable_isotopes.R` | Deuterium excess, snow fraction, d18O variability, Cedar City sensitivity | 2.5, 3.3; Tables S2, S6 |
| `07_major_ions.R` | Major-ion composition by flow path | 3.2, 4.1; Table S7 |
| `08_seasonal.R` | Seasonal means and within-campaign trends | 3.1; Table S5 |
| `09_site_geology_watershed.R` | Watershed, geology, site setting, side of river and elevation | 2.1; Table 1; Fig. 1 |
| `10_response_time.R` | Recharge zones, flow paths, confinement, Markagunt scenario, Monte Carlo response times | 2.6, 3.4, 4.2 |
| `11_water_levels_flowpaths.R` | Navajo water levels (smoothing spline through wells, spring emergence points and the river at the Navajo-Kayenta contact), interpreted flow, young fraction on geology | 4.2; Fig. 7a-b; ESM S4, Figs. S5, S6 |
| `12_figures_main.R` | Figs. 3, 4 and 6 (conceptual model) and ESM Fig. S4 (dual-isotope plot) | |
| `13_esm_figures.R` | ESM Figs. S1-S3; Table 3; Table S6 rows | |
| `14_baseflow_climate.R` | Baseflow of the North Fork Virgin River (USGS 09405500; river-only series excluding the Springdale Canal) vs GridMET climate: separation, recession, memory, trends | 2.8, 3.6, 4.2; ESM S7, Fig. S8; Tables S10, S11 |
| `15_spring_sensitivity_index.R` | First-order discharge sensitivity of each spring to a 10-20% recharge decline | 2.7, 3.5, 4.4; Table 4; Table S9 |
| `16_spring_catchments.R` | Spring outlet height above the river and the Navajo-Kayenta contact; local surface catchments (Whitebox) | 2.9, 4.1; ESM S5, Fig. S7 |
| `17_figure_timescales.R` | Fig. 5: response time vs water age, and baseflow memory | 3.4, 3.6, 4.2 |
| `18_figure_tritium_precip.R` | Fig. 2: GNIP tritium in precipitation, the moving modern/premodern thresholds (Lindsey et al. 2019) and the 1988-89 and 2025 samples | 2.3, 3.1 |
| `19_spring_emergence_band.R` | Spring emergence points from Google Earth imagery: elevation, geologic unit, height above the Navajo-Kayenta contact and the river; emergence elevation band | 4.1; Fig. 7 |
| `20_figure_3d_emergence.R` | 3D view of Zion Canyon (rayshader): emergence points, the emergence elevation band and where the river climbs into it. Not in `run_all.R` (full render takes a few minutes) | |

Script 20 runs after 09 and 19. A quick preview is the default; for the full-size figure set `PREVIEW <- FALSE` in the script or run

```bash
ZION_3D_FINAL=1 Rscript paper/scripts/20_figure_3d_emergence.R
```

Script 20 uses the 3D helper functions in `paper/R/map_3d.R` (copied from the CEBR_WSVA repository).

Helper functions are in `paper/R/functions.R`.

## Data layout

All data live in the symlinked data directory under `data/park/ZION/paper/data/`:

- `raw/`: copies of the original inputs (never edited), with checksums in `raw_manifest.csv`. Sources:
  - `Tritium/Park_Samples/Zion_Tritium.csv`: tritium and conductivity. 1988-89 from Kimball and Christensen (1996); 2025 from the University of Utah Dissolved Gas Laboratory (`0235 - Data Summary - Tritium.xlsm`).
  - `Tritium/GNIP WISER/*.xlsx`: IAEA/WMO GNIP precipitation (WISER; free account required, so check the IAEA terms before redistributing).
  - `Carbon_14/`: Kimball and Christensen (1996) transcriptions, Beta Analytic reports for the 2025 samples, and the correction method.
  - `isotope.csv`: stable isotopes of springs, wells and precipitation, including Cedar City event precipitation (Friedman et al. 2002).
  - `site_visit/`: 2025 field summary (GPS fixes, field parameters, sampling notes).
  - `Spring Emergence Points from google earth.kml`: spring emergence points picked on Google Earth imagery (K. Cognac, November 2026), with Google Earth elevations.
  - `geology/zion_geology_gdb_pro/`: NPS Geologic Resources Inventory digital geologic map.
- `external/`: cached downloads, which are the record of the web data the paper used:
  - NWIS gage and Navajo wells
  - NLDI basin, WBD hydrologic units, NHDPlus flowlines
  - NPS boundaries
  - AWS terrain tiles (DEMs and site elevations)
  - base maps
  - USGS tritium in precipitation by quadrangle (Jurgens 2018, doi:10.5066/P92CEFXN)
  - USGS 09405500 daily streamflow (WY1926-2025) and its NLDI basin; 09405501 (river only) and 09405499 (Springdale Canal), 1968-1991; GridMET daily climate averaged over the basin (WY1980-2025)
  - Zoom-14 terrain tiles over the canyon springs (`dem_canyon_z14.tif`) for spring catchments
  - Cuthbert et al. (2019) water table ratio, groundwater response time and perennial-stream spacing rasters (doi:10.6084/m9.figshare.7393304, CC BY 4.0), clipped on the native 1-km grid to the North Fork Virgin gage basin + 20 km and the Zion NP boundary (`cuthbert2019/`)
- `derived/`: intermediate results written by the scripts.

Delete a file in `external/` to download it again. Values can change when services update (for example new NWIS water levels). Terrain-tile point elevations on cliff faces can also differ by about 10 m between requests.

## Outputs

`paper/outputs/figures/` and `paper/outputs/tables/` hold the figures and tables (CSV) used in the manuscript and the electronic supplementary material.
Some figure files keep earlier figure numbers in their names:

| Manuscript | File | Script |
|---|---|---|
| Fig. 1 | `fig1_site_map.png` | 09 |
| Fig. 2 | `fig_tritium_precip_thresholds.png` | 18 |
| Fig. 3 | `fig2_tritium_young_fraction_by_site.png` | 12 |
| Fig. 4 | `fig3_conductivity_vs_young_fraction.png` | 12 |
| Fig. 5 | `fig4_timescales.png` | 17 |
| Fig. 6 | `fig5_conceptual_model.png` | 12 |
| Fig. 7 | `fig7_flow_emergence.png` | 19 |
| 3D figure | `fig_3d_zion_emergence_south*.png` (final, callouts, raw and preview versions) | 20 |
Monte Carlo steps use a fixed seed (`SEED` in `00_setup.R`).

## Notes

- Taylor Creek Well (Kolob Canyons) lies outside Zion Canyon. It is reported with the tritium results but not assigned a flow path.
- The primary tritium young-fraction estimate is the reported estimate. The bomb-peak endmember and the age models are sensitivity analyses (paper section 4.3).
- The manuscript itself (LaTeX/Word) is built separately in `Reports/paper_workspace/`.
