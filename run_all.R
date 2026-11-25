# Runs the full analysis from the raw data. Run from this folder (or open the
# .Rproj file first).

scripts <- c(
  "R/01_build_dataset.R",
  "R/02_equipment_sam.R",
  "R/03_equipment_spectra_figures.R",
  "R/04_equipment_permanova.R",
  "R/05_site_comparison.R",
  "R/06_summary_tables.R"
)

for (s in scripts) {
  message("Running ", s)
  source(s, local = new.env())
}

writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")
