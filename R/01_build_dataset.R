# ==============================================================================
# 1. BUILD THE SPECTRAL LIBRARY FROM RAW FILES
# ==============================================================================
# Reads the spectrometer exports (data/raw/spectra) and field sheets
# (data/raw/metadata) and saves one long table of reflectance by wavelength,
# with field metadata and classification levels attached.

source("R/functions.R")

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

library_long <- read_spectral_library("data/raw")

saveRDS(library_long, "data/processed/library_long.rds")

cat("Spectra:", n_distinct(library_long$spectra_id), "\n")
cat("Wavelengths:", min(library_long$wavelength_nm), "to", max(library_long$wavelength_nm), "nm\n")
