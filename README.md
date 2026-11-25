# Coral reef spectral library: equipment and site comparisons

Data and R code for Smart et al., *Building Transferable Spectral Libraries for Coral Reef Benthic Habitat Mapping: Equipment and Geographic Considerations*.

In situ reflectance spectra of coral reef and seagrass benthic features were collected with an ASD HandHeld 2 spectrometer at Heron Reef (Queensland), Coral Bay (Western Australia) and Eastern Banks, Moreton Bay (Queensland). The analysis compares:

1. **Equipment configuration**: ambient vs artificial (LED) light, and unit optics vs a fibre optic cable, in a fully crossed design (Heron Reef and Eastern Banks).
2. **Site**: the same benthic classes measured at Heron Reef and Coral Bay with unit optics under ambient light, at four thematic levels (coarse to very fine).

Differences are quantified with the spectral angle (SAM) between mean spectra.

## Running the analysis

Open `final_analysis_github_upload.Rproj` (or set the working directory to this folder) and run:

```r
source("run_all.R")
```

This rebuilds everything from the raw files in about two minutes. Scripts can also be run one at a time, in order.

Required packages: `tidyverse`, `readxl`, `signal`, `zoo`, `vegan`, `permute`. The versions used are recorded in `outputs/session_info.txt`.

## Folder structure

```
data/raw/spectra/       spectrometer exports, one file per field day (YYYYMMDD_spectra.csv)
data/raw/metadata/      field sheets, one per field day (YYYYMMDD_fielddata.xlsx)
data/processed/         intermediate files written by the scripts (not tracked)
R/functions.R           functions shared by all scripts
R/01_build_dataset.R    raw files to one long table with metadata and class names
R/02_equipment_sam.R    SAM between configurations; Figure 5, Supplementary Table 3
R/03_equipment_spectra_figures.R   Figures 3 and 4
R/04_equipment_permanova.R         Supplementary Table 5
R/05_site_comparison.R  Figure 6, Supplementary Tables 2 and 4
R/06_summary_tables.R   Table 1, Supplementary Table 1, values quoted in the Results
outputs/figures/        figures as they appear in the paper
outputs/tables/         tables as CSV
```

## Data

Each spectrometer export has one column per spectrum (`SpectrumNNNNN.asd`). The field sheet for the same day gives the first and last spectrum number for each target, its location, depth, equipment configuration and classification at five levels. Spectra are renamed `<target id>_<replicate>`. In the equipment experiment the target ID carries a letter for the configuration (for example `20241105_01a` to `20241105_01d`); `target` in the processed data drops this letter.

Reflectance was computed on the instrument from white reference and dark current readings taken immediately before each target.

## Processing

- Wavelengths restricted to 400 to 680 nm.
- Each spectrum smoothed with a Savitzky-Golay filter (3rd order polynomial, 11-point window). No normalisation is applied.
- Equipment experiment only: reflectance values outside 0 to 3 are treated as missing and linearly interpolated before smoothing. Wavelengths still missing in any spectrum (at the ends of the range) are dropped before computing SAM, so equipment SAM uses 405 to 657 nm.
- Measurements on the equipment test days that were not part of the experiment are excluded (`equipment_exclusions` in `R/functions.R`). These include three of the four configurations for target 20241105_04, so this target contributes only to the unit optics, artificial light mean and is left out of the PERMANOVA, which requires all four configurations.
- Site comparison: unit optics and ambient light only; the Clean Rock class is excluded. Classes are compared at a thematic level only where they were measured at both sites.

## Notes on the outputs

- `outputs/tables/equipment_sam_values.csv` holds every SAM value plotted in Figure 5.
- `outputs/tables/results_values.csv` lists the summary values quoted in the Results.
- Within-class SAM is the mean pairwise angle between individual spectra in a class. Where a class is represented by a single target it reflects the repeatability of replicate measurements rather than variation between individuals; target counts are given in Supplementary Tables 1 and 2.
- PERMANOVA uses SAM distances between target-mean spectra, with permutations restricted within targets. R² values are reported descriptively.
