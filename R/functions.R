# Shared functions for the spectral library analysis.
# Sourced by every script in R/.

library(tidyverse)
library(readxl)
library(signal)
library(zoo)

wl_min <- 400
wl_max <- 680

# ------------------------------------------------------------------------------
# Spectral angle
# ------------------------------------------------------------------------------

# Spectral angle (degrees) between two spectra
spectral_angle <- function(spectrum1, spectrum2) {
  if (length(spectrum1) != length(spectrum2)) stop("Spectra must have the same length")
  if (anyNA(spectrum1) || anyNA(spectrum2)) return(NA_real_)
  cos_angle <- sum(spectrum1 * spectrum2) /
    (sqrt(sum(spectrum1^2)) * sqrt(sum(spectrum2^2)))
  acos(pmax(-1, pmin(1, cos_angle))) * 180 / pi
}

# Matrix of spectral angles (degrees) between all rows of X
sam_matrix <- function(X) {
  X_unit <- X / sqrt(rowSums(X^2))
  cos_angle <- X_unit %*% t(X_unit)
  cos_angle[cos_angle > 1] <- 1
  cos_angle[cos_angle < -1] <- -1
  acos(cos_angle) * 180 / pi
}

# All pairwise angles among the rows of X
pairwise_angles <- function(X) {
  X <- X[complete.cases(X), , drop = FALSE]
  if (nrow(X) < 2) return(numeric(0))
  A <- sam_matrix(X)
  A[upper.tri(A)]
}

# ------------------------------------------------------------------------------
# Raw data to a single long table
# ------------------------------------------------------------------------------

# Reads the spectrometer exports and field sheets, renames each spectrum to
# <target id>_<replicate>, joins the field metadata and standardises class names.
read_spectral_library <- function(raw_dir = "data/raw") {
  metadata_files <- list.files(file.path(raw_dir, "metadata"), pattern = "_fielddata\\.xlsx$",
                               full.names = TRUE)
  spectra_files <- list.files(file.path(raw_dir, "spectra"), pattern = "_spectra\\.csv$",
                              full.names = TRUE)
  metadata <- setNames(lapply(metadata_files, read_excel), basename(metadata_files))
  spectra <- setNames(lapply(spectra_files, read.csv), basename(spectra_files))

  # Rename SpectrumNNNNN.asd columns using the start/end numbers in the field sheet
  renamed <- list()
  for (f in names(spectra)) {
    date_key <- str_extract(f, "\\d{8}")
    md <- metadata[[paste0(date_key, "_fielddata.xlsx")]]
    if (is.null(md)) {
      warning("No field sheet for ", f)
      next
    }
    sp <- spectra[[f]]
    new_names <- names(sp)
    for (i in seq_len(nrow(md))) {
      for (rep in 1:(md$end_number[i] - md$start_number[i] + 1)) {
        col <- which(names(sp) == sprintf("Spectrum%05d.asd", md$start_number[i] + rep - 1))
        if (length(col) > 0) new_names[col] <- sprintf("%s_%02d", md$id_no[i], rep)
      }
    }
    names(sp) <- new_names
    renamed[[date_key]] <- sp
  }

  wavelength <- renamed[[1]]$wavelength_nm
  stopifnot(all(sapply(renamed, function(d) identical(d$wavelength_nm, wavelength))))
  spectra_wide <- cbind(
    data.frame(wavelength_nm = wavelength),
    do.call(cbind, unname(lapply(renamed, function(d) d[, setdiff(names(d), c("X", "wavelength_nm")), drop = FALSE])))
  ) %>%
    dplyr::select(-starts_with("Spectrum")) %>%   # spectra not listed in a field sheet
    dplyr::filter(wavelength_nm >= 300, wavelength_nm <= 700)

  all_metadata <- as.data.frame(do.call(rbind, metadata))

  spectra_wide %>%
    pivot_longer(-wavelength_nm, names_to = "spectra_id", values_to = "reflectance") %>%
    mutate(base_id = str_remove(spectra_id, "_\\d+$"),
           replicate_no = str_extract(spectra_id, "\\d+$")) %>%
    left_join(all_metadata, by = c("base_id" = "id_no"), relationship = "many-to-many") %>%
    mutate(date = as.Date(date)) %>%
    standardise_library()
}

# Filters, equipment configuration labels and class names
standardise_library <- function(d) {
  d %>%
    mutate(level_2 = case_when(
      level_2 %in% c("Soft substrate", "Soft Subtrate") ~ "Soft Substrate",
      TRUE ~ level_2
    )) %>%
    dplyr::filter(lighting != "shaded",
                  level_4 != "Macroalgae - Red",
                  level_5 != "Macroalgae - Brown",
                  wavelength_nm >= 400,
                  wavelength_nm <= 700) %>%
    mutate(
      lighting_clean = case_when(
        str_detect(tolower(lighting), "divelantern") ~ "artificial",
        tolower(trimws(lighting)) == "ambient" ~ "ambient",
        TRUE ~ tolower(trimws(lighting))
      ),
      fore_optic_clean = str_replace(tolower(`fore-optic`), "fibre[\\s-]*optic[\\s-]*cable",
                                     "fibre-optic cable"),
      equipment_config = paste(
        ifelse(fore_optic_clean %in% c("bare fibre", "bare"), "Bare", "Fibre"),
        case_when(lighting_clean == "artificial" ~ "Artificial",
                  lighting_clean == "ambient" ~ "Ambient",
                  TRUE ~ paste0("ERROR_", lighting_clean)),
        sep = "_"
      ),
      # target ID without the configuration letter (e.g. 20241105_01a -> 20241105_01)
      target = str_remove(base_id, "[a-z]$")
    ) %>%
    mutate(level_5 = case_when(
      level_5 == "Croymbose Hard Coral - Blue" ~ "Corymbose Hard Coral - Blue",
      level_5 == "Fluroescing Hard Coral" ~ "Fluorescing Hard Coral",
      level_5 == "Massive Hard Coral  - Red" ~ "Massive Hard Coral - Red",
      level_5 == "Tabular /Plating Hard Coral - Red" ~ "Tabular/Plating Hard Coral - Red",
      level_5 %in% c("Sea Cucumber", "Invertebrate - Sea cucumber",
                     "Invertebrate - sea cucumber") ~ "Invertebrate - Sea Cucumber",
      level_5 %in% c("Anemone", "Invertebrate - anemone") ~ "Invertebrate - Anemone",
      level_5 == "Clean seagrass - halophila spinulosa" ~ "Clean seagrass - Halophila spinulosa",
      level_5 == "Clean seagrass - halophila ovalis" ~ "Clean seagrass - Halophila ovalis",
      level_5 == "Clean seagrass - halophila uninervis" ~ "Clean seagrass - Halophila uninervis",
      level_5 == "Clean seagrass - oceana serrulata" ~ "Clean seagrass - Oceana serrulata",
      level_5 == "Clean seagrass - zostera muelleri" ~ "Clean seagrass - Zostera muelleri",
      level_5 == "Fouled seagrass - oceana serrulata" ~ "Fouled seagrass - Oceana serrulata",
      level_5 == "Fouled seagrass - zostera muelleri" ~ "Fouled seagrass - Zostera muelleri",
      level_5 == "lyngbya majuscula" ~ "Lyngbya majuscula",
      level_5 %in% c("Branching Hard Coral - Brown", "Branching Hard Coral - Yellow",
                     "Digitate Hard Coral - Brown/Yellow") ~ "Branching Hard Coral - Brown/Yellow",
      level_5 == "Tabular/Plating Hard Coral - Brown" ~ "Tabular/Plating Hard Coral - Brown/Yellow",
      level_5 == "Encrusting Hard Coral - Brown" ~ "Encrusting Hard Coral - Brown/Yellow",
      level_5 == "Foliose Hard Coral - Brown" ~ "Foliose Hard Coral - Brown/Yellow",
      level_5 == "Corymbose Hard Coral - Brown" ~ "Corymbose Hard Coral - Brown/Yellow",
      TRUE ~ level_5
    )) %>%
    mutate(level_5 = case_when(
      level_5 == "Octocoral - Brown" ~ "Brown Octocoral",
      level_5 == "Octocoral - Green" ~ "Green Octocoral",
      level_5 == "Octocoral - Red" ~ "Red Octocoral",
      level_5 == "Octocoral - Yellow" ~ "Yellow Octocoral",
      level_5 == "Macroalgae - Padina" ~ "Padina",
      level_5 == "Macroalgae - Sargassum" ~ "Sargassum",
      level_5 == "Macroalgae - Turbinaria" ~ "Turbinaria",
      level_5 == "Macroalgae - Chlorodesmis" ~ "Chlorodesmis",
      level_5 == "Macroalgae - Halimeda" ~ "Halimeda",
      level_5 == "Macroalgae - Red" ~ "Red Macroalgae",
      level_5 == "Macroalgae - Lobophora" ~ "Lobophora",
      level_5 == "Macroalgae - Caulerpa" ~ "Caulerpa",
      level_5 == "Macroalgae - Blue/Green" ~ "Blue/Green Macroalgae",
      level_5 == "Macroalgae - Green" ~ "Green Macroalgae",
      level_5 == "Macroalgae - Mixed" ~ "Mixed Macroalgae",
      level_5 %in% c("Clean dense seagrass - Oceana serrulata",
                     "Clean seagrass - Oceana serrulata") ~ "Clean Dense Strap-Like Seagrass - Oceana serrulata",
      level_5 %in% c("Fouled dense seagrass - Oceana serrulata",
                     "Fouled seagrass - Oceana serrulata") ~ "Fouled Dense Strap-Like Seagrass - Oceana serrulata",
      level_5 == "Clean dense seagrass - Zostera muelleri" ~ "Clean Dense Strap-Like Seagrass - Zostera muelleri",
      level_5 == "Fouled dense seagrass - Zostera muelleri" ~ "Fouled Dense Strap-Like Seagrass - Zostera muelleri",
      level_5 == "Clean sparse seagrass - Halophila uninervis" ~ "Clean Sparse Leaf-Like Seagrass - Halophila uninervis",
      level_5 %in% c("Clean sparse seagrass - Halophila ovalis",
                     "Clean seagrass - Halophila ovalis") ~ "Clean Sparse Leaf-Like Seagrass - Halophila ovalis",
      level_5 == "Clean sparse seagrass - Halophila spinulosa" ~ "Clean Sparse Leaf-Like Seagrass - Halophila spinulosa",
      level_5 == "Clean dense seagrass - Halophila spinulosa" ~ "Clean Dense Strap-Like Seagrass - Halophila spinulosa",
      level_5 %in% c("Non Carbonate Sand", "Non Carbonate Sand - Silt Detritus") ~ "Non-Carbonate Sand",
      level_5 == "Crustose coralline algae" ~ "Crustose Coralline Algae",
      level_5 == "Tabular Hard Coral - White" ~ "Tabular/Plating Hard Coral - White",
      level_5 == "Tabular Hard Coral - Blue" ~ "Tabular/Plating Hard Coral - Blue",
      level_5 == "Branching Hard Coral - Red" ~ "Branching Hard Coral - Red/Pink",
      level_5 == "Corymbose Hard Coral - Red" ~ "Corymbose Hard Coral - Red/Pink",
      level_5 == "Encrusting Hard Coral - Red" ~ "Encrusting Hard Coral - Red/Pink",
      level_5 == "Tabular/Plating Hard Coral - Red" ~ "Tabular/Plating Hard Coral - Red/Pink",
      level_5 == "Massive Hard Coral - Red" ~ "Massive Hard Coral - Red/Pink",
      TRUE ~ level_5
    )) %>%
    mutate(level_4 = case_when(
      level_4 == "Hard Coral - Brown/Yellow" ~ "Brown/Yellow Hard Coral",
      level_4 == "Hard Coral - Blue" ~ "Blue Hard Coral",
      level_4 == "Hard Coral - White" ~ "White Hard Coral",
      level_4 == "Hard Coral - Green" ~ "Green Hard Coral",
      level_4 == "Hard Coral - Red" ~ "Red/Pink Hard Coral",
      level_4 == "Hard Coral - Fluorescing" ~ "Fluorescing Hard Coral",
      level_4 == "Octocoral - Brown" ~ "Brown Octocoral",
      level_4 == "Octocoral - Green" ~ "Green Octocoral",
      level_4 %in% c("Octocoal - Red", "Octocoral - Red") ~ "Red Octocoral",
      level_4 == "Octocoral - Yellow" ~ "Yellow Octocoral",
      level_4 == "Macroalgae - Brown" ~ "Brown Macroalgae",
      level_4 == "Macroalgae - Green" ~ "Green Macroalgae",
      level_4 == "Macroalgae - Red" ~ "Red Macroalgae",
      level_4 == "Macroalgae - Blue/Green" ~ "Blue/Green Macroalgae",
      level_4 == "Macroalgae - Mixed" ~ "Mixed Macroalgae",
      level_4 %in% c("Clean Dense Strap Like Seagrass", "Clean Dense Leaf Like Seagrass") ~ "Clean Dense Seagrass",
      level_4 == "Fouled Dense Strap Like Seagrass" ~ "Fouled Dense Seagrass",
      level_4 %in% c("Clean Sparse Leaf Like Seagrass", "Clean Sparse Strap Like Seagrass") ~ "Clean Sparse Seagrass",
      level_4 == "Non Carbonate Sand" ~ "Non-Carbonate Sand",
      level_4 == "Crustose coralline algae" ~ "Crustose Coralline Algae",
      level_4 == "lyngbya majuscula" ~ "Lyngbya majuscula",
      TRUE ~ level_4
    ))
}

# ------------------------------------------------------------------------------
# Equipment comparison
# ------------------------------------------------------------------------------

equipment_dates <- as.Date(c("2024-11-05", "2024-11-06", "2024-12-05"))

# Measurements taken on the equipment-test days that were not part of the
# equipment experiment
equipment_exclusions <- c(sprintf("20241205_%02d", 8:16),
                          "20241105_04a", "20241105_04b", "20241105_04c")

equipment_comparisons <- tribble(
  ~comparison_id, ~comparison,                                  ~config1,          ~config2,
  "a",            "Artificial vs ambient light, unit optics",   "Bare_Ambient",    "Bare_Artificial",
  "b",            "Artificial vs ambient light, cable optics",  "Fibre_Ambient",   "Fibre_Artificial",
  "c",            "Cable vs unit optics, ambient light",        "Bare_Ambient",    "Fibre_Ambient",
  "d",            "Cable vs unit optics, artificial light",     "Bare_Artificial", "Fibre_Artificial"
)

# Equipment spectra: values outside 0-3 are treated as missing and linearly
# interpolated, then each spectrum is Savitzky-Golay smoothed (3rd order, 11 points).
# Wavelengths with a missing value in any spectrum are dropped from the wide table.
prepare_equipment <- function(library_long) {
  long <- library_long %>%
    dplyr::filter(date %in% equipment_dates,
                  wavelength_nm >= wl_min, wavelength_nm <= wl_max,
                  !base_id %in% equipment_exclusions,
                  level_5 != "Echinoderm - Sea Star") %>%
    arrange(spectra_id, wavelength_nm) %>%
    group_by(spectra_id) %>%
    mutate(reflectance = ifelse(reflectance < 0 | reflectance > 3, NA, reflectance),
           reflectance = na.approx(reflectance, x = wavelength_nm, na.rm = FALSE),
           reflectance_smooth = sgolayfilt(reflectance, p = 3, n = 11)) %>%
    ungroup()

  wide <- long %>%
    dplyr::select(spectra_id, base_id, target, level_5, level_4, level_3,
                  equipment_config, wavelength_nm, reflectance_smooth) %>%
    pivot_wider(names_from = wavelength_nm, values_from = reflectance_smooth,
                names_prefix = "wl_") %>%
    dplyr::select(where(~ !anyNA(.)))

  list(long = long, wide = wide, wl_cols = grep("^wl_", names(wide), value = TRUE))
}

# Mean spectrum of each category under each configuration (all replicate scans)
category_means <- function(eq) {
  eq$wide %>%
    group_by(level_5, equipment_config) %>%
    summarise(across(all_of(eq$wl_cols), mean), .groups = "drop")
}

# SAM between configurations for each category and comparison
equipment_sam <- function(eq) {
  means <- category_means(eq)
  spectrum <- function(category, config) {
    as.numeric(means[means$level_5 == category & means$equipment_config == config, eq$wl_cols])
  }
  crossing(level_5 = unique(means$level_5), equipment_comparisons) %>%
    rowwise() %>%
    mutate(sam = spectral_angle(spectrum(level_5, config1), spectrum(level_5, config2))) %>%
    ungroup()
}

# Within-class SAM: mean pairwise angle among scans of a category under one configuration
equipment_within_sam <- function(eq) {
  eq$wide %>%
    group_by(level_5, equipment_config) %>%
    summarise(within_sam = mean(pairwise_angles(as.matrix(pick(all_of(eq$wl_cols))))),
              n_scans = n(), .groups = "drop")
}

# ------------------------------------------------------------------------------
# Site comparison
# ------------------------------------------------------------------------------

thematic_levels <- c(level_2 = "Coarse", level_3 = "Intermediate",
                     level_4 = "Fine", level_5 = "Very fine")

# Unit optics, ambient light spectra from Heron Reef and Coral Bay, smoothed as above.
prepare_sites <- function(library_long) {
  long <- library_long %>%
    dplyr::filter(equipment_config == "Bare_Ambient",
                  !level_1 %in% c("Biotic", "Animalia"),
                  location %in% c("Heron Island", "Coral Bay"),
                  level_3 != "Clean Rock",
                  wavelength_nm >= wl_min, wavelength_nm <= wl_max) %>%
    arrange(spectra_id, wavelength_nm) %>%
    group_by(spectra_id) %>%
    mutate(reflectance_smooth = sgolayfilt(reflectance, p = 3, n = 11)) %>%
    ungroup()

  wide <- long %>%
    dplyr::select(spectra_id, base_id, location, level_5, level_4, level_3, level_2,
                  wavelength_nm, reflectance_smooth) %>%
    pivot_wider(names_from = wavelength_nm, values_from = reflectance_smooth,
                names_prefix = "wl_")

  list(long = long, wide = wide, wl_cols = grep("^wl_", names(wide), value = TRUE))
}

# For each thematic level and each category found at both sites: number of targets
# and measurements per site, between-site SAM (angle between site mean spectra) and
# within-class SAM (mean pairwise angle within each site, and pooled over both sites).
site_sam <- function(st) {
  map_df(names(thematic_levels), function(lv) {
    w <- st$wide %>% mutate(category = .data[[lv]])
    shared <- w %>%
      group_by(category) %>%
      dplyr::filter(n_distinct(location) == 2) %>%
      pull(category) %>%
      unique()
    map_df(shared, function(cat) {
      wc <- w %>% dplyr::filter(category == cat)
      X_hr <- as.matrix(wc[wc$location == "Heron Island", st$wl_cols])
      X_cb <- as.matrix(wc[wc$location == "Coral Bay", st$wl_cols])
      a_hr <- pairwise_angles(X_hr)
      a_cb <- pairwise_angles(X_cb)
      tibble(
        level = thematic_levels[[lv]],
        category = cat,
        targets_heron = n_distinct(wc$base_id[wc$location == "Heron Island"]),
        measurements_heron = nrow(X_hr),
        targets_coral_bay = n_distinct(wc$base_id[wc$location == "Coral Bay"]),
        measurements_coral_bay = nrow(X_cb),
        within_sam_heron = if (length(a_hr)) mean(a_hr) else NA_real_,
        within_sam_coral_bay = if (length(a_cb)) mean(a_cb) else NA_real_,
        within_sam = mean(c(a_hr, a_cb)),
        between_site_sam = spectral_angle(colMeans(X_hr), colMeans(X_cb))
      )
    })
  }) %>%
    mutate(level = factor(level, levels = thematic_levels))
}
