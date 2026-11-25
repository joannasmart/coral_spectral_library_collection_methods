# ==============================================================================
# 4. EQUIPMENT CONFIGURATION: VARIANCE PARTITIONING
# ==============================================================================
# PERMANOVA on SAM distances between target-mean spectra (Supplementary Table 5)
# and a test of multivariate dispersion among configurations.
#
# Replicate scans are averaged to one spectrum per target and configuration.
# Only targets measured under all four configurations are used. Permutations are
# restricted within targets, so the category term cannot be tested; R2 values are
# reported descriptively.

source("R/functions.R")
library(vegan)
library(permute)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
set.seed(20261013)
n_perm <- 9999

eq <- readRDS("data/processed/equipment.rds")

target_means <- eq$wide %>%
  group_by(target, category = level_5, equipment_config) %>%
  summarise(across(all_of(eq$wl_cols), mean), .groups = "drop") %>%
  group_by(target) %>%
  dplyr::filter(n_distinct(equipment_config) == 4) %>%
  ungroup() %>%
  mutate(light = word(equipment_config, 2, sep = "_"),
         optics = ifelse(word(equipment_config, 1, sep = "_") == "Bare", "Unit", "Cable"))

design <- target_means %>%
  dplyr::select(target, category, light, optics, equipment_config) %>%
  mutate(across(everything(), factor))

cat("Targets:", n_distinct(design$target), " Observations:", nrow(design), "\n")
cat("Targets per category:\n")
print(table(table(distinct(design, target, category)$category)))

D <- as.dist(sam_matrix(as.matrix(target_means[, eq$wl_cols])))

# ------------------------------------------------------------------------------
# PERMANOVA (Supplementary Table 5)
# ------------------------------------------------------------------------------

permanova <- adonis2(
  D ~ category + light + optics + light:optics + category:light + category:optics,
  data = design,
  permutations = how(nperm = n_perm, blocks = design$target),
  by = "terms"
)
print(permanova)

term_labels <- c(
  category = "Benthic category",
  light = "Light",
  optics = "Optics",
  `light:optics` = "Light x optics",
  `category:light` = "Benthic category x light",
  `category:optics` = "Benthic category x optics",
  Residual = "Residual",
  Total = "Total"
)

supp_table_5 <- tibble(
  term = term_labels[rownames(permanova)],
  df = permanova$Df,
  sum_of_squares = round(permanova$SumOfSqs, 1),
  r_squared = round(permanova$R2, 3)
)
write.csv(supp_table_5, "outputs/tables/supplementary_table_5.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Multivariate dispersion by configuration
# ------------------------------------------------------------------------------

# Distances are measured to each configuration's centroid
dispersion <- betadisper(D, design$equipment_config, type = "centroid")
dispersion_test <- permutest(dispersion, permutations = how(nperm = n_perm, blocks = design$target))
print(dispersion_test)

dispersion_table <- tibble(
  configuration = names(tapply(dispersion$distances, design$equipment_config, mean)),
  mean_distance_to_centroid = round(as.numeric(tapply(dispersion$distances, design$equipment_config, mean)), 2),
  F = round(dispersion_test$tab["Groups", "F"], 2),
  df = paste(dispersion_test$tab["Groups", "Df"], dispersion_test$tab["Residuals", "Df"], sep = ", "),
  p = dispersion_test$tab["Groups", "Pr(>F)"]
)
write.csv(dispersion_table, "outputs/tables/supplementary_table_5_dispersion.csv", row.names = FALSE)
