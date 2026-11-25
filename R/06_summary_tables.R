# ==============================================================================
# 6. SUMMARY TABLES AND VALUES REPORTED IN THE TEXT
# ==============================================================================
# Table 1 (within-class vs between-factor SAM by comparison and thematic level),
# Supplementary Table 1 (equipment experiment targets) and the values quoted in
# the Results.

source("R/functions.R")

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

library_long <- readRDS("data/processed/library_long.rds")
eq <- readRDS("data/processed/equipment.rds")
within <- readRDS("data/processed/equipment_within.rds")
sites <- readRDS("data/processed/site_sam.rds")
sam <- equipment_sam(eq)

n_pct <- function(x) sprintf("%d (%.0f%%)", sum(x), 100 * mean(x))

# ------------------------------------------------------------------------------
# Table 1
# ------------------------------------------------------------------------------
# Equipment: within-class SAM for a comparison is the mean of the within-class SAM
# of the two configurations compared. Sites: within-class SAM pooled over both sites.

equipment_rows <- sam %>%
  left_join(within %>% dplyr::select(level_5, config1 = equipment_config, w1 = within_sam),
            by = c("level_5", "config1")) %>%
  left_join(within %>% dplyr::select(level_5, config2 = equipment_config, w2 = within_sam),
            by = c("level_5", "config2")) %>%
  mutate(within_sam = (w1 + w2) / 2) %>%
  group_by(comparison_id, comparison) %>%
  summarise(n_categories = n(),
            median_within_class_sam = median(within_sam),
            median_between_factor_sam = median(sam),
            n_exceeding_within_class = n_pct(sam > within_sam),
            .groups = "drop") %>%
  transmute(panel = "Equipment configuration",
            comparison = paste0("(", comparison_id, ") ", comparison),
            n_categories, median_within_class_sam, median_between_factor_sam,
            n_exceeding_within_class)

site_rows <- sites %>%
  group_by(level) %>%
  summarise(n_categories = n(),
            median_within_class_sam = median(within_sam),
            median_between_factor_sam = median(between_site_sam),
            n_exceeding_within_class = n_pct(between_site_sam > within_sam),
            .groups = "drop") %>%
  transmute(panel = "Site (Heron Reef vs Coral Bay)", comparison = as.character(level),
            n_categories, median_within_class_sam, median_between_factor_sam,
            n_exceeding_within_class)

table_1 <- bind_rows(equipment_rows, site_rows) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

write.csv(table_1, "outputs/tables/table_1.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Supplementary Table 1: targets in the equipment experiment
# ------------------------------------------------------------------------------
# Measurements collected per configuration, and how many of those were excluded
# from the analysis (see equipment_exclusions in functions.R).

config_names <- c(Fibre_Artificial = "artificial_light_cable",
                  Fibre_Ambient = "ambient_light_cable",
                  Bare_Ambient = "ambient_light_unit_optics",
                  Bare_Artificial = "artificial_light_unit_optics")

collected <- library_long %>%
  dplyr::filter(date %in% equipment_dates,
                target %in% unique(eq$wide$target)) %>%
  distinct(spectra_id, target, equipment_config)

analysed <- eq$wide %>% distinct(spectra_id, target)

target_info <- library_long %>%
  dplyr::filter(target %in% unique(eq$wide$target)) %>%
  distinct(target, location, site, depth, level_5, level_4, level_3, level_2) %>%
  group_by(target) %>%
  slice(1) %>%
  ungroup()

supp_table_1 <- collected %>%
  count(target, equipment_config) %>%
  mutate(equipment_config = config_names[equipment_config]) %>%
  pivot_wider(names_from = equipment_config, values_from = n, values_fill = 0) %>%
  left_join(collected %>%
              mutate(analysed = spectra_id %in% analysed$spectra_id) %>%
              group_by(target) %>%
              summarise(measurements_excluded = sum(!analysed), .groups = "drop"),
            by = "target") %>%
  left_join(target_info, by = "target") %>%
  mutate(location = recode(location, "Heron Island" = "Heron Reef")) %>%
  dplyr::select(target, location, sampling_area = site, depth_m = depth,
                very_fine = level_5, fine = level_4, intermediate = level_3, coarse = level_2,
                any_of(unname(config_names)), measurements_excluded) %>%
  arrange(target)

write.csv(supp_table_1, "outputs/tables/supplementary_table_1.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Values reported in the Results
# ------------------------------------------------------------------------------

value <- function(section, quantity, x, detail = "") {
  tibble(section, quantity, value = as.character(x), detail)
}
r2 <- function(x) sprintf("%.2f", x)

eq_means <- sam %>% group_by(comparison_id, comparison) %>% summarise(m = mean(sam), .groups = "drop")
over_20 <- sam %>% dplyr::filter(sam > 20) %>% arrange(comparison_id, desc(sam))
supp3 <- read.csv("outputs/tables/supplementary_table_3.csv")

site_range <- function(l) {
  x <- sites %>% dplyr::filter(level == l)
  c(min = r2(min(x$between_site_sam)), min_cat = x$category[which.min(x$between_site_sam)],
    max = r2(max(x$between_site_sam)), max_cat = x$category[which.max(x$between_site_sam)])
}
very_fine <- sites %>% dplyr::filter(level == "Very fine") %>% arrange(desc(between_site_sam))
vf_over_20 <- very_fine %>% dplyr::filter(between_site_sam > 20)

results_values <- bind_rows(
  map_df(seq_len(nrow(eq_means)), function(i)
    value("Equipment", paste0("Mean SAM (", eq_means$comparison_id[i], ") ", eq_means$comparison[i]),
          r2(eq_means$m[i]))),
  value("Equipment", "Minimum SAM", r2(min(sam$sam)),
        with(sam[which.min(sam$sam), ], paste0(level_5, ", comparison ", comparison_id))),
  value("Equipment", "Maximum SAM", r2(max(sam$sam)),
        with(sam[which.max(sam$sam), ], paste0(level_5, ", comparison ", comparison_id))),
  value("Equipment", "Median SAM", r2(median(sam$sam))),
  value("Equipment", "Categories with SAM > 20 degrees", n_distinct(over_20$level_5),
        paste0(over_20$level_5, " (", over_20$comparison_id, ", ", r2(over_20$sam), ")", collapse = "; ")),
  value("Equipment", "Categories where between-configuration SAM exceeds within-class SAM",
        paste(sum(supp3$between_exceeds_within), "of", nrow(supp3)),
        paste("Exceptions:", paste(supp3$benthic_category[!supp3$between_exceeds_within], collapse = "; "))),
  map_df(levels(sites$level), function(l) {
    r <- site_range(l)
    value("Sites", paste(l, "level between-site SAM range"), paste(r["min"], "to", r["max"]),
          paste(r["min_cat"], "to", r["max_cat"]))
  }),
  value("Sites", "Very fine categories with SAM > 20 degrees",
        sprintf("%d of %d (%.0f%%)", nrow(vf_over_20), nrow(very_fine), 100 * nrow(vf_over_20) / nrow(very_fine)),
        paste0(vf_over_20$category, " (", r2(vf_over_20$between_site_sam), ")", collapse = "; ")),
  map_df(levels(sites$level), function(l) {
    x <- sites %>% dplyr::filter(level == l)
    value("Sites", paste(l, "level categories exceeding within-class SAM"),
          n_pct(x$between_site_sam > x$within_sam))
  })
)

write.csv(results_values, "outputs/tables/results_values.csv", row.names = FALSE)
