# ==============================================================================
# 5. SITE COMPARISON: HERON REEF AND CORAL BAY
# ==============================================================================
# SAM between site mean spectra at four thematic levels (Figure 6), within-class
# SAM (Supplementary Table 4) and replication at each site (Supplementary Table 2).

source("R/functions.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

library_long <- readRDS("data/processed/library_long.rds")

st <- prepare_sites(library_long)
sites <- site_sam(st)
saveRDS(sites, "data/processed/site_sam.rds")

# ------------------------------------------------------------------------------
# Supplementary Table 2: very fine classes, their classification at each level,
# and the number of targets and measurements at each site
# ------------------------------------------------------------------------------

classes <- st$wide %>% distinct(level_5, level_4, level_3, level_2)

supp_table_2 <- sites %>%
  dplyr::filter(level == "Very fine") %>%
  left_join(classes, by = c("category" = "level_5")) %>%
  transmute(very_fine = category, fine = level_4, intermediate = level_3, coarse = level_2,
            coral_bay_targets = targets_coral_bay, coral_bay_measurements = measurements_coral_bay,
            heron_reef_targets = targets_heron, heron_reef_measurements = measurements_heron) %>%
  arrange(very_fine)

write.csv(supp_table_2, "outputs/tables/supplementary_table_2.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Supplementary Table 4
# ------------------------------------------------------------------------------

supp_table_4 <- sites %>%
  transmute(category, level,
            between_site_sam = between_site_sam,
            mean_within_class_sam = within_sam,
            ratio_between_to_within = between_site_sam / within_sam,
            between_exceeds_within = between_site_sam > within_sam) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>%
  arrange(level, category)

write.csv(supp_table_4, "outputs/tables/supplementary_table_4.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Figure 6
# ------------------------------------------------------------------------------

figure_6 <- ggplot(sites, aes(x = between_site_sam, y = reorder(category, between_site_sam),
                              colour = between_site_sam)) +
  geom_segment(aes(x = 0, xend = between_site_sam, y = category, yend = category),
               colour = "grey70", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "#2ECC71", mid = "#F39C12", high = "#E74C3C",
                         midpoint = 30, breaks = c(0, 10, 20, 30, 40)) +
  facet_wrap(~ level, scales = "free_y", ncol = 1) +
  labs(x = "SAM (degrees)", y = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

ggsave("outputs/figures/figure_6.png", figure_6, width = 8, height = 12, dpi = 300, bg = "white")
