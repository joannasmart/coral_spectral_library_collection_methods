# ==============================================================================
# 2. EQUIPMENT CONFIGURATION: SPECTRAL ANGLES
# ==============================================================================
# SAM between equipment configurations for each benthic category (Figure 5) and
# comparison with within-class SAM (Supplementary Table 3).

source("R/functions.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

library_long <- readRDS("data/processed/library_long.rds")

eq <- prepare_equipment(library_long)
saveRDS(eq, "data/processed/equipment.rds")

# ------------------------------------------------------------------------------
# Between-configuration SAM
# ------------------------------------------------------------------------------

sam <- equipment_sam(eq)

write.csv(sam %>% dplyr::select(level_5, comparison_id, comparison, sam),
          "outputs/tables/equipment_sam_values.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Within-class vs between-configuration SAM (Supplementary Table 3)
# ------------------------------------------------------------------------------

within <- equipment_within_sam(eq)
saveRDS(within, "data/processed/equipment_within.rds")

supp_table_3 <- within %>%
  group_by(level_5) %>%
  summarise(mean_within_class_sam = mean(within_sam), .groups = "drop") %>%
  left_join(sam %>%
              group_by(level_5) %>%
              summarise(mean_between_configuration_sam = mean(sam), .groups = "drop"),
            by = "level_5") %>%
  mutate(ratio_between_to_within = mean_between_configuration_sam / mean_within_class_sam,
         between_exceeds_within = mean_between_configuration_sam > mean_within_class_sam) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>%
  rename(benthic_category = level_5)

write.csv(supp_table_3, "outputs/tables/supplementary_table_3.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Figure 5
# ------------------------------------------------------------------------------

panel_titles <- c(
  a = "(a) Impact of artificial light\n(unit optics)",
  b = "(b) Impact of artificial light\n(cable optics)",
  c = "(c) Impact of cable optics\n(ambient light)",
  d = "(d) Impact of cable optics\n(artificial light)"
)

figure_5 <- sam %>%
  mutate(panel = factor(panel_titles[comparison_id], levels = panel_titles)) %>%
  ggplot(aes(x = sam, y = reorder(level_5, sam), colour = sam)) +
  geom_segment(aes(x = 0, xend = sam, y = level_5, yend = level_5), colour = "grey70") +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "#2ECC71", mid = "#F39C12", high = "#E74C3C", midpoint = 25) +
  facet_wrap(~ panel, ncol = 2) +
  labs(x = "SAM (degrees)", y = "Habitat Class") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA))

ggsave("outputs/figures/figure_5.png", figure_5, width = 10, height = 8, dpi = 300)
