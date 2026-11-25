# ==============================================================================
# 3. EQUIPMENT CONFIGURATION: SPECTRA AND REFLECTANCE DIFFERENCES
# ==============================================================================
# Mean spectra under each configuration (Figure 3) and reflectance differences
# between configurations at each wavelength (Figure 4).

source("R/functions.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

eq <- readRDS("data/processed/equipment.rds")

config_labels <- c(
  Bare_Ambient = "Unit optics, Ambient light",
  Bare_Artificial = "Unit optics, Artificial light",
  Fibre_Ambient = "Fiber optic cable, Ambient light",
  Fibre_Artificial = "Fiber optic cable, Artificial light"
)

category_labels <- c(
  "Branching Hard Coral - Blue" = "Branching Coral\n(Blue)",
  "Branching Hard Coral - Brown/Yellow" = "Branching Coral\n(Brown/Yellow)",
  "Carbonate Sand" = "Carbonate Sand",
  "Clean Dense Strap-Like Seagrass - Oceana serrulata" = "Clean Seagrass\n(O. serrulata)",
  "Clean Dense Strap-Like Seagrass - Zostera muelleri" = "Clean Seagrass\n(Z. muelleri)",
  "Encrusting Hard Coral - Brown/Yellow" = "Encrusting Coral\n(Brown/Yellow)",
  "Fouled Dense Strap-Like Seagrass - Oceana serrulata" = "Fouled Seagrass\n(O. serrulata)",
  "Fouled Dense Strap-Like Seagrass - Zostera muelleri" = "Fouled Seagrass\n(Z. muelleri)",
  "Mixed Macroalgae" = "Macroalgae\n(Mixed)",
  "Massive Hard Coral - Blue" = "Massive Coral\n(Blue)",
  "Massive Hard Coral - Brown/Yellow" = "Massive Coral\n(Brown/Yellow)",
  "Microphytobenthos" = "Microphytobenthos",
  "Mushroom/Free Living Hard Coral - Brown/Yellow" = "Mushroom Coral\n(Brown/Yellow)",
  "Non-Carbonate Sand" = "Non-Carbonate Sand",
  "Brown Octocoral" = "Octocoral\n(Brown)",
  "Green Octocoral" = "Octocoral\n(Green)",
  "Turf Algae" = "Turf Algae",
  "Lyngbya majuscula" = "Cyanobacteria\n(Lyngbya majuscula)"
)

# Mean and SD of smoothed reflectance across replicate scans
spectra_summary <- eq$long %>%
  group_by(wavelength_nm, level_5, equipment_config) %>%
  summarise(mean_reflectance = mean(reflectance_smooth, na.rm = TRUE),
            sd_reflectance = sd(reflectance_smooth, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(equipment_config = factor(config_labels[equipment_config], levels = config_labels),
         category = factor(category_labels[level_5], levels = category_labels))

stopifnot(!anyNA(spectra_summary$category))

# ------------------------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------------------------

figure_3 <- ggplot(spectra_summary,
                   aes(x = wavelength_nm, y = mean_reflectance,
                       group = equipment_config, colour = equipment_config)) +
  geom_ribbon(aes(ymin = mean_reflectance - sd_reflectance,
                  ymax = mean_reflectance + sd_reflectance,
                  fill = equipment_config),
              alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(limits = c(wl_min, wl_max), breaks = seq(400, 650, 50),
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_colour_viridis_d(option = "turbo", end = 0.9) +
  scale_fill_viridis_d(option = "turbo", end = 0.9) +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       colour = "Equipment\nConfiguration", fill = "Equipment\nConfiguration") +
  facet_wrap(~ category, nrow = 7, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.5),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.5),
    strip.text = element_text(size = 9, face = "bold", margin = margin(2, 2, 2, 2)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.5, "cm"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4)
  )

ggsave("outputs/figures/figure_3.png", figure_3, width = 9, height = 12, dpi = 600)

# ------------------------------------------------------------------------------
# Figure 4
# ------------------------------------------------------------------------------
# Treatment minus reference, so positive values mean artificial light (or cable
# optics) increased reflectance.

differences <- spectra_summary %>%
  dplyr::select(category, wavelength_nm, equipment_config, mean_reflectance) %>%
  pivot_wider(names_from = equipment_config, values_from = mean_reflectance) %>%
  transmute(
    category, wavelength_nm,
    `Impact of artificial light (unit optics)` =
      `Unit optics, Artificial light` - `Unit optics, Ambient light`,
    `Impact of artificial light (cable optics)` =
      `Fiber optic cable, Artificial light` - `Fiber optic cable, Ambient light`,
    `Impact of cable optics (ambient light)` =
      `Fiber optic cable, Ambient light` - `Unit optics, Ambient light`,
    `Impact of cable optics (artificial light)` =
      `Fiber optic cable, Artificial light` - `Unit optics, Artificial light`
  ) %>%
  pivot_longer(starts_with("Impact"), names_to = "comparison", values_to = "difference") %>%
  mutate(comparison = factor(comparison, levels = c(
    "Impact of artificial light (unit optics)",
    "Impact of artificial light (cable optics)",
    "Impact of cable optics (ambient light)",
    "Impact of cable optics (artificial light)"
  )))

max_abs_diff <- max(abs(differences$difference), na.rm = TRUE)

figure_4 <- ggplot(differences, aes(x = wavelength_nm, y = comparison, fill = difference)) +
  geom_tile() +
  facet_wrap(~ category, ncol = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-max_abs_diff, max_abs_diff),
                       name = "Reflectance\ndifference") +
  scale_x_continuous(breaks = seq(400, 650, 50), expand = c(0, 0)) +
  coord_cartesian(xlim = c(wl_min, wl_max)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Wavelength (nm)", y = "") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    strip.text = element_text(face = "bold", size = 12, margin = margin(2, 2, 2, 2)),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave("outputs/figures/figure_4.png", figure_4, width = 10, height = 12, dpi = 600)
