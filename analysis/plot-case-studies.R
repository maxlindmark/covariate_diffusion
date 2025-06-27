# Load libraries
library(tidyr)
library(dplyr)
library(ggtext)
library(ggstats)
library(patchwork)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())
library(viridis)
library(sf)
library(readr)
library(forcats)
library(ggrepel)
library(egg)
library(stringr)

root_dir <- here::here(".")

### Plot results from the case studies
# Read data
ebs_names <- read_csv(paste0(root_dir, "/data/clean_EBS_species.csv")) |>
  dplyr::select(X, abb_name, sc_name, common_name2)

ebs <-
  # read.csv(paste0(root_dir, "/results/2024-08-12_identity_LNP/Results_cz.csv")) |>
  #read.csv(paste0(root_dir, "/results/2024-09-24_identity_LNP/Results_cz.csv")) |>
  read.csv(paste0(root_dir, "/results/2025-06-27_identity_LNP/Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  rename(covar_corr = corr_depth) |>
  left_join(ebs_names, by = "X")

ebs |> distinct(abb_name) |> arrange()

# Fixme: this species is now NA...
#ebs |> filter(abb_name == "<i>P. camtschaticus</i> (Bb)")
ebs <- ebs |> filter(!abb_name == "<i>P. camtschaticus</i> (Bb)")

# Breeding bird case
bb <-
  #read.csv(paste0(root_dir, "/results/2024-08-07_LNP/Results_bb_cz.csv")) |>
  read.csv(paste0(root_dir, "/results/2025-06-27_LNP/Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  rename(covar_corr = corr_pop_dens) |>
  separate(X, "_", into = c("family", "species")) |>
  mutate(
    abb_name = substring(family, 1, 1),
    abb_name = paste0("<i>", paste(paste0(abb_name, "."), species), "</i>")
  )

dd <- bind_rows(
  ebs |> mutate(case = "Eastern Bering sea fishes"),
  bb |> mutate(case = "Breeding bird survey")
)

ebs |> nrow()
ebs |> filter(deltaAIC > 2) |> nrow()
ebs |> filter(deltaCAIC > 2) |> nrow()

p <- ggplot(dd, aes(deltaAIC, reorder(abb_name, desc(deltaAIC)), fill = covar_corr)) +
  geom_rect(aes(xmin = 2, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "grey95"
  ) +
  geom_vline(xintercept = c(-2, 2), alpha = 0.5, linetype = 2, linewidth = 0.35) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 1, linewidth = 0.35) +
  geom_point(shape = 21, color = "grey10", stroke = 0.01, size = 2.3) +
  labs(x = "ΔmAIC", y = "Species", fill = "Correlation between raw and diffused covariate") +
  scale_x_continuous(
    trans = "fourth_root_power",
    breaks = c(-2, 0, 2, 40, 80, 120)
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~case, scales = "free", ncol = 2) +
  theme(
    axis.text.y = element_markdown(),
    axis.text.x = element_text(size = 6.5),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.direction = "horizontal"
  ) +
  scale_fill_viridis(option = "rocket")

tag_facet(p, fontface = 1, size = 3.5, hjust = -9.5)

ggsave(paste0(here::here(), "/results/figures/case_summary.pdf"), width = 17, height = 19, unit = "cm", device = cairo_pdf)


# Plot maps from the case studies
# BB
source(file.path(root_dir, "functions/bb-map-plot.R"))

bb_stuff <-
  readRDS(paste0(here::here(), "/results/2024-08-07_LNP/stuff_gz_df.rds")) |>
  separate(species, "_", into = c("family", "sp"), remove = FALSE) |>
  mutate(
    abb_name = substring(family, 1, 1),
    abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>")
  )

sub <- bb_stuff |> filter(species %in% c("Sturnus_vulgaris", "Pheucticus_melanocephalus"))

diff_label <- bb |>
  filter(abb_name %in% unique(sub$abb_name)) |>
  mutate(label = paste0("log(\u03BA)=", round(ln_kappa2, digits = 2)))

# Main plot
p1 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = orig), color = NA) +
  labs(
    tag = "(a) log(population density)",
    fill = NULL
  ) +
  facet_wrap(~NA) +
  guides(fill = guide_colorbar(position = "inside", title.position = "top", title.hjust = 0.5)) +
  theme(
    legend.position.inside = c(0.16, 0.04),
    plot.tag.position = c(0.59, 0.76),
    plot.tag = element_text(color = "grey30")
  ) +
  annotate("text", label = "Canada", x = -113, y = 49.9, color = "gray50", size = 2.8) +
  annotate("text", label = "Mexico", x = -110, y = 30.5, color = "gray50", size = 2.8) +
  annotate("text", label = "Pacific\nOcean", x = -123.8, y = 35, color = "gray50", size = 2.8)

p2 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = diffused), color = NA) +
  labs(
    tag = "(b) Diffused log(population density)",
    fill = NULL
  ) +
  facet_wrap(~ factor(abb_name, levels = c(
    "<i>S. vulgaris</i>",
    "<i>P. melanocephalus</i>"
  )), ncol = 1) +
  geom_text(
    data = diff_label, aes(label = label),
    x = -109.3, y = 30.6, color = "gray20", size = 2.8
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag.position = c(0.48, 1),
    plot.tag = element_text(color = "grey30")
  )

p3 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = `log count diffused`), color = NA) +
  labs(
    tag = "(c) log(predicted count)",
    fill = NULL
  ) +
  scale_fill_viridis(option = "magma") +
  facet_wrap(~ factor(abb_name, levels = c(
    "<i>S. vulgaris</i>",
    "<i>P. melanocephalus</i>"
  )), ncol = 1) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(color = "grey30")
  )

(p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/bb_maps.pdf"), width = 20, height = 16, unit = "cm", device = cairo_pdf)


# Supporting info plot
mp_bb_fc +
  geom_sf(data = bb_stuff, aes(fill = diffused - orig), color = NA) +
  labs(fill = "Difference between diffused and original covariate")

ggsave(paste0(here::here(), "/results/figures/supporting/bb_diffused_original.pdf"), width = 22, height = 23, unit = "cm")


# EBS map plots
source(file.path(root_dir, "functions/ebs-map-plot.R"))

ebs_stuff <-
  readRDS(paste0(here::here(), "/results/2024-09-24_identity_LNP/stuff_gz_df.rds")) |>
  left_join(ebs_names, by = c("species" = "X"))

sub <- ebs_stuff |> dplyr::filter(species %in% c("a_poll", "cap"))

diff_label <- ebs |>
  mutate(species = abb_name) |> # FIXME: use latin names
  filter(abb_name %in% unique(sub$species)) |>
  mutate(label = paste0("log(\u03BA)=", round(ln_kappa2, digits = 2)))

# Main plot
p1 <- mp_ebs_s +
  facet_wrap(~NA) +
  geom_sf(data = sub, aes(fill = orig, color = orig), linewidth = 0.01) +
  labs(
    fill = NULL,
    tag = "(a) Depth"
  ) +
  theme(
    legend.position.inside = c(0.16, 0.07),
    plot.tag.position = c(0.59, 0.76),
    plot.tag = element_text(color = "grey30")
  ) +
  scale_fill_viridis(
    option = "mako", na.value = NA, # trans = "fourth_root_power",
    # breaks = c(0, 1, 4),
    direction = -1
  ) +
  scale_color_viridis(
    option = "mako", na.value = NA, # trans = "fourth_root_power",
    direction = -1
  ) +
  geom_sf(linewidth = 0.4, color = "gray40") +
  annotate("text", label = "USA", x = -160, y = 62, color = "gray50", size = 2.8) +
  annotate("text", label = "Bering\nSea", x = -176, y = 65, color = "gray50", size = 2.8)

p2a <- mp_ebs_s +
  geom_sf(data = sub |> filter(species == "a_poll"), aes(fill = diffused, color = diffused), linewidth = 0.01) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag.position = c(0.48, 1),
    legend.position.inside = c(0.17, 0.07),
    plot.tag = element_text(color = "grey30")
  ) +
  geom_text(
    data = diff_label |> filter(species == "a_poll"), aes(label = label),
    x = -158, y = 53.8, color = "grey20", size = 2.8
  ) +
  labs(
    tag = "(b) Diffused depth\n",
    fill = NULL
  ) +
  scale_fill_viridis(option = "mako", na.value = NA, direction = -1) +
  scale_color_viridis(option = "mako", na.value = NA, direction = -1) +
  geom_sf(linewidth = 0.4, color = "gray40")

p2b <- mp_ebs_s +
  geom_sf(data = sub |> filter(species == "cap"), aes(fill = diffused, color = diffused), linewidth = 0.01) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position.inside = c(0.17, 0.07)
  ) +
  geom_text(
    data = diff_label |> filter(species == "cap"), aes(label = label),
    x = -158, y = 53.8, color = "grey20", size = 2.8
  ) +
  labs(fill = NULL) +
  scale_fill_viridis(option = "mako", na.value = NA, direction = -1, breaks = c(-0.3, 0.2, 0.7)) +
  scale_color_viridis(option = "mako", na.value = NA, direction = -1) +
  geom_sf(linewidth = 0.4, color = "gray40")

p2 <- p2a / p2b

p3a <- mp_ebs_s +
  geom_sf(
    data = sub |> filter(species == "a_poll"),
    aes(fill = `log(dens) diffused`, color = `log(dens) diffused`),
    linewidth = 0.01
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag.position = c(0.5, 1),
    legend.position.inside = c(0.17, 0.07),
    plot.tag = element_text(color = "grey30")
  ) +
  labs(
    tag = "(c) log(predicted density)\n",
    fill = NULL
  ) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA) +
  geom_sf(linewidth = 0.4, color = "gray40")

p3b <- mp_ebs_s +
  geom_sf(
    data = sub |> filter(species == "cap"),
    aes(fill = `log(dens) diffused`, color = `log(dens) diffused`),
    linewidth = 0.01
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position.inside = c(0.17, 0.07)
  ) +
  labs(fill = NULL) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA) +
  geom_sf(linewidth = 0.4, color = "gray40")

p3 <- p3a / p3b

(p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/ebs_maps.pdf"), width = 20, height = 16, unit = "cm", device = cairo_pdf)

# Plotting correlations between omega
king <- ebs_stuff |> dplyr::filter(species %in% c("rking"))

p1a <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `orig`,
        Diffused = `diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Original"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  labs(
    fill = NULL,
    tag = "(a) Depth"
  ) +
  theme(
    plot.tag.position = c(0.59, 1.02),
    plot.tag = element_text(color = "grey30"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position.inside = c(0.17, 0.06)
  ) +
  scale_fill_viridis(option = "mako", na.value = NA, direction = -1) +
  scale_color_viridis(option = "mako", na.value = NA, direction = -1) +
  geom_sf(linewidth = 0.4, color = "gray40") +
  annotate("text", label = "USA", x = -160, y = 62, color = "gray50", size = 2.8) +
  annotate("text", label = "Bering\nSea", x = -176, y = 65, color = "gray50", size = 2.8)

p1b <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `orig`,
        Diffused = `diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Diffused"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  labs(fill = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position.inside = c(0.17, 0.06)
  ) +
  scale_fill_viridis(option = "mako", na.value = NA, direction = -1, breaks = c(-0.5, 0.3, 1.1)) +
  scale_color_viridis(option = "mako", na.value = NA, direction = -1) +
  geom_sf(linewidth = 0.4, color = "gray40")

p1 <- p1a / p1b

p2a <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `omega orig`,
        Diffused = `omega diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Original"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.tag.position = c(0.48, 1.02),
    plot.tag = element_text(color = "grey30"),
    legend.position.inside = c(0.17, 0.06)
  ) +
  labs(tag = "(b) Spatial random effects", fill = NULL) +
  scale_fill_gradient2() +
  scale_color_gradient2() +
  geom_sf(linewidth = 0.4, color = "gray40")

p2b <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `omega orig`,
        Diffused = `omega diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Diffused"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position.inside = c(0.17, 0.06)
  ) +
  labs(fill = NULL) +
  scale_fill_gradient2() +
  scale_color_gradient2() +
  geom_sf(linewidth = 0.4, color = "gray40")

p2 <- p2a / p2b

p3a <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `log(dens) orig`,
        Diffused = `log(dens) diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Original"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag.position = c(0.5, 1.02),
    plot.tag = element_text(color = "grey30"),
    legend.position.inside = c(0.17, 0.06)
  ) +
  labs(tag = "(c) log(predicted density)", fill = NULL) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA) +
  geom_sf(linewidth = 0.4, color = "gray40")

p3b <- mp_ebs_s +
  geom_sf(
    data = king |>
      rename(
        Original = `log(dens) orig`,
        Diffused = `log(dens) diffused`
      ) |>
      pivot_longer(c(Original, Diffused), names_to = "var") |>
      filter(var == "Diffused"),
    aes(fill = value, color = value), linewidth = 0.01
  ) +
  facet_wrap(~ factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position.inside = c(0.17, 0.06)
  ) +
  labs(fill = NULL) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA) +
  geom_sf(linewidth = 0.4, color = "gray40")

p3 <- p3a / p3b

(p1 | p2 | p3) +
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/supporting/king_map.pdf"), width = 20, height = 16, unit = "cm")

# Plot correlation
bb_cor <-
  readRDS(paste0(here::here(), "/results/2024-08-07_LNP/stuff_gz_df.rds")) |>
  separate(species, "_", into = c("family", "sp"), remove = FALSE) |>
  mutate(
    abb_name = substring(family, 1, 1),
    abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>")
  ) |>
  as_tibble() |>
  summarise(
    cor_pred = cor(`log count orig`, `log count diffused`),
    cor_omega = cor(`omega diffused`, `omega orig`),
    cor_covar = cor(`orig`, `diffused`),
    .by = abb_name
  )

ebs_cor <-
  readRDS(paste0(here::here(), "/results/2024-09-24_identity_LNP/stuff_gz_df.rds")) |>
  left_join(ebs_names, by = c("species" = "X")) |>
  filter(species %in%
    filter(ebs, `Diffusion\nfavoured` == "Y")$X) |>
  drop_na() |>
  as_tibble() |>
  summarise(
    cor_pred = cor(`log(dens) orig`, `log(dens) diffused`),
    cor_omega = cor(`omega diffused`, `omega orig`),
    cor_covar = cor(`diffused`, `orig`),
    .by = abb_name
  )

cor <- bind_rows(
  ebs_cor |> mutate(case = "Fish"),
  bb_cor |> mutate(case = "Bird")
) |>
  mutate(abb_name = str_remove_all(abb_name, "<i>|</i>"))

p1 <- cor |>
  ggplot(aes(cor_covar, cor_pred, label = abb_name, shape = case)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, fontface = "italic") +
  theme(aspect.ratio = 1) +
  labs(
    x = "Correlation between\ncovariates",
    y = "Correlation between\n predictions",
    shape = ""
  ) +
  xlim(0.15, 1) +
  ylim(0.15, 1)

p2 <- cor |>
  ggplot(aes(cor_covar, cor_omega, label = abb_name, color = cor_pred, shape = case)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, fontface = "italic") +
  theme(aspect.ratio = 1) +
  xlim(0.15, 1) +
  ylim(0.15, 1) +
  labs(
    x = "Correlation between\ncovariates",
    y = "Correlation between\nspatial random effects",
    color = "Correlation\nbetween\npredictions",
    shape = ""
  ) +
  NULL

(p1 | p2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")", tag_prefix = "(") +
  plot_layout(guides = "collect")

ggsave(paste0(here::here(), "/results/figures/supporting/correlations.pdf"), width = 20, height = 9, unit = "cm")


# Test plotting a few cases where AIC favours diffusion but the correlation is 1
test_sp <- ebs |>
  filter(deltaAIC > 2) |>
  arrange(covar_corr)

test_df <- ebs |>
  filter(abb_name %in% c(
    head(test_sp$abb_name, 3),
    tail(test_sp$abb_name, 3)
  ))

ebs_stuff_test <- ebs_stuff |>
  #rename(abb_name = species) |>
  filter(abb_name %in% unique(test_df$abb_name))

ebs_stuff_test <- ebs_stuff_test |>
  left_join(test_df |> dplyr::select(abb_name, deltaAIC, covar_corr))


# Plot partial effects
diff_sp <- ebs_stuff |>
  dplyr::filter(abb_name %in% test_sp$abb_name) |>
  dplyr::filter(!species == "a_yil") |>
  rename(
    Original = `pdepth orig`,
    Diffused = `pdepth diffused`
  ) |>
  mutate(diff = Diffused - Original) |>
  mutate(
    diff_scaled = 2 * (diff - min(diff, na.rm = TRUE)) / (max(diff, na.rm = TRUE) - min(diff, na.rm = TRUE)) - 1,
    .by = species
  )

mp_ebs_fc +
  geom_sf(
    data = diff_sp |>
      filter(species %in% c("a_kam", "a_akplce", "cap", "j_akplce", "j_great", "rking")),
    aes(fill = diff_scaled, color = diff_scaled), linewidth = 0.01
  ) +
  facet_wrap(~abb_name) +
  scale_fill_gradient2(na.value = NA, name = "Scaled difference between Diffusion and Original partial effect") +
  scale_color_gradient2(na.value = NA) +
  geom_sf(linewidth = 0.4, color = "gray40") +
  guides(color = "none")

ggsave(paste0(here::here(), "/results/figures/supporting/ebs_partial.pdf"), width = 20, height = 18, unit = "cm")
