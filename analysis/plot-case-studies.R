# Load libraries
library(tidyr)
library(dplyr)
library(ggtext)
library(ggstats)
library(patchwork)
library(ggplot2)
library(ggsidekick); theme_set(theme_sleek())
library(viridis)
library(sf)
library(forcats)

root_dir <- here::here(".")

### Plot results from the case studies
# Read data
ebs <- read.csv(paste0(root_dir, "/results/2024-08-12_identity_LNP/Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y"),
         abs_delta_aic = ifelse(abs(deltaAIC) > 2, "<2", ">2")) |>
  rename(covar_corr = corr_depth) |>
  mutate(abb_name = X)

# Breeding bird case
bb <- read.csv(paste0(root_dir, "/results/2024-08-07_LNP/Results_bb_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y"),
         abs_delta_aic = ifelse(abs(deltaAIC) > 2, "<2", ">2")) |>
  rename(covar_corr = corr_pop_dens) |>
  separate(X, "_", into = c("family", "species")) |>
  mutate(abb_name = substring(family, 1, 1),
         abb_name = paste0("<i>", paste(paste0(abb_name, "."), species), "</i>"))

dd <- bind_rows(ebs |> mutate(case = "Eastern Bering sea fishes"),
                bb |> mutate(case = "Breeding bird survey"))

rect <- data.frame()

ggplot(dd, aes(deltaAIC, reorder(abb_name, desc(deltaAIC)), fill = covar_corr)) +
  geom_rect(aes(xmin = 2, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey95") +
  geom_vline(xintercept = c(-2, 2), alpha = 0.5, linetype = 2, linewidth = 0.35) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 1, linewidth = 0.35) +
  geom_point(shape = 21, color = "grey10", stroke = 0.01, size = 2.3) +
  labs(x = "ΔAIC", y = "Species", fill = "Correlation between raw and diffused covariate") +
  scale_x_continuous(trans = "fourth_root_power") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~case, scales = "free", ncol = 2) +
  theme(axis.text.y = element_markdown(),
        axis.text.x = element_text(size = 6.5),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.direction = "horizontal") +
  scale_fill_viridis(option = "rocket")

ggsave(paste0(here::here(), "/results/figures/case_summary.pdf"), width = 17, height = 16, unit = "cm", device = cairo_pdf)



# Plot maps from the case studies
# BB
source(file.path(root_dir, "functions/bb-map-plot.R"))

bb_stuff <-
  readRDS(paste0(here::here(), "/results/2024-08-07_LNP/stuff_gz_df.rds")) |>
  separate(species, "_", into = c("family", "sp"), remove = FALSE) |>
  mutate(abb_name = substring(family, 1, 1),
         abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>"))

sub <- bb_stuff |> filter(species %in% c("Sturnus_vulgaris", "Pheucticus_melanocephalus"))

# Main plot?!
p1 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = orig), color = NA) +
  labs(tag = "(a) log(population density)",
       fill = NULL) +
  facet_wrap(~NA) +
  guides(fill = guide_colorbar(position = "inside", title.position = "top", title.hjust = 0.5)) +
  theme(legend.position.inside = c(0.16, 0.07),
        plot.tag.position = c(0.59, 0.76),
        plot.tag = element_text(color = "grey30"))

p2 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = diffused), color = NA) +
  labs(tag = "(b) Diffused log(population density)",
       fill = NULL) +
  facet_wrap(~factor(abb_name, levels = c("<i>S. vulgaris</i>",
                                          "<i>P. melanocephalus</i>")), ncol = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.tag.position = c(0.48, 1),
        plot.tag = element_text(color = "grey30"))

p3 <- mp_bb_s +
  geom_sf(data = sub, aes(fill = `log count diffused`), color = NA) +
  labs(tag = "(c) log(predicted count)",
       fill = NULL) +
  scale_fill_viridis(option = "magma") +
  facet_wrap(~factor(abb_name, levels = c("<i>S. vulgaris</i>",
                                          "<i>P. melanocephalus</i>")), ncol = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.tag.position = c(0.5, 1),
        plot.tag = element_text(color = "grey30"))

(p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/bb_maps.pdf"), width = 20, height = 16, unit = "cm")


# Supporting info plot
mp_bb_fc +
  geom_sf(data = bb_stuff, aes(fill = diffused - orig), color = NA) +
  labs(fill = "Difference between diffused and original covariate")

ggsave(paste0(here::here(), "/results/figures/supporting/bb_diffused_original.pdf"), width = 22, height = 23, unit = "cm")


# EBS map plots
source(file.path(root_dir, "functions/ebs-map-plot.R"))

ebs_stuff <-
  readRDS(paste0(here::here(), "/results/2024-08-07_identity_LNP/stuff_gz_df.rds")) #|>
# separate(species, "_", into = c("family", "sp"), remove = FALSE) |>
# mutate(abb_name = substring(family, 1, 1),
#        abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>"))

sub <- ebs_stuff |> dplyr::filter(species %in% c("a_poll", "cap"))

# Main plot?!
p1 <- mp_ebs_s +
  facet_wrap(~NA) +
  geom_sf(data = sub, aes(fill = orig, color = orig), linewidth = 0.01) +
  labs(fill = NULL,
       tag = "(a) log(depth)") +
  theme(legend.position.inside = c(0.16, 0.07),
        plot.tag.position = c(0.59, 0.76),
        plot.tag = element_text(color = "grey30")) +
  scale_fill_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                     breaks = c(1, 3, 6, 9), direction = -1) +
  scale_color_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                      direction = -1)

p2 <- mp_ebs_s +
  geom_sf(data = sub, aes(fill = diffused, color = diffused), linewidth = 0.01) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.tag.position = c(0.48, 1),
        plot.tag = element_text(color = "grey30")) +
  labs(tag = "(b) Diffused log(depth)",
       fill = NULL) +
  scale_fill_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                     breaks = c(1, 3, 6, 9), direction = -1) +
  scale_color_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                      direction = -1)

p3 <- mp_ebs_s +
  geom_sf(data = sub, aes(fill = `log(dens) diffused`, color = `log(dens) diffused`),
         linewidth = 0.01) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.tag.position = c(0.5, 1),
        plot.tag = element_text(color = "grey30")) +
  labs(tag = "(c) log(predicted density)",
       fill = NULL) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA)

(p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/ebs_maps.pdf"), width = 20, height = 16, unit = "cm")


# Supporting info plot
sp <- unique(ebs_stuff$species)

# ebs_stuff |>
#   as_tibble() |>
#   mutate(diff = diffused - orig) |>
#   summarise(min = min(diff),
#             max = max(diff),
#             .by = species) |>
#   arrange(min) |>
#   as.data.frame()
#
# ggplot(ebs_stuff |> filter(species == "a_yil"), aes(diff)) +
#   geom_histogram()

# Trim some quantiles to make the plots clearer
ebs_stuff <- ebs_stuff |>
  mutate(diff = diffused - orig) |>
  mutate(diff_lwr = quantile(diff, probs = 0.01),
         diff_upr = quantile(diff, probs = 0.99),
         .by = species) |>
  mutate(diff = ifelse(diff < diff_lwr, diff_lwr, diff),
         diff = ifelse(diff > diff_upr, diff_upr, diff))


cors <- ebs |>
  rename(species = X) |>
  dplyr::select(covar_corr, species)

mp_ebs_fc +
  geom_sf(data = ebs_stuff |> filter(species %in% sp[1:22]), aes(fill = diff),
          color = NA) +
  geom_text(data = cors |>
              filter(species %in% sp[1:22]),
            color = "grey20",
            aes(x = -Inf, y = -Inf, label = round(covar_corr, digits = 3)),
            hjust = -0.1, vjust = -14, size = 3) +
  facet_wrap(~species, ncol = 6) +
  labs(fill = "Difference between diffused and original covariate")

ggsave(paste0(here::here(), "/results/figures/supporting/ebs_diffused_original_1.pdf"), width = 22, height = 23, unit = "cm")

mp_ebs_fc +
  geom_sf(data = ebs_stuff |> filter(species %in% sp[23:length(sp)]), aes(fill = diff),
          color = NA) +
  geom_text(data = cors |>
              filter(species %in% sp[23:length(sp)]),
            color = "grey20",
            aes(x = -Inf, y = -Inf, label = round(covar_corr, digits = 3)),
            hjust = -0.1, vjust = -14, size = 3) +
  facet_wrap(~species, ncol = 6) +
  labs(fill = "Difference between diffused and original covariate")

ggsave(paste0(here::here(), "/results/figures/supporting/ebs_diffused_original_2.pdf"), width = 22, height = 23, unit = "cm")


# FIXME: In the EBS data, I think because the grid is so high-res, I see cell boundaries from geom_sf, even with color = null.
# The only way to remove this is to set color also to a variable (the fill). I wonder if this distorts the pattern, because the line has a width...


# Plotting correlations between omega
king <- ebs_stuff |> dplyr::filter(species %in% c("rking"))

p1 <- mp_ebs_s +
  geom_sf(data = king |>
            rename(Original = `orig`,
                   Diffused = `diffused`) |>
            pivot_longer(c(Original, Diffused), names_to = "var"),
          aes(fill = value, color = value), linewidth = 0.01) +
  facet_wrap(~factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  labs(fill = NULL,
       tag = "(a) log(depth)") +
  theme(plot.tag.position = c(0.59, 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),) +
  scale_fill_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                     breaks = c(1, 3, 6, 9), direction = -1) +
  scale_color_viridis(option = "mako", na.value = NA, trans = "fourth_root_power",
                      direction = -1)

p2 <- mp_ebs_s +
  geom_sf(data = king |>
            rename(Original = `omega orig`,
                   Diffused = `omega diffused`) |>
            pivot_longer(c(Original, Diffused), names_to = "var"),
          aes(fill = value, color = value), linewidth = 0.01) +
  facet_wrap(~factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.tag.position = c(0.48, 1)) +
  labs(tag = "(b) Spatial random effects",
       fill = NULL) +
  scale_fill_gradient2() +
  scale_color_gradient2()

p3 <- mp_ebs_s +
  geom_sf(data = king |>
            rename(Original = `log(dens) orig`,
                   Diffused = `log(dens) diffused`) |>
            pivot_longer(c(Original, Diffused), names_to = "var"),
          aes(fill = value, color = value), linewidth = 0.01) +
  facet_wrap(~factor(var, levels = c("Original", "Diffused")), ncol = 1) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.tag.position = c(0.5, 1)) +
  labs(tag = "(c) log(predicted density)",
       fill = NULL) +
  scale_fill_viridis(option = "magma", na.value = NA) +
  scale_color_viridis(option = "magma", na.value = NA)

(p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 11))

ggsave(paste0(here::here(), "/results/figures/supporting/king_map.pdf"), width = 20, height = 16, unit = "cm")

# Plot correlation
bb_cor <-
  readRDS(paste0(here::here(), "/results/2024-08-07_LNP/stuff_gz_df.rds")) |>
  separate(species, "_", into = c("family", "sp"), remove = FALSE) |>
  mutate(abb_name = substring(family, 1, 1),
         abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>")) |>
  as_tibble() |>
  summarise(cor_pred = cor(`log count orig`, `log count diffused`),
            cor_omega = cor(`omega diffused`, `omega orig`),
            cor_covar = cor(`orig`, `diffused`),
            .by = abb_name)

ebs_cor <-
  readRDS(paste0(here::here(), "/results/2024-08-07_identity_LNP/stuff_gz_df.rds")) |>
  filter(species %in%
           filter(ebs, `Diffusion\nfavoured` == "Y")$abb_name) |>
  drop_na() |>
  as_tibble() |>
  summarise(cor_pred = cor(`log(dens) orig`, `log(dens) diffused`),
            cor_omega = cor(`omega diffused`, `omega orig`),
            cor_covar = cor(`diffused`, `orig`),
            .by = species)

cor <- bind_rows(ebs_cor |> mutate(case = "Fish"),
                 bb_cor |> mutate(case = "Bird") |> rename(species = abb_name))

p1 <- cor |>
  ggplot(aes(cor_covar, cor_pred, label = species, shape = case)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5) +
  theme(aspect.ratio = 1) +
  labs(x = "Correlation between\ncovariates",
       y = "Correlation between\n predictions",
       shape = "") +
  xlim(0.15, 1) +
  ylim(0.15, 1)

p2 <- cor |>
  ggplot(aes(cor_covar, cor_omega, label = species, color = cor_pred, shape = case)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5) +
  theme(aspect.ratio = 1) +
  xlim(0.15, 1) +
  ylim(0.15, 1) +
  labs(x = "Correlation between\ncovariates",
       y = "Correlation between\nspatial random effects",
       color = "Correlation\nbetween\npredictions",
       shape = "") +
  NULL

(p1 | p2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")", tag_prefix = "(") +
  plot_layout(guides = "collect")

ggsave(paste0(here::here(), "/results/figures/supporting/correlations.pdf"), width = 20, height = 9, unit = "cm")


