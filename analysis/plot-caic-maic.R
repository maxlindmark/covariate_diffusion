# Script to produce a supporting figure showing difference between marginal and conditional AIC

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

root_dir <- here::here(".")

### Plot results from the case studies
# Read data
ebs_names <- read_csv(paste0(root_dir, "/data/clean_EBS_species.csv")) |>
  dplyr::select(X, abb_name, sc_name, common_name2)

ebs <-
  read.csv(paste0(root_dir, "/results/2025-07-14_identity_LNP/Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  rename(covar_corr = corr_depth) |>
  left_join(ebs_names, by = "X")

ebs <- ebs |> filter(!abb_name == "<i>P. camtschaticus</i> (Bb)")

options(scipen=999)
ebs |>
  arrange(deltaCAIC) |>
  dplyr::select(X, deltaCAIC, deltaAIC, range, ln_kappa2)

ebs |>
  arrange(deltaAIC) |>
  dplyr::select(X, deltaCAIC, deltaAIC, range, ln_kappa2)

# Breeding bird case
bb <-
  #read.csv(paste0(root_dir, "/results/2025-06-28_LNP/Results_cz.csv")) |>
  read.csv(paste0(root_dir, "/results/2025-07-14_LNP/Results_cz.csv")) |>
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
) |>
  mutate(cAIC_mAIC_diff = deltaAIC - deltaCAIC)

head(dd)

dd2 <- dd |>
  mutate(test = "No difference in favoured model",
         test = ifelse(deltaAIC > 2 & deltaCAIC < 2, "ΔmAIC > 2 & ΔcAIC < 2", test),
         test = ifelse(deltaAIC < 2 & deltaCAIC > 2, "ΔmAIC < 2 & ΔcAIC > 2", test))

dd2 |> arrange(deltaCAIC)

# dd2 <- dd2 |>
#   filter(cAIC_mAIC_diff > -500 & cAIC_mAIC_diff < 500)

dd2 |>
  summarise(n = n(), .by = c(case, test))

dd2 |> filter(case == "Eastern Bering sea fishes") |> filter(deltaAIC > 2) |> nrow()
dd2 |> filter(case == "Breeding bird survey") |> filter(deltaAIC > 2) |> nrow()

dd2 |> filter(case == "Eastern Bering sea fishes") |> filter(deltaCAIC > 2) |> nrow()
dd2 |> filter(case == "Breeding bird survey") |> filter(deltaCAIC > 2) |> nrow()

dd2 |>
  summarise(n = length(unique(abb_name)), .by = case)

dd2 |> filter(is.na(cAIC_mAIC_diff))

dd2 |>
  drop_na(cAIC_mAIC_diff) |>
  ggplot(aes(cAIC_mAIC_diff, reorder(abb_name, desc(cAIC_mAIC_diff)), color = test)) +
  geom_vline(xintercept = c(-2, 2), alpha = 0.5, linetype = 2, linewidth = 0.35) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 1, linewidth = 0.35) +
  geom_point(size = 2.3) +
  labs(x = "ΔmAIC - ΔcAIC", y = "Species") +
  guides(color = guide_legend(ncol = 1)) +
  facet_wrap(~case, scales = "free", ncol = 2) +
  scale_x_continuous(
    trans = "fourth_root_power",
    breaks = c(-2, 0, 2, 50, 250)
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    axis.text.y = element_markdown(),
    axis.text.x = element_text(size = 6.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.direction = "horizontal"
  )

ggsave(paste0(here::here(), "/results/figures/supporting/maic_caic.pdf"), width = 17, height = 19, unit = "cm", device = cairo_pdf)
