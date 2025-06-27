# Load libraries
library(tidyr)
library(dplyr)
library(ggtext)
library(patchwork)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())
library(viridis)
library(readr)

root_dir <- here::here(".")

### Plot results from the case studies
ebs_names <-
  read_csv(paste0(root_dir, "/data/clean_EBS_species.csv")) |>
  dplyr::select(X, abb_name, sc_name, common_name2)

# Read data
ebs <-
  read.csv(paste0(root_dir, "/results/2025-06-27_identity_LNP_main_cutoff/main_mesh_Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  rename(covar_corr = corr_depth) |>
  left_join(ebs_names, by = "X")

ebs_high <-
  read.csv(paste0(root_dir, "/results/2025-06-27_identity_LNP_high_cutoff/high_mesh_Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  rename(covar_corr = corr_depth) |>
  left_join(ebs_names, by = "X")

ebs_low <-
  read.csv(paste0(root_dir, "/results/2025-06-27_identity_LNP_low_cutoff/low_mesh_Results_cz.csv")) |>
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC < 0, "N", "Y")) |>
  left_join(ebs_names, by = "X")

dd <- bind_rows(
  ebs_low |> mutate(`Mesh cutoff` = "0.11"),
  ebs |> mutate(`Mesh cutoff` = "0.1"),
  ebs_high |> mutate(`Mesh cutoff` = "0.09")
)

d_wide <- left_join(
  ebs_low |> dplyr::select(abb_name, deltaAIC) |> rename(low_deltaAIC = deltaAIC),
  ebs_high |> dplyr::select(abb_name, deltaAIC) |> rename(high_deltaAIC = deltaAIC)
) |>
  left_join(ebs |> dplyr::select(abb_name, deltaAIC) |> rename(mid_deltaAIC = deltaAIC)) |>
  mutate(low_deltaAIC = ifelse(low_deltaAIC <= -2, NA, low_deltaAIC),
         mid_deltaAIC = ifelse(mid_deltaAIC <= -2, NA, mid_deltaAIC),
         high_deltaAIC = ifelse(high_deltaAIC <= -2, NA, high_deltaAIC)) |>
  mutate(xmin = pmin(low_deltaAIC, mid_deltaAIC, high_deltaAIC, na.rm = TRUE),
         xmax = pmax(low_deltaAIC, mid_deltaAIC, high_deltaAIC, na.rm = TRUE))

d_wide |> filter(abb_name == "<i>H. stenolepis</i>")

order <- dd |> summarise(deltaAIC = mean(deltaAIC), .by = abb_name) |> arrange(desc(deltaAIC)) |> pull(abb_name)

sort(unique(dd$abb_name))

# dd <- dd |>
#   tidylog::filter(deltaAIC > -100)

df <- dd |>
  dplyr::select(X, deltaAIC, `Mesh cutoff`) |>
  pivot_wider(names_from = `Mesh cutoff`, values_from = deltaAIC) |>
  rowwise() |>
  mutate(
    test = ifelse(2 >= min(c(`0.11`, `0.1`, `0.09`), na.rm = TRUE) &
                    2 <= max(c(`0.11`, `0.1`, `0.09`), na.rm = TRUE),
                  "No", "Yes")
  ) |>
  ungroup()

dd <- left_join(dd, df, by = "X") |> drop_na(test)

ddlong <- dd |>
  pivot_longer(c(`0.11`, `0.1`, `0.09`)) |>
  distinct(abb_name, name, .keep_all = TRUE) |>
  filter(value >-2)

ddlong |>
  distinct(abb_name, .keep_all = TRUE) |>
  summarise(n = n(), .by = test)

ggplot(ddlong, aes(value, factor(abb_name, levels = order), shape = name, color = test)) +
  geom_vline(xintercept = c(-2, 2), alpha = 0.5, linetype = 2, linewidth = 0.35) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 1, linewidth = 0.35) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 1, linewidth = 0.35) +
  geom_segment(data = d_wide,
               aes(xmin, xend = xmax, factor(abb_name, levels = order)),
               inherit.aes = FALSE, alpha = 0.3) +
  scale_shape_manual(values = c(21, 22, 23)) +
  geom_point(size = 2.3) +
  labs(x = "ΔAIC", y = "Species", color = "mAIC consistency", shape = "Mesh cutoff") +
  scale_x_continuous(
    trans = "fourth_root_power",
    breaks = c(-2, 0, 2, 40, 80, 120)
  ) +
  theme(
    axis.text.y = element_markdown(),
    axis.text.x = element_text(size = 6.5),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.direction = "horizontal"
  ) +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  guides(
    color = guide_legend(ncol = 1,
                         title.position = "top"),
    shape = guide_legend(
      title.position = "top",
      ncol = 1,
      title.hjust = 0.5,
      override.aes = list(alpha = 1, fill = "black", stroke = 0.5)
    )
  ) +
  NULL

ggsave(paste0(here::here(), "/results/figures/supporting/mesh-sensi-ebs.pdf"), width = 17, height = 19, unit = "cm", device = cairo_pdf)

