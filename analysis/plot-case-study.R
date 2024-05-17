# Load libraries
library(tidyr)
library(dplyr)
library(ggtext)
library(ggstats)
library(patchwork)
library(ggplot2)
library(ggsidekick); theme_set(theme_sleek())
library(terra) # Load terra to use extract... needed for prep-data

### Explore results from the case study

# Read data
bering <- read.csv(paste0(here::here(), "/results/2024-05-16_identity_LNP/Results_cz.csv")) %>%
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC > 0, "N", "Y"),
         abs_delta_aic = ifelse(abs(deltaAIC) > 2, "<2", ">2"))

bb <- read.csv(paste0(here::here(), "/results/2024-05-16__Poisson/Results_bb_cz.csv")) %>%
  mutate("Diffusion\nfavoured" = ifelse(deltaAIC > 0, "N", "Y"),
         abs_delta_aic = ifelse(abs(deltaAIC) > 2, "<2", ">2"))

# Plot $\Delta$AIC, correlation between depth and diffused depth, the range estimate and ln_kappa2
# The range of range estimates is pretty long so I log them for now!

# Bering Sea case
p1 <- ggplot(bering, aes(deltaAIC, reorder(X, desc(deltaAIC)),
                         color = `Diffusion\nfavoured`, shape = abs_delta_aic)) +
  geom_point(fill = NA) +
  labs(x = "ΔAIC", y = "Species") +
  guides(color = guide_legend(position = "inside"),
         shape = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.77, 0.75)) +
  scale_shape_manual(values = c(19, 21), name = "abs(ΔAIC)") +
  scale_color_brewer(palette = "Dark2", direction = -1)

p2 <- ggplot(bering, aes(corr_depth, reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(x = "Correlation between depth\nand diffused depth",
       y = "Species")

p3 <- ggplot(bering, aes(log(range), reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(y = "Species")

p4 <- ggplot(bering, aes(ln_kappa2, reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(y = "Species")

p1 + p2 + p3 + p4 +
  plot_layout(axes = "collect") +
  plot_annotation("Bering Sea case study", tag_levels = "A") &
  geom_stripped_rows(aes(y = X), inherit.aes = FALSE)

ggsave(paste0(here::here(), "/results/bering_case_overview.pdf"), width = 20, height = 28, unit = "cm", device = cairo_pdf)

# Looks like in adult yellowfin sole (in red axis text), the diffusion model is parsimonious,
# the correlation between diffused depth is intermediate (low correlations are very clear in the map,
# showing as very smoothed/diffused), and the range is not extreme either. This makes it a great example.
# It also seems in about half the species, the diffusion model is not parsimonious, whereas for the other half,
# either the diffusion model is parsimonious or they are similar. There doesn't seem to be an relationship
# between the $\Delta$AIC and the correlation coefficient, nor the range.

# Breeding bird case
p1 <- ggplot(bb, aes(deltaAIC, reorder(X, desc(deltaAIC)),
                     color = `Diffusion\nfavoured`, shape = abs_delta_aic)) +
  geom_point(fill = NA) +
  labs(x = "ΔAIC", y = "Species") +
  guides(color = guide_legend(position = "inside"),
         shape = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.77, 0.75)) +
  scale_shape_manual(values = c(19, 21), name = "abs(ΔAIC)") +
  scale_color_brewer(palette = "Dark2", direction = -1)

p2 <- ggplot(bb, aes(corr_pop_dens, reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(x = "Correlation between log pop dens\nand diffused pop dens",
       y = "Species")

p3 <- ggplot(bb, aes(log(range), reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(y = "Species")

p4 <- ggplot(bb, aes(ln_kappa2, reorder(X, desc(deltaAIC)))) +
  geom_point() +
  labs(y = "Species")

p1 + p2 + p3 + p4 +
  plot_layout(axes = "collect") +
  plot_annotation("Breeding bird case study", tag_levels = "A") &
  geom_stripped_rows(aes(y = X), inherit.aes = FALSE)

ggsave(paste0(here::here(), "/results/bb_case_overview.pdf"), width = 20, height = 28, unit = "cm", device = cairo_pdf)
