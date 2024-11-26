# FROM: C:\Users\James.Thorson\Desktop\Work files\Collaborations\2024 -- movement kernel for covariates\Old code\covariate_diffusion_2024-02-05.R

library(fmesher)
library(Matrix)
library(sf)
library(ggplot2)
library(viridis)
library(ggsidekick)
library(tidyr)
library(dplyr)
library(patchwork)

set.seed(194301)

#
loc <- matrix(runif(1000), ncol = 2)
mesh <- fm_mesh_2d(loc)
spde <- fm_fem(mesh)
invM0 <- spde$c0
diag(invM0) <- 1 / diag(spde$c0)

#
sf_mesh <- st_union(fm_as_sfc(mesh))
sf_grid <- st_make_grid(sf_mesh, n = 100)
sf_grid <- st_intersection(sf_grid, sf_mesh)
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(sf_grid)))$proj$A

# Loop through a couple of kappa values to visualize diffusion

stuff_list <- list()
sum_list <- list()

# Using IID normal deviates. Hold this constant across kappas
vec2 <- rnorm(mesh$n)

for (i in c("0", "2.5", "5")) {
  ln_kappa <- as.numeric(i)
  ln_tau <- log(1 / (1 + exp(2 * ln_kappa)))

  invD <- exp(ln_tau) * invM0 %*% (spde$c0 + exp(2 * ln_kappa) * spde$c0 + spde$g1)
  D <- solve(invD)

  which_mid <- which.min(rowSums(scale(mesh$loc[, 1:2])^2))

  # Diffusion from a specific point
  vec1 <- rep(0, mesh$n)
  vec1[which_mid] <- 1

  # Using IID normal deviates
  # vec2 = rnorm(mesh$n)

  stuff <- st_sf(sf_grid,
    "orig" = as.numeric(A_gs %*% vec1),
    "orig2" = as.numeric(A_gs %*% vec2),
    "proj" = as.numeric(A_gs %*% solve(invD, vec1)),
    "proj2" = as.numeric(A_gs %*% solve(invD, vec2))
  ) # , log(as.numeric(vec2b)))

  plot(stuff, cex = 2, pch = 19, border = NA)

  stuff_list[[i]] <- stuff |> mutate(kappa = i)

  sum_list[[i]] <- tibble(
    kappa = i,
    sum_orig = sum(stuff$orig),
    sum_proj = sum(stuff$proj),
    sum_orig2 = sum(stuff$orig2),
    sum_proj2 = sum(stuff$proj2)
  )
}

stuff <- bind_rows(stuff_list)
sums <- bind_rows(sum_list)

sums

# Reorganize data for plotting origin as a kappa-case
sums2 <- sums |>
  mutate(kappa = paste0("log(\u03BA)=", kappa)) |>
  bind_rows(sums |>
    filter(kappa == min(kappa)) |>
    dplyr::select(-sum_proj, -sum_proj2) |>
    rename(
      sum_proj = sum_orig,
      sum_proj2 = sum_orig2
    ) |>
    mutate(kappa = "Original")) |>
  dplyr::select(-sum_orig, -sum_orig2) |>
  mutate() |>
  pivot_longer(-kappa, names_to = "var", values_to = "sum") |>
  mutate(
    type = ifelse(var == "sum_proj", "point", "gmrf"),
    sum_proj = paste0("Sum=", round(sum, digits = 2))
  )

scale_values <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

stuff2 <- stuff |>
  mutate(kappa = paste0("log(\u03BA)=", kappa)) |>
  bind_rows(stuff |>
    filter(kappa == min(kappa)) |>
    dplyr::select(-proj, -proj2) |>
    rename(
      proj = orig,
      proj2 = orig2
    ) |>
    mutate(kappa = "Original")) |>
  dplyr::select(-orig, -orig2) |>
  pivot_longer(cols = c("proj", "proj2")) |>
  mutate(type = ifelse(name %in% c("orig2", "proj2"), "gmrf", "point")) |>
  mutate(
    value_sc = scale_values(value),
    .by = c(name, kappa)
  )

p <- ggplot() +
  facet_wrap(~ factor(kappa, levels = c("Original", "log(κ)=5", "log(κ)=2.5", "log(κ)=0")),
    ncol = 4
  ) +
  facet_wrap(~ factor(kappa, levels = c("Original", "log(κ)=5", "log(κ)=2.5", "log(κ)=0")),
    ncol = 4
  ) +
  scale_fill_viridis(name = "Scaled\ncovariate", option = "G") +
  scale_color_viridis(option = "G") +
  theme_void() +
  guides(color = "none") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )

t <- p +
  geom_sf(
    data = stuff2 |> filter(type == "point"),
    aes(fill = value_sc, color = value_sc)
  ) +
  geom_text(
    data = sums2 |> filter(type == "point"),
    aes(x = -Inf, y = Inf, label = sum_proj),
    hjust = -0.3, vjust = 2.5, color = "white", size = 2, fontface = "plain"
  ) +
  labs(y = "Point") +
  theme(
    axis.title.y = element_text(size = 8, angle = 90),
    strip.text.x.top = element_text(size = 10, margin = unit(rep(0.08, 4), "cm"))
  )

b <- p +
  geom_sf(
    data = stuff2 |> filter(type == "gmrf"),
    aes(fill = value_sc, color = value_sc)
  ) +
  geom_text(
    data = sums2 |> filter(type == "gmrf"),
    aes(x = -Inf, y = Inf, label = sum_proj),
    hjust = -0.3, vjust = 2.5, color = "white", size = 2, fontface = "plain"
  ) +
  labs(y = "IID standard normal distributions") +
  theme(
    axis.title.y = element_text(size = 8, angle = 90),
    strip.text.x.top = element_blank()
  )

# b +
#   scale_fill_viridis(name = "Scaled\ncovariate", option = "G") +
#   scale_color_viridis(option = "G")

(t / b) + plot_layout(guides = "collect") &
  coord_sf(xlim = c(0.05, 0.99), ylim = c(0.05, 0.99)) &
  theme(
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 8)
  )

ggsave(paste0(here::here(), "/results/figures/theory_check.pdf"),
  width = 18, height = 9, unit = "cm", device = cairo_pdf
)

################
# Time-scaling
################
run <- function(n) {
  cat(n, "\n")
  # Simulate locations
  loc <- matrix(rnorm(n * 2), ncol = 2)
  # Make SPDE objects
  mesh <- fm_mesh_2d(loc)
  spde <- fm_fem(mesh)
  invM0 <- spde$c0
  diag(invM0) <- 1 / diag(spde$c0)
  # Create inverse-D matrix
  ln_kappa <- log(5)
  ln_tau <- log(1 / (1 + exp(2 * ln_kappa)))
  invD <- exp(ln_tau) * invM0 %*% (spde$c0 + exp(2 * ln_kappa) * spde$c0 + spde$g1)
  x <- rnorm(mesh$n)
  # Do benchmark
  f1 <- function() solve(invD, x)
  f2 <- function() xx <- solve(invD) %*% x
  bench::mark(f1(), f2(), check = FALSE)
}

#size <- c(10, 100, 200, 500, 1000, 2000)
size <- c(10, 100, 200, 500, 1000, 1500, 2000)

b <- lapply(size, run)

times <- purrr::map_dfr(seq_along(b), \(i)

  data.frame(
    approach = c("Using inverse diffusion matrix directly", "Computing diffusion matrix"),
    time = as.numeric(b[[i]]$median),
    size = size[i]
))

g <- ggplot(times, aes(size, time, colour = approach)) +
  geom_line() +
  labs(x = "N", y = "Seconds per calculation", colour = "Approach") +
  scale_colour_brewer(palette = "Set2") +
  ggsidekick::theme_sleek() +
  scale_y_continuous(lim = c(-0.01, NA), expand = expansion(mult = c(0, 0.05))) +
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom")

g

ggsave(paste0(here::here(), "/results/figures/supporting/time_scaling.pdf"),
  width = 9, height = 9, unit = "cm", device = cairo_pdf)

g + scale_y_log10()
