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
  # OLD
  #invD <- exp(ln_tau) * invM0 %*% (spde$c0 + exp(2 * ln_kappa) * spde$c0 + spde$g1)
  # UPDATED
  invD = Diagonal(n=mesh$n) + exp(-2 * ln_kappa) * invM0 %*% spde$g1
  x <- rnorm(mesh$n)
  # Do benchmark
  f1 <- function() solve(invD, x)
  f2 <- function() xx <- solve(invD) %*% x
  bench::mark(f1(), f2(), check = FALSE)
}

#size <- c(10, 100, 200, 500, 1000, 2000)
size <- seq(1, 2000, by = 10)

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

g + scale_y_log10()

g

ggsave(paste0(here::here(), "/results/figures/supporting/time_scaling.pdf"),
  width = 9, height = 9, unit = "cm", device = cairo_pdf)
