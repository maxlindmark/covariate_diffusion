# Timing figure -----------------------------------------------------------

# This is a trimmed copy of Sean's sdmTMB paper code for only sdmTMB. Next I will modify this to work with TMB...

library(sf)
library(INLA)
library(fmesher)
library(Matrix)
library(TMB)
library(dplyr)
library(ggplot2)
library(ggsidekick); theme_set(theme_sleek())

# Code to simulate spatially structured data based on https://github.com/seananderson/sdmTMB-paper/blob/main/analysis/timing.R

# Function to simulate data
simulate_dat <- function(n_obs = 100,
                         cutoff = 0.1,
                         range = 0.5,
                         phi = 0.1,
                         tweedie_p = 1.5,
                         max.edge = 0.1,
                         family = poisson(),
                         sigma_O = 0.2,
                         iter = 1,
                         plot = FALSE,
                         seed = sample.int(.Machine$integer.max, 1L)) {
  set.seed(seed)
  predictor_dat <- data.frame(
    X = runif(n_obs), Y = runif(n_obs),
    a1 = rnorm(n_obs)
  )

  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- inla.mesh.segment(loc.bnd)

  me <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )

  mesh <- sdmTMB::make_mesh(predictor_dat, xy_cols = c("X", "Y"), mesh = me)

  sim_dat <- sdmTMB::sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    mesh = mesh,
    family = family,
    range = range,
    phi = phi,
    tweedie_p = tweedie_p,
    sigma_O = sigma_O,
    seed = seed,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  if (plot) {
    g <- ggplot(sim_dat, aes(X, Y, colour = observed)) +
      geom_point() +
      scale_color_gradient2()
    print(g)
  }
  list(mesh = mesh, dat = sim_dat)
}

# Simulate a simple data set and try and fit a model
s <- simulate_dat()
hist(s$dat$observed)
head(s$dat)

#### ML: try and fit a LNP diffusion model to data generated with above function
# root_dir <- here::here(".")
# tmb_dir <- file.path(root_dir, "tmb")
# setwd(tmb_dir)
# ML: replace these lines with the path to where movement_kernel_sim.cpp is
compile("movement_kernel_sim.cpp", framework = "TMBad")
dyn.load(dynlib("movement_kernel_sim"))

# Prepare data for MakeADFun
distribution <- c("Tweedie", "Poisson", "LNP")[2]
loc <- st_multipoint(as.matrix(s$dat[, c("X", "Y")]), dim = "XY")
spde <- fm_fem(s$mesh$mesh)
A_is <- fm_evaluator(s$mesh$mesh, loc = st_coordinates(loc))$proj$A # ML: sf_loc -> loc

# Set up a grid... This differs from the previous simulation which was based on the EBS domain shape
grid_spacing <- 0.025
grid_coords <- expand.grid(
  X = seq(0, 1, by = grid_spacing),
  Y = seq(0, 1, by = grid_spacing)
)
# grid_df <- as.data.frame(grid_coords)
# names(grid_df) <- c("X", "Y")
# ggplot(grid_df, aes(X, Y)) +
#   geom_point(alpha = 0.5) +
#   geom_point(data = s$dat, aes(X, Y, color = "tomato"), size = 2) +
#   coord_fixed()

grid_coords <- as.matrix(grid_coords)
A_gs <- fm_evaluator(s$mesh$mesh, loc = grid_coords)$proj$A
invM0 <- invsqrtM0 <- spde$c0
depthprime_s <- s$dat$a1 # This is our predictor

Params <- list(
  "beta0" = 0,
  "beta_j" = 0.1,
  "ln_tau" = 0,
  "ln_kappa" = 0,
  "omega_s" = rep(0, nrow(spde$c0)),
  "ln_kappa2" = exp(-1)
)
Params

Random <- "omega_s"

# This is for LNP
# Params$ln_sigma_eta <- log(0.1)
# Params$eta_i <- rnorm(nrow(s$dat))
# Random <- c(Random, "eta_i")

Data <- list(
  "method" = "diffusion", #"method" = "diffusion",
  "dist" = distribution,
  "c_i" = as.integer(s$dat$observed),
  "A_is" = A_is,
  "A_gs" = A_gs,
  "M0" = spde$c0,
  "M1" = spde$g1,
  "M2" = spde$g2,
  "invsqrtM0" = invsqrtM0,
  "invM0" = invM0,
  "depth_s" = depthprime_s,
  "sim_gmrf" = 1L
)

# # This crashes!
# Obj_diff <- MakeADFun(
#   data = Data,
#   parameters = Params,
#   random = Random,
#   checkParameterOrder = TRUE,
#   DLL = "movement_kernel_sim"
# )

