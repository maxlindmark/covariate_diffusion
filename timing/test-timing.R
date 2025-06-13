# Timing figure -----------------------------------------------------------

# This is a trimmed copy of Sean's sdmTMB paper code for only sdmTMB. Next I will modify this to work with TMB...

library(sf)
library(fmesher)
library(Matrix)
library(TMB)
library(dplyr)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())

# Code to simulate spatially structured data based on https://github.com/seananderson/sdmTMB-paper/blob/main/analysis/timing.R

# Function to simulate data
simulate_dat <- function(n_obs = 100,
                         cutoff = 0.1,
                         range = 0.2,
                         phi = 0.1,
                         tweedie_p = 1.5,
                         max.edge = 0.1,
                         family = poisson(),
                         sigma_O = 0.2,
                         iter = 1,
                         plot = FALSE,
                         seed = sample.int(.Machine$integer.max, 1L)) {
  set.seed(seed)

  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- fmesher::fm_segm(loc.bnd)

  me <- fmesher::fm_mesh_2d_inla(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )

  loc <- me$loc[, 1:2] # mesh vertices
  predictor_dat <- data.frame(
    X = c(runif(n_obs), loc[, 1]),
    Y = c(runif(n_obs), loc[, 2]),
    a1 = rnorm(n_obs + nrow(loc))
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
    B = c(0.2, -0.5) # B0 = intercept, B1 = a1 slope
  )
  if (plot) {
    g <- ggplot(sim_dat, aes(X, Y, colour = observed)) +
      geom_point() +
      scale_color_gradient2()
    print(g)
  }
  sim_dat_knots <- sim_dat[(n_obs + 1):nrow(sim_dat), ]
  sim_dat <- sim_dat[1:n_obs, ]
  list(mesh = mesh, dat = sim_dat, dat_knots = sim_dat_knots)
}


# Simulate a simple data set and try and fit a model
s <- simulate_dat(n_obs = 5000)
hist(s$dat$observed)
head(s$dat)

library(sdmTMB)
d <- s$dat
m <- sdmTMB(observed ~ 1 + a1, data = d, family = poisson(), mesh = make_mesh(d, c("X", "Y"), mesh = s$mesh$mesh))
m$model$par
# matches

#### ML: try and fit a LNP diffusion model to data generated with above function
# root_dir <- here::here(".")
# tmb_dir <- file.path(root_dir, "tmb")
# setwd(tmb_dir)
# ML: replace these lines with the path to where movement_kernel_sim.cpp is
compile("timing/covariate_diffusion.cpp", framework = "TMBad")
dyn.load(dynlib("timing/covariate_diffusion"))

compile("timing/simple.cpp", framework = "TMBad")
dyn.load(dynlib("timing/simple"))

# compile("tmb/movement_kernel_sim.cpp", framework = "TMBad")
# dyn.load(dynlib("tmb/movement_kernel_sim"))

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
diag(invsqrtM0) <- 1 / sqrt(diag(spde$c0))
diag(invM0) <- 1 / diag(spde$c0)
depthprime_s <- s$dat_knots$a1 # This is our predictor *at the knots!*

# simple non-difussion version:

Params <- list(
  "beta0" = 0,
  "beta_j" = c(0),
  "ln_tau" = 0,
  "ln_kappa" = 0,
  "omega_s" = rep(0, nrow(spde$c0))
)
Random <- "omega_s"

get_spde_matrices <- function(x) {
  x <- x$spde[c("c0", "g1", "g2")]
  names(x) <- c("M0", "M1", "M2") # legacy INLA names needed!
  x
}
spde_tmb <- get_spde_matrices(list(spde = spde))

Data <- list(
  "dist" = distribution,
  "c_i" = as.integer(s$dat$observed),
  "depth_i_obs" = s$dat$a1,
  "M0" = spde$c0,
  "M1" = spde$g1,
  "M2" = spde$g2,
  # spde = spde_tmb,
  "A_is" = A_is
)
Obj <- MakeADFun(
  data = Data,
  parameters = Params,
  random = Random,
  checkParameterOrder = TRUE,
  DLL = "simple"
)
Obj$env$beSilent()

tictoc::tic()
Opt <- nlminb(
  start = Obj$par,
  obj = Obj$fn,
  grad = Obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)
tictoc::toc()

Opt$par
m$model$par

# original/diffusion .cpp version:

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
  "method" = "null", # "method" = "diffusion",
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
  "depth_i_obs" = s$dat$a1,
  "covariate_type" = "at observation",
  "sim_gmrf" = 0L
)

Obj <- MakeADFun(
  data = Data,
  parameters = Params,
  random = Random,
  checkParameterOrder = TRUE,
  DLL = "covariate_diffusion"
)
Obj$env$beSilent()

tictoc::tic()
Opt <- nlminb(
  start = Obj$par,
  obj = Obj$fn,
  grad = Obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)
tictoc::toc()

tictoc::tic()
sdr <- sdreport(Obj)
tictoc::toc()

Opt$par
m$model$par

Data <- list(
  "method" = "diffusion", # "method" = "diffusion",
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
  "depth_i_obs" = s$dat$a1,
  "covariate_type" = "null",
  "sim_gmrf" = 0L
)

Params <- list(
  "beta0" = 0,
  "beta_j" = -0.2,
  "ln_tau" = -2.5,
  "ln_kappa" = 2.5,
  "omega_s" = rep(0, nrow(spde$c0)),
  "ln_kappa2" = exp(-1)
)
Obj_diff <- MakeADFun(
  data = Data,
  parameters = Params,
  random = Random,
  checkParameterOrder = TRUE,
  DLL = "covariate_diffusion"
)
# Obj_diff$env$beSilent()
Obj_diff$fn()
Obj_diff$gr()

r <- Obj_diff$report()
names(r)
r$invD
r$Range2
r$ln_kappa2
r$depth_s
r$depth_i
r$depth_g
r$mu_i

tictoc::tic()
Opt <- nlminb(
  start = Obj_diff$par,
  obj = Obj_diff$fn,
  grad = Obj_diff$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)
tictoc::toc()

tictoc::tic()
sdr <- sdreport(Obj_diff)
tictoc::toc()

Opt$par
