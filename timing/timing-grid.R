# Code to run comparison of time in seconds to estimate parameters and to calculate standard errors for a diffused model, and a non-diffused model

library(sf)
library(fmesher)
library(Matrix)
library(TMB)
library(dplyr)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())

# Code to simulate spatially structured data based on sdmTMB paper: https://github.com/seananderson/sdmTMB-paper/blob/main/analysis/timing.R

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
set.seed(1)
s <- simulate_dat(n_obs = 5000)

compile("timing/covariate_diffusion.cpp", framework = "TMBad")
dyn.load(dynlib("timing/covariate_diffusion"))

fit_models <- function(sim_object) {
  # compile("timing/simple.cpp", framework = "TMBad")
  # dyn.load(dynlib("timing/simple"))
  s <- sim_object

  distribution <- c("Tweedie", "Poisson", "LNP")[2]
  loc <- st_multipoint(as.matrix(s$dat[, c("X", "Y")]), dim = "XY")
  spde <- fm_fem(s$mesh$mesh)
  A_is <- fm_evaluator(s$mesh$mesh, loc = st_coordinates(loc))$proj$A # ML: sf_loc -> loc

  grid_spacing <- 0.025
  grid_coords <- expand.grid(
    X = seq(0, 1, by = grid_spacing),
    Y = seq(0, 1, by = grid_spacing)
  )

  grid_coords <- as.matrix(grid_coords)
  A_gs <- fm_evaluator(s$mesh$mesh, loc = grid_coords)$proj$A
  invM0 <- invsqrtM0 <- spde$c0
  diag(invsqrtM0) <- 1 / sqrt(diag(spde$c0))
  diag(invM0) <- 1 / diag(spde$c0)
  depthprime_s <- s$dat_knots$a1 # This is our predictor *at the knots!*

  Params <- list(
    "beta0" = 0,
    "beta_j" = 0.1,
    "ln_tau" = 0,
    "ln_kappa" = 0,
    "omega_s" = rep(0, nrow(spde$c0)),
    "ln_kappa2" = exp(-1)
  )
  Random <- "omega_s"
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

  times <- list()

  out <- system.time({
    Opt <- nlminb(
      start = Obj$par,
      obj = Obj$fn,
      grad = Obj$gr,
      control = list(eval.max = 1e4, iter.max = 1e4, trace = 0)
    )
  })
  times$standard_nlminb <- out[["elapsed"]]

  out <- system.time({
    sdr <- sdreport(Obj)
  })
  times$standard_sdreport <- out[["elapsed"]]

  Data$method <- "diffusion"
  Data$covariate_type <- "null"

  Obj_diff <- MakeADFun(
    data = Data,
    parameters = Params,
    random = Random,
    checkParameterOrder = TRUE,
    DLL = "covariate_diffusion"
  )
  Obj_diff$env$beSilent()

  out <- system.time({
    Opt_diff <- nlminb(
      start = Obj_diff$par,
      obj = Obj_diff$fn,
      grad = Obj_diff$gr,
      control = list(eval.max = 1e4, iter.max = 1e4, trace = 0)
    )
  })
  times$diffused_nlminb <- out[["elapsed"]]

  out <- system.time({
    sdr <- sdreport(Obj)
  })
  times$diffused_sdreport <- out[["elapsed"]]

  times$n <- nrow(s$dat)
  times$knots <- nrow(s$dat_knots)
  as.data.frame(times)
}

# fit_models(sim_object = s)

to_run <- expand.grid(
  n_obs = c(1e3, 1e4, 1e5),
  max.edge = c(0.05, 0.1, 0.15, 0.2),
  iter = c(1:50)
)
to_run$seed <- to_run$iter * 29212
nrow(to_run)

# don't run in parallel to ensure all timing is comparable:
system.time({
  sim_out <- purrr::pmap(
    to_run,
    simulate_dat,
    .progress = "timing"
  )
})

system.time({
  fit_out <- purrr::map_dfr(
    sim_out,
    fit_models,
    .progress = "timing"
  )
})
fit_out$seed <- to_run$seed

fit_out |>
  tidyr::pivot_longer(cols = c(-n, -knots, -seed)) |>
  summarise(median = median(value),
            upr = quantile(value, probs = 0.95),
            lwr = quantile(value, probs = 0.05),
            .by = c(n, knots, name)) |>
  mutate(Diffused = grepl("diffused", name)) |>
  mutate(type = ifelse(grepl("nlminb()", name), "nlminb", "sdreport()")) |>
  mutate(n_text = paste0("n = ", n)) |>
  ggplot(aes(knots, median, colour = Diffused)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Diffused), alpha = 0.3, color = NA) +
  geom_line() +
  facet_grid(type ~ n_text, scales = "free_y") +
  ylab("Time (s)") +
  xlab("Mesh vertices") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")

ggsave(paste0(here::here(), "/results/figures/supporting/timing.pdf"),
      width = 16, height = 9, unit = "cm")

# Relative plot
fit_out |>
  tidyr::pivot_longer(cols = c(-n, -knots, -seed)) |>
  summarise(median = median(value),
            upr = quantile(value, probs = 0.95),
            lwr = quantile(value, probs = 0.05),
            .by = c(n, knots, name)) |>
  dplyr::select(-upr, -lwr) |>
  tidyr::pivot_wider(names_from = name, values_from = median) |>
  mutate(ratio_nlminb = diffused_nlminb / standard_nlminb,
         ratio_sdreport = diffused_sdreport / standard_sdreport) |>
  tidyr::pivot_longer(c(ratio_nlminb, ratio_sdreport)) |>
  mutate(mean = mean(value), .by = c(name, n)) |>
  mutate(type = ifelse(grepl("nlminb()", name), "nlminb", "sdreport()")) |>
  mutate(n_text = paste0("n = ", n)) |>
  ggplot(aes(knots, value)) +
  geom_hline(aes(yintercept = mean, color = "mean")) +
  geom_line() +
  facet_grid(type ~ n_text, scales = "free_y") +
  ylab("Time (s)") +
  xlab("Mesh vertices") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
