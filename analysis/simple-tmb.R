source(here::here("analysis/prep-data.R"))

setwd(tmb_dir)
compile("movement_kernel.cpp", framework = "TMBad")
dyn.load(dynlib("movement_kernel"))
# dyn.unload(dynlib("movement_kernel"))

species_set <- colnames(region_data_all)[30:ncol(region_data_all)]
N_c <- colSums(ifelse(region_data_all[, species_set] > 0, 1, 0))
species_set <- species_set[N_c > 1000]

param_set <- c("obj_diffusion", "obj_null", "deltaAIC", "range", "ln_kappa2", "corr_depth")

cI <- match("a_yfs", species_set)
species <- species_set[cI]

if (f_depth == "identity") {
  depthprime_s <- (depth_s / 100)
} else if (f_depth == "log") {
  depthprime_s <- log(depth_s) - 4.5
}

Data <- list(
  "method" = "diffusion",
  "dist" = distribution,
  "c_i" = region_data_all[, species],
  "A_is" = A_is,
  "A_gs" = A_gs,
  "M0" = spde$c0,
  "M1" = spde$g1,
  "M2" = spde$g2,
  "invsqrtM0" = invsqrtM0,
  "invM0" = invM0,
  "depth_s" = depthprime_s
)

Params <- list(
  "beta0" = 0,
  "beta_j" = c(0.1, 0.1),
  "ln_tau" = 0,
  "ln_kappa" = 0,
  "omega_s" = rnorm(nrow(spde$c0)),
  "ln_kappa2" = exp(-1)
)

Random <- "omega_s"

# Special stuff
if (Data$dist == "Tweedie") {
  Params$ln_phi <- log(2)
  Params$finv_power <- 0
}
if (Data$dist == "LNP") {
  Params$ln_sigma_eta <- log(0.1)
  Params$eta_i <- rnorm(nrow(region_data_all))
  Random <- c(Random, "eta_i")
}

Obj <- MakeADFun(
  data = Data,
  parameters = Params,
  # profile = c("beta0","beta_j"),    # ,"ln_kappa2"
  random = Random,
  checkParameterOrder = TRUE
)
Obj$env$beSilent()

# Optimize
Opt <- nlminb(
  start = Obj$par,
  obj = Obj$fn,
  grad = Obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)

# sdr <- sdreport(Obj)

setwd(here::here())
