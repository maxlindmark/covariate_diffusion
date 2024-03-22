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
  "depth_s" = depthprime_s,
  "sim_gmrf" = 0L #< new
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
  checkParameterOrder = TRUE,
  DLL = "movement_kernel"
)
Obj$env$beSilent()

# Optimize
Opt <- nlminb(
  start = Obj$par,
  obj = Obj$fn,
  grad = Obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)
Opt

r <- Obj$report()

# sdr <- sdreport(Obj)

# Simulation demo:

set.seed(1)
s <- Obj$simulate()

plot(log(s$c_i+1), log(Data$c_i+1))
plot(s$omega_s, r$omega_s) # not simulated

# now demonstrate also simulating the GMRF:
Data$sim_gmrf <- 1L

obj2 <- MakeADFun(
  data = Data,
  random = Obj$env$random,
  parameters = Obj$env$parList(),
  DLL = "movement_kernel"
)

s2 <- obj2$simulate(complete = FALSE)
plot(log(s2$c_i+1), log(Data$c_i+1))
plot(s2$omega_s, r$omega_s) # simulated

# demonstrate changing some parameters:
p <- obj2$env$last.par.best
table(names(p))
p[names(p) == "beta0"] <- 2
p[names(p) == "beta_j"] <- c(-1, 0.3)
p[names(p) == "ln_sigma_eta"] <- 0.3
p[names(p) == "ln_kappa"] <- -0.2
p[names(p) == "ln_kappa2"] <- 1.1
# ignore omega_s, it's getting simulated because Data$sim_gmrf is 1

s3 <- obj2$simulate(par = p, complete = FALSE)
plot(s3$omega_s, s2$omega_s)

# now fit back to our simulated data:

Data_sim <- Data
Data_sim$c_i <- s3$c_i

Obj <- MakeADFun(
  data = Data_sim,
  parameters = Params,
  random = Random,
  checkParameterOrder = TRUE,
  DLL = "movement_kernel"
)
Opt <- nlminb(
  start = Obj$par,
  obj = Obj$fn,
  grad = Obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4)
)
Opt
r <- Obj$report()
# sdr <- sdreport(Obj)
# pl <- as.list(sdr, "Estimate")
