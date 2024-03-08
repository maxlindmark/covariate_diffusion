library(RTMB)
library(Matrix)
source(here::here("analysis/prep-data.R"))

setwd(tmb_dir)

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

# -------------------------------------------------------------------


nll <- function(par) {
  getAll(par, Data)
  jnll <- 0
  Q <- (exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2) * exp(2 * ln_tau)
  jnll <- jnll - dgmrf(omega_s, 0, Q)

  # Probability of random effects
  if (method == "diffusion") {
    # FIXME: RTMB error on MakeADFun:
    # Error in SparseArith2(e1, e2, .Generic) : '*' not implemented
    # from the `invM0 * ...` part
    invD <- invM0 * (M0 + exp(2 * ln_kappa2) * M0 + M1) / (1 + exp(2 * ln_kappa2))
    invD <- Matrix::lu(invD)
    depth_sz <- Matrix::qr(as(depth_s, "sparseMatrix")) # FIXME looking for equivalent of lu.solve( depth_s.matrix()
    # depth_s <- depth_sz.col(0) FIXME...
    # depth_s <- depth_sz
    Range2 <- sqrt(8) / exp(ln_kappa2)
  }
  depth_i <- A_is %*% depth_s
  depth_g <- A_gs %*% depth_s
  omega_i <- A_is %*% omega_s
  omega_g <- A_gs %*% omega_s
  pdepth_i <- depth_i * beta_j[1] + depth_i^2 * beta_j[2]
  pdepth_g <- depth_g * beta_j[1] + depth_g^2 * beta_j[2]

  # Probability of data conditional on random effects
  mu_i <- as.vector(exp(beta0 + omega_i + pdepth_i))
  if (dist == "Poisson") {
    jnll <- jnll - sum(dpois(c_i, mu_i))
  } else if (dist == "Tweedie") {
    jnll <- jnll - sum(dtweedie(c_i, mu_i, exp(ln_phi), 1 + plogis(finv_power)))
  } else if (dist == "LNP") {
    jnll <- jnll - sum(dnorm(eta_i, 0, exp(ln_sigma_eta)))
    jnll <- jnll - sum(dpois(c_i, mu_i * exp(eta_i)))
  } else {
    stop("`dist` invalid")
  }
  mu_g <- as.vector(exp(beta0 + omega_g + pdepth_g))

  REPORT(invD)
  REPORT(Range2)
  REPORT(ln_kappa2)
  REPORT(depth_s)
  REPORT(depth_i)
  REPORT(depth_g)
  REPORT(Q)
  REPORT(omega_s)
  REPORT(mu_i)
  REPORT(mu_g)
  REPORT(pdepth_i)
  REPORT(pdepth_g)

  jnll
}

nll(Params)

Obj <- RTMB::MakeADFun(
  nll,
  parameters = Params,
  # profile = c("beta0","beta_j"),    # ,"ln_kappa2"
  random = Random
)
# Obj$env$beSilent()
