# SA: quick untested conversion

library(RTMB)

nll <- function(par) {
  getAll(par, tmb_data)

  # Derived quantities
  Q <- (exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2) *
    exp(2 * ln_tau)
  omega_s %~% dgmrf(0, Q)

  # Probability of random effects
  if (method == "diffusion") {
    invKd <- exp(ln_tau2) * invsqrtM0 * (exp(2 * ln_kappa2) * M0 + M1)
    invD <- invM0 * (M0 + exp(2 * ln_kappa2) * M0 + M1) /
      (1 + exp(2 * ln_kappa2))
    invD <- Matrix::lu(invD) # CHECK
    REPORT(invD)

    depth_sz <- Matrix::solve(depth_s)

    # Correct for proportionality constant
    depth_s <- depth_sz[, 1] # .col(0);

    Range2 <- sqrt(8) / exp(ln_kappa2)
    REPORT(Range2)
    REPORT(ln_kappa2)
  }
  depth_i <- A_is * depth_s # FIXME %*%?
  depth_g <- A_gs * depth_s # FIXME %*%?
  REPORT(depth_s)
  REPORT(depth_i)
  REPORT(depth_g)
  omega_i <- A_is * omega_s # FIXME %*%?
  omega_g <- A_gs * omega_s # FIXME %*%?
  pdepth_i <- depth_i * beta_j[1] + depth_i^2 * beta_j[2] # FIXME %*%?
  pdepth_g <- depth_g * beta_j[1] + depth_g^2 * beta_j[2] # FIXME %*%?

  # Probability of data conditional on random effects
  mu_i <- exp(beta0 + omega_i + pdepth_i)
  if (dist == "Poisson") {
    c_i %~% dpois(mu_i)
  }
  if (dist == "Tweedie") {
    c_i %~% dtweedie(mu_i, exp(ln_phi), 1 + plogis(finv_power))
  }
  if (dist == "LNP") {
    eta_i %~% dnorm(0, exp(ln_sigma_eta))
    c_i %~% dpois(c_i, mu_i * exp(eta_i))
  }
  mu_g <- exp(beta0 + omega_g + pdepth_g)

  REPORT(Q)
  REPORT(omega_s)
  REPORT(mu_i)
  REPORT(mu_g)
  REPORT(pdepth_i)
  REPORT(pdepth_g)
}
