root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

# Custom function for adding legend in multi-panel figure
# https://stackoverflow.com/questions/52975447/reorganize-sf-multi-plot-and-add-a-legend
source(file.path(root_dir, "functions/add-legend.R"))
source(here::here("analysis/prep-ebs-data.R"))
# source conditional AIC function
source(here::here("analysis/mod-cAIC.R"))

# Compile
setwd(tmb_dir)
compile("movement_kernel_ebs.cpp", framework = "TMBad")
dyn.load(dynlib("movement_kernel_ebs"))

# saving stuff
species_set <- colnames(region_data_all)[30:ncol(region_data_all)]
N_c <- colSums(ifelse(region_data_all[, species_set] > 0, 1, 0))
species_set <- species_set[N_c > 1000]

param_set <- c("obj_diffusion", "obj_null", "deltaCAIC", "deltaAIC", "range", "ln_kappa2", "corr_depth")
#param_set <- c("obj_diffusion", "obj_null", "deltaAIC", "range", "ln_kappa2", "corr_depth")
Results_cz <- array(NA,
  dim = c(length(species_set), length(param_set)),
  dimnames = list(species_set, param_set)
)

# Store output for plotting
stuff_gz_list <- list()

# For scaling and scaling back covariate
mean_depth_s <- mean(depth_s)
sd_depth_s <- sd(depth_s)

# Loop
# cI <- 26
# This species has a non positive definite Hessian, so it throws an error when doing sdreport!
species_set <- species_set[-26]

for (cI in seq_along(species_set)) {
  species <- species_set[cI]

  # Build object
  if (f_depth == "identity") {
    depthprime_s <- (depth_s - mean_depth_s) / sd_depth_s # (depth_s / 100)
  } else if (f_depth == "log") {
    depthprime_s <- log(depth_s) - 4.5
  }
  Data <- list(
    "method" = "diffusion",
    "dist" = distribution,
    "c_i" = region_data_all[, species],
    "weights_i" = rep(1, times = length(region_data_all[, species])),
    "A_is" = A_is,
    "A_gs" = A_gs,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2,
    "invsqrtM0" = invsqrtM0,
    "invM0" = invM0,
    "depth_s" = depthprime_s,
    "sim_gmrf" = 0L
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

  # profile = c("beta0","beta_j") is unstable for j_berfl and depthprime_s = (depth_s / 100)
  Obj <- MakeADFun(
    data = Data,
    parameters = Params,
    # profile = c("beta0","beta_j"),    # ,"ln_kappa2"
    random = Random,
    checkParameterOrder = TRUE
  )
  Obj$env$beSilent()

  if (species %in% c("lhdab", "a_yfs", "j_berfl")) {
    Obj$par["ln_kappa2"] <- 4
  }

  # Optimize
  Opt <- nlminb(
    start = Obj$par,
    obj = Obj$fn,
    grad = Obj$gr,
    control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
  )
  Report <- Obj$report()

  # Stationary distribution
  Results_cz[cI, "obj_diffusion"] <- Opt$objective
  Results_cz[cI, "range"] <- Report$Range2
  Results_cz[cI, "ln_kappa2"] <- Opt$par["ln_kappa2"]
  Results_cz[cI, "corr_depth"] <- cor(Report$depth_g, as.vector(A_gs %*% Data$depth_s))

  # Re-run without covariate diffusion
  Data2 <- Data
  Data2$method <- "null"
  Params2 <- Params[setdiff(names(Params), c("ln_tau2", "ln_kappa2"))]
  Obj2 <- MakeADFun(
    data = Data2,
    parameters = Params2,
    # profile = c("beta0","beta_j"),
    # profile = c("beta_j"),
    random = Random,
    checkParameterOrder = TRUE
  )
  Obj2$env$beSilent()

  # Optimize
  Opt2 <- nlminb(
    start = Obj2$par,
    obj = Obj2$fn,
    grad = Obj2$gr,
    control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
  )
  Report2 <- Obj2$report()

  # Calculate delta conditional AIC
  diff_cAIC <- cAIC.TMB(
    obj = Obj,
    tmb_data = Data,
    parlist = Params,
    random = Random,
    p = length(setdiff(names(Params), Random)),
    what = "cAIC"  # or "EDF"
  )

  null_cAIC <- cAIC.TMB(
    obj = Obj2,
    tmb_data = Data2,
    parlist = Params2,
    random = Random,
    p = length(setdiff(names(Params2), Random)),
    what = "cAIC"  # or "EDF"
  )

  Results_cz[cI, "deltaCAIC"] <- null_cAIC - diff_cAIC
  Results_cz[cI, "deltaAIC"] <- (2 * Opt2$objective + 2 * length(Opt2$par)) - (2 * Opt$objective + 2 * length(Opt$par))
  Results_cz[cI, "obj_null"] <- Opt2$objective
  write.csv(Results_cz, file = file.path(date_dir, "Results_cz.csv"))

  # Plot results
  trunc_dens <- function(logdens, ratio = 0.01) ifelse(logdens < (max(logdens, na.rm = TRUE) + log(ratio)), NA, logdens)
  png(file = file.path(date_dir, paste0(species, ".png")), width = 6, height = 6, res = 200, units = "in")
  stuff_gz <- st_sf(sf_grid,
    "orig" = as.vector(Data$A_gs %*% Data$depth_s),
    "log(dens) orig" = trunc_dens(log(Report2$mu_g)),
    "pdepth orig" = trunc_dens(Report2$pdepth_g),
    "omega orig" = Report2$omega_g,
    "diffused" = Report$depth_g,
    "log(dens) diffused" = trunc_dens(log(Report$mu_g)),
    "pdepth diffused" = trunc_dens(Report$pdepth_g),
    "omega diffused" = Report$omega_g
  )
  # plot(stuff_gz, border=NA )
  par(mfcol = c(4, 2))
  for (i in 1:8) {
    plot(stuff_gz[i], border = NA, key.pos = NULL, reset = FALSE)
    plot(sf_usa, add = TRUE, col = "grey")
    add_legend(round(range(stuff_gz[[i]], na.rm = TRUE), 1))
  }
  dev.off()

  stuff_gz_list[[cI]] <- stuff_gz |>
    dplyr::mutate(
      species = species_set[cI]
      # Scale back covariate for plotting? How about the diffused covariate?
      # ,
      # orig = (orig * sd_depth_s) + mean_depth_s,
      # diffused = (diffused * sd_depth_s) + mean_depth_s
    )
}

stuff_gz_df <- dplyr::bind_rows(stuff_gz_list)

saveRDS(stuff_gz_df, file = file.path(date_dir, "stuff_gz_df.rds"))
