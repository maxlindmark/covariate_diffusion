root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

# Custom function for adding legend in multi-panel figure
# https://stackoverflow.com/questions/52975447/reorganize-sf-multi-plot-and-add-a-legend
source(file.path(root_dir, "functions/add-legend.R"))
source(here::here("analysis/prep-bb-data.R"))

# Compile
setwd(tmb_dir)
compile("movement_kernel_bb.cpp", framework = "TMBad")
dyn.load(dynlib("movement_kernel_bb"))

# saving stuff
species_set <- colnames(sf_DF)[3:22]

sf_DF <- as.data.frame(sf_DF)

#
param_set <- c("obj_diffusion", "obj_null", "deltaAIC", "range", "ln_kappa2", "corr_pop_dens")
Results_cz <- array(NA,
  dim = c(length(species_set), length(param_set)),
  dimnames = list(species_set, param_set)
)

# Store output for plotting
stuff_gz_list <- list()

# Loop
for (cI in seq_along(species_set)) {
  species <- species_set[cI]

  # Build object
  pop_dens_prime_s <- scale(pop_dens_s)

  Data <- list(
    "method" = "diffusion",
    "dist" = distribution,
    "c_i" = sf_DF[, species],
    "A_is" = A_is,
    "A_gs" = A_gs,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2,
    "invsqrtM0" = invsqrtM0,
    "invM0" = invM0,
    "pop_dens_s" = pop_dens_prime_s,
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
    Params$eta_i <- rnorm(nrow(sf_DF))
    Random <- c(Random, "eta_i")
  }

  Obj <- MakeADFun(
    data = Data,
    parameters = Params,
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
  Report <- Obj$report()

  # Stationary distribution
  Results_cz[cI, "obj_diffusion"] <- Opt$objective
  Results_cz[cI, "range"] <- Report$Range2
  Results_cz[cI, "ln_kappa2"] <- Opt$par["ln_kappa2"]
  Results_cz[cI, "corr_pop_dens"] <- cor(Report$pop_dens_g, as.vector(A_gs %*% Data$pop_dens_s))

  # Re-run without covariate diffusion
  Data2 <- Data
  Data2$method <- "null"
  Params2 <- Params[setdiff(names(Params), c("ln_tau2", "ln_kappa2"))]
  Obj2 <- MakeADFun(
    data = Data2,
    parameters = Params2,
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

  #
  Results_cz[cI, "deltaAIC"] <- (2 * Opt2$objective + 2 * length(Opt2$par)) - (2 * Opt$objective + 2 * length(Opt$par))
  Results_cz[cI, "obj_null"] <- Opt2$objective
  write.csv(Results_cz, file = file.path(date_dir, "Results_bb_cz.csv"))

  # Plot results
  png(file = file.path(date_dir, paste0(species, ".png")), width = 6, height = 6, res = 200, units = "in")
  stuff_gz <- st_sf(sf_grid,
    "orig" = as.vector(Data$A_gs %*% Data$pop_dens_s),
    "log count orig" = log(Report2$mu_g),
    "ppop_count orig" = Report2$ppop_dens_g,
    "omega orig" = Report2$omega_g,
    "diffused" = Report$pop_dens_g,
    "log count diffused" = log(Report$mu_g),
    "ppop_count diffused" = Report$ppop_dens_g,
    "omega diffused" = Report$omega_g
  )
  par(mfcol = c(4, 2))
  for (i in 1:8) {
    plot(stuff_gz[i], border = NA, key.pos = NULL, reset = FALSE)
    add_legend(round(range(stuff_gz[[i]], na.rm = TRUE), 1))
  }
  dev.off()

  stuff_gz_list[[cI]] <- stuff_gz |> dplyr::mutate(species = species_set[cI])
}

stuff_gz_df <- dplyr::bind_rows(stuff_gz_list)

saveRDS(stuff_gz_df, file = file.path(date_dir, "stuff_gz_df.rds"))
