# Mesh sensitivity - default cutoff
# Note for fixed effects (not eta_i and omega_s) we used optimized parameters from HIGH mesh resolution (LOW CUTOFF)
library(sf)
library(Matrix)
library(TMB)
library(rnaturalearth)
library(terra)
library(fmesher)
library(ggsidekick)
library(ggplot2)


region <- c("EBS", "GOA")[1]

root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

# We'll fit this lower res model with optimized pars
optimized_pars <-
  #read.csv(paste0(root_dir, "/results/2025-06-27_identity_LNP_low_cutoff/optimized_pars.csv")) |>
  read.csv(paste0(root_dir, "/results/2025-07-16_identity_LNP_low_cutoff/optimized_pars.csv")) |>
  dplyr::rename(species_par = species)

region_data_all <- read.csv(file.path(data_dir, "all_EBS_data_2021.csv"))
shape1 <- st_read(file.path(data_dir, "EBSThorson")) # , promote_to_multi=FALSE
shape2 <- st_read(file.path(data_dir, "NBSThorson")) # , promote_to_multi=FALSE
domain_shape <- c(st_geometry(shape1), st_geometry(shape2))

bdepth <- rast(file.path(data_dir, "bdepth.tif"))

# Make domain from shapefile
# See: https://github.com/mstrimas/smoothr/blob/main/R/fill-holes.r to fill holes
domain_shape <- st_geometry(domain_shape)
domain_shape <- st_cast(domain_shape, to = "POLYGON")
if (region == "GOA") domain_shape <- domain_shape[which.max(st_area(domain_shape))]
domain_shape <- st_transform(domain_shape, crs = 4326)
domain_shape <- st_make_valid(domain_shape)
plot(domain_shape)

loc <- st_multipoint(as.matrix(region_data_all[, c("lon", "lat")]), dim = "XY")
sf_loc <- st_sfc(loc, crs = st_crs(bdepth))
sf_loc <- st_transform(sf_loc, crs = 4326)

# domain from data
domain <- st_concave_hull(sf_loc, ratio = 0.1) # 0.1 works for both GOA and EBS
domain <- st_sfc(domain, crs = st_crs(sf_loc))

# domain from mesh
sf_grid <- st_make_grid(domain_shape, cellsize = c(0.1, 0.1))
sf_grid <- st_intersection(sf_grid, domain_shape)
grid_loc <- st_coordinates(st_centroid(sf_grid))
# plot(sf_grid)
mesh <- fmesher::fm_mesh_2d(st_coordinates(sf_loc)[, 1:2], cutoff = 0.1)
mesh$n

# ggplot() +
#   geom_fm(data = mesh, fill = NA) +
#   theme_sleek() +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(aspect.ratio = 1)

# Other objects
spde <- fm_fem(mesh)
A_is <- fm_evaluator(mesh, loc = st_coordinates(sf_loc))$proj$A
A_gs <- fm_evaluator(mesh, loc = grid_loc)$proj$A
invM0 <- invsqrtM0 <- spde$c0
diag(invsqrtM0) <- 1 / sqrt(diag(spde$c0))
diag(invM0) <- 1 / diag(spde$c0)

# Get depth at vertices of SPDE mesh
mesh_points <- sf_project(mesh$loc[, 1:2], from = st_crs(4326), to = st_crs(bdepth))
depth_s <- terra::extract(bdepth, mesh_points)[, 1]

# FIll in missing covariates
# FIXME:  Could be done using 0 for land and positive values otherwise
depth_s <- ifelse(is.na(depth_s), mean(depth_s, na.rm = TRUE), depth_s)

f_depth <- c("identity", "log")[1]
distribution <- c("Tweedie", "Poisson", "LNP")[3] # Poisson is numerically unstable for j_nrs

Date <- Sys.Date()
date_dir <- file.path(root_dir, "results", paste0(Date, "_", f_depth, "_", distribution, "_", "main_cutoff"))
dir.create(date_dir, recursive = TRUE, showWarnings = FALSE)
sf_usa <- ne_countries(country = "united states of america", return = "sf")




###### Run case studies and save data
root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

# Compile
setwd(tmb_dir)
compile("movement_kernel_ebs_uw.cpp", framework = "TMBad")
dyn.load(dynlib("movement_kernel_ebs_uw"))

# saving stuff
species_set <- colnames(region_data_all)[30:ncol(region_data_all)]
N_c <- colSums(ifelse(region_data_all[, species_set] > 0, 1, 0))
species_set <- species_set[N_c > 1000]

#
param_set <- c("obj_diffusion", "obj_null", "deltaAIC", "range", "ln_kappa2", "corr_depth")
Results_cz <- array(NA,
                    dim = c(length(species_set), length(param_set)),
                    dimnames = list(species_set, param_set)
)

# For scaling and scaling back covariate
mean_depth_s <- mean(depth_s)
sd_depth_s <- sd(depth_s)

# Loop
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
  # replace with optimized params (note not omega_s!)
  Params$beta0 <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "diffusion", "beta0"]
  Params$beta_j <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "diffusion", c("beta_j_2", "beta_j_1")] |> as.numeric()
  Params$ln_tau <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "diffusion", "ln_tau"]
  Params$ln_kappa <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "diffusion", "ln_kappa"]
  Params$ln_kappa2 <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "diffusion", "ln_kappa2"]

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

  if (species %in% c("lhdab", "a_yfs", "j_berfl", "a_bigmouth", "lhdab")) {
    Obj$par["ln_kappa2"] <- 10
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

  # replace with optimized params (note not omega_s!)
  Params2$beta0 <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "null", "beta0"]
  Params2$beta_j <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "null", c("beta_j_2", "beta_j_1")] |> as.numeric()
  Params2$ln_tau <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "null", "ln_tau"]
  Params2$ln_kappa <- optimized_pars[optimized_pars$species_par == species & optimized_pars$model == "null", "ln_kappa"]

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

  #
  Results_cz[cI, "deltaAIC"] <- (2 * Opt2$objective + 2 * length(Opt2$par)) - (2 * Opt$objective + 2 * length(Opt$par))
  Results_cz[cI, "obj_null"] <- Opt2$objective
  write.csv(Results_cz, file = file.path(date_dir, "main_mesh_Results_cz.csv"))

}



