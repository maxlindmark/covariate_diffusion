library(sf)
library(Matrix)
library(TMB)
library(rnaturalearth)
library(terra)
library(fmesher)

region <- c("EBS", "GOA")[1]

root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

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
date_dir <- file.path(root_dir, "results", paste0(Date, "_", f_depth, "_", distribution))
dir.create(date_dir, recursive = TRUE, showWarnings = FALSE)
sf_usa <- ne_countries(country = "united states of america", return = "sf")
