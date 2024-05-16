# Max Lindmark 2024.04.23 - Modified version of this script: https://github.com/James-Thorson/Spatio-temporal-models-for-ecologists/blob/main/Chap_11/BBS.R

library(sf)
library(TMB)
library(ape)
library(rnaturalearth)
library(elevatr)
library(viridisLite)
library(raster)
library(terra)
library(ggplot2)
library(fmesher)
library(tidyterra) # to make pop_dens a spatraster so that I can extract values at mesh nodes
library(tidyr)
library(Matrix)


root_dir <- here::here(".")
data_dir <- file.path(root_dir, "data")
tmb_dir <- file.path(root_dir, "tmb")

# Load data
DF = read.csv(file.path(data_dir, "BB/Top20_Samples.csv"))
#trait_set = read.csv( "Top20_traits.csv" )

# Pivot wider so that 1 species = 1 column
DF <- DF %>%
  pivot_wider(names_from = Genus_species,
              values_from = SpeciesTotal,
              values_fill = 0)

# Load population density
pop_dens = st_read( file.path(data_dir, "BB/population_density.csv"), options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs=st_crs("+proj=longlat +datum=WGS84") )
pop_dens$Dens2020 = as.numeric(pop_dens$Dens2020)

# Load copNDVI (saved from rasterdiv)
copNDVI = raster( file.path(data_dir, "BB/NDVI.tif") )
#library(rasterdiv)

# Get spatial domain
sf_states = ne_states( c("United States of America"), return="sf")
sf_states = sf_states[pmatch(c("Cal", "Oregon", "Washington", "Idaho", "Montana", "Utah", "New Mex", "Arizona", "Wyoming", "Colorad", "Nevada"), sf_states$name_en),]
sf_states = st_union(sf_states)

# Create data-frame
sf_DF = st_as_sf( DF, coords=c("Longitude","Latitude"), crs="+proj=longlat +datum=WGS84")

# TODO: grid or continuous locations?
sf_fullgrid = st_make_grid( sf_DF, cellsize=1, square=FALSE )
sf_grid = st_make_valid(st_intersection( sf_fullgrid, sf_states ))
sf_grid = sf_grid[ st_area(sf_grid)>(0.01*max(st_area(sf_grid))) ]    # or 0.01

# make data frame of covariates
df_grid = st_centroid(sf_grid)
# st_as_sf to convert from sfc point type, else error
df_grid = get_elev_point( st_as_sf(df_grid), src = "aws" )
df_grid$log_elevation_km = log( ifelse(df_grid$elevation<1, 1, df_grid$elevation) / 1000 )
df_grid$NDVI = raster::extract( x=copNDVI, y=as(df_grid,"Spatial") )
df_grid$scale_NDVI = scale( df_grid$NDVI )[,1]
df_grid$pop_dens = pop_dens$Dens2020[ st_nearest_feature( sf_grid, pop_dens ) ]
df_grid$log_pop_dens = log(df_grid$pop_dens)
df_grid = data.frame(df_grid)


sf_DF = st_intersection( sf_DF, st_union(sf_grid) )
# Since I made the data wide to better match the Bering Sea case study script I skip this
#sf_DF$Genus_species = factor(sf_DF$Genus_species)
temp_DF = get_elev_point( sf_DF, src = "aws" )
sf_DF$log_elevation_km = log( ifelse(temp_DF$elevation<1, 1, temp_DF$elevation) / 1000 )
sf_DF$NDVI = raster::extract( x=copNDVI, y=as(sf_DF,"Spatial") )
sf_DF$scale_NDVI = scale( sf_DF$NDVI )[,1]
sf_DF$pop_dens = pop_dens$Dens2020[ st_nearest_feature( sf_DF, pop_dens ) ]
sf_DF$log_pop_dens = log(sf_DF$pop_dens)

# Note, only the log_pop_dens covariate is used for now

#
#taxa = levels( sf_DF$Genus_species )

# create mesh
loc_DF = st_coordinates(sf_DF)
loc_grid = st_coordinates(st_centroid(sf_grid))
mesh = fm_mesh_2d( loc_grid, refine=TRUE, cutoff=1)
# Create matrices in INLA
spde <- fm_fem(mesh, order=2)
#plot(mesh)

# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=loc_DF )$proj$A
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=loc_grid )$proj$A

# Here the script from https://github.com/James-Thorson/Spatio-temporal-models-for-ecologists/blob/main/Chap_11/BBS.R ends
# and I instead follow the Bering Sea case study

# Other objects
invM0 = invsqrtM0 = spde$c0
diag(invsqrtM0) = 1 / sqrt(diag(spde$c0))
diag(invM0) = 1 / diag(spde$c0)

# Get pop_dens at vertices of SPDE mesh
mesh_points = mesh$loc[,1:2]
pop_dens_s = raster::extract( as_spatraster(pop_dens), mesh_points )[,2] # Dens2020

# FIll in missing covariates
# FIXME:  Could be done using 0 for land and positive values otherwise
pop_dens_s = ifelse( is.na(pop_dens_s), mean(pop_dens_s,na.rm=TRUE), pop_dens_s )

pop_dens_s <- log(pop_dens_s)

# Explore
# sf_DF %>%
#   pivot_longer(4:23) %>%
#   ggplot(aes(log_pop_dens, value)) +
#   geom_smooth(alpha = 0.2) +
#   facet_wrap(~name, scales = "free")


################
# TMB
################

distribution = c("Tweedie", "Poisson", "LNP")[2] # Trying a Poisson first...

#
Date = Sys.Date()
date_dir = file.path(root_dir, "results", paste0(Date,"_","_",distribution) )
dir.create(date_dir, recursive=TRUE)
sf_usa = ne_countries( country="united states of america", return="sf" )
