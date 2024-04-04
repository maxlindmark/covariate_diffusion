
# Custom function for adding legend in multi-panel figure
# https://stackoverflow.com/questions/52975447/reorganize-sf-multi-plot-and-add-a-legend
add_legend <-
function( legend,
          col = sf.colors(),
          legend_x = c(0.9,  1.0),
          legend_y = c(0.05, 0.45),
          text_col = "black",
          ...){

    # Get the axis limits and calculate size
    axisLimits <- par()$usr
    xLength <- axisLimits[2] - axisLimits[1]
    yLength <- axisLimits[4] - axisLimits[3]

    xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
    xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
    yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
    yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
    if( diff(legend_y) > diff(legend_x) ){
      align = c("lt","rb")[2]
      gradient = c("x","y")[2]
    }else{
      align = c("lt","rb")[1]
      gradient = c("x","y")[1]
    }

    # Add the legend
    plotrix::color.legend( xl = xl,
                           xr = xr,
                           yb = yb,
                           yt = yt,
                           legend = legend,
                           rect.col = col,
                           gradient="y",
                           col = text_col,
                           ... )
}

###################
# Using GAP-EFH data
#
# NOTE:  Must connect to VPN and then open "Y:\RACE_EFH_variables\Variables\Variables_EBS_1km" in explorer
#
###################

library(sf)
library(Matrix)
library(TMB)
library(rnaturalearth)
library(terra)

region = c("EBS", "GOA")[1]

root_dir = here::here(".")
data_dir = file.path( root_dir, "data" )
tmb_dir = file.path( root_dir, "tmb" )

if( region=="GOA" ){
  load( file.path(data_dir,"region_data_all.rda") )
  load( file.path(data_dir,"raster_stack.rda") )
  domain_shape = st_read( file.path(data_dir,"GOAThorson") )  # , promote_to_multi=FALSE
}
if( region=="EBS" ){
  region_data_all = read.csv( file.path(data_dir,"all_EBS_data_2021.csv") )
  # Shapefile
  shape1 = st_read( file.path(data_dir,"EBSThorson") )  # , promote_to_multi=FALSE
  shape2 = st_read( file.path(data_dir,"NBSThorson") )  # , promote_to_multi=FALSE
  domain_shape = c( st_geometry(shape1), st_geometry(shape2) )

  # NOTE:  Must connect to VPN and then open "Y:\RACE_EFH_variables\Variables\Variables_EBS_1km" in explorer
  #load( file.path(data_dir,"EBS_rasterstack.rds") )
  #raster_stack = raster.stack
  #bdepth = raster_stack$bdepth

  #
  bdepth = rast( file.path(data_dir,"bdepth.tif") )
}

# Make domain from shapefile
# See: https://github.com/mstrimas/smoothr/blob/main/R/fill-holes.r to fill holes
domain_shape = st_geometry(domain_shape)
domain_shape = st_cast( domain_shape, to="POLYGON" )
if(region=="GOA") domain_shape = domain_shape[ which.max(st_area(domain_shape)) ]
domain_shape = st_transform( domain_shape, crs = 4326 )
domain_shape = st_make_valid( domain_shape )
plot(domain_shape)

#
loc = st_multipoint( as.matrix(region_data_all[,c('lon','lat')]), dim="XY" )
sf_loc = st_sfc( loc, crs=st_crs(bdepth) )
sf_loc = st_transform( sf_loc, crs = 4326 )

# domain from data
domain = st_concave_hull(sf_loc, ratio=0.1)  # 0.1 works for both GOA and EBS
domain = st_sfc( domain, crs=st_crs(sf_loc) )

# domain from mesh
sf_grid = st_make_grid( domain_shape, cellsize=c(0.1,0.1) )
sf_grid = st_intersection( sf_grid, domain_shape )
grid_loc = st_coordinates( st_centroid(sf_grid) )
#plot(sf_grid)

library(fmesher)
mesh = fm_mesh_2d( st_coordinates(sf_loc)[,1:2],
                   #loc.domain = domain,
                   #offset = c(-0.02),
                   #interior = fm_segm( domain, boundary=FALSE),
                   cutoff = 0.1 );

# Other objects
spde = fm_fem( mesh )
A_is = fm_evaluator( mesh, loc=st_coordinates(sf_loc) )$proj$A
A_gs = fm_evaluator( mesh, loc=grid_loc )$proj$A
invM0 = invsqrtM0 = spde$c0
diag(invsqrtM0) = 1 / sqrt(diag(spde$c0))
diag(invM0) = 1 / diag(spde$c0)

#sf_mesh = st_union(fm_as_sfc(mesh))
#sf_meshgrid = st_make_grid( sf_mesh, n=50 )
#sf_grid = st_intersection( sf_grid, sf_mesh )

# Get depth at vertices of SPDE mesh
mesh_points = sf_project( mesh$loc[,1:2], from=st_crs(4326), to=st_crs(bdepth) )
depth_s = extract( bdepth, mesh_points )[,1]

# FIll in missing covariates
# FIXME:  Could be done using 0 for land and positive values otherwise
depth_s = ifelse( is.na(depth_s), mean(depth_s,na.rm=TRUE), depth_s )


################
# TMB
################

#
f_depth = c("identity", "log")[1]
distribution = c("Tweedie", "Poisson", "LNP")[3]    # Poisson is numerically unstable for j_nrs

#
Date = Sys.Date()
date_dir = file.path(root_dir, "results", paste0(Date,"_",f_depth,"_",distribution) )
dir.create(date_dir, recursive=TRUE)
sf_usa = ne_countries( country="united states of america", return="sf" )

# Compile
setwd( tmb_dir )
compile( "movement_kernel.cpp", framework="TMBad" )
dyn.load( dynlib("movement_kernel") )
#dyn.unload( dynlib("movement_kernel") )

# saving stuff
species_set = colnames(region_data_all)[30:ncol(region_data_all)]
N_c = colSums( ifelse(region_data_all[,species_set]>0,1,0) )
species_set = species_set[ N_c > 1000 ]

#
param_set = c("obj_diffusion", "obj_null", "deltaAIC", "range", "ln_kappa2", "corr_depth")
Results_cz = array( NA,
                    dim = c(length(species_set),length(param_set)),
                    dimnames = list(species_set,param_set) )

# Loop
cI = match( "a_yfs", species_set )
for( cI in seq_along(species_set) ){
  species = species_set[cI]

  # Build object
  if(f_depth=="identity"){
    depthprime_s = (depth_s / 100)
  }else if(f_depth=="log"){
    depthprime_s = log(depth_s) - 4.5
  }
  Data = list( "method" = "diffusion",
               "dist" = distribution,
               "c_i" = region_data_all[,species],
               "A_is" = A_is,
               "A_gs" = A_gs,
               "M0" = spde$c0,
               "M1" = spde$g1,
               "M2" = spde$g2,
               "invsqrtM0" = invsqrtM0,
               "invM0" = invM0,
               "depth_s" = depthprime_s,
               "sim_gmrf" = 0L )
  Params = list( "beta0"=0,
                 "beta_j" = c(0.1 ,0.1),
                 "ln_tau"=0,
                 "ln_kappa"=0,
                 "omega_s"=rnorm(nrow(spde$c0)),
                 "ln_kappa2"=exp(-1) )
  Random = "omega_s"

  # Special stuff
  if( Data$dist == "Tweedie" ){
    Params$ln_phi = log(2)
    Params$finv_power = 0
  }
  if( Data$dist == "LNP" ){
    Params$ln_sigma_eta = log(0.1)
    Params$eta_i = rnorm(nrow(region_data_all))
    Random = c(Random, "eta_i")
  }

  # profile = c("beta0","beta_j") is unstable for j_berfl and depthprime_s = (depth_s / 100)
  Obj = MakeADFun( data = Data,
                   parameters = Params,
                   #profile = c("beta0","beta_j"),    # ,"ln_kappa2"
                   random = Random,
                   checkParameterOrder = TRUE )
  Obj$env$beSilent()

  # Optimize
  Opt = nlminb( start = Obj$par,
                obj = Obj$fn,
                grad = Obj$gr,
                control = list(eval.max=1e4, iter.max=1e4, trace=1) )
  Report = Obj$report()

  # Stationary distribution
  Results_cz[cI,'obj_diffusion'] = Opt$objective
  Results_cz[cI,'range'] = Report$Range2
  Results_cz[cI,'ln_kappa2'] = Opt$par['ln_kappa2']
  Results_cz[cI,'corr_depth'] = cor( Report$depth_g, as.vector(A_gs%*%Data$depth_s) )

  # Re-run without covariate diffusion
  Data2 = Data
  Data2$method = "null"
  Params2 = Params[setdiff(names(Params),c("ln_tau2","ln_kappa2"))]
  Obj2 = MakeADFun( data = Data2,
                    parameters = Params2,
                    #profile = c("beta0","beta_j"),
                    #profile = c("beta_j"),
                    random = Random,
                    checkParameterOrder = TRUE )
  Obj2$env$beSilent()

  # Optimize
  Opt2 = nlminb( start = Obj2$par,
                 obj = Obj2$fn,
                 grad = Obj2$gr,
                 control = list(eval.max=1e4, iter.max=1e4, trace=1) )
  Report2 = Obj2$report()

  #
  Results_cz[cI,'deltaAIC'] = (2*Opt2$objective+2*length(Opt2$par)) - (2*Opt$objective+2*length(Opt$par))
  Results_cz[cI,'obj_null'] = Opt2$objective
  write.csv( Results_cz, file=file.path(date_dir,"Results_cz.csv") )

  # Plot results
  trunc_dens = function(logdens, ratio=0.01) ifelse( logdens<(max(logdens,na.rm=TRUE)+log(ratio)), NA, logdens )
  png( file=file.path(date_dir,paste0(species,".png")), width=6, height=6, res=200, units="in" )
    stuff_gz = st_sf( sf_grid,
                       "orig"=as.vector(Data$A_gs%*%Data$depth_s),
                       "log(dens) orig" = trunc_dens(log(Report2$mu_g)),
                       "pdepth orig" = trunc_dens(Report2$pdepth_g),
                       "diffused"=Report$depth_g,
                       "log(dens) diffused" = trunc_dens(log(Report$mu_g)),
                       "pdepth diffused" = trunc_dens(Report$pdepth_g) )
    #plot(stuff_gz, border=NA )
    par( mfcol=c(3,2) )
    for( i in 1:6 ){
      plot( stuff_gz[i], border=NA, key.pos=NULL, reset=FALSE )
      plot( sf_usa, add=TRUE, col="grey" )
      add_legend( round(range(stuff_gz[[i]],na.rm=TRUE),1) )
    }
  dev.off()
}
