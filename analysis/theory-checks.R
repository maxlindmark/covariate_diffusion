
# FROM: C:\Users\James.Thorson\Desktop\Work files\Collaborations\2024 -- movement kernel for covariates\Old code\covariate_diffusion_2024-02-05.R

library(fmesher)
library(Matrix)
library(pracma)
library(sf)

#
#loc = matrix( runif(1000), ncol=2 )
#loc = as.matrix( expand.grid(seq(0,1,len=20),seq(0,1,len=20)) )
loc = poisson2disk(n=1000)
mesh = fm_mesh_2d( loc )
spde = fm_fem( mesh )
invM0 = spde$c0
diag(invM0) = 1/diag(spde$c0)

#
sf_mesh = st_union(fm_as_sfc(mesh))
sf_grid = st_make_grid( sf_mesh, n=100 )
sf_grid = st_intersection( sf_grid, sf_mesh )
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(sf_grid)) )$proj$A

#
ln_kappa = log(10)
#ln_tau = log(0.4)
ln_tau = log( 1 / (1 + exp(2 * ln_kappa)) )
#Q = exp(ln_tau*2) * (exp(ln_kappa*4) * spde$c0 + 2 * exp(ln_kappa*2) * spde$g1 + spde$g2)

#
invD = exp(ln_tau) * invM0 %*% (spde$c0 + exp(2*ln_kappa)*spde$c0 + spde$g1)
D = solve(invD)
# rowSums(Dt) = constant
# colSums(D) = constant ... conserves mass

#
#invDt = exp(ln_tau) * (spde$c0 + exp(2*ln_kappa)*spde$c0 + spde$g1) %*% invM0

# Sanity checks
f = function(kappa, tau) 1 / (1 + kappa^2) / tau
D = solve(invD)

# colSums( A_gs * x ) should equal colSums( A_gs * D * x )
# Check with midpoint of domain
which_mid = which.min( rowSums(scale(mesh$loc[,1:2])^2) )
vec1 = rep(0,mesh$n)
vec1[which_mid] = 1
# Cmopare the two
t(rep(1,length(sf_grid))) %*% A_gs %*% vec1
t(rep(1,length(sf_grid))) %*% A_gs %*% solve(invD,vec1)

# plot diffusion from a specified point
# plotting code
stuff = st_sf( sf_grid,
               "orig"=as.numeric(A_gs%*%vec1),
               "proj"=as.numeric(A_gs%*%solve(invD,vec1)) )   # , log(as.numeric(vec2b)))
plot( stuff, cex=2, pch=19, border=NA )
sum( stuff$orig )
sum( stuff$proj )

# Using IID normal deviates
vec1 = rnorm(mesh$n)
# Cmopare the two
t(rep(1,length(sf_grid))) %*% A_gs %*% vec1
t(rep(1,length(sf_grid))) %*% A_gs %*% solve(invD,vec1)

# plot diffusion from a specified point
# plotting code
stuff = st_sf( sf_grid,
               "orig"=as.numeric(A_gs%*%vec1),
               "proj"=as.numeric(A_gs%*%solve(invD,vec1)) )   # , log(as.numeric(vec2b)))
plot( stuff, cex=2, pch=19, border=NA )
sum( stuff$orig )
sum( stuff$proj )


###################
#
# Work out proportionality constant
#  Visualize 
#
###################

f = function( kappa, tau ){
  ln_kappa = log(kappa)
  ln_tau = log(tau)
  invD = exp(ln_tau) * invM0 %*% (spde$c0 + exp(2*ln_kappa)*spde$c0 + spde$g1)
  vec = solve(invD, rep(1,mesh$n))
  mean(vec)
}
X = seq(0,10,length=21)
Y = sapply( X, FUN=\(x)f(x,1) )

# Formula I worked out by pattern matching
f2 = \(kappa, tau) 1 / (1 + kappa^2) / tau
Y2 = sapply( X, FUN=\(x)f2(x,1) )

# Compare the two
matplot( x=X, y=cbind(1/Y,1/Y2) )


