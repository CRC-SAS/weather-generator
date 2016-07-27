# this method is 4x faster than grf
# we can't use sim.rf because the range is too large for circulant embedding approximation
# we will use Cholesky decomposition just one time for each month, each climate variable
# then
# we would want 12 "cSigma" for each climate variable
# here is for January precipitation

library(fields)

grid     <- list(x = unique(predloc.GK[,1]), y = unique(predloc.GK[,2]))
xg       <- make.surface.grid(grid)
bigSigma <- stationary.cov(xg, theta=params[1,3])
# bigSigma <- stationary.cov(xg, theta=298)
cSigma   <- t(chol(bigSigma))
# NOTE:  bigSigma == cSigma%*%t(cSigma)
# only need to find cholesky once.

system.time(w.occ <- cSigma%*%rnorm(nrow(xg)) * sqrt(params[1,2]))

# NOTE: cSigma%*%rnorm(nrow(xg)) has nugget 0, sill=1 by default.
# we want nugget zero for everything, so that's fine
# but we need to multiply by the sqrt of process variance
# for probit regression this is unnecessary (i.e., variance = 1)
# but for temperatures, this will alter the field by a constant factor

quilt.plot(predloc.GK, w.occ)
