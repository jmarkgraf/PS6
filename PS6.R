##############################################################
# Statistical Programming 
# Problem Set 6
# Author: Jonas Markgraf
# Date: April 17, 2017

library(testthat)
library(microbenchmark)
library()


# Function that can deal with multiple dimensions and parallelize estimation
# ===========================================================================

sg.int<-function(g,...,lower,upper, parallel = F){
  require("SparseGrid")
  require("dplyr")
  require("parallel")
  require("doMC")
  
  # identify and register number of cores for parallelization
  cores <- detectCores()
  registerDoMC(cores)
  
  # identify number of dimensions
  no.dim <- length(upper)
  
  # adjusting upper and lower bound downward/upward 
  lower<-floor(lower)
  upper<-ceiling(upper)
  
  # stop function if window for integration is wrongly specified
  if (any(lower>upper)) stop("lower must be smaller than upper")
  
  # function that creates matrix of all possible combinations of values between the lower 
  # and upper bound values, depending on number of dimensions
  gridss<-as.matrix(expand.grid(lapply(1:no.dim, function(x) seq(lower[x], upper[x]-1, by=1))))
  
  # function creates nodes and weights for integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension=no.dim, k=5 )
  # creating nodes for integration
  nodes<-gridss[1,]+sp.grid$nodes
  # creating weights for integration
  weights<-sp.grid$weights
  
  for (i in 2:nrow(gridss)) {
    nodes<-rbind(nodes,gridss[i,]+sp.grid$nodes)  
    weights<-c(weights,sp.grid$weights)
  }
  
  # parallelize 'apply' function
  gx.sp <- aaply(nodes, 1, g, .parallel = parallel)
  val.sp <- gx.sp %*%weights
  
  # return value
  val.sp
}


# Measure speed for parallelized and non-parallelized function
# =============================================================

# create test functions
test2dim <- function(x) x[1]^4 + 2*x[2]  # 2 dimensions
  
test3dim <- function(x) 4*x[1] + 3*x[2]^2 + 2*x[3]^3  # 3 dimensions

test4dim <- function(x) 4*x[1] + 3*x[2]^2 + 2*x[3]^3 + x[4]^4  # 4 dimensions

# measure speed for 2 dimensions
microbenchmark(
  "solo"=sg.int(g=test2dim, lower = c(0, 0), upper = c(4.2, 4.2)),
  "parallel"=sg.int(g=test2dim, lower = c(0, 0), upper = c(4.2, 4.2), parallel = T),
  times =20
)

# measure speed for 3 dimensions
microbenchmark(
  "solo"=sg.int(g=test3dim, lower = c(0, 0, 0), upper = c(4.2, 4.2, 4.2)),
  "parallel"=sg.int(g=test3dim, lower = c(0, 0, 0), upper = c(4.2, 4.2, 4.2), parallel = T),
  times =20
)

# measure speed for 3 dimensions
microbenchmark(
  "solo"=sg.int(g=test4dim, lower = c(0, 0, 0, 0), upper = c(4.2, 4.2, 4.2, 4.2)),
  "parallel"=sg.int(g=test4dim, lower = c(0, 0, 0, 0), upper = c(4.2, 4.2, 4.2, 4.2), parallel = T),
  times =20
)

