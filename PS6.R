##############################################################
# Statistical Programming 
# Problem Set 6
# Author: Jonas Markgraf
# Date: April 17-19, 2017
###########################

# load packages
library(testthat)
library(microbenchmark) 
library(cubature)

rm(list = ls())

# Function that can deal with multiple dimensions and parallelize estimation
# ===========================================================================

sg.int<-function(g,...,lower, upper, parallel = F){
  # load packages that are required to run function
  require("SparseGrid")
  require("dplyr")  # for parallelized function
  require("parallel")  # to detect number of cores
  require("doMC")  # to register number of cores
  
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
  # and upper bound values in steps of 1, depending on number of dimensions
  gridss<-as.matrix(expand.grid(lapply(1:no.dim, function(x) seq(lower[x], upper[x]-1, by=1))))
  
  # function creates nodes and weights for integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension=no.dim, k=5 )  # 'KPU'for unweighted integral; 'k' for accuracy
  # creating nodes for integration by summing gridds and previosuly defined nodes
  nodes<-gridss[1,]+sp.grid$nodes
  # taking weights values from createIntegrationGrid function result
  weights<-sp.grid$weights
  # looping over gridss object to create node for each point and weight for respective point
  for (i in 2:nrow(gridss)) {
    nodes<-rbind(nodes,gridss[i,]+sp.grid$nodes)  
    weights<-c(weights,sp.grid$weights)
  }
  
  # apply function 'g' for all nodes, option to parallelize function
  gx.sp <- aaply(nodes, 1, g, .parallel = parallel)
  # weigh the resulting value
  val.sp <- gx.sp %*%weights
  
  # return value
  val.sp
}

# Generate test functions ###########################
# ===================================================

# create test functions for unit testing and speed benchmarking

test2dim <- function(x) 4*x[1] + 3*x[2]^2  # 2 dimensions

test3dim <- function(x) 4*x[1] + 3*x[2]^2 + 2*x[3]  # 3 dimensions

test4dim <- function(x) 4*x[1] + 3*x[2]^2 + 2*x[3] + x[4]  # 4 dimensions


# Unit Testing ######################################
# ===================================================

# test that function estimates correct values
test_that("Test that sparsegrid calculates same values than adaptIntegrate",
  expect_equal(as.numeric(sg.int(test2dim, lower = c(0,0), upper = c(2,2))),
    as.numeric(adaptIntegrate(test2dim, lower = c(0,0), upper = c(2,2))[1])))
  
# test that function throws an error if upper bound is smaller than lower bound
test_that("Test that sparsegrid throws error if lower > upper bound",  
  expect_error(sg.int(test2dim, lower = c(2,2), upper = c(0,0))))

# Measure speed: ####################################
# for parallelized and non-parallelized function 
# & between functions sparsegrid vs adaptIntegrate
# ===================================================

# measure speed for 2 dimensions --------
microbenchmark(
  "solo"=sg.int(g=test2dim, lower = c(0, 0), upper = c(4, 4)),
  "parallel"=sg.int(g=test2dim, lower = c(0, 0), upper = c(4, 4), parallel = T),
  "adaptInt" = adaptIntegrate(f = test2dim, lower = c(0, 0, 0), upper = c(4, 4, 4)),
  times =20
)
## solo vs. parallel: function without parallelization is faster
## adaptInt vs. sparse grid: adaptInt is much faster

# measure speed for 3 dimensions --------
microbenchmark(
  "solo"=sg.int(g=test3dim, lower = c(0, 0, 0), upper = c(4, 4, 4)),
  "parallel"=sg.int(g=test3dim, lower = c(0, 0, 0), upper = c(4, 4, 4), parallel = T),
  "adaptInt" = adaptIntegrate(f = test3dim, lower = c(0, 0, 0), upper = c(4, 4, 4)),
  times =5
)
## solo vs. parallel: function without parallelization is faster
## adaptInt vs. sparse grid: adaptInt is much faster

# measure speed for 4 dimensions --------
microbenchmark(
  "solo"=sg.int(g=test4dim, lower = c(0, 0, 0, 0), upper = c(4, 4, 4, 4)),
  "parallel"=sg.int(g=test4dim, lower = c(0, 0, 0, 0), upper = c(4, 4, 4, 4), parallel = T),
  "adaptInt" = adaptIntegrate(f = test4dim, lower = c(0, 0, 0), upper = c(4, 4, 4)),
  times=5
)
## solo vs. parallel: function without parallelization is faster
## adaptInt vs. sparse grid: adaptInt is much faster

#################################### MAXIMIZATION IN R ################################
#######################################################################################

# generate function from slide
func <- function(x, y) {
  sin((x^2) / 2 - (y^2) / 4) * cos((2*x) - exp(y))
}

# use 'optimize' function to find local maxima
optimize(f = func, lower = c(-1,-1), upper = c(3,3), maximum = T)

# use 'optim' function to find local maxima
optim(par = c(2,1), fn = func, lower = c(-1,-1), upper = c(3,3), method = "L-BFGS-B")
