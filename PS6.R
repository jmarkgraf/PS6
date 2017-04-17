sg.int<-function(g,...,lower,upper, dimension = no.dim){
  require("SparseGrid")
  
  # adjusting upper and lower bound downward
  lower<-floor(lower)
  upper<-ceiling(upper)
  
  # stop function if window for integration is wrongly specified
  if (any(lower>upper)) stop("lower must be smaller than upper")
  
  # function that creates matrix of all possible combinations of values between the lower 
  # and upper bound values, depending on number of dimensions
  gridss<-as.matrix(expand.grid(lapply(function(x) seq(lower[x], upper[x]-1, by=1))))
  
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
  
  gx.sp <- apply(nodes, 1, g,...)
  val.sp <- gx.sp %*%weights
  val.sp
}