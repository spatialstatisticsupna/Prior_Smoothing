f.inits <- function(prior = "iid",
                    carto = NULL,
                    min.unif = NULL,
                    max.unif = NULL){
  
  if(!(prior %in% c("iid","iCAR","pCAR","LCAR","BYM","BYM2", "GP")))
    stop("invalid prior spatial argument")
  
  #number of areas
  S.area <- length(unique(carto$ID))
  
  if(prior =="iid"){
    inits = function(){
      list(alpha = 0, sigma = runif(1, min.unif, max.unif), 
           iid = rnorm(S.area, sd = 0.1))}
  }
  else if(prior == "iCAR"){
    inits = function(){
      list(alpha = 0, sigma = runif(1, min.unif, max.unif), 
           theta = rnorm(S.area, sd = 0.1))}
  }
  else if(prior =="pCAR"){
    inits = function(){
      list(alpha = 0, sigma = runif(1, min.unif, max.unif),
           theta = rnorm(S.area, sd = 0.1), 
           gamma = runif(1, -1, 1))}
  }
  else if(prior == "LCAR"){
    inits = function(){
      list(alpha = 0, sd.theta = runif(1, min.unif, max.unif), 
           rho = runif(1), 
           theta = rnorm(S.area, sd = 0.1))}
  }
  else if(prior == "BYM"){
    inits = function(){
      list(alpha = 0, 
           sigma.phi = runif(1, min.unif, max.unif), 
           sigma.iid = runif(1, min.unif, max.unif),
           iid = rnorm(S.area, sd = 0.1), 
           phi = rnorm(S.area, sd=0.1))}
  }
  else if(prior == "BYM2"){
    inits = function(){
      list(alpha = 0, 
           sigma.b = runif(1, min.unif, max.unif), 
           theta = rnorm(S.area, sd = 0.1), 
           phi = rnorm(S.area, sd = 0.1), 
           rho = runif(1))}
  }
  else if(prior == "GP"){
    g <- st_as_sf(carto)
    g_centroid <- st_point_on_surface(x = g)
    distMatrix <- st_distance(g_centroid, g_centroid)
    distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid)) 
    a <- max(distMatrix)
    
    inits = function(){
      list(alpha = 0, sigma = runif(1, min.unif, max.unif), 
           rho = runif(1, 0, a), 
           theta = rnorm(S.area, sd = 0.1))}
  }
  
  
  return(inits)
  
}
