f.constants <- function(prior = "iid",
                        data = NULL,
                        min_u = NULL,
                        max_u = NULL,
                        carto = NULL){
  
  if(!(prior %in% c("iid","iCAR","pCAR","LCAR","BYM","BYM2", "GP")))
    stop("invalid prior spatial argument")
  
  
  #number of areas
  S.area <- length(unique(carto$ID))
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  
  ##Constants for all priors
  constants <- list(pop = data$Pop,
                    rate = data$crude.rate/10^5,
                    N = S.area,
                    min_unif = min_u,
                    max_unif = max_u)
  
  
  
  if(prior == "GP"){
    g <- st_as_sf(carto)
    g_centroid <- st_point_on_surface(x = g)
    distMatrix <- st_distance(g_centroid, g_centroid)
    distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid)) 
    a <- max(distMatrix)
    
    constants <- append(constants,
                        list(dists = distMatrix,
                             a = a))
  }
  else if(prior%in%c("iCAR", "BYM")){
    constants <- append(constants,
                        list(L = length(nbInfo$adj),
                             num = nbInfo$num,
                             weights = nbInfo$weights,
                             adj = nbInfo$adj))
  }
  if(prior=="pCAR"){
    CM <- as.carCM(num = nbInfo$num, weights = nbInfo$weights, adj = nbInfo$adj)
    
    constants <- append(constants,
                        list(L = length(nbInfo$adj),
                             num = nbInfo$num,
                             adj = nbInfo$adj,
                             C = CM$C, M = CM$M))
    
  }
  else if(prior == "LCAR"){
    # Number of neighbors of each municipality
    nadj <- card(carto.nb)
    # Neighbors of each municipality
    map <- unlist(carto.nb)
    
    # Diagonal matrix with the number of neighbors of each area
    D <- diag(nadj)
    # Adjacency matrix
    W <- nb2mat(carto.nb, style = "B", zero.policy = TRUE)
    # Eigenvalues of D-W
    Lambda <- eigen(D - W)$values
    # Identity matrix
    I <- diag(rep(1, S.area))
    
    # All the neighborhoods j ~ k where k < j
    from.to <- cbind(rep(1:S.area, times = nadj), map); colnames(from.to) <- c("from", "to")
    from.to <- from.to[which(from.to[, 1] < from.to[, 2]), ]
    NDist <- nrow(from.to)
    
    constants <- append(constants,
                        list(NDist = NDist, 
                             Lambda = Lambda,
                             from.to = from.to))
  }
  else if(prior== "BYM2"){
    W <- nb2mat(carto.nb, style = "B", zero.policy = TRUE)
    W.scale <- -W
    diag(W.scale) <- abs(apply(W.scale, 1, sum))
    Q <- W.scale*exp(mean(log(diag(MASS::ginv(as(W.scale,"matrix"))))))
    scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))
    
    
    constants <- append(constants,
                        list(L = length(nbInfo$adj),
                             num = nbInfo$num,
                             weights = nbInfo$weights,
                             adj = nbInfo$adj,
                             scale = scale))
  }
  
  
  
  return(constants)
}