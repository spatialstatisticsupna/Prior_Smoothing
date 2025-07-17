################################################################################
################################################################################
##########              Theoretically expected smoothing              ##########
################################################################################
################################################################################
library(spdep)
# library(Matrix)
library(MASS)
library(sf)

load("../../Data_Spain/Carto_LUNG_LOCP.RData")
carto <- Carto_LUNG_LOCP
carto <- carto[order(carto$district),]

sf::sf_use_s2(FALSE)
carto.nb <- poly2nb(carto)
W.nb <- nb2mat(carto.nb, style="B")
nbInfo <- nb2WB(carto.nb)

neigh <- nbInfo$num

sigma <- c(0.001, 0.01, 0.03, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.5, 0.9)
lambda <- c(0.1, 0.5, 0.9)
rho <- c(1, 5, 9)

################################################################################
##########                            iid                             ##########
################################################################################
tab.iid <- matrix(NA, ncol = 1, nrow = length(sigma))
for (sg in 1:length(sigma)) {
  tab.iid[sg, 1] <- length(neigh)*sigma[sg]^2
}

tab.iid <- cbind(sigma, tab.iid)


################################################################################
##########                            GP                              ##########
################################################################################
g <- st_as_sf(carto)
g_centroid <- st_point_on_surface(x = g)
distMatrix <- st_distance(g_centroid, g_centroid)
distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid))
n <- 47

tab.GP <- matrix(NA, ncol = length(rho), nrow = length(sigma))
for (sg in 1:length(sigma)) {
  for (rh in 1:length(rho)) {
    R <- matrix(nrow = n, ncol = n)
    for(i in 1:n){
      for(j in 1:n){
        R[i, j] <- exp(-distMatrix[i,j]/rho[rh]) 
      }
    }
    inv.R <- MASS::ginv(as(R, "matrix"))
    
    tab.GP[sg, rh] <- round(sum(sigma[sg]^2/diag(inv.R)),3)
  }
}

tab.GP <- cbind(sigma, tab.GP)
tab.GP <- rbind(c(" ", rho), tab.GP)


################################################################################
##########                           iCAR                             ##########
################################################################################
tab.iCAR <- matrix(NA, ncol = 1, nrow = length(sigma))
for (sg in 1:length(sigma)) {
    tab.iCAR[sg, 1] <- round(sum(sigma[sg]^2/neigh),3)
}

tab.iCAR <- cbind(sigma, tab.iCAR)



################################################################################
##########                            BYM                             ##########
################################################################################
W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
W.scale <- -W
diag(W.scale) <- abs(apply(W.scale, 1, sum))
inv.R <- MASS::ginv(as(W.scale, "matrix"))


tab.BYM2 <- matrix(NA, ncol = length(sigma), nrow = length(sigma))
for (sg in 1:length(sigma)) {
  for (sg2 in 1:length(sigma)) {
    # tab.BYM[sg, sg2] <- round(sum(sigma[sg]^2/neigh + sigma[sg2]^2),3)
    
    # inv.Q <- MASS::ginv(sigma[sg]^2*inv.R + sigma[sg2]^2*diag(1, length(nbInfo$num)))
    # tab.BYM[sg, sg2] <- round(sum(1/diag(inv.Q)),3)
    
    inv.Q2 <- MASS::ginv(inv.R + (sigma[sg2]/sigma[sg])^2*diag(1, length(nbInfo$num)))
    tab.BYM2[sg, sg2] <- round(sum(sigma[sg]^2/diag(inv.Q2)),3)
  }
}

tab.BYM2 <- cbind(sigma, tab.BYM2)
tab.BYM2 <- rbind(c(" ", sigma), tab.BYM2)



################################################################################
##########                           pCAR                             ##########
################################################################################
tab.pCAR <- matrix(NA, ncol = 1, nrow = length(sigma))
for (sg in 1:length(sigma)) {
  tab.pCAR[sg, 1] <- round(sum(sigma[sg]^2/neigh),3)
}

tab.pCAR <- cbind(sigma, tab.pCAR)



################################################################################
##########                           LCAR                             ##########
################################################################################
tab.LCAR <- matrix(NA, ncol = length(lambda), nrow = length(sigma))
for (sg in 1:length(sigma)) {
  for (lb in 1:length(lambda)) {
    tab.LCAR[sg, lb] <- round(sum(sigma[sg]^2/(lambda[lb]*neigh + 1- lambda[lb])),3)
  }
}

tab.LCAR <- cbind(sigma, tab.LCAR)
tab.LCAR <- rbind(c(" ", lambda), tab.LCAR)



################################################################################
##########                           BYM2                             ##########
################################################################################
W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
nbInfo <- nb2WB(carto.nb)
# ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
# carto_use <- sf::st_as_sf(carto)
# W.nb <- nb2listw(poly2nb(carto_use))$neighbours
# nbInfo <- nb2WB(W.nb)



W.scale <- -W
diag(W.scale) <- abs(apply(W.scale, 1, sum))
# solve(W.scale) # this should not work since by definition the matrix is singular
library(INLA)
Q = inla.scale.model(W.scale, constr=list(A=matrix(1, nrow=1, ncol=nrow(W.scale)), e=0))
scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))

# Qs.scaled <- W.scale*exp(mean(log(diag(INLA:::inla.ginv(W.scale)))))
# scale1 = exp((1/nrow(W.scale))*sum(log(1/diag(Qs.scaled))))
# 
# Qs.scaled2 <- W.scale*exp(mean(log(diag(MASS::ginv(as(W.scale,"matrix"))))))
# scale = exp((1/nrow(W.scale))*sum(log(1/diag(Qs.scaled2))))
# nbInfo$weights <- 0.41179*nbInfo$weights
# nbInfo$weights <- rep(1,222)
# all.equal(Qs.scaled,Qs.scaled2)
library(matlib)

inv.R <- MASS::ginv(as(Q, "matrix"))

tab.BYM2 <- matrix(NA, ncol = length(lambda), nrow = length(sigma))
for (sg in 1:length(sigma)) {
  for (lb in 1:length(lambda)) {
    # tab.BYM2[sg, lb] <- round(sum(sigma[sg]^2*(lambda[lb]/diag(Q) + 1- lambda[lb])),3)
    inv.Q <- MASS::ginv(lambda[lb]*inv.R + (1-lambda[lb])*diag(1, length(nbInfo$num)))
    tab.BYM2[sg, lb] <- round(sum(sigma[sg]^2/diag(inv.Q)),3)
  }
}

tab.BYM2 <- cbind(sigma, tab.BYM2)
tab.BYM2 <- rbind(c(" ", lambda), tab.BYM2)

