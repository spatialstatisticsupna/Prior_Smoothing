rm(list=ls())
library(spdep)
library(MASS)
library(sf)


#################################################
## Load the cartography file  ##
#################################################
load("../../Data_Spain/Carto_LUNG_LOCP.RData")
carto <- Carto_LUNG_LOCP
carto <- carto[order(carto$district),]
S.area <- length(unique(carto$district))

################################################
## Load data  ##
#################################################
load("./Lung_cancer_Spain2021.Rdata")
##Order data as in carto
data <- data[order(data$Group.1),]
data$Group.1 == carto$ID
# data$Group.1 <- as.numeric(data$Group.1)
data$crude.rate <- data$O/data$POB*10^5

## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
sf::sf_use_s2(FALSE)
cv.nb <- poly2nb(carto)
W.nb <- nb2mat(cv.nb, style="B")
W <- nb2mat(cv.nb, zero.policy = TRUE, style = "B")
nbInfo <- nb2WB(cv.nb)


library(nimble)

################################################################################
######iid
################################################################################
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  # iid[1:N] ~ dmnorm(mean = rep(0, N), prec = tau * diag(N))
  sigma2 ~ dunif (0, 100)
  tau <- 1 / sigma2
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + iid[i]
    
    iid[i] ~ dnorm(0, tau) 
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    cond.var[i] <- sigma2
  }
})


set.seed(09052024)
constants <- list(N = S.area)
inits <- list(list(alpha = 0, sigma2 = runif(1), iid = rnorm(S.area)),
              list(alpha = 0, sigma2 = runif(1), iid = rnorm(S.area)),
              list(alpha = 0, sigma2 = runif(1), iid = rnorm(S.area)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5)

iid.res <- nimbleMCMC(code = code,
                      constants = constants,
                      data = data.nimble,
                      inits = inits,
                      nchains = 3,
                      niter = 300000,
                      nburnin = 50000,
                      thin = 75,
                      summary = TRUE,
                      samples = TRUE,
                      monitors = c('alpha', 'sigma2', 'tau',
                                   'r', 'MSS.r', 'iid',
                                   'RMSS.r', 'cond.var'),
                      samplesAsCodaMCMC = TRUE,
                      setSeed = TRUE,
                      WAIC = TRUE)




################################################################################
######iCAR
################################################################################

code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  theta[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 1)
  sigma2 ~ dunif (0, 100)
  tau <- 1 / sigma2
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    cond.var[i] <- sigma2/num[i]
  }
})


set.seed(09052024)
constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                  weights = nbInfo$weights, adj = nbInfo$adj)

inits <- list(list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5)

iCAR.res <- nimbleMCMC(code = code,
                       constants = constants,
                       data = data.nimble,
                       inits = inits,
                       nchains = 3,
                       niter = 300000,
                       nburnin = 50000,
                       thin = 75,
                       summary = TRUE,
                       samples = TRUE,
                       monitors = c('alpha', 'sigma2', 'tau',
                                    'r', 'MSS.r', 'theta',
                                    'RMSS.r', 'cond.var'),
                       samplesAsCodaMCMC = TRUE,
                       setSeed = TRUE,
                       WAIC = TRUE)


################################################################################
######pCAR
################################################################################
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  theta[1:N] ~ dcar_proper(mu[1:N], C[1:L], adj[1:L], num[1:N], M[1:N], tau, gamma)
  #M: is a diagonal matrix of conditional variances
  #C: Spatial proximity matrix. ??normalized weights!!  (each row sums to one)
  sigma2 ~ dunif (0, 100)
  tau <- 1 / sigma2
  gamma ~ dunif(-1,1)
  ## gamma ~ dunif(min.bound,max.bound)
  mu[1:N] <- rep(0, N)
  
  
  #Constraint
  zero.theta ~ dnorm(theta.mean, 10000)
  theta.mean <- mean(theta[1:N])
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    cond.var[i] <- sigma2/(num[i])
  }
})


set.seed(09052024)
CM <- as.carCM(num = nbInfo$num, weights = nbInfo$weights, adj = nbInfo$adj) #To create C and M.
constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                  adj = nbInfo$adj, C = CM$C, M = CM$M)

inits <- list(list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)),
              list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)),
              list(alpha = 0, sigma2 = runif(1), theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5,  zero.theta = 0)


pCAR.res <- nimbleMCMC(code = code,
                       constants = constants,
                       data = data.nimble,
                       inits = inits,
                       nchains = 3,
                       niter = 300000,
                       nburnin = 50000,
                       thin = 75,
                       summary = TRUE,
                       samples = TRUE,
                       monitors = c('alpha', 'sigma2', 'tau', 'gamma',
                                    'r', 'MSS.r', 'theta',
                                    'RMSS.r', 'cond.var'),
                       samplesAsCodaMCMC = TRUE,
                       setSeed = TRUE,
                       WAIC = TRUE)






################################################################################
######LCAR
################################################################################
# Number of neighbors of each municipality
nadj <- card(cv.nb)
# Neighbors of each municipality
map <- unlist(cv.nb)

# Diagonal matrix with the number of neighbors of each area
D <- diag(nadj)
# Adjacency matrix
W <- nb2mat(cv.nb, style = "B", zero.policy = TRUE)
# Eigenvalues of D-W
Lambda <- eigen(D - W)$values
# Identity matrix
I <- diag(rep(1, S.area))

# All the neighborhoods j ~ k where k < j
from.to <- cbind(rep(1:S.area, times = nadj), map); colnames(from.to) <- c("from", "to")
from.to <- from.to[which(from.to[, 1] < from.to[, 2]), ]
NDist <- nrow(from.to)


dcar_leroux <- nimbleFunction(
  name = 'dcar_leroux',
  run = function(x = double(1),        # Spatial random effect (vector)
                 rho = double(0),      # Amount of spatial dependence (scalar)
                 sd.theta = double(0), # Standard deviation (scalar)
                 Lambda = double(1),   # Eigenvalues of matrix D - W
                 from.to = double(2),  # Matrix of distinct pairs of neighbors from.to[, 1] < from.to[, 2]
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Number of small areas
    N <- dim(x)[1]
    # Number of distinct pairs of neighbors
    NDist <- dim(from.to)[1]
    # Required vectors
    x.from <- nimNumeric(NDist)
    x.to <- nimNumeric(NDist)
    for (Dist in 1:NDist) {
      x.from[Dist] <- x[from.to[Dist, 1]]
      x.to[Dist] <- x[from.to[Dist, 2]]
    }
    
    # Log-density
    logDens <- sum(dnorm(x[1:N], mean = 0, sd = sd.theta * pow(1 - rho, -1/2), log = TRUE)) -
      N/2 * log(1 - rho) + 1/2 * sum(log(rho * (Lambda[1:N] - 1) + 1)) -
      1/2 * pow(sd.theta, -2) * rho * sum(pow(x.from[1:NDist] - x.to[1:NDist], 2))
    if(log) return(logDens)
    else return(exp(logDens))
  }
)

# Another function of type rcar_leroux must be defined, otherwise it will cause an error
rcar_leroux <- nimbleFunction(
  name = 'rcar_leroux',
  run = function(n = integer(0),
                 rho = double(0),
                 sd.theta = double(0),
                 Lambda = double(1),
                 from.to = double(2)) {
    returnType(double(1))
    
    nimStop("user-defined distribution dcar_leroux provided without random generation function.")
    x <- nimNumeric(542)
    return(x)
  }
)

assign('dcar_leroux', dcar_leroux, envir = .GlobalEnv)
assign('rcar_leroux', rcar_leroux, envir = .GlobalEnv)


modelCode <- nimbleCode(
  {
    # Likelihood
    for (i in 1:N){
      O[i] ~ dpois(pop[i]*r[i])
      logit(r[i]) <- alpha + theta[i]
      
      MSS.r[i] <- (rate[i]-r[i])^2
      RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
      cond.var[i] <- sd.theta2/(rho*num[i] + (1-rho))
    }
    
    # Prior distributions
    alpha ~ dflat()
    
    # theta[1:NMuni] spatial random effect
    # LCAR distribution
    theta[1:N] ~ dcar_leroux(rho = rho,
                             sd.theta = sd.theta,
                             Lambda = Lambda[1:N],
                             from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    # sd.theta ~ dhalfflat()
    sd.theta2 ~ dunif (0, 100)
    sd.theta <- sqrt(sd.theta2)
    
    # Zero-mean constraint for theta[1:NMuni]
    zero.theta ~ dnorm(mean.thetas, 10000)
    mean.thetas <- mean(theta[1:N])
    
  }
)


set.seed(09052024)
constants <- list(N = S.area, NDist = NDist, Lambda = Lambda, 
                  from.to = from.to, num = nbInfo$num)


inits <- list(list(alpha = 0, sd.theta2 = runif(1), rho = runif(1), 
                   theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sd.theta2 = runif(1), rho = runif(1), 
                   theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sd.theta2 = runif(1), rho = runif(1), 
                   theta = rnorm(S.area, sd = 0.1)))


data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5,  zero.theta = 0)


LCAR.res <- nimbleMCMC(code = modelCode,
                       constants = constants,
                       data = data.nimble,
                       inits = inits,
                       nchains = 3,
                       niter = 300000,
                       nburnin = 50000,
                       thin = 75,
                       summary = TRUE,
                       samples = TRUE,
                       monitors = c('alpha', 'sd.theta2', 'rho',
                                    'r', 'MSS.r', 'theta',
                                    'RMSS.r', 'cond.var'),
                       samplesAsCodaMCMC = TRUE,
                       setSeed = TRUE,
                       WAIC = TRUE)



################################################################################
######BYM
################################################################################
# library(MASS)
# inv.condvar <- nimbleFunction(     
#   run = function(sigma.phi = double(0), sigma.iid = double(0), 
#                  adj = double(1), num = double(1), N = N) {
#     returnType(double(2))
#     Q.nb <- nimbleList(double(1))
#     pos <- 1
#     for (j in 1:N) {
#       # Q.nb[[j]] <- adj[(pos):(pos + num[j] - 1)]
#       pos <- pos + num[j]
#     }
#     class(Q.nb) <- "nb"
#     Q <- nb2mat(Q.nb, zero.policy = TRUE, style = "B")
#     Q <- -Q
#     diag(Q) <- abs(apply(Q, 1, sum))
#     inv.condvar <- ginv(ginv(as(Q, "matrix"))+ (sigma.iid/sigma.phi)^2*diag(N))
#     return(inv.condvar)
#   })
# inv.condvar <- compileNimble(inv.condvar)



code <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau.phi, zero_mean = 1) # its scaled so tau = 1
  # iid[1:N] ~ dmnorm(mean = rep(0, N), prec = tau.iid * diag(N))
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of theta
  sigma.phi2 ~ dunif (0, 100)
  tau.phi <- 1 / sigma.phi2                  # the variance of theta
  
  sigma.iid2 ~ dunif (0, 100)
  tau.iid <- 1 / sigma.iid2                  # the variance of theta
  
  # inv.condvar[1:N, 1:N] <- inv.condvar(sigma.phi = sigma.phi, sigma.iid = sigma.iid)
  # #Constraint
  # zero.phi ~ dnorm(phi.mean, 10000)
  # phi.mean <- mean(phi[1:N])
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + phi[i] + iid[i]
    
    iid[i] ~ dnorm(0, tau.iid) 
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    
    # cond.var[i] <- sigma.phi^2/diag(inv.condvar)[i]
  }
})

set.seed(09052024)
constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                  weights = nbInfo$weights, adj = nbInfo$adj)
inits <- list(list(alpha = 0, sigma.phi2 = runif(1), sigma.iid2 = runif(1),
                   iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)),
              list(alpha = 0, sigma.phi2 = runif(1), sigma.iid2 = runif(1),
                   iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)),
              list(alpha = 0, sigma.phi2 = runif(1), sigma.iid2 = runif(1),
                   iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5)

BYM.res <- nimbleMCMC(code = code,
                      constants = constants,
                      data = data.nimble,
                      inits = inits,
                      nchains = 3,
                      niter = 300000,
                      nburnin = 50000,
                      thin = 75,
                      summary = TRUE,
                      samples = TRUE,
                      monitors = c('alpha', 'sigma.phi2', 'tau.phi',
                                   'sigma.iid2', 'tau.iid', 'phi',
                                   'r', 'MSS.r', 'RMSS.r'),
                      samplesAsCodaMCMC = TRUE,
                      setSeed = TRUE,
                      WAIC = TRUE)



################################################################################
######BYM2
################################################################################
W.scale <- -W
diag(W.scale) <- abs(apply(W.scale, 1, sum))
# solve(W.scale) # this should not work since by definition the matrix is singular
library(INLA)
Q = inla.scale.model(W.scale, constr=list(A=matrix(1, nrow=1, ncol=nrow(W.scale)), e=0))
scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))

code <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau = 1, zero_mean = 1) # its scaled so tau = 1
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of the reparametrization
  sigma.b2 ~ dunif (0, 100)
  tau.b <- 1 / sigma.b2
  
  # precision parameter of theta
  # sigma.theta ~ dunif (0, 100)
  # tau.theta <- 1 / sigma.theta^2                       # the variance of theta
  
  # mixing parameter
  rho ~ dbeta(1, 1)                                    # prior for the mixing parameter
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + b[i]
    
    b[i] <- (1/sqrt(tau.b))*(sqrt((1-rho))*theta[i] + sqrt(rho/scale)*phi[i])
    # b[i] <- (1/sqrt(tau.b))*(sqrt((1-rho))*theta[i] + sqrt(rho)*phi[i])
    theta[i] ~ dnorm(0, tau = 1)               # area-specific RE
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    # cond.var[i] <- sigma.b^2*(sqrt(rho/scale)/num[i] + sqrt((1-rho)))
    
  }
})

set.seed(09052024)
constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                  weights = nbInfo$weights, adj = nbInfo$adj,
                  scale = scale)

inits <- list(list(alpha = 0, sigma.b2 = runif(1), theta = rnorm(S.area, sd = 0.1), 
                   phi = rnorm(S.area, sd = 0.1), rho = runif(1)),
              list(alpha = 0, sigma.b2 = runif(1), theta = rnorm(S.area, sd = 0.1), 
                   phi = rnorm(S.area, sd = 0.1), rho = runif(1)),
              list(alpha = 0, sigma.b2 = runif(1), theta = rnorm(S.area, sd = 0.1), 
                   phi = rnorm(S.area, sd = 0.1), rho = runif(1)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5)

BYM2.res <- nimbleMCMC(code = code,
                       constants = constants,
                       data = data.nimble,
                       inits = inits,
                       nchains = 3,
                       niter = 300000,
                       nburnin = 50000,
                       thin = 75,
                       summary = TRUE,
                       samples = TRUE,
                       monitors = c('alpha', 'sigma.b2', 'tau.b',
                                    'b', 'r', 'MSS.r', 'rho',
                                    'RMSS.r'),
                       samplesAsCodaMCMC = TRUE,
                       setSeed = TRUE,
                       WAIC = TRUE)




################################################################################
######GP
################################################################################
g <- st_as_sf(carto)
g_centroid <- st_point_on_surface(x = g)
# plot(st_geometry(g))
# plot(st_geometry(g_centroid), add = TRUE)
distMatrix <- st_distance(g_centroid, g_centroid)
distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid))
a <- max(distMatrix)

expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })
cExpcov <- compileNimble(expcov)


##GP
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  sigma2 ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  tau <- 1 / sigma2
  sigma <- sqrt(sigma2)
  rho ~ dunif(0, a)
  
  
  mu[1:N] <- rep(0, N)
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  theta[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    # cond.var[i] <- (sigma^2)*(exp(1/rho))
  }
})


set.seed(09052024)
constants <- list(N = S.area, dists = distMatrix, a = a)

inits <- list(list(alpha = 0, sigma2 = runif(1), rho = runif(1, 0, a), 
                   theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sigma2 = runif(1), rho = runif(1, 0, a), 
                   theta = rnorm(S.area, sd = 0.1)),
              list(alpha = 0, sigma2 = runif(1), rho = runif(1, 0, a), 
                   theta = rnorm(S.area, sd = 0.1)))

data.nimble <- list(O = data$O, pop = data$POB,
                    rate = data$crude.rate/10^5)

GP.res <- nimbleMCMC(code = code,
                     constants = constants,
                     data = data.nimble,
                     inits = inits,
                     nchains = 3,
                     niter = 300000,
                     nburnin = 50000,
                     thin = 75,
                     summary = TRUE,
                     samples = TRUE,
                     monitors = c('alpha', 'sigma2', 'rho',
                                  'r', 'MSS.r', 'theta',
                                  'RMSS.r'),
                     samplesAsCodaMCMC = TRUE,
                     setSeed = TRUE,
                     WAIC = TRUE)



results <- list(iid.res = iid.res,
                iCAR.res = iCAR.res,
                pCAR.res = pCAR.res,
                LCAR.res = LCAR.res,
                BYM.res = BYM.res,
                BYM2.res = BYM2.res,
                GP.res = GP.res)

save(results, file = paste0("./Results_Lung_cancer_Spain2021_new.Rdata"))


