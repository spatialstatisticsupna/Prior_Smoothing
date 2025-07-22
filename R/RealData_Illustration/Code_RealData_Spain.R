################################################################################
################################################################################
##########                    5. Data Illustration                    ##########
##########                   Code to fit Spain data                   ##########
################################################################################
################################################################################

#################################################
###    Required libraries   ###
#################################################
library(spdep)
library(MASS)
library(sf)
library(nimble)

#################################################
###    Define folder to save results   ###
#################################################
if (!file.exists("./Results")){
  dir.create("./Results")
}


#################################################
###    Required Constants   ###
#################################################
### Number of areas
n.areas <- c("47", "100", "300")
### Priors
prior <- c("iid", "GP", "iCAR", "BYM", "pCAR", "LCAR","BYM2")

### name for each case
name.s <- c("small", "medium", "large", "non")
### min values for the uniform distribution
min.unif <- c(0, 0.01, 0.16, 0)
### max values for the uniform distribution
max.unif <- c(0.01, 0.16, 100, 100)


#################################################
###    NIMBLE Codes   ###
#################################################
########################
###iid
########################
code.iid <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  
  sigma ~ dunif (min_unif, max_unif)
  tau <- 1 / sigma
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + iid[i]
    
    iid[i] ~ dnorm(0, tau) 
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})


########################
###GP
########################
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
code.GP <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  
  sigma ~ dunif(min_unif, max_unif) 
  tau <- 1 / sigma
  rho ~ dunif(0, a)
  
  
  mu[1:N] <- rep(0, N)
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sqrt(sigma))
  theta[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})



########################
###iCAR
########################
code.iCAR <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  
  theta[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 1)
  sigma ~ dunif (min_unif, max_unif)
  tau <- 1 / sigma
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})


########################
###pCAR
########################
code.pCAR <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  
  theta[1:N] ~ dcar_proper(mu[1:N], C[1:L], adj[1:L], num[1:N], M[1:N], tau, gamma)
  #M: is a diagonal matrix of conditional variances
  #C: Spatial proximity matrix. ??normalized weights!!  (each row sums to one)
  sigma ~ dunif (min_unif, max_unif)
  tau <- 1 / sigma
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
  }
})


########################
###BYM
########################
code.BYM <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau.phi, zero_mean = 1) 
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of theta
  sigma.phi ~ dunif (min_unif, max_unif)
  tau.phi <- 1 / sigma.phi                  
  
  sigma.iid ~ dunif (min_unif, max_unif)
  tau.iid <- 1 / sigma.iid                  
  
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + phi[i] + iid[i]
    
    iid[i] ~ dnorm(0, tau.iid) 
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})



########################
###LCAR
########################
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


code.LCAR <- nimbleCode(
  {
    # Likelihood
    for (i in 1:N){
      O[i] ~ dpois(pop[i]*r[i])
      logit(r[i]) <- alpha + theta[i]
      
      MSS.r[i] <- (rate[i]-r[i])^2
      RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    }
    
    # Prior distributions
    alpha ~ dflat()
    
    # theta[1:NMuni] spatial random effect
    # LCAR distribution
    theta[1:N] ~ dcar_leroux(rho = rho,
                             sd.theta = sqrt(sd.theta),
                             Lambda = Lambda[1:N],
                             from.to = from.to[1:NDist, 1:2])
    
    # Hyperparameters of the spatial random effect
    rho ~ dunif(0, 1)
    # sd.theta ~ dhalfflat()
    sd.theta ~ dunif (min_unif, max_unif)
    
    # Zero-mean constraint for theta[1:NMuni]
    zero.theta ~ dnorm(mean.thetas, 10000)
    mean.thetas <- mean(theta[1:N])
    
  }
)


########################
###BYM2
########################
code.BYM2 <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau = 1, zero_mean = 1) # its scaled so tau = 1
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of the reparametrization
  sigma.b ~ dunif (min_unif, max_unif)
  tau.b <- 1 / sigma.b
  
  
  # mixing parameter
  rho ~ dbeta(1, 1)                                  # prior for the mixing parameter
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + b[i]
    
    b[i] <- (1/sqrt(tau.b))*(sqrt((1-rho))*theta[i] + sqrt(rho/scale)*phi[i])
    theta[i] ~ dnorm(0, tau = 1)               # area-specific RE
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})






for (mu in 1:length(name.s)) {
  ###
  source("./Function_constants.R")
  ###
  source("./Function_inits.R")
  
  
  for (na in 1:length(n.areas)) {
    #################################################
    ###    Load the cartography   ###
    #################################################
    load(paste0("../../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
    carto <- Carto.areas
  
    #################################################
    ###    Load the data   ###
    #################################################
    load(paste0("../../Data/Data_Spain_",n.areas[na],"areas.Rdata"))
    Data.areas$crude.rate <- Data.areas$O/Data.areas$Pop*10^5
    
    
    for (pr in 1:length(prior)) {
      
      constants <- f.constants(prior = paste(prior[pr]),
                               data = Data.areas,
                               min_u = min.unif[mu],
                               max_u = max.unif[mu],
                               carto = carto)
      
      
      data.nimble <- list(O = Data.areas$O)
      if(prior[pr]%in%c("LCAR", "pCAR")){
        data.nimble <- append(data.nimble, list(zero.theta = 0))
      }
      
      
      
      inits <- f.inits(prior = paste(prior[pr]),
                       carto = carto,
                       min.unif =  min.unif[mu],
                       max.unif = max.unif[mu])
      mon <- names(inits())
      
      eval(parse(text = paste0(prior[pr],".res <- nimbleMCMC(code = code.",prior[pr],",
                             constants = constants,
                             data = data.nimble,
                             inits = inits(),
                             nchains = 3,
                             niter = 30000,
                             nburnin = 5000,
                             thin = 75,
                             summary = TRUE,
                             samples = TRUE,
                             monitors = c('r', 'MSS.r', 'RMSS.r', mon),
                             samplesAsCodaMCMC = TRUE,
                             setSeed = c(20112023, 54782021, 04062025),
                             WAIC = TRUE)")))
      
      rm(list = c("constants", "data.nimble", "inits"))
      
    }
    
    results <- list(iid.res = iid.res,
                    iCAR.res = iCAR.res,
                    pCAR.res = pCAR.res,
                    LCAR.res = LCAR.res,
                    BYM.res = BYM.res,
                    BYM2.res = BYM2.res,
                    GP.res = GP.res)
    
    save(results, file = paste0("./Results/Results_Spain_",n.areas[na],"_",name.s[mu],".Rdata"))
    
    rm(list = c("results", "iid.res", "iCAR.res",
                "pCAR.res", "LCAR.res",
                "BYM.res", "BYM2.res",
                "GP.res"))
  }
}





