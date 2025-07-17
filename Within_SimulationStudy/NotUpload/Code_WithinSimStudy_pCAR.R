################################################################################
################################################################################
##########                     Simulation Study                       ##########
################################################################################
################################################################################
library(spdep)
library(MASS)
library(sf)

if (!file.exists("./pCAR")){
  dir.create("./pCAR")
}

#################################################
## Load the cartography file  ##
#################################################
hotspot <- c("Scenario2")
n.areas <- c("47", "100", "300")
sd.value <- c(0.01, 0.05, 0.09, 0.2, 0.5)

library(nimble)
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##p-CAR
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  theta[1:N] ~ dcar_proper(mu[1:N], C[1:L], adj[1:L], num[1:N], M[1:N], tau, gamma)
  sigma <- fix.sd
  tau <- 1 / sigma^2
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
    TCV[i] <- sigma^2/(num[i])
  }
})

for (na in 1:length(n.areas)) {
  load(paste0("../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  
  # sf::sf_use_s2(FALSE)
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  load(paste0("../Data/Data_SimulationStudy_",hotspot,"_",n.areas[na],"areas.Rdata"))
  S.area <- length(unique(DataSIM$code))
  
  for (sd in 1:length(sd.value)) {
    n.sim <- 1000
    act.sim <- 1
    c.save <- 0
    
    pCAR.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
      
      ##########################################################################
      ##########################################################################
      ####p-CAR
      CM <- as.carCM(num = nbInfo$num, weights = nbInfo$weights, adj = nbInfo$adj) #To create C and M.
      constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                        adj = nbInfo$adj, C = CM$C, M = CM$M, fix.sd = sd.value[sd])
      
      inits <- list(list(alpha = 0, theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)),
                    list(alpha = 0, theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)),
                    list(alpha = 0, theta = rnorm(S.area, sd = 0.1), gamma = runif(1, -1, 1)))
      
      data.nimble <- list(O = data$observed, pop = data$population,
                          rate = data$crude.rate/10^5,  zero.theta = 0)
      
      
      mcmc.out <- nimbleMCMC(code = code,
                             constants = constants,
                             data = data.nimble,
                             inits = inits,
                             nchains = 3,
                             niter = 30000,
                             nburnin = 5000,
                             thin = 75,
                             summary = TRUE,
                             samples = TRUE,
                             monitors = c('alpha', 'sigma', 'tau', 'gamma',
                                          'r', 'MSS.r', 'theta',
                                          'RMSS.r', 'TCV'),
                             samplesAsCodaMCMC = TRUE,
                             setSeed = c(20112023, 54782021, 04062025),
                             WAIC = TRUE)
      
      
      
      pCAR.res[act.sim-c.save] <- list(mcmc.out)
      
      if (act.sim%%25 == 0){
        save(list = c("pCAR.res"),
             file = paste0("./pCAR/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        pCAR.res <- list()
      }
      
      act.sim <- act.sim +1
      
      if(act.sim==n.sim+1){break}
    }
    
  }
  
}








  

