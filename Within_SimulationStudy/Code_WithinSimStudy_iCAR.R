################################################################################
################################################################################
##########              4.1 Within prior simulation study             ##########
##########                 Code to fit the iCAR prior                 ##########
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
if (!file.exists("./iCAR")){
  dir.create("./iCAR")
}

#################################################
###    Required Constants   ###
#################################################
### Used scenario
hotspot <- c("Scenario2")
### Number of areas
n.areas <- c("47", "100", "300")
### Fix values for hyperparameters
sd.value <- c(0.01, 0.05, 0.09, 0.2, 0.5)


#################################################
###    NIMBLE Code for iCAR  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##iCAR
code <- nimbleCode({
  
  # Priors:
  # intercept
  alpha ~ dflat()                               # vague uniform prior
  
  # ICAR
  theta[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 1)
  
  # precision parameter of theta
  sigma <- fix.sd                               # fixed
  tau <- 1 / sigma^2
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + theta[i]
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})



#####################################################################
###   Code for all nº of areas and fixed hyperparameters values   ###
####################################################################

for (na in 1:length(n.areas)) {
  
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  #################################################
  ###    Load the simulated data   ###
  #################################################
  load(paste0("../Data/Data_SimulationStudy_",hotspot,"_",n.areas[na],"areas.Rdata"))
  S.area <- length(unique(DataSIM$code))
  
  for (sd in 1:length(sd.value)) {
    n.sim <- 1000
    act.sim <- 1
    c.save <- 0
    
    iCAR.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
      
      ##########################################################################
      ##########################################################################
      ####iCAR
      constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                        weights = nbInfo$weights, adj = nbInfo$adj, fix.sd = sd.value[sd])
      
      inits <- list(list(alpha = 0, theta = rnorm(S.area, sd = 0.1)),
                    list(alpha = 0, theta = rnorm(S.area, sd = 0.1)),
                    list(alpha = 0, theta = rnorm(S.area, sd = 0.1)))
      
      data.nimble <- list(O = data$observed, pop = data$population,
                          rate = data$crude.rate/10^5)
      
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
                             monitors = c('alpha', 'sigma', 'tau',
                                          'r', 'MSS.r', 'theta',
                                          'RMSS.r'),
                             samplesAsCodaMCMC = TRUE,
                             setSeed = c(20112023, 54782021, 04062025),
                             WAIC = TRUE)
      
      
      iCAR.res[act.sim-c.save] <- list(mcmc.out)
      
      #################################################
      ###    Save results each 25 simulations   ###
      #################################################
      if (act.sim%%25 == 0){
        save(list = c("iCAR.res"),
             file = paste0("./iCAR/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        iCAR.res <- list()
      }
      
      act.sim <- act.sim +1
      
      if(act.sim==n.sim+1){break}
    }
  }
}






  

