################################################################################
################################################################################
##########             4.2 Across priors simulation study             ##########
##########                  Code to fit the iid prior                 ##########
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
if (!file.exists("./SimulationStudy_iid")){
  dir.create("./SimulationStudy_iid")
}



#################################################
###    NIMBLE Code for iid  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##iid
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  
  sigma ~ dunif (0, 100)
  tau <- 1 / sigma^2
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + iid[i]
    
    iid[i] ~ dnorm(0, tau)
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
  }
})




################################################################################
##########                      Peninsular Spain                      ##########
################################################################################
#################################################
###    Required Constants   ###
#################################################
### Scenarios
hotspot <- c("Scenario1", "Scenario2", "Scenario3")
### Number of areas
n.areas <- c("47", "100", "300")



#################################################
###   Code for all nº of areas   ###
#################################################

for (na in 1:length(n.areas)) {
  for (h in 1:length(hotspot)) {
    
    #################################################
    ###    Load the simulated data   ###
    #################################################
    load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_",n.areas[na],"areas.Rdata"))
    S.area <- length(unique(DataSIM$code))
    
    
    n.sim <- 1000
    act.sim <- 1
    c.save <- 0
    
    iid.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
      
      ##########################################################################
      ##########################################################################
      constants <- list(pop = data$population,
                        rate = data$crude.rate/10^5,
                        N = S.area)
      
      inits = function(){
        list(alpha = 0, sigma = runif(1), iid = rnorm(S.area, sd = 0.1))}
      
      data.nimble <- list(O = data$observed)
      
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
                             monitors = c('r', 'MSS.r', 'RMSS.r'),
                             samplesAsCodaMCMC = TRUE,
                             setSeed = c(20112023, 54782021, 04062025),
                             WAIC = TRUE)
      
      iid.res[act.sim-c.save] <- list(mcmc.out)
      
      #################################################
      ###    Save results each 25 simulations   ###
      #################################################
      if (act.sim%%25 == 0){
        save(list = c("iid.res"),
             file = paste0("./SimulationStudy_iid/Results_SimulationStudy_",n.areas[na],"_",hotspot[h],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        iid.res <- list()
      }
      
      act.sim <- act.sim +1
      
      if(act.sim==n.sim+1){break}
    }
  }
}




################################################################################
##########                           England                          ##########
################################################################################
#################################################
###    Required Constants   ###
#################################################
### Scenarios
hotspot <- c("Scenario1", "Scenario2", "Scenario3")


for (h in 1:length(hotspot)) {
    
  #################################################
  ###    Load the simulated data   ###
  #################################################
  load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_England.Rdata"))
  S.area <- length(unique(DataSIM$code))
    
    
  n.sim <- 1000
  act.sim <- 1
  c.save <- 0
    
  iid.res <- list()
  repeat{
    print(act.sim)
    data <- DataSIM[which(DataSIM$sim==act.sim),]
    
    ##########################################################################
    ##########################################################################
    constants <- list(pop = data$population,
                      rate = data$crude.rate/10^5,
                      N = S.area)
      
    inits = function(){
      list(alpha = 0, sigma = runif(1), iid = rnorm(S.area, sd = 0.1))}
    
    data.nimble <- list(O = data$observed)
    
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
                           monitors = c('r', 'MSS.r', 'RMSS.r'),
                           samplesAsCodaMCMC = TRUE,
                           setSeed = c(20112023, 54782021, 04062025),
                           WAIC = TRUE)
      
    iid.res[act.sim-c.save] <- list(mcmc.out)
    
    #################################################
    ###    Save results each 25 simulations   ###
    #################################################
    if (act.sim%%25 == 0){
      save(list = c("iid.res"),
           file = paste0("./SimulationStudy_iid/Results_SimulationStudy_England_",hotspot[h],"_",act.sim,".Rdata"))
      
      c.save<-c.save+25
      iid.res <- list()
    }
      
    act.sim <- act.sim +1
    
    if(act.sim==n.sim+1){break}
  }
}













  

