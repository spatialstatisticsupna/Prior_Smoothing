################################################################################
################################################################################
##########             4.2 Across priors simulation study             ##########
##########                 Code to fit the BYM prior                  ##########
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
if (!file.exists("./SimulationStudy_BYM")){
  dir.create("./SimulationStudy_BYM")
}

#################################################
###    NIMBLE Code for BYM  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##BYM
code <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau.phi, zero_mean = 1) 
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of theta
  sigma.phi ~ dunif (0, 100)
  tau.phi <- 1 / sigma.phi^2                  
  
  sigma.iid ~ dunif (0, 100)
  tau.iid <- 1 / sigma.iid^2                  
  
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + phi[i] + iid[i]
    
    iid[i] ~ dnorm(0, tau.iid) 
    
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




#####################################################################
###   Code for all nº of areas  ###
####################################################################

for (na in 1:length(n.areas)) {
  
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  for (h in 1:length(hotspot)) {
    #################################################
    ###    Load the simulated data   ###
    #################################################
    load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_",n.areas[na],"areas.Rdata"))
    S.area <- length(unique(DataSIM$code))
    
    n.sim <- 1000
    act.sim <- 1
    c.save <- 0
      
    BYM.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
        
      ##########################################################################
      ##########################################################################
      ####BYM
      constants <- list(pop = data$population,
                        rate = data$crude.rate/10^5,
                        N = S.area, 
                        L = length(nbInfo$adj), 
                        num = nbInfo$num,
                        weights = nbInfo$weights, 
                        adj = nbInfo$adj)
        
      inits = function(){
        list(alpha = 0, sigma.phi = runif(1), sigma.iid = runif(1),
             iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1))}
      
      data.nimble <- list(O = data$observed)
        
      mcmc.out <- nimbleMCMC(code = code,
                             constants = constants,
                             data = data.nimble,
                             inits = inits(),
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

        
      BYM.res[act.sim-c.save] <- list(mcmc.out)
      
      #################################################
      ###    Save results each 25 simulations   ###
      #################################################
      if (act.sim%%25 == 0){
        save(list = c("BYM.res"),
             file = paste0("./SimulationStudy_BYM/Results_SimulationStudy_",n.areas[na],"_",hotspot[h],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        BYM.res <- list()
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


#################################################
###    Load the cartography   ###
#################################################
carto <- st_read("../../Data/Carto_England/carto_england.shp")
carto <- carto[order(carto$Code), ]

carto.nb <- poly2nb(carto)
W.nb <- nb2mat(carto.nb, style="B")
nbInfo <- nb2WB(carto.nb)
  
for (h in 1:length(hotspot)) {
  #################################################
  ###    Load the simulated data   ###
  #################################################
  load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_England.Rdata"))
  S.area <- length(unique(DataSIM$code))
    
  n.sim <- 1000
  act.sim <- 1
  c.save <- 0
    
  BYM.res <- list()
  repeat{
    print(act.sim)
    data <- DataSIM[which(DataSIM$sim==act.sim),]
    
    ##########################################################################
    ##########################################################################
    ####BYM
    constants <- list(pop = data$population,
                      rate = data$crude.rate/10^5,
                      N = S.area, 
                      L = length(nbInfo$adj), 
                      num = nbInfo$num,
                      weights = nbInfo$weights, 
                      adj = nbInfo$adj)
    
    inits = function(){
      list(alpha = 0, sigma.phi = runif(1), sigma.iid = runif(1),
           iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1))}
    
    data.nimble <- list(O = data$observed)
    
    mcmc.out <- nimbleMCMC(code = code,
                           constants = constants,
                           data = data.nimble,
                           inits = inits(),
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
    
      
    BYM.res[act.sim-c.save] <- list(mcmc.out)
    
    #################################################
    ###    Save results each 25 simulations   ###
    #################################################
    if (act.sim%%25 == 0){
      save(list = c("BYM.res"),
           file = paste0("./SimulationStudy_BYM/Results_SimulationStudy_England_",hotspot[h],"_",act.sim,".Rdata"))
      
      c.save<-c.save+25
      BYM.res <- list()
    }
      
    act.sim <- act.sim +1
      
    if(act.sim==n.sim+1){break}
  }
}



