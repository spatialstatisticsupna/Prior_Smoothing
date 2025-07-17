#################################################################################
################################################################################
##########              4.1 Within prior simulation study             ##########
##########                  Code to fit the BYM prior                 ##########
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
if (!file.exists("./BYM")){
  dir.create("./BYM")
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
v.value <- c(0.25, 1, 4)

#################################################
###    NIMBLE Code for iCAR  ###
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
  sigma.phi <- fix.sd
  tau.phi <- 1 / sigma.phi^2                  # fixed
  
  sigma.iid <- fix.tau
  tau.iid <- 1 / sigma.iid^2                  # fixed
  
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + phi[i] + iid[i]
    
    iid[i] ~ dnorm(0, tau.iid) 
    
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
  
  for (v in 1:length(v.value)) {
    for (sd in 1:length(sd.value)) {
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
        constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                          weights = nbInfo$weights, adj = nbInfo$adj, fix.sd = sd.value[sd],
                          fix.tau = sqrt(v.value[v])*sd.value[sd])
        inits <- list(list(alpha = 0, 
                           iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)),
                      list(alpha = 0, 
                           iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)),
                      list(alpha = 0,
                           iid = rnorm(S.area, sd = 0.1), phi = rnorm(S.area, sd=0.1)))
        
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
                               monitors = c('alpha', 'sigma.phi', 'tau.phi',
                                            'sigma.iid', 'tau.iid', 'phi',
                                            'r', 'MSS.r', 'RMSS.r'),
                               samplesAsCodaMCMC = TRUE,
                               setSeed = c(20112023, 54782021, 04062025),
                               WAIC = TRUE)
        
        BYM.res[act.sim-c.save] <- list(mcmc.out)
        
        #################################################
        ###    Save results each 25 simulations   ###
        #################################################
        if (act.sim%%25 == 0){
          save(list = c("BYM.res"),
               file = paste0("./BYM/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",v.value[v],"_",act.sim,".Rdata"))
          
          c.save<-c.save+25
          BYM.res <- list()
        }
        
        act.sim <- act.sim +1
        
        if(act.sim==n.sim+1){break}
      }
    }
  }
}






  

