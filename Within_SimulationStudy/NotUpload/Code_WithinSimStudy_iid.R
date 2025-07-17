################################################################################
################################################################################
##########                     Simulation Study                       ##########
################################################################################
################################################################################
library(spdep)
library(MASS)
library(sf)

if (!file.exists("./iid")){
  dir.create("./iid")
}

#################################################
## Load the cartography file  ##
#################################################
hotspot <- c("Scenario2")
n.areas <- c("47", "100", "300")
sd.value <- c(0.01, 0.05, 0.09, 0.2, 0.5)

library(nimble)
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##iid
code <- nimbleCode({
  
  # priors:
  alpha ~ dflat()
  overallRR <- exp(alpha)
  
  # iid[1:N] ~ dmnorm(mean = rep(0, N), prec = tau * diag(N))
  sigma  <- fix.sd
  tau <- 1 / sigma^2
  
  # likelihood:
  for(i in 1:N){
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + iid[i]
    
    iid[i] ~ dnorm(0, tau)
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    TCV[i] <- sigma^2
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
    
    iCAR.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
      
      ##########################################################################
      ##########################################################################
      constants <- list(N = S.area, fix.sd = sd.value[sd])
      inits <- list(list(alpha = 0, iid = rnorm(S.area)),
                    list(alpha = 0, iid = rnorm(S.area)),
                    list(alpha = 0, iid = rnorm(S.area)))
      
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
                                          'r', 'MSS.r', 'iid',
                                          'RMSS.r', 'TCV'),
                             samplesAsCodaMCMC = TRUE,
                             setSeed = 20112023,
                             WAIC = TRUE)
      
      iid.res[act.sim-c.save] <- list(mcmc.out)
      
      if (act.sim%%25 == 0){
        save(list = c("iid.res"),
             file = paste0("./iid/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        iCAR.res <- list()
      }
      
      act.sim <- act.sim +1
      
      if(act.sim==n.sim+1){break}
    }
    
  }
  
}






  

