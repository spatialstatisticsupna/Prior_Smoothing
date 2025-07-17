################################################################################
################################################################################
##########             4.2 Across priors simulation study             ##########
##########                 Code to fit the BYM2 prior                 ##########
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
if (!file.exists("./SimulationStudy_BYM2")){
  dir.create("./SimulationStudy_BYM2")
}


#################################################
###    Required Constants   ###
#################################################
### Scenarios
hotspot <- c("Scenario1", "Scenario2", "Scenario3")
### Number of areas
n.areas <- c("47", "100", "300")


#################################################
###    NIMBLE Code for BYM2  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)


##BYM2
code <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau = 1, zero_mean = 1) # its scaled so tau = 1
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of the reparametrization
  sigma.b ~ dunif (0, 100)
  tau.b <- 1 / sigma.b^2
  
  
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


#################################################
###   Code for all nº of areas   ###
#################################################

for (na in 1:length(n.areas)) {
  
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  carto.nb <- poly2nb(carto)
  W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
  nbInfo <- nb2WB(carto.nb)
  
  W.scale <- -W
  diag(W.scale) <- abs(apply(W.scale, 1, sum))
  Q <- W.scale*exp(mean(log(diag(MASS::ginv(as(W.scale,"matrix"))))))
  scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))
  
  
  for (h in 1:length(hotspot)) {
    
    #################################################
    ###    Load the simulated data   ###
    #################################################
    load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_",n.areas[na],"areas.Rdata"))
    S.area <- length(unique(DataSIM$code))
    
    
    n.sim <- 1000
    act.sim <- 1
    c.save <- 0
    
  
    BYM2.res <- list()
    repeat{
      print(act.sim)
      data <- DataSIM[which(DataSIM$sim==act.sim),]
      
      ##########################################################################
      ##########################################################################
      ####BYM2
      constants <- list(pop = data$population,
                        rate = data$crude.rate/10^5,
                        N = S.area, 
                        L = length(nbInfo$adj), 
                        num = nbInfo$num,
                        weights = nbInfo$weights, 
                        adj = nbInfo$adj,
                        scale = scale)
        
      inits = function(){
        list(alpha = 0, sigma.b = runif(1), theta = rnorm(S.area, sd = 0.1), 
             phi = rnorm(S.area, sd = 0.1), rho = runif(1))}
        
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
        
      BYM2.res[act.sim-c.save] <- list(mcmc.out)
      
      if (act.sim%%25 == 0){
        save(list = c("BYM2.res"),
             file = paste0("./SimulationStudy_BYM2/Results_SimulationStudy_",n.areas[na],"_",hotspot[h],"_",act.sim,".Rdata"))
        
        c.save<-c.save+25
        BYM2.res <- list()
      }
        
      act.sim <- act.sim +1
      
      if(act.sim==n.sim+1){break}
    }
  }
}







  

