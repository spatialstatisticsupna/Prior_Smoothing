################################################################################
################################################################################
##########              4.1 Within prior simulation study             ##########
##########                  Code to fit the GP prior                  ##########
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
if (!file.exists("./GP")){
  dir.create("./GP")
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
rho.value <- c(1, 5, 9)

#################################################
###    NIMBLE Code for GP  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)


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
  
  sigma <- fix.sd 
  tau <- 1 / sigma^2
  rho <- fix.rho
  
  
  mu[1:N] <- rep(0, N)
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  theta[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  
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
  load(paste0("../../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  g <- st_as_sf(carto)
  g_centroid <- st_point_on_surface(x = g)
  distMatrix <- st_distance(g_centroid, g_centroid)
  ###Taking into account the distances between areas, we decided to scale the distance to units of 
  ###100km.Since the distMatrix is originally in meters, we divided all values by 100,000, 
  ###resulting in a maximum distance sligthly above 10 (i.e., just over 1,000 km).
  distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid)) 
   
  
  #################################################
  ###    Load the simulated data   ###
  #################################################
  load(paste0("../../Data/Data_SimulationStudy_",hotspot,"_",n.areas[na],"areas.Rdata"))
  S.area <- length(unique(DataSIM$code))
  
  
  for (l in 1:length(rho.value)) {
    for (sd in 1:length(sd.value)) {
      n.sim <- 1000
      act.sim <- 1
      c.save <- 0
      
      GP.res <- list()
      repeat{
        print(act.sim)
        data <- DataSIM[which(DataSIM$sim==act.sim),]
        
        ##########################################################################
        ##########################################################################
        ####GP
        constants <- list(pop = data$population,
                          rate = data$crude.rate/10^5,
                          N = S.area, 
                          dists = distMatrix,
                          fix.sd = sd.value[sd],
                          fix.rho = rho.value[l])
        
        inits = function(){
          list(alpha = 0, theta = rnorm(S.area, sd = 0.1))}
        
        data.nimble <- list(O = data$observed)
        
        # mcmc.out <- nimbleMCMC(code = code,
        #                        constants = constants,
        #                        data = data.nimble,
        #                        inits = inits,
        #                        nchains = 3,
        #                        niter = 30000,
        #                        nburnin = 5000,
        #                        thin = 75,
        #                        summary = TRUE,
        #                        samples = TRUE,
        #                        monitors = c('r', 'MSS.r', 'RMSS.r'),
        #                        samplesAsCodaMCMC = TRUE,
        #                        setSeed = c(20112023, 54782021, 04062025),
        #                        WAIC = TRUE)
        
        nimModel <- nimbleModel(code = code, data = data.nimble,
                                constants = constants, inits = inits())
        compileNimble(nimModel)
        modBuilt <- buildMCMC(nimModel, monitors = c('r', 'MSS.r', 'RMSS.r'),
                              setSeed = TRUE)
        modComp <- compileNimble(modBuilt, project = nimModel)
        mcmc.out <- runMCMC(modComp, nchains = 3, niter = 30000, thin = 75,
                            nburnin = 5000, summary = TRUE)
        
        GP.res[act.sim-c.save] <- list(mcmc.out)
        
        #################################################
        ###    Save results each 25 simulations   ###
        #################################################
        if (act.sim%%25 == 0){
          save(list = c("GP.res"),
               file = paste0("./GP/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",rho.value[l],"_",act.sim,".Rdata"))
          
          c.save<-c.save+25
          GP.res <- list()
        }
        
        act.sim <- act.sim +1
        
        if(act.sim==n.sim+1){break}
      }
    }
  }
}

