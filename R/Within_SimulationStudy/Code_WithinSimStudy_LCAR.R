################################################################################
################################################################################
##########              4.1 Within prior simulation study             ##########
##########                 Code to fit the LCAR prior                 ##########
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
if (!file.exists("./LCAR")){
  dir.create("./LCAR")
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
lambda.value <- c(0.1, 0.5, 0.9)


#################################################
###    NIMBLE Code for LCAR  ###
#################################################
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)


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
    rho <- fix.lambda
    # sd.theta ~ dhalfflat()
    sd.theta <- fix.sd
    
    # Zero-mean constraint for theta[1:NMuni]
    zero.theta ~ dnorm(mean.thetas, 10000)
    mean.thetas <- mean(theta[1:N])
    
  }
)


#####################################################################
###   Code for all nº of areas and fixed hyperparameters values   ###
####################################################################

for (na in 1:length(n.areas)) {
  
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  S.area <- length(unique(carto$ID))
  

  cv.nb <- poly2nb(carto)
  W.nb <- nb2mat(cv.nb, style="B")
  nbInfo <- nb2WB(cv.nb)
  # Number of neighbors of each municipality
  nadj <- card(cv.nb)
  # Neighbors of each municipality
  map <- unlist(cv.nb)
  
  
  ##LCAR
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
  
  
  #################################################
  ###    Load the simulated data   ###
  #################################################
  load(paste0("../../Data/Data_SimulationStudy_",hotspot,"_",n.areas[na],"areas.Rdata"))
  S.area <- length(unique(DataSIM$code))
  
  for (l in 1:length(lambda.value)) {
    for (sd in 1:length(sd.value)) {
      n.sim <- 1000
      act.sim <- 1
      c.save <- 0
      
      LCAR.res <- list()
      repeat{
        print(act.sim)
        data <- DataSIM[which(DataSIM$sim==act.sim),]
        
        ##########################################################################
        ##########################################################################
        ####L-CAR
        constants <- list(pop = data$population,
                          rate = data$crude.rate/10^5,
                          N = S.area, 
                          NDist = NDist, 
                          Lambda = Lambda, 
                          from.to = from.to, 
                          fix.sd = sd.value[sd],
                          fix.lambda = lambda.value[l])
        
        
        inits = function(){
          list(alpha = 0, theta = rnorm(S.area, sd = 0.1))}
        
        
        data.nimble <- list(O = data$observed, zero.theta = 0)
        
        # mcmc.out <- nimbleMCMC(code = modelCode,
        #                        constants = constants,
        #                        data = data.nimble,
        #                        inits = inits(),
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
        
        nimModel <- nimbleModel(code = modelCode, data = data.nimble,
                                constants = constants, inits = inits())
        compileNimble(nimModel)
        modBuilt <- buildMCMC(nimModel, monitors = c('r', 'MSS.r', 'RMSS.r'),
                              setSeed = TRUE)
        modComp <- compileNimble(modBuilt, project = nimModel)
        mcmc.out <- runMCMC(modComp, nchains = 3, niter = 30000, thin = 75,
                                                 nburnin = 5000, summary = TRUE)
        
        LCAR.res[act.sim-c.save] <- list(mcmc.out)
        
        #################################################
        ###    Save results each 25 simulations   ###
        #################################################
        if (act.sim%%25 == 0){
          save(list = c("LCAR.res"),
               file = paste0("./LCAR/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",lambda.value[l],"_",act.sim,".Rdata"))
          
          c.save<-c.save+25
          LCAR.res <- list()
        }
        
        act.sim <- act.sim +1
        
        if(act.sim==n.sim+1){break}
      }
    }
  }
}






  

