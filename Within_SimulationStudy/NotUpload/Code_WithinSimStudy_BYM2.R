################################################################################
################################################################################
##########                     Simulation Study                       ##########
################################################################################
################################################################################
library(spdep)
library(MASS)
library(sf)

if (!file.exists("./BYM2")){
  dir.create("./BYM2")
}

#################################################
## Load the cartography file  ##
#################################################
hotspot <- c("Scenario2")
n.areas <- c("47", "100", "300")
sd.value <- c(0.01, 0.05, 0.09, 0.2, 0.5)
lambda.value <- c(0.1, 0.5, 0.9)

library(nimble)
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)


##BYM2
code <- nimbleCode({
  
  # Priors:
  # ICAR
  phi[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau = 1, zero_mean = 1) # its scaled so tau = 1
  
  # intercept
  alpha ~ dflat()                                      # vague uniform prior
  
  # precision parameter of the reparametrization
  sigma.b <- fix.sd
  tau.b <- 1 / sigma.b^2
  
  
  # mixing parameter
  rho <- fix.lambda                                  # prior for the mixing parameter
  
  # likelihood:
  for (i in 1:N){
    
    O[i] ~ dpois(pop[i]*r[i])
    logit(r[i]) <- alpha + b[i]
    
    b[i] <- (1/sqrt(tau.b))*(sqrt((1-rho))*theta[i] + sqrt(rho/scale)*phi[i])
    # b[i] <- (1/sqrt(tau.b))*(sqrt((1-rho))*theta[i] + sqrt(rho)*phi[i])
    theta[i] ~ dnorm(0, tau = 1)               # area-specific RE
    
    MSS.r[i] <- (rate[i]-r[i])^2
    RMSS.r[i] <- (rate[i]-r[i])^2/r[i]
    TCV[i] <- sigma.b^2/(sqrt(rho/scale)*num[i] + sqrt((1-rho)))
    # cond.var2[i] <- sigma.b^2*(1/(sqrt(rho/scale)*num[i]) + 1/sqrt(1-rho))
  }
})


for (na in 1:length(n.areas)) {
  load(paste0("../Data/Carto_Spain_",n.areas[na],"areas.Rdata"))
  carto <- Carto.areas
  
  # sf::sf_use_s2(FALSE)
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  W.scale <- -W
  diag(W.scale) <- abs(apply(W.scale, 1, sum))
  # solve(W.scale) # this should not work since by definition the matrix is singular
  library(INLA)
  Q = inla.scale.model(W.scale, constr=list(A=matrix(1, nrow=1, ncol=nrow(W.scale)), e=0))
  scale = exp((1/nrow(W.scale))*sum(log(1/diag(Q))))
  
  
  load(paste0("../Data/Data_SimulationStudy_",hotspot,"_",n.areas[na],"areas.Rdata"))
  S.area <- length(unique(DataSIM$code))
  
  for (l in 1:length(lambda.value)) {
    for (sd in 1:length(sd.value)) {
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
        constants <- list(N = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                          weights = nbInfo$weights, adj = nbInfo$adj,
                          scale = scale, fix.sd = sd.value[sd], fix.lambda = lambda.value[l])
        
        inits <- list(list(alpha = 0, theta = rnorm(S.area, sd = 0.1), 
                           phi = rnorm(S.area, sd = 0.1)),
                      list(alpha = 0, theta = rnorm(S.area, sd = 0.1), 
                           phi = rnorm(S.area, sd = 0.1)),
                      list(alpha = 0, theta = rnorm(S.area, sd = 0.1), 
                           phi = rnorm(S.area, sd = 0.1)))
        
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
                               monitors = c('alpha', 'sigma.b', 'tau.b',
                                            'b', 'r', 'MSS.r',
                                            'RMSS.r', 'TCV'),
                               samplesAsCodaMCMC = TRUE,
                               setSeed = 20112023,
                               WAIC = TRUE)
        
        BYM2.res[act.sim-c.save] <- list(mcmc.out)
        
        if (act.sim%%25 == 0){
          save(list = c("BYM2.res"),
               file = paste0("./BYM2/Results_SimulationStudy_",n.areas[na],"_",hotspot,"_",sd.value[sd],"_",lambda.value[l],"_",act.sim,".Rdata"))
          
          c.save<-c.save+25
          BYM2.res <- list()
        }
        
        act.sim <- act.sim +1
        
        if(act.sim==n.sim+1){break}
      }
      
    }
  }
  
  
}






  

