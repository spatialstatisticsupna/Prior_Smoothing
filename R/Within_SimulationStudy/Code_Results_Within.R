################################################################################
################################################################################
##########              4.1 Within prior simulation study             ##########
##########                           Results                          ##########
################################################################################
################################################################################

#################################################
###    Required libraries   ###
#################################################
library(spdep)
library(MASS)
# library(sf)


#################################################
###    Required Constants   ###
#################################################

### Fitted priors
priors <- c("GP", "iCAR", "BYM", "LCAR")

### Fix values for hyperparameters
sigma <- c(0.01, 0.05, 0.09, 0.2, 0.5) #Common for all priors
v <- c(0.25, 1, 4) # For BYM
lambda <- c(0.1, 0.5, 0.9) #For LCAR
rho <- c(1, 5, 9) #For GP

### Number of areas
areas <- c(47, 100, 300)
### Used scenario
scenario <- "Scenario2" 

### Number of simulations
n.sim <- 25
### We save results each 25 simulations
save.value <- seq(25, n.sim, by =25)



#################################################
###    Code for all nº of areas and priors   ###
#################################################

for (sa in 1:length(areas)) {
  
  #################################################
  ###    Load the cartography to compute TCV   ###
  ################################################
  load(paste0("../../Data/Carto_Spain_",areas[sa],"areas.RData"))
  carto <- Carto.areas
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  
  S.area <- length(unique(carto$ID))
  
  for (pr in 1:length(priors)) {
    if (priors[pr]%in%c("iid", "iCAR", "pCAR")){
      tab.results <- matrix(NA, ncol = 7, nrow = length(sigma))
        
      for (sc in 1:length(sigma)) {
        
        ################################################
        ###                Compute TCV               ###
        ################################################
        neigh <- nbInfo$num
        
        if(priors[pr]=="iid"){TCV <- length(neigh)*sigma[sc]^2}
        else{TCV <- round(sum(sigma[sc]^2/neigh),3)}
        
        
        #################################################
        ###         Compute Empirical metrics        ###
        ################################################
        meanMSS <- 0
        meanRMSS <- 0
        
        maxMSS <- 0
        maxRMSS <- 0
        
        SP <- 0
          
        for (sv in 1:length(save.value)){
          load(paste0("./",priors[pr],"/Results_SimulationStudy_",areas[sa],"_",scenario,"_",sigma[sc],"_",save.value[sv],".Rdata"))
          load(paste0("../../Data/Data_SimulationStudy_",scenario,"_",areas[sa],"areas.Rdata"))
            
          for (ns in 1:25) {
            eval(parse(text = paste0("res <- ",priors[pr],".res[[ns]]")))
            
            act.sim <- ns + (save.value[sv]-25)
            data <- DataSIM[which(DataSIM$sim==act.sim),]
            
            meanMSS <- meanMSS + res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10
            meanRMSS <- meanRMSS + res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5
            
            maxMSS <- maxMSS + max(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)
            maxRMSS <- maxRMSS + max(res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5)
              
            SP <- SP + sum(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)/sum((mean(data$crude.rate) - data$crude.rate)^2)
          }
        }
        tab.results[sc, 1]<- paste("(",sigma[sc],")")
          
        tab.results[sc, 2]<- TCV 
        
        tab.results[sc, 3]<- round(SP/n.sim,3)
          
        tab.results[sc, 4]<- round(mean(meanMSS/n.sim),3)
        tab.results[sc, 5]<- round(mean(meanRMSS/n.sim),3)
          
        tab.results[sc, 6]<- round(maxMSS/n.sim,3)
        tab.results[sc, 7]<- round(maxRMSS/n.sim,3)

      }
      
  
    tab.results <- rbind(c("parameters", "TCV", "SP", "MSS", "RMSS", "max MSS", "max RMSS"),
                           tab.results)
    latex_table <- xtable::xtable(tab.results, caption=paste("Results ",priors[pr]), digits=3)
    xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
      
  }
    
    if (priors[pr]%in%c("LCAR", "BYM2")){
      tab.results <- matrix(NA, ncol = 7, nrow = length(sigma)*length(lambda))
      
      for (sc in 1:length(sigma)) {
        for (lb in 1:length(lambda)) {
          ################################################
          ###                Compute TCV               ###
          ################################################
          neigh <- nbInfo$num
          
          if(priors[pr]=="LCAR"){TCV <- round(sum(sigma[sc]^2/(lambda[lb]*neigh + 1- lambda[lb])),3)}
          else{
            W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
            nbInfo <- nb2WB(carto.nb)
            
            W.scale <- -W
            diag(W.scale) <- abs(apply(W.scale, 1, sum))
            Q <- W.scale*exp(mean(log(diag(MASS::ginv(as(W.scale,"matrix"))))))
            inv.R <- MASS::ginv(as(Q, "matrix"))
            inv.Q <- MASS::ginv(lambda[lb]*inv.R + (1-lambda[lb])*diag(1, length(nbInfo$num)))
            
            TCV <- round(sum(sigma[sc]^2/diag(inv.Q)),3)
          } 
          
          
          #################################################
          ###         Compute Empirical metrics        ###
          ################################################
          meanMSS <- 0
          meanRMSS <- 0
            
          maxMSS <- 0
          maxRMSS <- 0
          
          SP <- 0
            
          for (sv in 1:length(save.value)){
            load(paste0("./",priors[pr],"/Results_SimulationStudy_",areas[sa],"_",scenario,"_",sigma[sc],"_",lambda[lb],"_",save.value[sv],".Rdata"))
            load(paste0("../../Data/Data_SimulationStudy_",scenario,"_",areas[sa],"areas.Rdata"))
              
            for (ns in 1:25) {
              eval(parse(text = paste0("res <- ",priors[pr],".res[[ns]]")))
              
              act.sim <- ns + (save.value[sv]-25)
              data <- DataSIM[which(DataSIM$sim==act.sim),]
              
              meanMSS <- meanMSS + res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10
              meanRMSS <- meanRMSS + res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5
              
              maxMSS <- maxMSS + max(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)
              maxRMSS <- maxRMSS + max(res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5)
              
              SP <- SP + sum(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)/sum((mean(data$crude.rate) - data$crude.rate)^2)
                
            }
          }
          tab.results[lb + length(lambda)*(sc-1), 1]<- paste("(",sigma[sc]," - ",lambda[lb],")")
            
          tab.results[lb + length(lambda)*(sc-1), 2]<- TCV
            
          tab.results[lb + length(lambda)*(sc-1), 3]<- round(SP/n.sim,3)
            
          tab.results[lb + length(lambda)*(sc-1), 4]<- round(mean(meanMSS/n.sim),3)
          tab.results[lb + length(lambda)*(sc-1), 5]<- round(mean(meanRMSS/n.sim),3)
            
          tab.results[lb + length(lambda)*(sc-1), 6]<- round(maxMSS/n.sim,3)
          tab.results[lb + length(lambda)*(sc-1), 7]<- round(maxRMSS/n.sim,3)
            
      }
      }
      
      tab.results <- rbind(c("parameters", "TCV", "SP", "MSS", "RMSS", "max MSS", "max RMSS"),
                           tab.results)
      latex_table <- xtable::xtable(tab.results, caption=paste("Results ",priors[pr]), digits=3)
      xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
      
    }
    
    if (priors[pr]%in%c("BYM")){
      tab.results <- matrix(NA, ncol = 7, nrow = length(sigma)*length(v))
      
      for (sc in 1:length(sigma)) {
        for (lb in 1:length(v)) {
          sigma2 <- sigma[sc]*sqrt(v[lb])
          
          ################################################
          ###                Compute TCV               ###
          ################################################
          W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
          W.scale <- -W
          diag(W.scale) <- abs(apply(W.scale, 1, sum))
          inv.R <- MASS::ginv(as(W.scale, "matrix"))
          
          inv.Q2 <- MASS::ginv(inv.R + (sigma2/sigma[sc])^2*diag(1, length(nbInfo$num)))
          TCV <- round(sum(sigma[sc]^2/diag(inv.Q2)),3)
          
          #################################################
          ###         Compute Empirical metrics        ###
          ################################################
          meanMSS <- 0
          meanRMSS <- 0
            
          maxMSS <- 0
          maxRMSS <- 0
           
          SP <- 0
            
          for (sv in 1:length(save.value)){
            load(paste0("./",priors[pr],"/Results_SimulationStudy_",areas[sa],"_",scenario,"_",sigma[sc],"_",v[lb],"_",save.value[sv],".Rdata"))
            load(paste0("../../Data/Data_SimulationStudy_",scenario,"_",areas[sa],"areas.Rdata"))
              
            for (ns in 1:25) {
              eval(parse(text = paste0("res <- ",priors[pr],".res[[ns]]")))
              
              act.sim <- ns + (save.value[sv]-25)
              data <- DataSIM[which(DataSIM$sim==act.sim),]
                
              meanMSS <- meanMSS + res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10
              meanRMSS <- meanRMSS + res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5
              
              maxMSS <- maxMSS + max(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)
              maxRMSS <- maxRMSS + max(res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5)
              
              SP <- SP + sum(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)/sum((mean(data$crude.rate) - data$crude.rate)^2)
                
            }
          }
          tab.results[lb + length(v)*(sc-1), 1]<- paste("(",sigma[sc]," - ",v[lb],")")
            
          tab.results[lb + length(v)*(sc-1), 2]<- TCV
          
          tab.results[lb + length(v)*(sc-1), 3]<- round(SP/n.sim,3)
          
          tab.results[lb + length(v)*(sc-1), 4]<- round(mean(meanMSS/n.sim),3)
          tab.results[lb + length(v)*(sc-1), 5]<- round(mean(meanRMSS/n.sim),3)
          
          tab.results[lb + length(v)*(sc-1), 6]<- round(maxMSS/n.sim,3)
          tab.results[lb + length(v)*(sc-1), 7]<- round(maxRMSS/n.sim,3)
          
          tab.results[lb + length(v)*(sc-1), 3]<- round(SP/n.sim,3)
        
      }
    }
      
      tab.results <- rbind(c("parameters", "TCV", "SP", "MSS", "RMSS", "max MSS", "max RMSS"),
                           tab.results)
      latex_table <- xtable::xtable(tab.results, caption=paste("Results ",priors[pr]), digits=3)
      xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
      
    }
    
    if (priors[pr]%in%c("GP")){
      tab.results <- matrix(NA, ncol = 7, nrow = length(sigma)*length(rho))
      
      for (sc in 1:length(sigma)) {
        for (lb in 1:length(rho)) {
          ################################################
          ###                Compute TCV               ###
          ################################################
          g <- st_as_sf(carto)
          g_centroid <- st_point_on_surface(x = g)
          distMatrix <- st_distance(g_centroid, g_centroid)
          distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid))
          
          R <- matrix(nrow = S.area, ncol = S.area)
          for(i in 1:S.area){
            for(j in 1:S.area){
              R[i, j] <- exp(-distMatrix[i,j]/rho[lb]) 
            }
          }
          inv.R <- MASS::ginv(as(R, "matrix"))
          
          TCV <- round(sum(sigma[sc]^2/diag(inv.R)),3)
          
          #################################################
          ###         Compute Empirical metrics        ###
          ################################################
          meanMSS <- 0
          meanRMSS <- 0
            
          maxMSS <- 0
          maxRMSS <- 0
            
          SP <- 0
            
          for (sv in 1:length(save.value)){
            load(paste0("./",priors[pr],"/Results_SimulationStudy_",areas[sa],"_",scenario,"_",sigma[sc],"_",rho[lb],"_",save.value[sv],".Rdata"))
            load(paste0("../../Data/Data_SimulationStudy_",scenario,"_",areas[sa],"areas.Rdata"))
              
            for (ns in 1:25) {
              eval(parse(text = paste0("res <- ",priors[pr],".res[[ns]]")))
              
              meanMSS <- meanMSS + res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10
              meanRMSS <- meanRMSS + res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5
                
              maxMSS <- maxMSS + max(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)
              maxRMSS <- maxRMSS + max(res$summary$all.chains[paste0("RMSS.r[", 1:S.area, "]"), "Mean"]*10^5)
                
              act.sim <- ns + (save.value[sv]-25)
              data <- DataSIM[which(DataSIM$sim==act.sim),]
              SP <- SP + sum(res$summary$all.chains[paste0("MSS.r[", 1:S.area, "]"), "Mean"]*10^10)/sum((mean(data$crude.rate) - data$crude.rate)^2)
                
            }
          }
          tab.results[lb + length(rho)*(sc-1), 1]<- paste("(",sigma[sc]," - ",rho[lb],")")
            
          tab.results[lb + length(rho)*(sc-1), 2]<- TCV
            
          tab.results[lb + length(rho)*(sc-1), 3]<- round(SP/n.sim,3)
            
          tab.results[lb + length(rho)*(sc-1), 4]<- round(mean(meanMSS/n.sim),3)
          tab.results[lb + length(rho)*(sc-1), 5]<- round(mean(meanRMSS/n.sim),3)
            
          tab.results[lb + length(rho)*(sc-1), 6]<- round(maxMSS/n.sim,3)
          tab.results[lb + length(rho)*(sc-1), 7]<- round(maxRMSS/n.sim,3)
        }
      }
      
      tab.results <- rbind(c("parameters", "TCV", "SP", "MSS", "RMSS", "max MSS", "max RMSS"),
                           tab.results)
      latex_table <- xtable::xtable(tab.results, caption=paste("Results ",priors[pr]), digits=3)
      xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
      
    }
  }
}






