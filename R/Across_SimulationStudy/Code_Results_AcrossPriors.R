################################################################################
################################################################################
##########             4.2 Across priors simulation study             ##########
##########                           Results                          ##########
################################################################################
################################################################################

#################################################
###    Required Constants   ###
#################################################
### Fitted priors
priors <- c("iid", "GP", "iCAR", "BYM", "pCAR", "LCAR", "BYM2")
### Number of areas
areas <- c(47, 100, 300)
### Scenarios
scenarios <- c("Scenario1", "Scenario2", "Scenario3")

### Number of simulations
n.sim <- 1000
### We save results each 25 simulations
save.value <- seq(25, n.sim, by =25)



#################################################
###    Code for all nº of areas and priors   ###
#################################################

for (sa in 1:length(areas)) {
  S.area <- areas[sa]
  
  #################################################
  ###    Load the cartography    ###
  ################################################
  load(paste0("../../Data/Carto_Spain_",S.area,"areas.RData"))
  carto <- Carto.areas
  

  for (sc in 1:length(scenarios)) {
    tab.results <- matrix(NA, ncol = 5, nrow = length(priors))
    
    for (pr in 1:length(priors)) {
      meanMSS <- 0
      meanRMSS <- 0
      
      maxMSS <- 0
      maxRMSS <- 0
      
      SP <- 0
      
      for (sv in 1:length(save.value)){
        load(paste0("./SimulationStudy_",priors[pr],"/Results_SimulationStudy_",S.area,"_",scenarios[sc],"_",save.value[sv],".Rdata"))
        load(paste0("../../Data/Data_SimulationStudy_",scenarios[sc],"_",S.area,"areas.Rdata"))
        
        
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
      
      tab.results[pr, 1]<- round(mean(meanMSS/n.sim),3)
      tab.results[pr, 2]<- round(mean(meanRMSS/n.sim),3)
      
      tab.results[pr, 3]<- round(maxMSS/n.sim,3)
      tab.results[pr, 4]<- round(maxRMSS/n.sim,3)
      
      tab.results[pr, 5]<- round(SP/n.sim,3)
      
    }
    
    tab.results <- cbind(priors, tab.results)
    tab.results <- rbind(c(" ", "MSS", "RMSS", "max MSS", "max RMSS", "SP"),
                         tab.results)
    latex_table <- xtable::xtable(tab.results, caption=paste("Results Scenario",scenarios[sc]), digits=3)
    xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
    
    
  }
  
}


  