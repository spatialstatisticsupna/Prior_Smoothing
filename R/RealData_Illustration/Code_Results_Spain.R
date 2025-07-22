################################################################################
################################################################################
##########                    5. Data Illustration                    ##########
##########                       Results Spain                        ##########
################################################################################
################################################################################

#################################################
###    Required libraries   ###
#################################################
rm(list=ls())
library(spdep)
library(MASS)
library(sf)
library("ggplot2")


#################################################
###    Required Constants   ###
#################################################
### name for each case
unif <- c("small", "medium", "large", "non")
### Priors
priors <- c("iid", "GP", "iCAR", "BYM", "pCAR", "LCAR", "BYM2")
### Number of areas
S.area <- c(47, 100, 300)


################################################################################
########## Table
################################################################################

for (sa in 1:length(S.area)) {
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../../Data/Carto_Spain_",S.area[sa],"areas.Rdata"))
  carto <- Carto.areas
  n.area <- length(unique(carto$ID))
  
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  #################################################
  ###    Load the data   ###
  #################################################
  load(paste0("../../Data/Data_Spain_",S.area[sa],"areas.Rdata"))
  Data.areas$crude.rate <- Data.areas$O/Data.areas$Pop*10^5
  
  
  carto$rate <- Data.areas$crude.rate
  carto$pop <- Data.areas$Pop
  carto$count <- Data.areas$O
  
  for (sg in 1:length(unif)) {
    tab.results <- matrix(NA, ncol = 8, nrow = length(priors))
    load(paste0("./Results/Results_Spain_",S.area[sa],"_",unif[sg],".Rdata"))
    
    for (pr in 1:length(priors)) {
      eval(parse(text = paste0("res <- results$",priors[pr],".res")))
      
      tab.results[pr, 1] <- round(mean(res$summary$all.chains[paste0("MSS.r[", 1:S.area[sa], "]"), "Mean"]*10^10),3)
      tab.results[pr, 2] <- round(mean(res$summary$all.chains[paste0("RMSS.r[", 1:S.area[sa], "]"), "Mean"]*10^5),3)
      
      tab.results[pr, 3] <- round(max(res$summary$all.chains[paste0("MSS.r[", 1:S.area[sa], "]"), "Mean"]*10^10),3)
      tab.results[pr, 4] <- round(max(res$summary$all.chains[paste0("RMSS.r[", 1:S.area[sa], "]"), "Mean"]*10^5),3)
      
      tab.results[pr, 8] <- round(sum(res$summary$all.chains[paste0("MSS.r[", 1:S.area[sa], "]"), "Mean"]*10^10)/sum((mean(Data.areas$crude.rate) - Data.areas$crude.rate)^2),3)
      
      
      if(priors[pr]%in%c("iid")){
        tab.results[pr, 5] <- round(sum(length(nbInfo$num)*res$summary$all.chains["sigma", "Mean"]),3)
        tab.results[pr, 6] <- round(res$summary$all.chains["sigma", "Mean"],3)
      }
      if(priors[pr]%in%c("iCAR", "pCAR")){
        tab.results[pr, 5] <- round(sum(res$summary$all.chains["sigma", "Mean"]/nbInfo$num),3)
        tab.results[pr, 6] <- round(res$summary$all.chains["sigma", "Mean"],3)
      }
      else if(priors[pr]=="LCAR"){
        tab.results[pr, 5] <- round(sum(res$summary$all.chains["sd.theta", "Mean"]/(res$summary$all.chains["rho", "Mean"]*nbInfo$num + 1- res$summary$all.chains["rho", "Mean"])),3)
        tab.results[pr, 6] <- round(res$summary$all.chains["sd.theta", "Mean"],3)
        tab.results[pr, 7] <- round(res$summary$all.chains["rho", "Mean"],3)
      }
      else if(priors[pr]=="BYM"){
        tab.results[pr, 6] <- round(res$summary$all.chains["sigma.phi", "Mean"],3)
        tab.results[pr, 7] <- round(res$summary$all.chains["sigma.iid", "Mean"],3)
        
        carto.nb <- poly2nb(carto)
        W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
        W.scale <- -W
        diag(W.scale) <- abs(apply(W.scale, 1, sum))
        inv.R <- MASS::ginv(as(W.scale, "matrix"))
        inv.condvar <- ginv(inv.R + (tab.results[pr, 7]/tab.results[pr, 6])*diag(1, length(carto.nb)))
        tab.results[pr, 5] <- round(sum(tab.results[pr, 6]/diag(inv.condvar)),3)
        
      }
      else if(priors[pr]=="BYM2"){
        tab.results[pr, 6] <- round(res$summary$all.chains["sigma.b", "Mean"],3)
        tab.results[pr, 7] <- round(sum(res$summary$all.chains["rho", "Mean"]),3)
        
        W <- nb2mat(carto.nb, zero.policy = TRUE, style = "B")
        nbInfo <- nb2WB(carto.nb)
        
        W.scale <- -W
        diag(W.scale) <- abs(apply(W.scale, 1, sum))
        Q <- W.scale*exp(mean(log(diag(MASS::ginv(as(W.scale,"matrix"))))))
        inv.R <- MASS::ginv(as(Q, "matrix"))
        
        inv.condvar <- ginv(tab.results[pr, 7]*inv.R + (1-tab.results[pr, 7])*diag(1, length(carto.nb)))
        tab.results[pr, 5] <- round(sum(tab.results[pr, 6]/diag(inv.condvar)),3)
        
      }
      else if(priors[pr]=="GP"){
        tab.results[pr, 6] <- round(res$summary$all.chains["sigma", "Mean"],3)
        tab.results[pr, 7] <- round(res$summary$all.chains["rho", "Mean"],3)
        
        g <- st_as_sf(carto)
        g_centroid <- st_point_on_surface(x = g)
        distMatrix <- st_distance(g_centroid, g_centroid)
        distMatrix <- matrix(distMatrix/100000, ncol = nrow(g_centroid))
        R <- matrix(nrow = n.area, ncol = n.area)
        for(i in 1:n.area){
          for(j in 1:n.area){
            R[i, j] <- exp(-distMatrix[i,j]/tab.results[pr, 7]) 
          }
        }
        inv.R <- MASS::ginv(as(R, "matrix"))
        tab.results[pr, 5] <- round(sum(tab.results[pr, 6]/diag(inv.R)),3)
      }
      
    }
    
    tab.results <- cbind(priors, tab.results)
    tab.results <- rbind(c(" ", "MSS", "RMSS", "max MSS", "max RMSS", "TCV", "sigma", "param", "prop"),
                         tab.results)
    latex_table <- xtable::xtable(tab.results, caption=paste("Results Spain", S.area[sa],"areas with", unif[sg], "informative prior"),
                                  digits=3)
    xtable::print.xtable(latex_table, include.rownames = FALSE,include.colnames = FALSE, comment=FALSE, caption.placement = getOption("xtable.caption.placement", "top"))
    
  }
}




################################################################################
########## Maps
################################################################################
### name for figures
prior.dist <- c("Small-size informative prior distribution",
                "Medium-size informative prior distribution",
                "Large-size informative prior distribution",
                "Non-informative prior distribution")


for (sa in 1:length(S.area)) {
  #################################################
  ###    Load the cartography   ###
  #################################################
  load(paste0("../../Data/Carto_Spain_",S.area[sa],"areas.Rdata"))
  carto <- Carto.areas
  district <- carto$ID
  
  #################################################
  ###    Load the data   ###
  #################################################
  load(paste0("../../Data/Data_Spain_",S.area[sa],"areas.Rdata"))
  Data.areas$crude.rate <- Data.areas$O/Data.areas$Pop*10^5
  data <- Data.areas
  
  carto$rate <- data$crude.rate
  carto$pop <- data$Pop
  carto$count <- data$O

for (sg in 1:length(unif)) {
  data.results <- cbind(rep("crude rate", S.area[sa]), district, data$crude.rate)
  break_points <- round(quantile(data$crude.rate, probs = seq(0, 1, 1/5)),0)
  break_points <- as.numeric(break_points)[1:5]
  load(paste0("./Results/Results_Spain_",S.area[sa],"_",unif[sg],".Rdata"))

  for (pr in 1:length(priors)) {
    eval(parse(text = paste0("res <- results$",priors[pr],".res")))

    data.results <- rbind(data.results,
                          cbind(rep(priors[pr], S.area[sa]), district, 
                                res$summary$all.chains[paste0("r[", 1:S.area[sa], "]"), "Mean"]*10^5))

  }

  data.results <- data.frame(data.results)
  colnames(data.results) <- c("prior", "ID", "est.rate")
  data.results$prior <- as.factor(data.results$prior)
  levels(data.results$prior)
  data.results$prior<- factor(data.results$prior, levels = c("crude rate", priors))
  data.results$est.rate <- as.numeric(data.results$est.rate)

  data2 <- merge(carto, data.results, by = c("ID"))


  palet <- c("#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404")

  if(sg==2 | sg==4){
    usa_1 <- ggplot(data = data2, aes(group = ID)) +
      geom_sf(aes(fill = est.rate)) +
      facet_wrap(~prior, ncol = 4) +
      labs(title = paste(prior.dist[sg]),
           subtitle = "",
           x = NULL,
           y = NULL) +
      binned_scale(aesthetics = 'fill', scale_name = 'custom',
                   palette = ggplot2:::pal_binned(scales::manual_pal(values = palet)),
                   guide = 'bins',
                   breaks = break_points)+
      theme(
        plot.title = element_text(size = 33, face = 'bold',
                                  hjust=0.5,color='gray35'),
        plot.margin = margin(0, 0.5, 0.5, 0.5, 'cm'),
        panel.background = element_rect(fill = '#f5f5f2', color = NA),
        legend.position = 'bottom',
        legend.title = element_text(size = 25, color='gray35'),
        legend.text = element_text(size = 25, color='gray35'),
        strip.text = element_text(size=25, color='gray35'),
        legend.key = element_rect(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())  +
      guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                    title='rates per 100,000 inhabitants',
                                    title.position='top',
                                    title.hjust=0.5,
                                    ticks.colour='#f5f5f2',
                                    ticks.linewidth=5,
                                    barwidth = 75,
                                    barheight = 2))

  }
  else{
    usa_1 <- ggplot(data = data2, aes(group = ID)) +
      geom_sf(aes(fill = est.rate)) +
      facet_wrap(~prior, ncol = 4) +
      labs(title = paste(prior.dist[sg]),
           subtitle = "",
           x = NULL,
           y = NULL) +
      binned_scale(aesthetics = 'fill', scale_name = 'custom',
                   palette = ggplot2:::pal_binned(scales::manual_pal(values = palet)),
                   guide = 'bins',
                   breaks = break_points)+
      theme(
        plot.title = element_text(size = 33, face = 'bold',
                                  hjust=0.5,color='gray35'),
        plot.margin = margin(0, 0.5, 0.5, 0.5, 'cm'),
        panel.background = element_rect(fill = '#f5f5f2', color = NA),
        legend.position = 'none',
        legend.title = element_text(size = 25, color='gray35'),
        legend.text = element_text(size = 25, color='gray35'),
        strip.text = element_text(size=25, color='gray35'),
        legend.key = element_rect(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())  +
      guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                    title='rates per 100,000 inhabitants',
                                    title.position='top',
                                    title.hjust=0.5,
                                    ticks.colour='#f5f5f2',
                                    ticks.linewidth=5,
                                    barwidth = 75,
                                    barheight = 2))

  }


  eval(parse(text = paste0("p",sg," <- usa_1")))




}



#png file
library(ggpubr)
ga1 <- ggarrange(NULL, p1, NULL, p2,
                 ncol = 1, nrow = 4, heights = c(0.03, 0.43, 0.04, 0.5))

ga2 <- ggarrange( NULL, p3, NULL, p4,
                  ncol = 1, nrow = 4, heights = c(0.03, 0.43, 0.04, 0.5))


ggplot2::ggsave(filename = paste0("./Maps_rates_Spain_",S.area[sa],"_1.pdf"),
                plot = ga1,
                device = "pdf",
                dpi = 600,
                width = 20,
                height = 25,
                units = "in")

ggplot2::ggsave(filename = paste0("./Maps_rates_Spain_",S.area[sa],"_2.pdf"),
                plot = ga2,
                device = "pdf",
                dpi = 600,
                width = 20,
                height = 25,
                units = "in")
}


