################# Plot envrionmental varaibles along cine #####################
### Script to create plots to show change of envrionmental varaition along the
### shore for each hybrid zone.

### James Reeve - University of Gothenburg
### 14/08/2024

### Preparation
rm(list = ls())
dev.off()

### Packages
library(ggplot2)
library(ggpubr)
library(ape)


#### A: Access data ####
### Filepaths
PATH <- "path/to/environmental/data"

### Midpoint to separate each bay
CZA_split <- 206.7
CZB_split <- 110.3
CZD_split <- 133.325

### Envrionmental variables
env_var <- read.csv(paste0(PATH, "all_islands_habitat_20240216.csv"), header = TRUE)

### Major habitat transitions
# Using centre parameter from fixed cline models (Westram et al. 2021, supp. mat.)
ht_ANG <- 68.22
ht_CZA_L <- 62.755; ht_CZA_R <- 62.755
ht_CZB_L <- 20.45; ht_CZB_R <- 20.45
ht_CZD_L <- 42.495; ht_CZD_R <- 42.495



#### B: Data wrangelling ####
### Remove unneeded parameters
env_var <- env_var[,c(1, 6, 8:9, 12:14)]

### Split by site
# Note: nrow(env_var) != nrow(all sites) because it contains the CZC enviromental data
ANG <- env_var[grep("ANG", env_var$snail_ID),]
CZA_L <- env_var[grepl("CZA", env_var$snail_ID) & env_var$LCmeanDist < CZA_split,]
CZA_R <- env_var[grepl("CZA", env_var$snail_ID) & env_var$LCmeanDist >= CZA_split,]
CZB_L <- env_var[grepl("CZB", env_var$snail_ID) & env_var$LCmeanDist < CZB_split,]
CZB_R <- env_var[grepl("CZB", env_var$snail_ID) & env_var$LCmeanDist >= CZB_split,]
CZD_L <- env_var[grepl("CZD", env_var$snail_ID) & env_var$LCmeanDist < CZD_split,]
CZD_R <- env_var[grepl("CZD", env_var$snail_ID) & env_var$LCmeanDist >= CZD_split,]

### Subtract distances from the centre of the boulder field at each bay
# For left sites this reverses the order so plots go Crab-Wave habitat
CZA_L$LCmeanDist <- abs(CZA_L$LCmeanDist - CZA_split)
CZB_L$LCmeanDist <- abs(CZB_L$LCmeanDist - CZB_split)
CZD_L$LCmeanDist <- abs(CZD_L$LCmeanDist - CZD_split)
# For right sites this reorients the distances to start from the centre each bay
CZA_R$LCmeanDist <- abs(CZA_R$LCmeanDist - CZA_split)
CZB_R$LCmeanDist <- abs(CZB_R$LCmeanDist - CZB_split)
CZD_R$LCmeanDist <- abs(CZD_R$LCmeanDist - CZD_split)

### Adjust height based on maximum at each site
height_adj <- function(data){data$height - min(data$height)}
ANG$height <- height_adj(ANG)
CZA_L$height <- height_adj(CZA_L); CZA_R$height <- height_adj(CZA_R)
CZB_L$height <- height_adj(CZB_L); CZB_R$height <- height_adj(CZB_R)
CZD_L$height <- height_adj(CZD_L); CZD_R$height <- height_adj(CZD_R)

### Define ecotypes
ANG$Ecotype <- ifelse(ANG$LCmeanDist < ht_ANG, "Crab", "Wave")
CZA_L$Ecotype <- ifelse(CZA_L$LCmeanDist < ht_CZA_L, "Crab", "Wave")
CZA_R$Ecotype <- ifelse(CZA_R$LCmeanDist < ht_CZA_R, "Crab", "Wave")
CZB_L$Ecotype <- ifelse(CZB_L$LCmeanDist < ht_CZB_L, "Crab", "Wave")
CZB_R$Ecotype <- ifelse(CZB_R$LCmeanDist < ht_CZB_R, "Crab", "Wave")
CZD_L$Ecotype <- ifelse(CZD_L$LCmeanDist < ht_CZD_L, "Crab", "Wave")
CZD_R$Ecotype <- ifelse(CZD_R$LCmeanDist < ht_CZD_R, "Crab", "Wave")


#### C: Plot ####
### Function
env_plot <- function(data, variable, habitat.transition){
  #std.var <- (data[,variable] - mean(data[,variable]))/sqrt(var(data[,variable]))
  # Set maximum bounds for variables
  if(variable == "height")LIM <- c(0, 2.642) # Limits descided by manually inspecting data
  if(variable == "topo")LIM <- c(0, 2.6)
  if(variable %in% c("SN_rock", "low_barnacle", "low_fucus"))LIM <- c(-1, 1)

  # Scatterplot
  p1 <- ggplot(data)+
      geom_vline(xintercept = habitat.transition, lty = 2)+
      geom_point(aes(x = LCmeanDist, y = .data[[variable]]))+
      lims(y = LIM)+
      theme_classic()+
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 8))

  # Boxplot
  p2 <- ggplot(data)+
      geom_boxplot(aes(x = Ecotype,  y = .data[[variable]], fill = Ecotype))+
      scale_fill_manual(values = c("white", "grey80"))+
      labs(x = "", y = "")+
      lims(y = LIM)+
      theme_classic()+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = unit(c(0.1,0,0.1,0), "cm"),
            legend.position = "none")
  
  ggarrange(p1, p2, ncol = 2, widths = c(1, 0.2))
}

### Height
pANG.H <- env_plot(ANG, "height", ht_ANG)
pCZA_L.H <- env_plot(CZA_L, "height", ht_CZA_L)
pCZA_R.H <- env_plot(CZA_R, "height", ht_CZA_R)
pCZB_L.H <- env_plot(CZB_L, "height", ht_CZB_L)
pCZB_R.H <- env_plot(CZB_R, "height", ht_CZB_R)
pCZD_L.H <- env_plot(CZD_L, "height", ht_CZD_L)
pCZD_R.H <- env_plot(CZD_R, "height", ht_CZD_R)

### rock (substrate)
pANG.R <- env_plot(ANG, "SN_rock", ht_ANG)
pCZA_L.R <- env_plot(CZA_L, "SN_rock", ht_CZA_L)
pCZA_R.R <- env_plot(CZA_R, "SN_rock", ht_CZA_R)
pCZB_L.R <- env_plot(CZB_L, "SN_rock", ht_CZB_L)
pCZB_R.R <- env_plot(CZB_R, "SN_rock", ht_CZB_R)
pCZD_L.R <- env_plot(CZD_L, "SN_rock", ht_CZD_L)
pCZD_R.R <- env_plot(CZD_R, "SN_rock", ht_CZD_R)

### topo (shore slope)
pANG.T <- env_plot(ANG, "topo", ht_ANG)
pCZA_L.T <- env_plot(CZA_L, "topo", ht_CZA_L)
pCZA_R.T <- env_plot(CZA_R, "topo", ht_CZA_R)
pCZB_L.T <- env_plot(CZB_L, "topo", ht_CZB_L)
pCZB_R.T <- env_plot(CZB_R, "topo", ht_CZB_R)
pCZD_L.T <- env_plot(CZD_L, "topo", ht_CZD_L)
pCZD_R.T <- env_plot(CZD_R, "topo", ht_CZD_R)

### Fucus
pANG.F <- env_plot(ANG, "low_fucus", ht_ANG)
pCZA_L.F <- env_plot(CZA_L, "low_fucus", ht_CZA_L)
pCZA_R.F <- env_plot(CZA_R, "low_fucus", ht_CZA_R)
pCZB_L.F <- env_plot(CZB_L, "low_fucus", ht_CZB_L)
pCZB_R.F <- env_plot(CZB_R, "low_fucus", ht_CZB_R)
pCZD_L.F <- env_plot(CZD_L, "low_fucus", ht_CZD_L)
pCZD_R.F <- env_plot(CZD_R, "low_fucus", ht_CZD_R)

### Barnacles
pANG.B <- env_plot(ANG, "low_barnacle", ht_ANG)
pCZA_L.B <- env_plot(CZA_L, "low_barnacle", ht_CZA_L)
pCZA_R.B <- env_plot(CZA_R, "low_barnacle", ht_CZA_R)
pCZB_L.B <- env_plot(CZB_L, "low_barnacle", ht_CZB_L)
pCZB_R.B <- env_plot(CZB_R, "low_barnacle", ht_CZB_R)
pCZD_L.B <- env_plot(CZD_L, "low_barnacle", ht_CZD_L)
pCZD_R.B <- env_plot(CZD_R, "low_barnacle", ht_CZD_R)

### Plot figure
p <- ggarrange(
  ggarrange(pANG.H, pANG.T, pANG.R, pANG.B, pANG.F, nrow = 1),
  ggarrange(pCZA_L.H, pCZA_L.T, pCZA_L.R, pCZA_L.B, pCZA_L.F, nrow = 1),
  ggarrange(pCZA_R.H, pCZA_R.T, pCZA_R.R, pCZA_R.B, pCZA_R.F, nrow = 1),
  ggarrange(pCZB_L.H, pCZB_L.T, pCZB_L.R, pCZB_L.B, pCZB_L.F, nrow = 1),
  ggarrange(pCZB_R.H, pCZB_R.T, pCZB_R.R, pCZB_R.B, pCZB_R.F, nrow = 1),
  ggarrange(pCZD_L.H, pCZD_L.T, pCZD_L.R, pCZD_L.B, pCZD_L.F, nrow = 1),
  ggarrange(pCZD_R.H, pCZD_R.T, pCZD_R.R, pCZD_R.B, pCZD_R.F, nrow = 1),
  ncol = 1)

annotate_figure(p, bottom = "Least-cost distance (m)")
