##################### Figures S2: Correlations between environmental variables #######################
### This script calcaultes how much of the variance is shared between environmental variables. This is split
### into two anlyses: i) correlations among the five envrionmental variables, and ii) correlations of each 
### envrionmental variable with the cline predictions of each inversion and hybrid zone. These were plotted
### as two panels which were combined into Figure S2 using Inkscape.

### James Reeve
### 2024-11-12

### Preparation
rm(list = ls())
dev.off()

### Packages
library(ggplot2)
library(ggcorrplot)
library(ggpubr)

### Filepaths
PATH <- "path/to/envrionmental/data"

### Parameters
# List of inversion names
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", "LGC7.2", "LGC9.1", 
          "LGC10.1", "LGC10.2", "LGC11.1", "LGC14.1", "LGC14.3", "LGC17.1")
# List of hybrid zones
CZs <- c("CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R", "ANG")



#### 1. Access envrionmental data ####

### A: Access full records
env.data <- read.csv(paste0(PATH, "all_islands_habitat_20240216.csv"), header = TRUE)
env.data <- env.data[,c("snail_ID", "LCmeanDist", "height", "SN_rock", "low_barnacle", "low_fucus", "topo")]

### B: Loop through sites to standardise scores
tmp <- lapply(CZs, function(s){
  
  ## Site details
  site.prefix <- substr(s, 1, 3)
  site.suffix <- ifelse(grepl("L", s), "L", "R")
  
  ## Define hybrid zone boundary (Westram et al., 2021)
  if(site.prefix == "ANG") hzb <- 0
  if(site.prefix == "CZA") hzb <- 206.7
  if(site.prefix == "CZB") hzb <- 110.3
  if(site.prefix == "CZD") hzb <- 133.325
  
  ## Subset to bay
  dat <- env.data[grepl(site.prefix, substr(env.data$snail_ID, 1, 3)),]
  
  ## Subset to snails on left or right side of bay
  if(site.suffix == "L"){
    dat <- dat[dat$LCmeanDist < hzb,] # Left of bay 
  } else {
    dat <- dat[dat$LCmeanDist >= hzb,] # Right of bay
  }
  
  ## Subtract hybrid zone boundary from LCmeanDist 
  if(s != "ANG"){dat$LCmeanDist <- abs(dat$LCmeanDist - hzb)}
  
  ## Standardize variables
  stnd.var <- function(x){(x-mean(x))/sd(x)}
  vars <- c("SN_rock", "low_barnacle", "low_fucus", "topo", "height")
  dat[,vars] <- apply(dat[,vars], 2, stnd.var)
  
  ## Add a site ID column
  dat$site_ID <- s
  
  ## Rename variables to match other plots
  colnames(dat) <- c("snail_ID", "LCmeanDist", "height", "rock", "barn", "fucus", "topo", "site_ID")
  
  return(dat)
})

env.data2 <- do.call(rbind.data.frame, tmp)



#### 2. Cline predictions ####

tmp <- lapply(CZs, function(s){
  
  ## Site details
  site.prefix <- substr(s, 1, 3)
  site.suffix <- ifelse(grepl("L", s), "L", "R")
  
  ## Access least cost distances from envrionmental data (needed for cline fit estimation)
  env.dat <- env.data2[env.data2$site_ID == s, ]
  
  ## Access parameters for each bay (note there is a bit of redunancy here, but I decided to split the 
  ## job to run things more efficently).
  params <- read.table(paste0(PATH, site.prefix, "_inv_clines_free_201907", ifelse(grepl("A", s), "16", "17"), ".txt"), header = T)
  
  ## Subset by hybrid zone
  if(site.suffix == "L"){params <- params[params$Side == "left",] } else { params <- params[params$Side == "right",] }
  
  ## Reverse transform parameter estimates
  params$w <- exp(params$logWidth)
  params$pC <- exp(params$lp_crab)/(1+exp(params$lp_crab))
  params$pW <- exp(params$lp_wave)/(1+exp(params$lp_wave))
  
  ## Predict inversion frequencies from cline model looping over inversions
  cfs <- lapply(INVs, function(i){
    try({
      # Subset parameters
      inv.params <- params[params$Inv == i,]
      
      # Fit cline
      cf <- inv.params$pC+((inv.params$pW-inv.params$pC)/(1+exp(-4*((env.dat$LCmeanDist-inv.params$Centre)/inv.params$w))))
      
      # Pearson correlation between cline estimate and each environmental variable
      r.cline_barn = cor(env.dat$barn, cf, 
                           use = "complete.obs", method = "pearson")
      r.cline_fucus = cor(env.dat$fucus, cf, 
                            use = "complete.obs", method = "pearson")
      r.cline_height = cor(env.dat$height, cf, 
                             use = "complete.obs", method = "pearson")
      r.cline_rock = cor(env.dat$rock, cf, 
                           use = "complete.obs", method = "pearson")
      r.cline_topo = cor(env.dat$topo, cf, 
                           use = "complete.obs", method = "pearson")
      
      # Write as data frame
      res1 <- data.frame("site" = s, "inversion" = i,
                         "env.var" = c("barn", "fucus", "height", "rock", "topo"),
                         "r" = c(r.cline_barn, r.cline_fucus, r.cline_height, 
                                 r.cline_rock, r.cline_topo))
        return(res1)
    })
  })
  
  # Drop inversions missing cline estimates
  # These are flagged as "Errors" in the output
  cfs <- cfs[grep("Error", cfs, invert = TRUE)]
  
  ## Rowbind inversion data frames
  res2 <- do.call(rbind.data.frame, cfs)
  
  return(res2)
})

cline.env.corrs <- do.call(rbind.data.frame, tmp)



#### 3. Figure S2A: correlations among envrionmental variables ####

FigS2A <- ggarrange(plotlist = lapply(CZs, function(s){
  # Correlation matrix
  tmp <- cor(env.data2[env.data2$site_ID == s, c("height", "rock", "topo", "barn", "fucus")], method = "pearson")
  
  ggcorrplot(tmp, method = "square", type = "upper", title = s)+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 4),
          plot.margin = unit(c(0,0.2,0,1), 'lines'))
}), nrow = 4, ncol = 2, common.legend = TRUE)

ggsave("path/to/Plots/FigS2A.tiff", FigS2A, device = "tiff", dpi = 300, width = 15, height = 33.75, units = "cm")



#### 4. Figure S2B: Correlations between cline predictions and environmental variables ####

FigS2B <- ggplot(cline.env.corrs)+
  geom_tile(aes(x = "", y = env.var, fill = r), colour = "black")+
  facet_grid(rows = vars(factor(inversion, levels = INVs)), cols = vars(factor(site, levels = CZs)))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, breaks = seq(-1, 1, 0.5))+
  labs(x = "", y = "Envrionmental variables")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(family = "arial", face = "bold", size = 8, angle = 70),
        strip.text.y = element_text(family = "arial", face = "bold", size = 8, angle = 360))

ggsave("path/to/Plots/FigS2B.tiff", FigS2B, device = "tiff", dpi = 300, width = 12, height = 40, units = "cm")
