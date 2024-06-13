############################# Plot fitted values ############################
### Plotting the fitted models for each inversion and hybrid zone. Three data
### sets are plotted on the same axes, observed inversion frequencies, fitted
### values from a cline+envrionment model, and cline fit estimates.
 
### James Reeve - University of Gothenburg
### 2024-01-18

#### 0: Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Packages
library(tidyverse) 
library(ggpubr)  # Multi-panel plots

### Set file path to data
PATH1 <- "path/to/cline/fit/parameters"
PATH2 <-  "path/to/GLM_fitted_values"

### Create list of inversion names
inversions <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", 
                "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC14.1",
                "LGC14.3", "LGC17.1")



#### 1: Function to access data - warning this is a big function ####
CZ_fit_vals <- function(inversion, site){
  site.prefix <- substr(site, 1, 3)
  
  # Define hybrid zone boundary based on Westram et al. 2021
  if(site.prefix == "ANG") hzb <- 0
  if(site.prefix == "CZA") hzb <- 206.7
  if(site.prefix == "CZB") hzb <- 110.3
  if(site.prefix == "CZD") hzb <- 133.325
  
  ### A: cline fit paramters
  params <- read.table(paste0(PATH1, site.prefix, "_inv_clines_free_201907", 
                              ifelse(grepl("A", site), "16", "17"),
                              ".txt"), header = T)
  # Subset to inversion
  params <- params[params$Inv == inversion,]
  # Subset to site
  if(grepl("_L", site)){params <- params[params$Side == "left",] } else {
    params <- params[params$Side == "right",]}
  
  # Transform the cline width and edge frequency estimates
  params$w <- exp(params$logWidth)
  params$pC <- exp(params$lp_crab)/(1+exp(params$lp_crab))
  params$pW <- exp(params$lp_wave)/(1+exp(params$lp_wave))
  
  
  ### B: Envrionmental varaibles
  # Access data for all islands
  env_var <- read.csv(paste0(PATH1, "all_islands_habitat_20240216.csv"), header = TRUE)
  env_var <- env_var[,c("snail_ID", "LCmeanDist",
                        "height", "SN_rock", "low_barnacle", "low_fucus", "topo")]
  
  # Subset to bay by looking at first 3 characters from the site name
  env_var <- env_var[grepl(site.prefix, substr(env_var$snail_ID, 1, 3)),]
  if(grepl("_L", site)){
    env_var <- env_var[env_var$LCmeanDist < hzb,] # Left of bay 
  } else {
    env_var <- env_var[env_var$LCmeanDist >= hzb,] # Right of bay
  }
  
  # Subtract hybrid zone boundary from LCmeanDist 
  if(site != "ANG"){env_var$LCmeanDist <- abs(env_var$LCmeanDist - hzb)}
  
  
  # Standardize variables
  stnd.var <- function(x){(x-mean(x))/sd(x)}
  vars <- c("SN_rock", "low_barnacle", "low_fucus", "topo", "height")
  env_var[,vars] <- apply(env_var[,vars], 2, stnd.var)
  
  
  ### C: Inversion count per snail
  inv_count <- read.csv(paste0(PATH1, site.prefix, "_15july3.csv"),header = T)
  inv_count <- inv_count[, c("snail_ID", inversion)]
  # Flip the focal inversion to match the cline fits if params$Alleles == 0
  if(params$Allele == 0) inv_count[,inversion] <- 2 - inv_count[,inversion]
  
  
  ### D: Merge inversion counts and envrionmental variables
  # Note: more snails are in the inversion data than the envrionmental data
  #         ANG: 589 env and 367 inv
  #       CZA_L: 366 env and 354 inv
  #       CZA_R: 203 env and 354 inv
  #       CZB_L: 315 env and 381 inv
  #       CZB_R: 282 env and 381 inv
  #       CZD_L: 291 env and 369 inv
  #       CZD_R: 319 env and 369 inv
  dat <- inner_join(inv_count, env_var, by = "snail_ID")
  
  # Shortern column names: this prevents an annoying bug in glmulti due to long
  # filenames being ignored as interaction terms.
  dat <- dat[, c("snail_ID", inversion, "LCmeanDist",
                 "height", "SN_rock", "low_barnacle", "low_fucus", "topo")]
  colnames(dat) <- c(colnames(dat)[1], "Ninv", "LCmeanDist",
                     "height", "rock", "barn", "fucus", "topo")
  
  
  ### E: Cline fit
  # Cline parameters
  cnt <- params$Centre
  cw <- params$w
  fc <- params$pC
  fw <- params$pW
  
  # Cline fit function
  cf <- fc+((fw-fc)/(1+exp(-4*((dat$LCmeanDist-cnt)/cw))))
  
  # Addto data frame
  dat$cf <- cf
  
  
  ### F: Estimate fitted values from aggregate model
  
  # Aggregate model
  agg.fit <- read.csv(paste0(PATH2, "aggregate_fits/glmuti_multimodel_inference_", 
                             site, "_", inversion, "_m2.v4.csv"))
  # Order by estimate
  agg.fit <- agg.fit[order(agg.fit$X),]
  
  # Estimate fitted values
  fv <- agg.fit$Estimate[1]+   # Intercept
    dat$barn*agg.fit$Estimate[2] +   
    dat$barn*dat$fucus*agg.fit$Estimate[3] + dat$barn*dat$height*agg.fit$Estimate[4] +
    dat$barn*dat$rock*agg.fit$Estimate[5] + dat$barn*dat$topo*agg.fit$Estimate[6] +
    dat$cf*agg.fit$Estimate[7] +  
    dat$fucus*agg.fit$Estimate[8] +   
    dat$fucus*dat$height*agg.fit$Estimate[9] + dat$fucus*dat$rock*agg.fit$Estimate[10] + 
    dat$fucus*dat$topo*agg.fit$Estimate[11] + 
    dat$height*agg.fit$Estimate[12] +   
    dat$height*dat$rock*agg.fit$Estimate[13] + dat$height*dat$topo*agg.fit$Estimate[14] +
    dat$rock*agg.fit$Estimate[15] +  
    dat$rock*dat$topo*agg.fit$Estimate[16] +
    dat$topo*agg.fit$Estimate[17]
  
  
  # Backtransform from logit values
  dat$fv <- exp(fv)/(1+exp(fv))
  
  
  ### G: Return data
  res <- dat[,c("snail_ID", "LCmeanDist", "Ninv", "cf", "fv")]
  return(res)
}



#### 2: Window-based scores ####
CZ_wind_scores <- function(inversion, site, window.size){
  dat <- CZ_fit_vals(inversion, site)
  
  # Assign windows
  Winds <- seq(0, ceiling(max(dat$LCmeanDist))+window.size, by = window.size)
  
  dat2 <- array(0, dim = c(length(Winds) - 1, 8))
  for(i in 1:(length(Winds)-1)){
    
    # Subset data to current window
    tmp <- dat[dat$LCmeanDist > Winds[i] & dat$LCmeanDist <= Winds[i+1], ]
    
    # Find average position of snails in window
    Avg.pos <- mean(tmp$LCmeanDist)
    
    # Average inversion frequency
    Avg.inv <- mean(tmp$Ninv / 2)
    # Average fitted value
    Avg.fv <- mean(tmp$fv)
    
    # Number of snails in window
    N <- nrow(tmp)
    # Observed heterozygosity
    Hobs <- sum(tmp$Ninv == 1) / N
    # Expected heterozygosity
    Hexp <- sum(2*tmp$cf*(1 - tmp$cf)) / N
    
    # Enter values
    dat2[i,1] <- Winds[i]
    dat2[i,2] <- Winds[i+1]
    dat2[i,3] <- Avg.pos
    dat2[i,4] <- N
    dat2[i,5] <- Avg.inv
    dat2[i,6] <- Avg.fv
    dat2[i,7] <- Hobs
    dat2[i,8] <- Hexp
  }
  rm(i)
  
  dat2 <- as.data.frame(dat2)
  colnames(dat2) <- c("Swind", "Ewind", "avgPos", "N", "avgObs", "avgFits", "Hobs", "Hexp")
  
  # Remove windows with NAs
  dat2 <- dat2[dat2$avgObs != "NA",]
  
}



#### 3: Plotting function ####
CZ_fit_plot <- function(inversion, site, window.size){

# Using try() here to skip past inversions which are missing cline fit parameters
  
  try({
    ###A: Get cline centre and width from parameter
    params <- read.table(paste0(PATH1, substr(site, 1, 3), "_inv_clines_free_201907", 
                                ifelse(grepl("A", site), "16", "17"), ".txt"), header = T)
    params <- params[params$Inv == inversion & 
                       params$Side == ifelse(grepl("_L", site), "left", "right"),]
    # Cline centre
    cnt <- params$Centre 
    # Cline width
    cw <- exp(params$logWidth)
    
    ### B: Fitted values
    dat <- CZ_fit_vals(inversion, site)
    
    ### C: Window scores
    dat2 <- CZ_wind_scores(inversion, site, window.size)
    
    ### D: Fitted values plot
    p1 <- ggplot()+ 
      geom_rect(aes(xmin = cnt - cw/2, xmax = cnt + cw/2, ymin = 0, ymax = 1), 
                fill = "grey95")+
      geom_vline(xintercept = cnt, lty = 2, col = "grey50")+
      geom_line(data = dat, aes(x = LCmeanDist, y = cf), col = "blue")+
      geom_point(data = dat2, aes(x = avgPos, y = avgObs), col = "black")+
      geom_point(data = dat2, aes(x = avgPos, y = avgFits), pch = 5)+
      labs(y = "p",
           title = paste0(site, ": ", inversion, " - ", window.size, "m windows"))+
      coord_cartesian(xlim = c(0, max(dat$LCmeanDist)), ylim = c(0, 1))+
      theme_bw()+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    
    ### E: Heterozygosity plot
    p2 <- ggplot(dat2, aes(x = avgPos))+
      geom_rect(aes(xmin = cnt - cw/2, xmax = cnt + cw/2, ymin = 0, ymax = 1),
                fill = "grey95")+
      geom_vline(xintercept = cnt, lty = 2, col = "grey50")+
      geom_point(aes(y = Hexp), colour = "navy", pch = 5)+
      geom_point(aes(y = Hobs), colour = "blue")+
      labs(x = "Least-costs distance (m)", y = "H")+
      coord_cartesian(xlim = c(0, max(dat$LCmeanDist)), ylim = c(0, 1))+
      theme_classic()
    
    
    return(ggarrange(p1, p2, nrow = 2, heights = c(1.8, 1)))
    })
}



#### 4: Investigate plots ####

# Plot different window sizes for ANG:LGC6.1/2
ggarrange(
  plotlist = lapply(1:12, CZ_fit_plot, 
                    inversion = "LGC6.1", site = "ANG"),
  nrow = 3, ncol = 4)

# Plot ANG
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "ANG", window.size = 5),
  nrow = 3, ncol = 5)
# Plot CZA_L
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZA_L", window.size = 5),
          nrow = 3, ncol = 5)
# Plot CZA_R
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZA_R", window.size = 5),
          nrow = 3, ncol = 5)
# Plot CZB_L
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZB_L", window.size = 5),
          nrow = 3, ncol = 5)
# Plot CZB_R
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZB_R", window.size = 5),
          nrow = 3, ncol = 5)
# Plot CZD_L
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZD_L", window.size = 5),
          nrow = 3, ncol = 5)
# Plot CZD_R
ggarrange(plotlist = lapply(inversions, CZ_fit_plot, site = "CZD_R", window.size = 5),
          nrow = 3, ncol = 5)



#### 5: Plot and save each fit ####
lapply(c("ANG", "CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R"),
       function(SITE){
                 ggarrange(plotlist = lapply(inversions, function(INV){
                   try({
                     # Environmental model
                     p <- CZ_fit_plot(INV, SITE, 5)
                     
                     # Save plot
                     ggsave(paste0(PATH2, "Fitted_values_m2_", SITE, "_", INV, "_5m_windows.tiff"), 
                            p, device = "tiff", dpi = 300, width = 17.78, height = 16, units = "cm")
                     })
                   
                   }))
               })
