######################################### Variance Decomposition of models #####################################################
### This script parses out the variance in inversion arragement frequencies explained by local scale envrionments, from that 
### explained by the major habitat transitions (i.e., the clines). Given this is a logistic regression, the typical approach of 
### looking at variance in the residuals is not possible, due to the logit transformation. Instead the deviance in the models is 
### used as a proxy for the variance. As the cline+environment models are the aggregate of 100 top sub-models called my glmumti, 
### the median deviance of these sub-model is used.

### James Reeve
### 2024-12-05

### Preparation
rm(list = ls())
dev.off()

### Packages
library(tidyverse)

### Filepaths
PATH1 <- "path/to/cline_data/"
PATH2 <- "path/to/inversion_data/"

### List inversions
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", "LGC7.2", 
        "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC14.1", "LGC14.3", "LGC17.1")

### List of sites
CZs <- c("ANG", "CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R")


#### Custom function to run variance decomposition ####
var.decomp <- function(site, inversion){
  try({
    ### 0. Site details:
    site.prefix <- substr(site, 1, 3)
    site.suffix <- ifelse(grepl("L", site), "L", "R")
    
    # Define hybrid zone boundary (Westram et al., 2021)
    if(site.prefix == "ANG") hzb <- 0
    if(site.prefix == "CZA") hzb <- 206.7
    if(site.prefix == "CZB") hzb <- 110.3
    if(site.prefix == "CZD") hzb <- 133.325
    
    
    
    ### 1. Access data
    ## A: Cline fit parameters (Westram et al., 2021)
    params <- read.table(paste0(PATH2, site.prefix, "_inv_clines_free_201907", 
                                ifelse(grepl("A", site), "16", "17"), ".txt"), header = T)
    # Subset to inversion
    params <- params[params$Inv == inversion,]
    # Subset to site
    if(site.suffix == "L"){ params <- params[params$Side == "left",] } else { params <- params[params$Side == "right",] }
    
    # Transform the cline width and edge frequency estimates
    params$w <- exp(params$logWidth)
    params$pC <- exp(params$lp_crab)/(1+exp(params$lp_crab))
    params$pW <- exp(params$lp_wave)/(1+exp(params$lp_wave))
    
    
    ## B: Envrionmental data
    env.data <- read.csv(paste0(PATH2, "all_islands_habitat_20240216.csv"), header = TRUE)
    env.data <- env.data[,c("snail_ID", "LCmeanDist", "height", "SN_rock", "low_barnacle", "low_fucus", "topo")]
    # Select bay by looking at first 3 characters in the site name
    env.data <- env.data[grepl(site.prefix, substr(env.data$snail_ID, 1, 3)),]
    
    # Select hybrid zone on left or right side of bay
    if(site.suffix == "L"){
      env.data <- env.data[env.data$LCmeanDist < hzb,] # Left of bay 
    } else {
      env.data <- env.data[env.data$LCmeanDist >= hzb,] # Right of bay
    }
    
    # Subtract hybrid zone boundary from LCmeanDist 
    if(site != "ANG"){env.data$LCmeanDist <- abs(env.data$LCmeanDist - hzb)}
    
    # Standardize variables
    stnd.var <- function(x){(x-mean(x))/sd(x)}
    vars <- c("SN_rock", "low_barnacle", "low_fucus", "topo", "height")
    env.data[,vars] <- apply(env.data[,vars], 2, stnd.var)
    
    # Rename variables to match other plots
    colnames(env.data) <- c("snail_ID", "LCmeanDist", "height", "rock", "barn", "fucus", "topo")
    
    
    ## C: Inversion genotypes
    inv_count <- read.csv(paste0(PATH2, site.prefix, "_15july3.csv"),header = T)
    inv_count <- inv_count[, c("snail_ID", inversion)]
    # Flip the focal inversion to match the cline fits if params$Alleles == 0
    if(params$Allele == 0) inv_count[,inversion] <- 2 - inv_count[,inversion]
    
    
    ## D: Get deviance from glmulti
    dev <- read.csv(paste0(PATH1, "AICs/", site, "_", inversion, "_AIC_m2.v4.csv"))
    
    
    
    ### 2. Data wrangling
    ## A: Merge inversion genotypes with envrionmental data
    dat <- inner_join(inv_count, env.data, by = "snail_ID")
    
    
    ## B: Estimate cline fits
    # Cline fit function
    cf <- params$pC+((params$pW-params$pC)/(1+exp(-4*((dat$LCmeanDist-params$Centre)/params$w))))
    # Add cline estimate to envrionmental dataset
    dat$cf <- cf
    
    
    
    ### 3. Fit glm to get null and simple cline deviances
    ## A: Fit models
    # Convert observed inversions to count of both alleles. These counts will be encoded as the number of successes and failures 
    # in the response variable of a binomial model
    y <- cbind(dat[,inversion], 2 - dat[,inversion])
    
    # Model 0 - total variance from inversions
    mod0 <- glm(y ~ 1, family = binomial(link = "logit"))
    # Model 1 - variance explained by simple cline fit
    mod1 <- glm(y ~ dat$cf, family = binomial(link = "logit"))
    
    
    ## B: Calculate deviance from each model
    dev0 <- mod0$deviance
    dev1 <- mod1$deviance
    dev2 <- median(dev$dev)
    
    
    ## C: Write output
    res <- data.frame("site" = site,
                      "inversion" = inversion,
                      "model" = c("null", "cline", "cline_env"),
                      "deviance" = c(dev0, dev1, dev2))
    
    
    return(res)
  })
}


#### Run variance decomposition ####

### 1. Run across inversion and sites
vd <- lapply(CZs, function(s){
  
  # Loop across all inversion
  tmp <- lapply(INVs, var.decomp, site = s)
  
  # Drop inversions missing cline estimates
  # These are flagged as "Errors" in the output
  tmp <- tmp[grep("Error", tmp, invert = TRUE)]
  
  # Rowbind outputs together
  tmp <- do.call(rbind.data.frame, tmp)
  
  return(tmp)
})

vd <- do.call(rbind.data.frame, vd)


### 2. Save results
write.csv(vd, paste0(PATH2, "deviance_across_models_v2.csv"), row.names = FALSE, quote = FALSE)
