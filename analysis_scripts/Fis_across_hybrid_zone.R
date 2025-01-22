######################## Calculate Fis across hybrid zones ###########################
### This script calcuates Fis across each hybrid zones and at the cline edges to 
### idenfity if inversion frequencies are departing from expectations, either due to 
### adaptation to microhabitats withing the cline, or balancing selection.

### James Reeve - University of Gothenburg
### 2024-01-14

#### 0: Preparation ####
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(reshape2)  # Transforming data tables
library(ggpubr)    # Multi-panel plot

### Filepath
PATH <- "path/to/inversion/data"

### List of inversion
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", "LGC7.2", 
          "LGC9.1", "LGC10.1","LGC10.2", "LGC11.1", "LGC14.1", "LGC14.3", "LGC17.1")

### List of site names
CZs <- c("ANG", "CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R")



#### A: Load inversion frequency data ####
Load.Inv.data <- function(inversion, site){
  
  # Locate files
  tmp <- list.files(path = PATH, pattern = substr(site, 1, 3))
  
  Freqs <- read.csv(paste0(PATH, tmp[1]))
  ClineParams <- read.table(paste0(PATH, tmp[2]), header = TRUE)
  LCdist <- read.csv(paste0(PATH, tmp[3]))
  
  # Remove unessecary columns
  Freqs <- Freqs[, c("snail_ID", inversion)]
  LCdist <- LCdist[, c("snail_ID", "LCmeanDist")]
  
  # Add least-cost distances to inversion frequencies
  Freqs <- merge(Freqs, LCdist, by = "snail_ID", all = F)
  
  
  # Split by side of transect
  # Habitat transitions from Anja Westram
  if(site == "ANG") ht <- 0
  if(grepl("CZA", site)) ht <- 206.7
  if(grepl("CZB", site)) ht <- 110.3
  if(grepl("CZD", site)) ht <- 133.325
  
  side <- gsub(".*_", "", site)
  if(side == "L"){Freqs <- Freqs[Freqs$LCmeanDist <= ht, ]
  } else {Freqs <- Freqs[Freqs$LCmeanDist > ht, ]}
  
  # Subtract habitat transition distance from LC distances
  # Note: this is very important! If not done, it changes the results a lot!
  if(site != "ANG"){Freqs$LCmeanDist <- abs(Freqs$LCmeanDist - ht)}
  
  # Reverse the genotype order if the 'Allele' column in
  # the parametes file is '0'
  Allele <- ClineParams[ClineParams$Inv == inversion & 
                          ClineParams$Side == ifelse(side == "L", "left", "right"), "Allele"]
  
  if(Allele == 0) Freqs[, inversion] <- -Freqs[, inversion] + 2
  
  
  return(Freqs)
}



#### B: Calculating expected heterozygosity from a cline model ####
cline_fit <- function(data, inversion, site){
  if(class(inversion) != "character"){stop("Error: the inversion argument must be a string!")}
  if(class(site) != "character"){stop("Error: the site argument must be a string!")}
  if(class(data) != "data.frame"){stop("Error: data must be formatted as a data frame with three columns!")}
  
  
  # Access focal cline's parameters
  tmp <- list.files(path = PATH, pattern = substr(site, 1, 3))
  Params <- read.table(paste0(PATH, tmp[2]), header = TRUE)
  Params <- Params[Params$Inv == inversion, ]
  
  # Subset to just parameters for focal side of shore
  if(gsub(".*_", "", site) == "L"){
    Params <- Params[Params$Side == "left", ]
  } else {
    Params <- Params[Params$Side == "right", ]}
  
  # Set warning if no cline estimate was possible for focal inversion/site
  if(is.na(Params$Centre)) paste("Warning: no cline for", inversion)
  
  # Store parameters
  cnt <- Params$Centre
  w <- Params$logWidth
  p_c <- Params$lp_crab
  p_w <- Params$lp_wave
  
  # Back transform from log-fits
  w <- exp(w)
  p_c <- exp(p_c) / (1 + exp(p_c))
  p_w <- exp(p_w) / (1 + exp(p_w))
  
  # Calculate expected inversion frequency
  p_x <-1/(1+exp(-4*((data$LCmeanDist - cnt) / w)))
  
  # Limit fit by Pcrab and Pwave
  p_x <- p_c + (p_w - p_c)*p_x
  
  # Calculate expected frequencies of each haplotype
  P_i <- p_x^2 # E(AA)
  Q_i <- 2*p_x*(1-p_x) # E(Aa)
  R_i <- (1 - p_x)^2 # E(aa)
  
  # Save as dataset
  # Note: A = inversion, a = non-inverted allele
  res <- data.frame("snail_ID" = data$snail_ID, 
                    "LCmeanDist" = data$LCmeanDist,
                    "inv_freq" = p_x,
                    "inv" = rep(inversion, nrow(data)),
                    "Exp_AA" = P_i, 
                    "Exp_Aa" = Q_i, 
                    "Exp_aa" = R_i)
  
  return(res)
}



calculate.heterozygosity <- function(data, position = range(data[,3])){
  if(class(data) != "data.frame") stop("Error: input data must be a data frame!")
  if(ncol(data) != 3) stop("Error: data must have three columns: sample name, genotypes & distance along transect")
  if(length(position) != 2 | class(position) != "numeric") stop("Error: positions must be entered as min and max in a numeric vector!")
  
  tmp <- data[data[,3] >= position[1] & data[,3] <= position[2],]
  
  # Remove 'NAs'
  tmp <- tmp[complete.cases(tmp),]
  
  N <- nrow(tmp)
  
  AA <- sum(tmp[,2] == 2)
  Aa <- sum(tmp[,2] == 1)
  aa <- sum(tmp[,2] == 0)
  
  Het <- Aa / N
  
  res <- data.frame("Inv" = colnames(tmp)[2],
                    "sPos" = position[1],
                    "ePos" = position[2],
                    "N" = N, 
                    "Obs_AA" = AA, 
                    "Obs_Aa" = Aa, 
                    "Obs_aa" = aa,
                    "Ho" = Het) 
  
  return(res)
}



#### D: Calculate Fis ####
calculate.Fis <- function(data, inversion, site, position = range(data[,3])){
  
  # Calculate expected heterozygosity
  tmp1 <- cline_fit(data = data, inversion = inversion, site = site)
  tmp1 <- tmp1[tmp1$LCmeanDist >= position[1] & tmp1$LCmeanDist <= position[2],]
  He <- sum(tmp1$Exp_Aa)
  
  # Calculate observed heterozygosity
  tmp2 <- calculate.heterozygosity(data = data, position = position)
  Ho <- tmp2$Obs_Aa
  N <- tmp2$N
  
  # Calculate Fis
  Fis <- ifelse(is.na(He) | He < 1, NA, 1 - (Ho/He))
  #try({if(Fis < -1) Fis <- -1}) # Limit Fis to -1
  
  ToH.chisq <- function(N, He, Ho){
    (Ho-He)^2/He + ((N-Ho)-(N-He))^2/(N-He)
  }
  # Create matrix of heterozygosity counts
  if(is.na(He)){
    if(is.na(Ho)){
      
      # Matrix if observed heterozygosity is NA
      #tmp_matrix <- matrix(c(0, 0,  N,  nrow(tmp1)), nrow = 2)
      ChiSq <- ToH.chisq(N, He = 0, Ho = 0)
    } else {
      
      # Matrix if expected heterozygosity is NA
      #tmp_matrix <- matrix(c(Ho, 0, N-Ho,  nrow(tmp1)), nrow = 2)
      ChiSq <- ToH.chisq(N, He = 0, Ho)
    }
  } else {
    
    # Normal matrix
    #tmp_matrix <- matrix(c(Ho, He, N-Ho, nrow(tmp1)-He), nrow = 2)
    ChiSq <- ToH.chisq(N, He, Ho)
  }
  
  
  # Chi-squared test for heterozygosity differences
  #chi <- chisq.test(tmp_matrix)
  p.val <- pchisq(ChiSq, df = 1, lower.tail = FALSE)
  
  ### Output
  res <- data.frame("Inv" = inversion,
                    "Site" = site,
                    "sPos" = position[1],
                    "ePos" = position[2],
                    "N" = N,
                    "nHo" = Ho,
                    "pHo" = Ho / N,
                    "nHe" = He,
                    "pHe" = He / N,
                    "chisq" = ChiSq,#chi$statistic,
                    "p.value" = p.val,#chi$p.value,
                    "Fis" = Fis)
  
  return(res)
}



#### E: Calculate Fis across the whole hybrid zone ####

### Run functions in loop across all inversions and all hybrid zones
Fis_WCZ <- lapply(CZs, function(i){
  
  tmp <- lapply(INVs, function(ii){
    
    calculate.Fis(data =  Load.Inv.data(ii, i),
                  inversion = ii,
                  site = i) })
  
  tmp <- do.call(rbind.data.frame, tmp)
  
  return(tmp[complete.cases(tmp),]) })

Fis_WCZ <- do.call(rbind.data.frame, Fis_WCZ)

# Add site name prefix
Fis_WCZ$Site.prefix <- substr(Fis_WCZ$Site, 1, 3)

# Benjamini-Hochberg correction on P-values
Fis_WCZ <- Fis_WCZ %>% group_by(Site) %>% 
  mutate("adj.p.value" = p.adjust(p.value, method = "BH"))



#### F: Calculating Fis in the Crab habitat ####
Fis_Crab <- lapply(CZs, function(i){
  
  tmp <- lapply(INVs, function(ii){
    
    tryCatch({
      # Load cline fit parameters
      tmp <- list.files(path = PATH, pattern = substr(i, 1, 3))
      Params <- read.table(paste0(PATH, tmp[2]), header = TRUE)
      Params <- Params[Params$Inv == ii,]
      
      # Subset to just parameters for focal side of shore
      if(gsub(".*_", "", i) == "L"){
        Params <- Params[Params$Side == "left", ]
      } else {
        Params <- Params[Params$Side == "right", ]}
      
      # Define Crab habitat as all snails left of cline width
      cnt <- Params$Centre
      w <- Params$logWidth
      w <- exp(w)
      C.edge <- cnt - w/2
      
      # Calculate Fis
      calculate.Fis(data =  Load.Inv.data(ii, i),
                    inversion = ii,
                    site = i,
                    position = c(0, C.edge)) 
    }, error = function(e){cat(paste0("Error: ", i, ": ", ii,
                                   " does not have a cline fit estimate.")) }) 
  })
  tmp <- do.call(rbind.data.frame, tmp)
  
  return(tmp[complete.cases(tmp),]) })

Fis_Crab <- do.call(rbind.data.frame, Fis_Crab)
Fis_Crab$Site.prefix <- substr(Fis_Crab$Site, 1, 3)
Fis_Crab <- Fis_Crab %>% group_by(Site) %>% 
  mutate("adj.p.value" = p.adjust(p.value, method = "BH"))



#### G: Calcaulting Fis in the Wave habitat ####
Fis_Wave <- lapply(CZs, function(i){
  
  tmp <- lapply(INVs, function(ii){
    
    tryCatch({
      # Load cline fit parameters
      tmp <- list.files(path = PATH, pattern = substr(i, 1, 3))
      Params <- read.table(paste0(PATH, tmp[2]), header = TRUE)
      Params <- Params[Params$Inv == ii,]
      
      # Subset to just parameters for focal side of shore
      if(gsub(".*_", "", i) == "L"){
        Params <- Params[Params$Side == "left", ]
      } else {
        Params <- Params[Params$Side == "right", ]}
      
      # Define Crab habitat as all snails left of cline width
      cnt <- Params$Centre
      w <- Params$logWidth
      w <- exp(w)
      W.edge <- cnt + w/2
      
      # Calculate Fis
      calculate.Fis(data =  Load.Inv.data(ii, i),
                    inversion = ii,
                    site = i,
                    position = c(W.edge, 3e2)) 
    }, error = function(e){cat(paste0("Error: ", i, ": ", ii,
                                      " does not have a cline fit estimate.")) }) 
  })
  tmp <- do.call(rbind.data.frame, tmp)
  
  return(tmp[complete.cases(tmp),]) })

Fis_Wave <- do.call(rbind.data.frame, Fis_Wave)
Fis_Wave$Site.prefix <- substr(Fis_Wave$Site, 1, 3)
Fis_Wave <- Fis_Wave %>% group_by(Site) %>% 
  mutate("adj.p.value" = p.adjust(p.value, method = "BH"))



#### H: Save data ####
write.csv(Fis_WCZ, paste0(PATH, "CZ_Fis_across_hybrid_zones.csv"), row.names = FALSE, quote = FALSE)
write.csv(Fis_Crab, paste0(PATH, "CZ_Fis_Crab_habitat.csv"), row.names = FALSE, quote = FALSE)
write.csv(Fis_Wave, paste0(PATH, "CZ_Fis_Wave_habitat.csv"), row.names = FALSE, quote = FALSE)
