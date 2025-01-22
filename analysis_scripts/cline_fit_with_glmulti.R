############### Fitting a GLM of cline with envrionmental data ##################
### This script runs multiple genralized linear models for models of inversion
### frequencies and the envrionmental measurements in Littorina saxatilis.
### Measurements were taken from 7 hybrid zones between the 'crab' and 'wave'
### ecotype. This script has been designed to run of one inversion and site at a
### time. Details about sampling are found in Westram et al. 2021
### (DOI: 10.1111/mec.15861). Two models are tested:
###             m0 = a neutral clinal models, where inversion frequencies are
###                  determined by geneflow [cline model]
###             m1 = the clinal+envrionmental model, which test for diviations
###                  from the cline fit model, which might be caused by the
###                  envrionment
### m1 is run thousands of times using the model selection package 'glmulti' to 
### permute parameters. An aggregate model is then constructed from the weighted
### estimates of the top 100 models. This will take a lot of time on a laptop, so
### I recomend running this script on a server.

### James Reeve - University of Gothenburg
### 2024-01-25


#### 0. Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse) # Computing resources
library(glmulti)   # GLM permutations


#### 1. Parameters ####

### Filepaths
# Access data
PATH <- "path/to/data/dir"
# Model output directory
OUT <- "path/to/where/data/stored"

### Inversion name
INV <- "LGC1.1"
# Possible values: "LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1-2", "LGC7.1", 
#                  "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
#                  "LGC12.2", "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1-2", 
#                  "LGC14.3", or "LGC17.1"

### Hybrid zone name
CZ <- "ANG"
# Possible values: "ANG", "CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", or "CZD_R"
site.prefix <- substr(CZ, 1, 3)

### Position of site split.
# CZA, CZB and CZD are three separate bays which were split in the middle into
# two hybrid zones per site. See Fig. 1 in the manuscript for images of the bays.
if(site.prefix == "ANG") site.split <- 0
if(site.prefix == "CZA") site.split <- 206.7
if(site.prefix == "CZB") site.split <- 110.3
if(site.prefix == "CZD") site.split <- 133.325


#### 2. Data handeling ####

### A: Cline fit parameters
params <- read.table(paste0(PATH, "cline_fit_parameters/", site.prefix, "_inv_clines_free_201907",
                            ifelse(grepl("A", CZ), "16", "17"),
                            ".txt"), header = T)
params <- params[params$Inv == INV,]
if(grepl("L", CZ)){
  params <- params[params$Side == "left",]
  } else {
  params <- params[params$Side == "right",]}

### B: Environmental variables
# Access data for all islands
env_var <- read.csv(paste0(PATH, "environment/all_islands_habitat_20240216.csv"))
# Remove unnecessary columns
env_var <- env_var[,c("snail_ID", "height", "LCmeanDist", "SN_rock", "topo", "low_barnacle", "low_fucus")]#-c(1:2, 4:5, 7, 10:11, 15:16)]

# Subset to bay
env_var <- env_var[grepl(site.prefix, substr(env_var$snail_ID, 1, 3)),]

# Subset to site
if(grepl("_L", CZ)){
  env_var <- env_var[env_var$LCmeanDist < site.split,] # Left of bay 
} else {
  env_var <- env_var[env_var$LCmeanDist >= site.split,] # Right of bay
}

# Subtract distances from the centre of the boulder field at each bay
env_var$LCmeanDist <- abs(env_var$LCmeanDist - site.split)

# Standardize environmental variables
stnd.var <- function(x){(x-mean(x))/sd(x)}
vars <- c("SN_rock", "low_barnacle", "low_fucus", "topo", "height")
env_var[,vars] <- apply(env_var[,vars], 2, stnd.var)


### C: Inversion count per snail
inv_count <- read.csv(paste0(PATH, site.prefix, "_15july3.csv"),header = T)
# Remove unnessecary columns
inv_count <- inv_count[, c("snail_ID", INV)]

# Switch to alternative allele if it was more common in Wave habitat
if(params$Allele == 0) inv_count[,INV] <- 2 - inv_count[,INV]


### D: Merge inversion counts and envrionmental variables
dat <- inner_join(inv_count, env_var, by = "snail_ID")

# Rename columns: this prevents an annoying bug in glmulti that looses barnacle*fucus interactions
colnames(dat) <- c(colnames(dat)[1], "Ninv", colnames(dat)[3:4], "rock", "topo", "barn", "fucus")


#### 3. Cline model estimates ####

### Cline positions from environmental data
lcd <- dat$LCmeanDist

### Extract parameters required for the model
cnt <- params$Centre # Cline centre
cw <- params$logWidth # Cline width
fc <- params$lp_crab # Allele frequency Crab edge
fw <- params$lp_wave # Allele frequency Wave edge

### Back-transform parameters in logit and log scale
# Log
cw <- exp(cw)
# Logit
fc <- exp(fc)/(1+exp(fc))
fw <- exp(fw)/(1+exp(fw))

### Cline fit function
cf <- fc+(fw-fc)*(1/(1+exp(-4*((lcd-cnt)/cw))))

### Add cline fit estimates to data
dat$cline <- cf


#### 4. m0: Fitting cline model ####

### Run GLM on cline fit
# GLM function - some inversions are missing cline fit parameters. The if()
# statement sets a flat model for missing parameters.
if(sum(is.na(dat$cline)) == nrow(dat)){
  f <- formula(cbind(Ninv, 2-Ninv)~1)
} else {
  f <- formula(cbind(Ninv, 2-Ninv)~cline)
}

# Run the GLM
fit <- glm(f, family = binomial(link = "logit"), data = dat)

### Return a dataframe of AIC scores
m0.AIC <- data.frame("Inv" = INV,
                     "Site" = CZ,
                     "aic" = fit$aic,
                     "dev" = fit$deviance,
                     "model" = if(sum(is.na(dat$cline)) != nrow(dat)){
                       paste("Pinv", fit$formula[3], sep = "~")} else {
                         "No-cline"
                       })

write.csv(m0.AIC, paste0(OUT, "AICs/", CZ, "_", INV, "_AIC_m0.v4.csv"), row.names = FALSE, quote = c(4))


#### 5. m1: Environmental + Cline fit model ####

### A: Set up custom function for glmulti
# Note: due to a problem with the way functions are nested in R, I have to have
# call this function outside the wrapper function for the env+cline model.
myglm <- function(f, data) {
  if(missing(data)) data <- envrionment(f)

  # Null model
  nullos <- glm(cbind(Ninv, 2-Ninv)~1,
                family = binomial(link = "logit"), data = data)

  # Get a list of all possible terms in the model
  termz = terms(f,data=data)
  # Order of all terms
  orderz = attr(termz,"order")
  # All 1st level factors
  sfact <- which(orderz == 1)
  # All 2nd level pairwise interactions
  intz = which(orderz == 2)

  # Locate the rows which corresponding to 'cline'
  index=which(dimnames(attr(termz,"factors"))[[1]] == "cline")


  # Check that any model terms are present
  if(length(orderz > 1)){

    # If 'cline' is not in the model return Null model
    if(sum(attr(termz,"factors")[index,sfact]) != 1) return(nullos)
    # If 'cline' is in any interaction, return Null model
    if(sum(attr(termz,"factors")[index,intz]) != 0) return(nullos)

  } else return(nullos)

  # If 'cline' is present and not in any interactions, return a glm
  return(glm(formula=f, family = binomial(link = "logit"), data = data))
}


### B: Run GLM
# Function for full model including cline fit
f <- formula(cbind(Ninv, 2-Ninv)~height*rock*barn*fucus*topo*cline)
# Run the glm with the modified fit function
fits <- do.call("glmulti", list(f, data = dat, fitfunc = myglm))


### C: Save plots of fit behavior
tiff(file = paste0(OUT, "fit_plots/", CZ, "_", INV, "_m1fit.tiff"))
  par(mfrow = c(3,  1))
  plot(fits, type = "s")
  plot(fits, type = "p")
  plot(fits, type = "w")
dev.off()


### D: Save AIC of all models
m1.AIC <- data.frame("Inv" = rep(INV, length(fits@formulas)),
                     "Site" = rep(CZ, length(fits@formulas)),
                     "rank" = 1:length(fits@formulas),
                     "aic" = sapply(fits@objects, function(X){X$aic}),
                     "dev" = sapply(fits@objects, function(X){X$deviance}),
                     "model" = paste("Pinv", sapply(fits@formulas, function(X){X[[3]]}), sep = " ~ "))

write.csv(m1.AIC, paste0(OUT, "AICs/", CZ, "_", INV, "_AIC_m1.csv"), row.names = FALSE, quote = c(4))


### E: Aggregate model
# Estimates from all "top" 2000 models
m1.agg.fit <- coef.glmulti(fits)
m1.agg.fit <- as.data.frame(m1.agg.fit)
# Transform variance to get standard deviation
m1.agg.fit <- data.frame("Estimate" = m1.agg.fit$Est,
                      "sd" = sqrt(m1.agg.fit$Uncond),
                      "Importance" = m1.agg.fit$Importance,
                      row.names = row.names(m1.agg.fit))
# Empirical P-value
m1.agg.fit$z <- m1.agg.fit$Estimate / m1.agg.fit$sd
m1.agg.fit$p <- 2*pnorm(abs(m1.agg.fit$z), lower.tail=FALSE)
# Rename columns
names(m1.agg.fit) <- c("Estimate", "sd", "Importance", "z value", "Pr(>|z|)")
# Order factors by importance
m1.agg.fit <- m1.agg.fit[order(m1.agg.fit$Importance, decreasing=TRUE), c(1,2,4:5,3)]

### Save results of aggregate model
write.csv(m1.agg.fit, paste0(OUT, "aggregate_fits/glmuti_multimodel_inference_", CZ, "_", INV, "_m1.csv"), quote = FALSE)
