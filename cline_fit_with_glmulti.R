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
### 2023-06-02


#### 0. Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse) # Computing resources
library(glmulti)   # GLM permutations
library(ggpubr)    # Multi-panel plots


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

### Position where a site it split.
# CZA, CZB and CZD are three separate bays which were split in the middle into
# two hybrid zones per site. See Fig. 1 in the manuscript for images of the bays.
site.split <- 0
# Values for each site: ANG = 0
#                       CZA = XXX
#                       CZB = YYY
#                       CZD = ZZZ


#### 2. Data handeling ####

### A: Environmental variables
# Access data for all bays
env_var <- read.csv(paste0(PATH, "all_islands_habitat_20190530.csv"), header = TRUE)
# Remove unnessecary columns
env_var <- env_var[,-c(1:2, 4:5, 7, 10:11, 15:16)]

# Subset data to focal site
env_var <- env_var[grepl(site.prefix, substr(env_var$snail_ID, 1, 3)),]
if(grepl("L", CZ)){
  env_var <- env_var[env_var$LCmeanDist < site.split,] # Left of bay
  env_var$LCmeanDist <- abs(env_var$LCmeanDist - site.split)
} else {
  env_var <- env_var[env_var$LCmeanDist >= site.split,] # Right of bay
  env_var$LCmeanDist <- env_var$LCmeanDist - site.split}

# Standardize variables
stnd.var <- function(x){(x-mean(x))/sd(x)}
vars <- c("SN_rock", "low_barnacle", "low_fucus", "topo", "height")
env_var[,vars] <- apply(env_var[,vars], 2, stnd.var)


### B: Inversion count per snail
inv_count <- read.csv(paste0(PATH, site.prefix, "_15july3.csv"),header = T)
# Remove unnessecary columns
inv_count <- inv_count[, c("snail_ID", INV)]


### C: Merge inversion counts and envrionmental variables
dat <- inner_join(inv_count, env_var, by = "snail_ID")

# Rename columns: this prevents an annoying bug in glmulti that looses barnacle*fucus interactions
colnames(dat) <- c(colnames(dat)[1], "Ninv", colnames(dat)[3:4], "rock", "topo", "barn", "fucus")


#### 3. Cline model estimates ####

### Snail positions along hybrid zone
dists <- dat$LCmeanDist

### Get cline fit parameter estimates (from Westram et. al 2021)
params <- read.table(paste0(PATH, site.prefix, "_inv_clines_free_201907",
                            ifelse(grepl("A", CZ), "16", "17"),
                            ".txt"), header = T)
params <- params[params$Inv == INV,]
if(grepl("L", CZ)){params <- params[params$Side == "left",]} else {params <- params[params$Side == "right",]}

### Extrat parameters required for the model
# Cline centre
cnt <- params$Centre
# Cline width
cw <- params$logWidth
# Allele frequency crab end
fc <- params$lp_crab
fc <- exp(fc)/(1+exp(fc)) # Back-transform from logit
# Allele frequency wave end
fw <- params$lp_wave
fw <- exp(fw)/(1+exp(fw))
# Least cost distances
lcd <- dists

### Cline fit function
cf <- fc+(fw-fc)*(1/(1+exp(-4*((lcd-cnt)/exp(cw)))))
# Reverse the fit values if the '2' alleles was determined to be an inversion
cf <- if(params$Allele == 2) cf else 1-cf
# Standardise cline fit
cf <- (cf - mean(cf)) / sd(cf)


#### 4. m0: Fitting cline model ####

### Add cline model estimates to data
dat$cline <- cf

### Run GLM on cline fit
# GLM function - some inversions are missing cline fit parameters. The if satement
# will replace these with a flat null model (y ~ 1).
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

write.csv(m0.AIC, paste0(OUT, "AICs/", CZ, "_", INV, "_AIC_m0.csv"), row.names = FALSE, quote = c(4))


#### 5. m1: Environmental + Cline fit model ####

### Set up custom function for glmulti
# Note: due to a problem with the way functions are nested in R, I have to have
# call this function outside the wrapper function for the env+cline model.
myglm <- function(y, data) {
  if(missing(data)) data <- envrionment(y)

  # Null model if fit contains cline interactions
  nullos <- glm(cbind(Ninv, 2-Ninv)~1,
                family = binomial(link = "logit"), data = data)

  # Get a list of all possible terms in the model
  termz = terms(y,data=data)
  # Order of all terms
  orderz = attr(termz,"order")
  # All pairwise interactions, if any. Otherwise the formula is okay
  intz = which(orderz==2)
  # Locate the row corresponding to the undesired variable
  index=which(dimnames(attr(termz,"factors"))[[1]] == "cline")
  # we simply test that all interactions exclude the effect
  # otherwise we return the crappy null model
  if(length(orderz > 1)){
    if(sum(attr(termz,"factors")[index,intz])!=0) return(nullos)
  } else return(nullos)
  # if all is good we just call glm as usual
  return(glm(formula=y, family = binomial(link = "logit"), data = data))
}

### Run GLM
# Function for full model including cline fit
f <- formula(cbind(Ninv, 2-Ninv)~height*rock*barn*fucus*topo*cline)
# Run the glm with the modified fit function
fits <- do.call("glmulti", list(y = f, data = dat, fitfunc = myglm))

# Save plots of fit behavior
pdf(file = paste0(OUT, "fit_plots/", CZ, "_", INV, "_m1fit.png"))
par(mfrow = c(3,  1))
plot(fits, type = "s")
plot(fits, type = "p")
plot(fits, type = "w")
dev.off()

# Save AIC of all models
m1.AIC <- data.frame("Inv" = rep(INV, length(fits@formulas)),
                      "Site" = rep(CZ, length(fits@formulas)),
                      "rank" = 1:length(fits@formulas),
                      "aic" = sapply(fits@objects, function(X){X$aic}),
                      "dev" = sapply(fits@objects, function(X){X$deviance}),
                      "model" = paste("Pinv", sapply(fits@formulas, function(X){X[[3]]}), sep = " ~ "))

write.csv(m1.AIC, paste0(OUT, "AICs/", CZ, "_", INV, "_AIC_m1.csv"), row.names = FALSE, quote = c(4))

### Aggregate model
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
