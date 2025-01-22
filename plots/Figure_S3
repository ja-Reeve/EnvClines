################# Figure S3: Variance explained between habitats and model importance ##############
### This script is a test of a specific hypothesis in our manuscript. We assumed that model terms which
### varaied among habitats should tend to be more important (i.e., represented among the top 100 
### submodels selected by glmulti) in cline+environment models. There were a few hypothetical outcomes:
###       H0: no association between environmental variance and model importance
###       H1: variation between habitats is strongly correlated with importance
###       H2a: variance within just Crab habitat is correlated with importance
###       H2b: varaince within just the Wave habitat is correlated with importance
###       H3: variance within both habitats is correlated with importance
### These hypotheses were tested by comparing the varaince explained by environmental variable clines
### (i.e., clines estimated from just the envrionmental variables; data given in Table S5 of Westram 
### et al. 2021) with the importance scores from the cline+environment models. These were tested using 
### correlation tests.

### James Reeve
### 31/12/2024

#### 0: Preparation ####
rm(list = ls())
dev.off()

### Packages
library(ggplot2)

### Filepaths
PATH <- "path/to/model/fits"

### Parameters
# List of inversion names
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", "LGC7.2", "LGC9.1", 
          "LGC10.1", "LGC10.2", "LGC11.1", "LGC14.1", "LGC14.3", "LGC17.1")
# List of hybrid zones
CZs <- c("CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R", "ANG")



#### 1. Access glmulti results ####

### Access aggregate model fits
## Access function
CZ.access <- function(site, inversion.list, model){
  fits <- lapply(inversion.list, function(INV){
    f <- paste0(PATH, "Cline_fits/aggregate_fits/glmuti_multimodel_inference_", site, "_", INV, "_", model, ".v4.csv")
    if(file.exists(f)){
      tmp <- read.csv(f, stringsAsFactors = FALSE)
      tmp$Inv <- INV
      tmp$site <- site
      tmp$model <- model
      return(tmp)
    }
  })
  fits <- fits[fits != 'NULL']
}

# ANG
ANG.fits <- CZ.access("ANG", INVs, "m2")
# CZA
CZA_L.fits <- CZ.access("CZA_L", INVs, "m2")
CZA_R.fits <- CZ.access("CZA_R", INVs, "m2")
# CZB
CZB_L.fits <- CZ.access("CZB_L", INVs, "m2")
CZB_R.fits <- CZ.access("CZB_R", INVs, "m2")
# CZD
CZD_L.fits <- CZ.access("CZD_L", INVs, "m2")
CZD_R.fits <- CZ.access("CZD_R", INVs, "m2")

## Merge into one data frame
multifit <- do.call(rbind.data.frame, c(ANG.fits, CZA_L.fits, CZA_R.fits, CZB_L.fits, CZB_R.fits, CZD_L.fits, CZD_R.fits))
names(multifit) <- c("Env.Var", "Estimate", "sd", "z.value", "p.value", "Importance", "Inv", "site", "model")
multifit$Significant <- multifit$p.value < 0.05

#Remove the intercept
multifit <- multifit[multifit$Env.Var != "(Intercept)", ]

## Set order for variables
multifit$Env.Var <- factor(multifit$Env.Var,
                           levels = c("barn", "fucus", "height", "rock", "topo", 
                                      "barn:fucus", "barn:height", "barn:rock", "barn:topo",
                                      "fucus:height", "fucus:rock", "fucus:topo",
                                      "height:rock", "height:topo", "rock:topo"))


## Subset to key environmental varaibles
multifit <- multifit[grep(":", multifit$Env.Var, invert = TRUE),]
multifit <- multifit[complete.cases(multifit),]



#### 2. Access environmental clines ####
### This was coppied from Table S5 of Westram et al. 2021.
### Note that the order of variables was different for ANG.

d <- data.frame("site" = c("ANG", "ANG", "ANG", "ANG", "ANG",
                      "CZA_L", "CZA_L", "CZA_L", "CZA_L", "CZA_L",
                      "CZA_R", "CZA_R", "CZA_R", "CZA_R", "CZA_R",
                      "CZB_L", "CZB_L", "CZB_L", "CZB_L", "CZB_L",
                      "CZB_R", "CZB_R", "CZB_R", "CZB_R", "CZB_R",
                      "CZD_L", "CZD_L", "CZD_L", "CZD_L", "CZD_L",
                      "CZD_R", "CZD_R", "CZD_R", "CZD_R", "CZD_R"),
           "variable" = c("rock", "barn", "fucus", "topo", "height",
                          "barn", "fucus", "topo", "rock", "height",
                          "barn", "fucus", "topo", "rock", "height",
                          "barn", "fucus", "topo", "rock", "height",
                          "barn", "fucus", "topo", "rock", "height",
                          "barn", "fucus", "topo", "rock", "height",
                          "barn", "fucus", "topo", "rock", "height"),
           "V.exp" = c(0.860, 0.250, 0.110, 0.310, 0.162,
                       0.686, 0.025, 0.189, 0.739, 0.277,
                       0.922, 0.393, 0.007, 0.853, 0.037,
                       0.432, 0.182, 0.263, 0.898, 0.263,
                       0.707, 0.116, 0.059, 0.955, 0.059,
                       0.452, 0.100, 0.532, 0.864, 0.443,
                       0.639, 0.040, 0.026, 0.785, 0.006))
Env.exp <- read.csv(paste0(PATH, "Env_cline_varExp.csv"))


### Merge with aggregate model
tmp <- merge(multifit, Env.exp, by.x = c("site", "Env.Var"), by.y = c("site", "variable") )



#### 3. Figure S3A-B: Scatter plots ####

### Figure S3A: All points coloured by hybrid zone
p1 <- ggplot(tmp, aes(x = Importance, y = V.exp, group = Env.Var, colour = Env.Var))+
  geom_point(aes(shape = Env.Var), size = 2)+
  geom_smooth(aes(linetype = Env.Var), method = "glm", se = FALSE)+
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
  scale_colour_manual(values = c("#b4bdd4ff", "#aa9d62ff", "#5b6fa0ff", "grey80","grey10"))+
  scale_shape_manual(values = c(15, 16, 15, 16, 1))+
  scale_linetype_manual(values = c(1, 1, 1, 3, 3))+
  labs(x = "Importance", y = "Variance Explained")+
  theme_bw()

# Figure S3B: Scatterplots split into panels for each inversion
p2 <- ggplot(tmp, aes(x = V.exp, y = Importance))+
  geom_point()+
  geom_smooth(method = "glm", se = F)+
  facet_grid(rows = vars(factor(Inv, levels = INVs)), cols = vars(Env.Var))+
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
  labs(x = "Variance explained by cline model", y = "Importance")+
  theme_bw()

### Grid of all hybrid zones
# Exluding, as all plots are vertical lines. You're welcome to test it out.
#ggplot(tmp, aes(x = V.exp, y = Importance))+
#  geom_point()+
#  geom_smooth(method = "glm", se = F)+
#  facet_grid(rows = vars(factor(site, levels = CZs)), cols = vars(Env.Var))+
#  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
#  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))+
#  labs(x = "Variance explained by cline model", y = "Importance")+
#  theme_bw()



#### 4. Figure S3C: correlations ####
### To get an overall idea of the relationship between variance explained across hybrid zones and environmental 
### variable importance, correlation scores were calculated for each inversion and envrionmental variable. If there 
### was a balanced number of positive and negative correlations, then there was no clear relationship between 
### environmental clines and model importance.

### Loop through inversions
Inv.loop <- lapply(INVs, function(inv){
  
  # Loop through environmental variables
  Env.loop <- lapply(unique(tmp$Env.Var), function(env){
    tmp2 <- tmp[tmp$Inv == inv & tmp$Env.Var == env,]
    # Correlation
    Corr <- cor.test(tmp2$V.exp, tmp2$Importance)
    # Save as dataframe
    res <- data.frame("inversion" = inv, 
                      "variable" = env,
                      "r" = Corr$estimate, 
                      "t" = Corr$statistic,
                      "df" = Corr$parameter,
                      "p.value" = Corr$p.value)
    
    return(res)
  })
  
  # Merge environmental variable results
  Env.Corr <- do.call(rbind.data.frame, Env.loop)
  
  return(Env.Corr)
})

corr.dat <- do.call(rbind.data.frame, Inv.loop)

### Check how many correlations are +ve or -ve
sum(corr.dat$r < 0) # 35 negative
sum(corr.dat$r > 0) # 35 positive

# barn
sum(corr.dat$r[corr.dat$variable == "barn"] < 0) # 5 negative
sum(corr.dat$r[corr.dat$variable == "barn"] > 0) # 9 positive

# fucus
sum(corr.dat$r[corr.dat$variable == "fucus"] < 0) # 9 negative
sum(corr.dat$r[corr.dat$variable == "fucus"] > 0) # 5 positive

# height
sum(corr.dat$r[corr.dat$variable == "height"] < 0) # 7 negative
sum(corr.dat$r[corr.dat$variable == "height"] > 0) # 7 positive

# topo
sum(corr.dat$r[corr.dat$variable == "topo"] < 0) # 8 negative
sum(corr.dat$r[corr.dat$variable == "topo"] > 0) # 6 positive

# rock
sum(corr.dat$r[corr.dat$variable == "rock"] < 0) # 6 negative
sum(corr.dat$r[corr.dat$variable == "rock"] > 0) # 8 positive

### Tile plot of correlations
p3 <- ggplot(corr.dat)+
  geom_tile(aes(x = variable, y = factor(inversion, levels = rev(INVs)), fill = r, alpha = 1-p.value), colour = "grey30")+
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0)+
  labs(x = "Environmental variable", y = "")+
  theme_classic()



#### 5. Multipanel plot ####
FigS3 <- ggarrange(p1, ggarrange(p2, p3, widths = c(0.6, 0.4)), ncol = 1, heights = c(0.3, 0.7))

ggsave("path/to/Plots/FigS3.tiff", FigS3, device = "tiff", dpi = 300, width = 53.34, height = 14, units = "cm")

ggsave("/Users/james/Dropbox/PhD/CZ_env_analysis/Plots/supplementary_figures/Fig.S3.tiff",
       FigS3, width = 28, height = 18, units = "cm", dpi = 150)


