################# GLM model fit plots: Figure 2 #################
### The following code generates Figure 2 in the manuscript.

### James Reeve - University of Gothenburg
### 2023-04-02

#### 0: Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Packages
library(tidyverse)
library(ggpubr)  # Multi-panel plots

### Set file path to data
PATH1 <- "path/to/inversion/data"
PATH2 <-  "path/to/model/fits"

### Create list of inversion names
inversions <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", "LGC7.1", 
                "LGC7.2", "LGC9.1", "LGC10.1", "LGC10.2", "LGC11.1", "LGC14.1",
                "LGC14.3", "LGC17.1")



#### 1: Fig. 2A ####

### A: Access aggregate model fits
# Access function
CZ.access <- function(site, inversion.list, model){
  fits <- lapply(inversion.list, function(INV){
    f <- paste0(PATH2, "aggregate_fits/glmuti_multimodel_inference_", site, "_", INV, "_", model, ".csv")
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
ANG.fits <- CZ.access("ANG", inversions, "m2")
# CZA
CZA_L.fits <- CZ.access("CZA_L", inversions, "m2")
CZA_R.fits <- CZ.access("CZA_R", inversions, "m2")
# CZB
CZB_L.fits <- CZ.access("CZB_L", inversions, "m2")
CZB_R.fits <- CZ.access("CZB_R", inversions, "m2")
# CZD
CZD_L.fits <- CZ.access("CZD_L", inversions, "m2")
CZD_R.fits <- CZ.access("CZD_R", inversions, "m2")


### B: Merge into one data frame
multifit <- do.call(rbind.data.frame, c(ANG.fits, CZA_L.fits, CZA_R.fits, 
                                        CZB_L.fits, CZB_R.fits, CZD_L.fits, CZD_R.fits))
names(multifit) <- c("Env.Var", "Estimate", "sd", "z.value", "p.value", "Importance", "Inv", "site", "model")
multifit$Significant <- multifit$p.value < 0.05

#Remove the intercept
multifit <- multifit[multifit$Env.Var != "(Intercept)", ]


### C: Set order of variables
multifit$Env.Var <- factor(multifit$Env.Var,
                           levels = c("cline", "barn", "fucus", "height", "rock", "topo", 
                                      "barn:fucus", "barn:height", "barn:rock", "barn:topo",
                                      "fucus:height", "fucus:rock", "fucus:topo",
                                      "height:rock", "height:topo", "rock:topo"))


### D: Remove inversions which are part of larger complexes on LG6 & LG14
inversions2 <- inversions[!(inversions %in% c("LGC6.2a", "LGC6.2b", "LGC6.2c",
                                              "LGC14.2a",  "LGC14.2b",  "LGC14.2c"))]


### E: Tile plot of factor significance and importance
Fig2A <- ggplot(multifit[multifit$Inv %in% inversions,])+
  geom_tile(aes(y = site, x = Env.Var, fill = Significant, alpha = Importance))+
  geom_vline(xintercept = c(1.5), col = "black")+
  geom_vline(xintercept = c(6.5), col = "black")+
  facet_wrap(~ factor(Inv, levels = inversions), nrow = 3)+
  labs(fill = "P < 0.05", alpha = "Importance")+
  scale_fill_manual(values = c("#998ec3", "#f1a340"))+
  scale_x_discrete(limits = levels(multifit$Env.Var))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 60, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.position = "bottom")

Fig2A <- annotate_figure(Fig2A,
                left = text_grob("Environmental Factor", face = "bold", size = 24, rot = 90))

# Save panel
ggsave("path/to/Plots/Fig2A.tiff", 
       Fig2A, device = "tiff", dpi = 300, width = 53.34, height = 21, units = "cm")



#### Fig. 2B ####
### Upload AIC of top 100 models
### Note: try only best AIC vs AIC_cline & vs AIC_fix (i.e. AIC assuming centre is at habitat transition!)

### A: Function to access AIC scores
AIC.access <- function(site, inversion.list, model){
  aics <- lapply(inversions, function(INV){
    f <- paste0(PATH2, "AICs/", site, "_", INV, "_AIC_", model,".csv")
    if(file.exists(f)){
      tmp <- read.csv(f, stringsAsFactors = FALSE)
      return(tmp)
    }
  })
  aics <- aics[aics != 'NULL'] # Remove missing data
}


### B: Access AIC data in a loop
AIC.dat <- lapply(c("ANG", "CZA_L", "CZA_R", "CZB_L", "CZB_R", "CZD_L", "CZD_R"), function(Site){
  # Get AICs
  aic.m0 <- AIC.access(Site, inversions, "m0")
  aic.m1 <- AIC.access(Site, inversions, "m1")
  
  # Merge files for each inversion
  aic.m0 <- do.call(rbind.data.frame, aic.m0)
  aic.m1 <- do.call(rbind.data.frame, aic.m1)
  
  # Add flag for model and site
  aic.m1$analysis <- paste(Site, "m1", sep = "_")
  
  # Rename output
  dAIC <- aic.m1
  
  # Calculate difference in AIC m1 - m0
  dAIC$dAIC <- 0
  for(i in inversions) {
    nullaic <- aic.m0[aic.m0$Inv == i, "aic"]
    dAIC[dAIC$Inv == i,"dAIC"] <- dAIC[dAIC$Inv == i,"aic"] - nullaic
  }; rm(i)
  
  return(dAIC)
})

AIC.dat <- do.call(rbind.data.frame, AIC.dat)

# Remove sub inversions in complexes on LG6 & LG14
AIC.dat <- AIC.dat[!(AIC.dat$Inv %in% c("LGC6.2a", "LGC6.2b", "LGC6.2c",
                                        "LGC14.2a",  "LGC14.2b",  "LGC14.2c")),]

# Set a threshold for "good fits". This is needed to colour plot by ∆AIC
AIC.dat <- AIC.dat %>% group_by(Site, Inv) %>% mutate("Good_Fit" = min(dAIC) < -10)


### C: Viloin plots of ∆AIC
Fig3B <- ggplot(AIC.dat, aes(y = Site, x = dAIC))+
  geom_rect(aes(xmin = -10, xmax = 0, ymin = "ANG", ymax = "CZD_R"), 
            fill = "grey95", lty = 2)+
  geom_vline(xintercept = 0, col = "grey25")+
  geom_violin(aes(fill = Good_Fit, colour = Good_Fit), width = 0.8)+
  scale_fill_manual(values = c("#998ec3", "#f1a340"))+
  scale_colour_manual(values = c("#998ec3", "#f1a340"))+
  facet_wrap(vars(factor(Inv, levels = inversions)), ncol = 5, nrow = 3)+
  labs(x = expression(Delta~AIC~" (Env - Cline)"))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = "none")


# Save panel
ggsave("path/to/Plots/Fig2B.tiff", 
       Fig3B, device = "tiff", dpi = 300, width = 53.34, height = 14, units = "cm")
