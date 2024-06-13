################# Fis plots: Figure 3 #################
### Plotting Fis along hybrid zones. 

### James Reeve - University of Gothenburg
### 2024-01-14

#### 0: Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Packages
library(ggplot2)

### Filepath
PATH <- "path/to/inversions/data/"

### List of inversion
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC6.1", 
          "LGC7.1", "LGC7.2", "LGC9.1", "LGC10.1","LGC10.2", 
          "LGC11.1", "LGC14.1", "LGC14.3", "LGC17.1")



#### 1: Plots across whole hybrid zones ####

### Access data
datHZ <- read.csv(paste0(PATH, "CZ_Fis_across_hybrid_zones.csv"))

### Plot
pHZ <- ggplot(datHZ, aes(y = Inv, x = Fis, shape = Site, fill = Site))+
  geom_vline(xintercept = 0)+
  geom_point(aes(size = adj.p.value < 0.05), alpha = 0.6, col = "grey20", show.legend = FALSE)+ 
  facet_wrap(facets = vars(Site.prefix), nrow = 1)+
  scale_size_manual(values = c(3, 6))+
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a",
                               "#b2df8a", "#33a02c", "#33a02c"))+
  scale_shape_manual(values = c(21, 24, 21, 24, 21, 24, 21))+
  scale_y_discrete(limits = INVs)+
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1))+
  labs(title = "Whole hybrid zone", x = expression("F"[IS]))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

### Save plot
ggsave("path/to/Plots/Fig3A.tiff",
       pHZ, device = "tiff", width = 320, height = 120, units = "mm")



#### 2: Plots of Crab habitat ####
datC <- read.csv(paste0(PATH, "CZ_Fis_Crab_habitat.csv"))

pC <- ggplot(datC, aes(y = Inv, x = Fis, shape = Site, fill = Site))+
  geom_vline(xintercept = 0)+
  geom_point(aes(size = adj.p.value < 0.05), alpha = 0.6, col = "grey20", show.legend = FALSE)+ 
  facet_wrap(facets = vars(Site.prefix), nrow = 2, ncol = 2)+
  scale_size_manual(values = c(3, 6))+
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a",
                               "#b2df8a", "#33a02c", "#33a02c"))+
  scale_shape_manual(values = c(21, 24, 21, 24, 21, 24, 21))+
  scale_y_discrete(limits = INVs)+
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1))+
  labs(title = "Crab habitat", x = expression("F"[IS]))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave("path/to/Plots/Fig3B.tiff",
       pC, device = "tiff", width = 180, height = 240, units = "mm")



#### 3: Plots of Wave habitat ####
datW <- read.csv(paste0(PATH, "CZ_Fis_Wave_habitat.v2.csv"))
# Adjustment for an extremly negative Fis
datW[datW$Fis <-1, "Fis"] <- -1

pW <- ggplot(datW, aes(y = Inv, x = Fis, shape = Site, fill = Site))+
  geom_vline(xintercept = 0)+
  geom_point(aes(size = adj.p.value < 0.05), alpha = 0.6, col = "grey20", show.legend = FALSE)+ 
  facet_wrap(facets = vars(Site.prefix), nrow = 2, ncol = 2)+
  scale_size_manual(values = c(3, 6))+
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a",
                               "#b2df8a", "#33a02c", "#33a02c"))+
  scale_shape_manual(values = c(21, 24, 21, 24, 21, 24, 21))+
  scale_y_discrete(limits = INVs)+
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1))+
  labs(title = "Wave habitat", x = expression("F"[IS]))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave("path/to/Plots/Fig3C.tiff",
       pW, device = "tiff", width = 180, height = 240, units = "mm")
