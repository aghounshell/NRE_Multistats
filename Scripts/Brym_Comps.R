### Script to calculate yearly average and SD for BEPOM parameters for comparisons to Brym et al. 2014
### Data from Brym et al., 2014
### A Hounshell, 01 June 2020

# Load in libraries
pacman::p_load(tidyverse,ggpubr)

# Load in data: from Brym et al. 2014 and from this study
brym <- read.csv('C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Brym_Comps.csv')
database <- read.csv('C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Database_DOSat.csv')

# Separate database by station (and all samples) then calculate average and stdev
# SURFACE SAMPLES ONLY!!!!!
sta0 <- database %>% filter(Station == "0" & Depth == "S")
sta20 <- database %>% filter(Station == "20" & Depth == "S")
sta30 <- database %>% filter(Station == "30" & Depth == "S")
sta50 <- database %>% filter(Station == "50" & Depth == "S")
sta60 <- database %>% filter(Station == "60" & Depth == "S")
sta70 <- database %>% filter(Station == "70" & Depth == "S")
sta100 <- database %>% filter(Station == "100" & Depth == "S")
sta120 <- database %>% filter(Station == "120" & Depth == "S")
sta140 <- database %>% filter(Station == "140"& Depth == "S")
sta160 <- database %>% filter(Station == "160" & Depth == "S")
sta180 <- database %>% filter(Station == "180" & Depth == "S")

# Create matrix of avg and sd for each station and each parameter
data15 <- matrix(NA,nrow = 11, ncol=7)
colnames(data15) <- c("Station","Mean_Sal","SD_Sal","Mean_Chla","SD_Chla","Mean_a254p","SD_a254p")
Station <-  c("0","20","30","50","60","70","100","120","140","160","180")
# then cbind station