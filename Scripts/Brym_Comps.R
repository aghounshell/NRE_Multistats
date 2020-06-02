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
data15 <- matrix(NA,nrow = 11, ncol=6)
colnames(data15) <- c("Mean_Sal","SD_Sal","Mean_Chla","SD_Chla","Mean_a254p","SD_a254p")
Station <-  as.numeric(c("0","20","30","50","60","70","100","120","140","160","180"))

data15 <- cbind.data.frame(Station,data15)

data15$Mean_Sal[1] <- mean(sta0$Sal,na.rm=TRUE)
data15$SD_Sal[1] <- sd(sta0$Sal,na.rm=TRUE)
data15$Mean_Chla[1] <- mean(sta0$Chla,na.rm=TRUE)
data15$SD_Chla[1] <- sd(sta0$Chla,na.rm=TRUE)
data15$Mean_a254p[1] <- mean(sta0$a254_POM,na.rm=TRUE)
data15$SD_a254p[1] <- sd(sta0$a254_POM,na.rm=TRUE)

data15$Mean_Sal[2] <- mean(sta20$Sal,na.rm=TRUE)
data15$SD_Sal[2] <- sd(sta20$Sal,na.rm=TRUE)
data15$Mean_Chla[2] <- mean(sta20$Chla,na.rm=TRUE)
data15$SD_Chla[2] <- sd(sta20$Chla,na.rm=TRUE)
data15$Mean_a254p[2] <- mean(sta20$a254_POM,na.rm=TRUE)
data15$SD_a254p[2] <- sd(sta20$a254_POM,na.rm=TRUE)

data15$Mean_Sal[3] <- mean(sta30$Sal,na.rm=TRUE)
data15$SD_Sal[3] <- sd(sta30$Sal,na.rm=TRUE)
data15$Mean_Chla[3] <- mean(sta30$Chla,na.rm=TRUE)
data15$SD_Chla[3] <- sd(sta30$Chla,na.rm=TRUE)
data15$Mean_a254p[3] <- mean(sta30$a254_POM,na.rm=TRUE)
data15$SD_a254p[3] <- sd(sta30$a254_POM,na.rm=TRUE)

data15$Mean_Sal[4] <- mean(sta50$Sal,na.rm=TRUE)
data15$SD_Sal[4] <- sd(sta50$Sal,na.rm=TRUE)
data15$Mean_Chla[4] <- mean(sta50$Chla,na.rm=TRUE)
data15$SD_Chla[4] <- sd(sta50$Chla,na.rm=TRUE)
data15$Mean_a254p[4] <- mean(sta50$a254_POM,na.rm=TRUE)
data15$SD_a254p[4] <- sd(sta50$a254_POM,na.rm=TRUE)

data15$Mean_Sal[5] <- mean(sta60$Sal,na.rm=TRUE)
data15$SD_Sal[5] <- sd(sta60$Sal,na.rm=TRUE)
data15$Mean_Chla[5] <- mean(sta60$Chla,na.rm=TRUE)
data15$SD_Chla[5] <- sd(sta60$Chla,na.rm=TRUE)
data15$Mean_a254p[5] <- mean(sta60$a254_POM,na.rm=TRUE)
data15$SD_a254p[5] <- sd(sta60$a254_POM,na.rm=TRUE)

data15$Mean_Sal[6] <- mean(sta70$Sal,na.rm=TRUE)
data15$SD_Sal[6] <- sd(sta70$Sal,na.rm=TRUE)
data15$Mean_Chla[6] <- mean(sta70$Chla,na.rm=TRUE)
data15$SD_Chla[6] <- sd(sta70$Chla,na.rm=TRUE)
data15$Mean_a254p[6] <- mean(sta70$a254_POM,na.rm=TRUE)
data15$SD_a254p[6] <- sd(sta70$a254_POM,na.rm=TRUE)

data15$Mean_Sal[7] <- mean(sta100$Sal,na.rm=TRUE)
data15$SD_Sal[7] <- sd(sta100$Sal,na.rm=TRUE)
data15$Mean_Chla[7] <- mean(sta100$Chla,na.rm=TRUE)
data15$SD_Chla[7] <- sd(sta100$Chla,na.rm=TRUE)
data15$Mean_a254p[7] <- mean(sta100$a254_POM,na.rm=TRUE)
data15$SD_a254p[7] <- sd(sta100$a254_POM,na.rm=TRUE)

data15$Mean_Sal[8] <- mean(sta120$Sal,na.rm=TRUE)
data15$SD_Sal[8] <- sd(sta120$Sal,na.rm=TRUE)
data15$Mean_Chla[8] <- mean(sta120$Chla,na.rm=TRUE)
data15$SD_Chla[8] <- sd(sta120$Chla,na.rm=TRUE)
data15$Mean_a254p[8] <- mean(sta120$a254_POM,na.rm=TRUE)
data15$SD_a254p[8] <- sd(sta120$a254_POM,na.rm=TRUE)

data15$Mean_Sal[9] <- mean(sta140$Sal,na.rm=TRUE)
data15$SD_Sal[9] <- sd(sta140$Sal,na.rm=TRUE)
data15$Mean_Chla[9] <- mean(sta140$Chla,na.rm=TRUE)
data15$SD_Chla[9] <- sd(sta140$Chla,na.rm=TRUE)
data15$Mean_a254p[9] <- mean(sta140$a254_POM,na.rm=TRUE)
data15$SD_a254p[9] <- sd(sta140$a254_POM,na.rm=TRUE)

data15$Mean_Sal[10] <- mean(sta160$Sal,na.rm=TRUE)
data15$SD_Sal[10] <- sd(sta160$Sal,na.rm=TRUE)
data15$Mean_Chla[10] <- mean(sta160$Chla,na.rm=TRUE)
data15$SD_Chla[10] <- sd(sta160$Chla,na.rm=TRUE)
data15$Mean_a254p[10] <- mean(sta160$a254_POM,na.rm=TRUE)
data15$SD_a254p[10] <- sd(sta160$a254_POM,na.rm=TRUE)

data15$Mean_Sal[11] <- mean(sta180$Sal,na.rm=TRUE)
data15$SD_Sal[11] <- sd(sta180$Sal,na.rm=TRUE)
data15$Mean_Chla[11] <- mean(sta180$Chla,na.rm=TRUE)
data15$SD_Chla[11] <- sd(sta180$Chla,na.rm=TRUE)
data15$Mean_a254p[11] <- mean(sta180$a254_POM,na.rm=TRUE)
data15$SD_a254p[11] <- sd(sta180$a254_POM,na.rm=TRUE)

## Plot Brym and 2015 data on same plot as points
data15 <- as.data.frame(data15)
data15 <- data15 %>% arrange(Station)

sal <- ggplot()+
  geom_point(data15,mapping=aes(x=Station,y=Mean_Sal,color="2015-2016"),size=3)+
  geom_line(data15,mapping=aes(x=Station,y=Mean_Sal,color="2015-2016"),size=1)+
  geom_errorbar(data15,mapping=aes(x=Station,y=Mean_Sal,ymin=Mean_Sal-SD_Sal,ymax=Mean_Sal+SD_Sal,
                                   color="2015-2016"),size=1)+
  geom_point(brym,mapping=aes(x=Station,y=Mean_Sal,color="2011-2012"),size=3)+
  geom_line(brym,mapping=aes(x=Station,y=Mean_Sal,color="2011-2012"),size=1)+
  geom_errorbar(brym,mapping=aes(x=Station,y=Mean_Sal,ymin=Mean_Sal-SD_Sal,ymax=Mean_Sal+SD_Sal,
                                 color="2011-2012"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2011-2012"), values=c('#173f5f','#3caea3'))+
  ylab("Salinity")+
  theme_classic(base_size=15)+
  theme(legend.title = element_blank())

chla <- ggplot()+
  geom_point(data15,mapping=aes(x=Station,y=Mean_Chla,color="2015-2016"),size=3)+
  geom_line(data15,mapping=aes(x=Station,y=Mean_Chla,color="2015-2016"),size=1)+
  geom_errorbar(data15,mapping=aes(x=Station,y=Mean_Chla,ymin=Mean_Chla-SD_Chla,ymax=Mean_Chla+SD_Chla,
                                   color="2015-2016"),size=1)+
  geom_point(brym,mapping=aes(x=Station,y=Mean_Chla,color="2011-2012"),size=3)+
  geom_line(brym,mapping=aes(x=Station,y=Mean_Chla,color="2011-2012"),size=1)+
  geom_errorbar(brym,mapping=aes(x=Station,y=Mean_Chla,ymin=Mean_Chla-SD_Chla,ymax=Mean_Chla+SD_Chla,
                                 color="2011-2012"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2011-2012"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("Chla (", mu,"g L"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.title = element_blank())

a254 <- ggplot()+
  geom_point(data15,mapping=aes(x=Station,y=Mean_a254p,color="2015-2016"),size=3)+
  geom_line(data15,mapping=aes(x=Station,y=Mean_a254p,color="2015-2016"),size=1)+
  geom_errorbar(data15,mapping=aes(x=Station,y=Mean_a254p,ymin=Mean_a254p-SD_Sal,ymax=Mean_a254p+SD_a254p,
                                   color="2015-2016"),size=1)+
  geom_point(brym,mapping=aes(x=Station,y=Mean_a254p,color="2011-2012"),size=3)+
  geom_line(brym,mapping=aes(x=Station,y=Mean_a254p,color="2011-2012"),size=1)+
  geom_errorbar(brym,mapping=aes(x=Station,y=Mean_a254p,ymin=Mean_a254p-SD_a254p,ymax=Mean_a254p+SD_a254p,
                                 color="2011-2012"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2011-2012"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("a254 (m"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.title = element_blank())

ggarrange(sal,chla,a254,common.legend=TRUE,ncol=2,nrow=2)

ggplot()+
  geom_point(data15,mapping=aes(x=Station,y=Mean_a254p,color="2015-2016"),size=3)+
  geom_line(data15,mapping=aes(x=Station,y=Mean_a254p,color="2015-2016"),size=1)+
  geom_errorbar(data15,mapping=aes(x=Station,y=Mean_a254p,ymin=Mean_a254p-SD_Sal,ymax=Mean_a254p+SD_a254p,
                                   color="2015-2016"),size=1)+
  geom_point(brym,mapping=aes(x=Station,y=Mean_a254p,color="2011-2012"),size=3)+
  geom_line(brym,mapping=aes(x=Station,y=Mean_a254p,color="2011-2012"),size=1)+
  geom_errorbar(brym,mapping=aes(x=Station,y=Mean_a254p,ymin=Mean_a254p-SD_a254p,ymax=Mean_a254p+SD_a254p,
                                 color="2011-2012"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2011-2012"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("a254 (m"^-1*")")))+
  ylim(-10,20)+
  theme_classic(base_size=15)+
  theme(legend.title = element_blank())


