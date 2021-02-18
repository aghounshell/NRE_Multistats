# Script to plot box plots for 2015-2016 and calculate long term median (2000-2019)
# A Hounshell, 11 Mar 2020
# Updated: 09 Feb 2021 in response to revisions
# Removed long-term data (13 Feb 2021) - doesn't add!

# Load in libraries
pacman::p_load(tidyverse,PerformanceAnalytics,GGally,dplyr,ggpubr,ggplot2,akima,lubridate,colorRamps,
               RColorBrewer,rLakeAnalyzer)

# Load in data (Database_DOSat.csv)
my_data <- read.csv('C:/Users/ahoun/Desktop/NRE_Multistats/Data/Database_DOSat.csv')
my_data$Date <- as.POSIXct(strptime(my_data$Date, "%m/%d/%Y", tz="EST"))
# Separate into surface and bottom
mydata_all_s <- my_data %>% 
  filter(Depth == "S")
mydata_all_b <- my_data %>% 
  filter(Depth == "B")
# Calculate water density
mydata_all_s <- mydata_all_s %>% 
  mutate(density = water.density(mydata_all_s$Temp,sal = mydata_all_s$Sal))
mydata_all_b <- mydata_all_b %>% 
  mutate(density = water.density(mydata_all_b$Temp,sal = mydata_all_b$Sal))

Strat = mydata_all_b$density - mydata_all_s$density
Strat <- cbind.data.frame(mydata_all_s$Date,mydata_all_s$Season,mydata_all_s$Station,Strat)
Strat <- Strat %>% 
  rename("Date" = "mydata_all_s$Date","Season" = "mydata_all_s$Season","Station" = "mydata_all_s$Station")
Strat <- Strat[complete.cases(Strat),]
Strat <- Strat %>% 
  mutate(Season = ifelse(Date < '2015-09-01' & Season == "Summer","Summer15",
                         ifelse(Date > '2016-06-01' & Season == "Summer","Summer16",
                                Season)))
Strat$Season<-factor(Strat$Season, levels=c("Summer15", "Fall", "Winter", "Spring","Summer16"))

# Plot be season and by station
sum15_strat <- Strat %>% 
  filter(Season == "Summer15")
fall_strat <- Strat %>% 
  filter(Season == "Fall")
winter_strat <- Strat %>% 
  filter(Season == "Winter")
spring_strat <- Strat %>% 
  filter(Season == "Spring")
sum16_strat <- Strat %>% 
  filter(Season == "Summer16")

jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS10.jpg",width=400,height=350,units="mm",res=800)

sum15 <- ggplot(sum15_strat,mapping=aes(as.factor(Station),Strat,fill=Season))+
  geom_boxplot()+
  xlab("Distance down estuary (km)")+
  ylab(expression(paste("Stratification (kg m"^-3*")")))+
  scale_fill_manual(values=c("#D81B60"))+
  ylim(c(0,10))+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

fall <- ggplot(fall_strat,mapping=aes(as.factor(Station),Strat,fill=Season))+
  geom_boxplot()+
  xlab("Distance down estuary (km)")+
  ylab(expression(paste("Stratification (kg m"^-3*")")))+
  scale_fill_manual(values=c("#FFC107"))+
  ylim(c(0,10))+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

winter <- ggplot(winter_strat,mapping=aes(as.factor(Station),Strat,fill=Season))+
  geom_boxplot()+
  xlab("Distance down estuary (km)")+
  ylab(expression(paste("Stratification (kg m"^-3*")")))+
  scale_fill_manual(values=c("#1E88E5"))+
  ylim(c(0,10))+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

spring <- ggplot(spring_strat,mapping=aes(as.factor(Station),Strat,fill=Season))+
  geom_boxplot()+
  xlab("Distance down estuary (km)")+
  ylab(expression(paste("Stratification (kg m"^-3*")")))+
  scale_fill_manual(values=c("#004D40"))+
  ylim(c(0,10))+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

sum16 <- ggplot(sum16_strat,mapping=aes(as.factor(Station),Strat,fill=Season))+
  geom_boxplot()+
  xlab("Distance down estuary (km)")+
  ylab(expression(paste("Stratification (kg m"^-3*")")))+
  scale_fill_manual(values=c("#F7C1BB"))+
  ylim(c(0,10))+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

ggarrange(sum15,fall,winter,spring,sum16,
          nrow=3,ncol=2)

dev.off()

# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]
#my_data2$Date <- as.POSIXct(strptime(my_data2$Date, "%m/%d/%Y", tz="EST"))

# Updated season to reflect summer 2015 and summer 2016
my_data2 <- my_data2 %>% 
  mutate(Season = ifelse(Date < '2015-09-01' & Season == "Summer","Summer15",
                         ifelse(Date > '2016-06-01' & Season == "Summer","Summer16",
                                Season)))

# Plot salinity and chla
# Separate into S and B
mydata_s <- my_data2 %>% 
  filter(Depth == "S") %>% 
  mutate(Date = as.Date(Date))
mydata_b <- my_data2 %>% 
  filter(Depth == "B") %>% 
  mutate(Date=as.Date(Date))

## Calculate median for the yearly data (2015-2016) for each season
med_s <- mydata_s %>% select(Season,Sal,Chla,DOC_mg,POC_mg,,HIX_DOM,HIX_POM,Flushing_Time) %>% group_by(Season) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE)
med_b <- mydata_b %>% select(Season,Sal,Chla,DOC_mg,POC_mg,,HIX_DOM,HIX_POM,Flushing_Time) %>% group_by(Season) %>% 
  summarise_if(is.numeric,median,na.rm=TRUE)

# Define seasons
mydata_s$Season<-factor(mydata_s$Season, levels=c("Summer15", "Fall", "Winter", "Spring","Summer16"))
mydata_b$Season<-factor(mydata_b$Season, levels=c("Summer15", "Fall", "Winter", "Spring","Summer16"))
my_data2$Season<-factor(my_data2$Season, levels=c("Summer15","Fall","Winter","Spring","Summer16"))

# Define depth
my_data2$Depth <- factor(my_data2$Depth, levels=c("S","B"))

############ Heatmaps for spatial/temporal visualizations of key parameters? #######################
## In response to paper revisions: 26 Jan 2021
interp_sal <- interp(x=mydata_s$Date, y=mydata_s$Station, z=mydata_s$Sal,
                     xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                     yo = seq(0,180,by = 1),
                     extrap = F, linear = T, duplicate = "strip")
interp_sal <- interp2xyz(interp_sal,data.frame=T)
interp_sal$Date <- as.Date(interp_sal$x)

interp_sal_b <- interp(x=mydata_b$Date, y=mydata_b$Station, z=mydata_b$Sal,
                       xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                       yo = seq(0,180,by = 1),
                       extrap = F, linear = T, duplicate = "strip")
interp_sal_b <- interp2xyz(interp_sal_b,data.frame=T)
interp_sal_b$Date <- as.Date(interp_sal_b$x)

interp_chla <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$Chla,
                      xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_chla <- interp2xyz(interp_chla,data.frame=T)
interp_chla$Date <- as.Date(interp_chla$x)

interp_chla_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$Chla,
                        xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                        yo = seq(0,180,by = 1),
                        extrap = F, linear = T, duplicate = "strip")
interp_chla_b <- interp2xyz(interp_chla_b,data.frame=T)
interp_chla_b$Date <- as.Date(interp_chla_b$x)

# Plot salinity heatmap + boxplots from above - Figure 3 in MS

jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure3.jpg",width=400,height=440,units="mm",res=800)

# Salinty surface heatmap
sal_s_heat <- ggplot()+
  geom_tile(interp_sal,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,20))+
  labs(x = "",y="Distance (km)",fill="Sal")+
  theme_classic(base_size=25)

# Salinity surface boxplot
# Convert to ggboxplot
sal_s_box <- ggplot(data = mydata_s,aes(Season,Sal))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,20))+
  ylab("Surface Salinity")+
  theme_classic(base_size = 21)

# Salinity bottom heatmap
sal_b_heat <- ggplot()+
  geom_tile(interp_sal_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,20))+
  labs(x = "",y="Distance (km)",fill="Sal")+
  theme_classic(base_size=25)

# Salinity bottom boxplot
sal_b_box <- ggplot(data = mydata_b,aes(Season,Sal))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,20))+
  ylab("Bottom Salinity")+
  theme_classic(base_size = 21)

# Chla surface heatmap
chla_s_heat <- ggplot()+
  geom_tile(interp_chla,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,130))+
  labs(x = "",y="Distance (km)",fill="Chla")+
  theme_classic(base_size=25)

# Chla surface boxplot
ylab.text=expression(paste("Surface Chla (",mu,"g L"^"-1"*")"))
chla_s_box <- ggplot(data = mydata_s,aes(Season,Chla))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,135))+
  ylab(ylab.text)+
  theme_classic(base_size = 21)

# Chla bottom heatmap
chla_b_heat <- ggplot()+
  geom_tile(interp_chla_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,130))+
  labs(x = "",y="Distance (km)",fill="Chla")+
  theme_classic(base_size=25)

# Chla bottom boxplot
ylab.text=expression(paste("Bottom Chla (",mu,"g L"^"-1"*")"))
chla_b_box <- ggplot(data = mydata_b,aes(Season,Chla))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,135))+
  ylab(ylab.text)+
  theme_classic(base_size = 21)

ggarrange(sal_s_heat,sal_s_box,sal_b_heat,sal_b_box,chla_s_heat,chla_s_box,chla_b_heat,chla_b_box,
          nrow=4,ncol=2,widths=c(2,1))

dev.off()

## Update Figure 4 in MS with heatmaps and boxplots of DOC and POC concentration
interp_doc <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$DOC_mg,
                     xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                     yo = seq(0,180,by = 1),
                     extrap = F, linear = T, duplicate = "strip")
interp_doc <- interp2xyz(interp_doc,data.frame=T)
interp_doc$Date <- as.Date(interp_doc$x)

interp_doc_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$DOC_mg,
                       xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                       yo = seq(0,180,by = 1),
                       extrap = F, linear = T, duplicate = "strip")
interp_doc_b <- interp2xyz(interp_doc_b,data.frame=T)
interp_doc_b$Date <- as.Date(interp_doc_b$x)

interp_poc <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$POC_mg,
                     xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                     yo = seq(0,180,by = 1),
                     extrap = F, linear = T, duplicate = "strip")
interp_poc <- interp2xyz(interp_poc,data.frame=T)
interp_poc$Date <- as.Date(interp_poc$x)

interp_poc_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$POC_mg,
                       xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                       yo = seq(0,180,by = 1),
                       extrap = F, linear = T, duplicate = "strip")
interp_poc_b <- interp2xyz(interp_poc_b,data.frame=T)
interp_poc_b$Date <- as.Date(interp_poc_b$x)

## Create Figure 4 - DOC and POC heatmaps + boxplots
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure4.jpg",width=400,height=440,units="mm",res=800)

# DOC heatmap S
docs_heat <- ggplot()+
  geom_tile(interp_doc,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(4,15))+
  labs(x = "",y="Distance (km)",fill="DOC")+
  theme_classic(base_size=25)

# DOC S boxplot
docs_box <- ggplot(data = mydata_s,aes(Season,DOC_mg))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,15))+
  ylab(expression("Surface DOC (mg L"^-1*")"))+
  theme_classic(base_size = 21)

# DOC B heat map
docb_heat <- ggplot()+
  geom_tile(interp_doc_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(4,15))+
  labs(x = "",y="Distance (km)",fill="DOC")+
  theme_classic(base_size=25)

# DOC B boxplot
docb_box <- ggplot(data = mydata_b,aes(Season,DOC_mg))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,15))+
  ylab(expression("Bottom DOC (mg L"^-1*")"))+
  theme_classic(base_size = 21)

# POC heatmap S
pocs_heat <- ggplot()+
  geom_tile(interp_poc,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,6))+
  labs(x = "",y="Distance (km)",fill="POC")+
  theme_classic(base_size=25)

# POC boxplot S
pocs_box <- ggplot(data = mydata_s,aes(Season,POC_mg))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,6))+
  ylab(expression("Surface POC (mg L"^-1*")"))+
  theme_classic(base_size = 21)

# POC heatmap B
pocb_heat <- ggplot()+
  geom_tile(interp_poc_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_vline(xintercept=as.Date("2015-09-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2015-12-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-03-01"),color="grey",size=1)+
  geom_vline(xintercept=as.Date("2016-06-01"),color="grey",size=1)+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.9)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,6))+
  labs(x = "",y="Distance (km)",fill="POC")+
  theme_classic(base_size=25)

# POC boxplot B
pocb_box <- ggplot(data = mydata_b,aes(Season,POC_mg))+
  geom_boxplot()+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  ylim(c(0,6))+
  ylab(expression("Bottom POC (mg L"^-1*")"))+
  theme_classic(base_size = 21)

ggarrange(docs_heat,docs_box,docb_heat,docb_box,pocs_heat,pocs_box,pocb_heat,pocb_box,
          nrow=4,ncol=2,widths=c(2,1))

dev.off()

############ Figure 5 - boxplots of OM indicators with station (C:N, SUVA, BIX, HIX)
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure5.jpg",width=400,height=440,units="mm",res=800)

doctodon <- ggplot()+
  geom_boxplot(data = my_data2,aes(as.factor(Station),DOC_DON,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("DOC:DON")+
  ylim(0,50)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

poctopn <- ggplot()+
  geom_boxplot(data = my_data2,aes(as.factor(Station),POCtoPN,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("POC:PN")+
  theme_classic(base_size = 21)+
  theme(legend.position=c(0.8,0.8))

suva_dom <- ggplot()+
  geom_boxplot(data = my_data2,aes(as.factor(Station),SUVA_DOM,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab(expression("DOM SUVA (L mg"^-1*"C m"^-1*")"))+
  ylim(0,5)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

suva_pom <- ggplot()+
  geom_boxplot(data = my_data2,aes(as.factor(Station),SUVA_POC,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab(expression("POM SUVA (L mg"^-1*"C m"^-1*")"))+
  ylim(0,5)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

bix_dom <- ggplot()+
  geom_hline(yintercept = 0.6,linetype="dashed")+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(as.factor(Station),BIX_DOM,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("DOM BIX")+
  ylim(0,1.3)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

bix_pom <- ggplot()+
  geom_hline(yintercept = 0.6,linetype="dashed")+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(as.factor(Station),BIX_POM,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("POM BIX")+
  ylim(0,1.3)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

hix_dom <- ggplot()+
  geom_hline(yintercept = 6,linetype="dashed")+
  geom_hline(yintercept = 10,linetype="dashed")+
  geom_hline(yintercept = 16,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(as.factor(Station),HIX_DOM,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("DOM HIX")+
  ylim(0,26)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

hix_pom <- ggplot()+
  geom_hline(yintercept = 6,linetype="dashed")+
  geom_hline(yintercept = 10,linetype="dashed")+
  geom_hline(yintercept = 16,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(as.factor(Station),HIX_POM,color=Depth))+
  scale_color_manual(values=c("black","grey53"))+
  xlab("Distance down estuary (km)")+
  ylab("POM HIX")+
  ylim(0,26)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

ggarrange(doctodon,poctopn,suva_dom,suva_pom,bix_dom,bix_pom,hix_dom,hix_pom,
          nrow=4,ncol=2)

dev.off()

############ Figure 5 - boxplots of OM indicators with season: C:N, SUVA, BIX, HIX
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure6_new.jpg",width=400,height=440,units="mm",res=800)

doctodon <- ggplot()+
  geom_boxplot(data = my_data2,aes(Season,DOC_DON,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("DOC:DON")+
  ylim(0,50)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

poctopn <- ggplot()+
  geom_boxplot(data = my_data2,aes(Season,POCtoPN,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("POC:PN")+
  ylim(0,50)+
  theme_classic(base_size = 21)+
  theme(legend.position=c(0.8,0.8))

suva_dom <- ggplot()+
  geom_boxplot(data = my_data2,aes(Season,SUVA_DOM,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab(expression("DOM SUVA (L mg"^-1*"C m"^-1*")"))+
  ylim(0,5)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

suva_pom <- ggplot()+
  geom_boxplot(data = my_data2,aes(Season,SUVA_POC,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab(expression("POM SUVA (L mg"^-1*"C m"^-1*")"))+
  ylim(0,5)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

bix_dom <- ggplot()+
  geom_hline(yintercept = 0.6,linetype="dashed")+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(Season,BIX_DOM,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("DOM BIX")+
  ylim(0,1.3)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

bix_pom <- ggplot()+
  geom_hline(yintercept = 0.6,linetype="dashed")+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(Season,BIX_POM,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("POM BIX")+
  ylim(0,1.3)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

hix_dom <- ggplot()+
  geom_hline(yintercept = 6,linetype="dashed")+
  geom_hline(yintercept = 10,linetype="dashed")+
  geom_hline(yintercept = 16,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(Season,HIX_DOM,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("DOM HIX")+
  ylim(0,26)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

hix_pom <- ggplot()+
  geom_hline(yintercept = 6,linetype="dashed")+
  geom_hline(yintercept = 10,linetype="dashed")+
  geom_hline(yintercept = 16,linetype="dashed")+
  geom_boxplot(data = my_data2,aes(Season,HIX_POM,color=Depth))+
  scale_x_discrete(labels=c("Summer15" = "Sum '15","Fall","Winter","Spring","Summer16" = "Sum '16"))+
  scale_color_manual(values=c("black","grey53"))+
  ylab("POM HIX")+
  ylim(0,26)+
  theme_classic(base_size = 21)+
  theme(legend.position="none")

ggarrange(doctodon,poctopn,suva_dom,suva_pom,bix_dom,bix_pom,hix_dom,hix_pom,
          nrow=4,ncol=2)

dev.off()

#### Plot HIX vs. BIX Biplots for Supplementary (Fig S3)
# Updated with Summer 2015/Summer 2016 seasons
# Updated to include surface and bottom
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS3.jpg",width=400,height=250,units="mm",res=800)

hixbix_dom_s <- ggplot(mydata_s,mapping=aes(x=BIX_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season))+
  geom_point(size=2.5)+
  scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_shape_manual(values=c(24,21,22,23,25))+
  ylab("Surface DOM HIX")+
  xlab("Surface DOM BIX")+
  xlim(c(0,1.3))+
  ylim(c(0,30))+
  theme_classic(base_size=20)

hixbix_dom_b <- ggplot(mydata_b,mapping=aes(x=BIX_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season))+
  geom_point(size=2.5)+
  scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_shape_manual(values=c(24,21,22,23,25))+
  ylab("Bottom DOM HIX")+
  xlab("Bottom DOM BIX")+
  xlim(c(0,1.3))+
  ylim(c(0,30))+
  theme_classic(base_size=20)

hixbix_pom_s <- ggplot(mydata_s,mapping=aes(x=BIX_POM,y=HIX_POM,color=Season,shape=Season,fill=Season))+
  geom_point(size=2.5)+
  scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_shape_manual(values=c(24,21,22,23,25))+
  ylab("Surface POM HIX")+
  xlab("Surface POM BIX")+
  xlim(c(0,1.3))+
  ylim(c(0,30))+
  theme_classic(base_size=20)

hixbix_pom_b <- ggplot(mydata_b,mapping=aes(x=BIX_POM,y=HIX_POM,color=Season,shape=Season,fill=Season))+
  geom_point(size=2.5)+
  scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
  scale_shape_manual(values=c(24,21,22,23,25))+
  ylab("Bottom POM HIX")+
  xlab("Bottom POM BIX")+
  xlim(c(0,1.3))+
  ylim(c(0,30))+
  theme_classic(base_size=20)

ggarrange(hixbix_dom_s,hixbix_dom_b,hixbix_pom_s,hixbix_pom_b,
          nrow=2,ncol=2)

dev.off()

######################## OLD CODE ##############################
# Interp DO
interp_do <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$DO_Sat,
                    xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                    yo = seq(0,180,by = 1),
                    extrap = F, linear = T, duplicate = "strip")
interp_do <- interp2xyz(interp_do,data.frame=T)
interp_do$Date <- as.Date(interp_do$x)

interp_do_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$DO_Sat,
                      xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_do_b <- interp2xyz(interp_do_b,data.frame=T)
interp_do_b$Date <- as.Date(interp_do_b$x)

# Plot DO heatmaps - probably not going to use these?
ggplot()+
  geom_tile(interp_do,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray")+
  labs(x = "",y="Distance down estuary (km)",fill="% DO")+
  theme_classic(base_size=15)

ggplot()+
  geom_tile(interp_do_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray")+
  labs(x = "",y="Distance down estuary (km)",fill="% DO")+
  theme_classic(base_size=15)

### Then interpret and plot various Fl and Abs parameters? What would be best here?
# Start with BIX and HIX
interp_hixd <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$HIX_DOM,
                      xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_hixd <- interp2xyz(interp_hixd,data.frame=T)
interp_hixd$Date <- as.Date(interp_hixd$x)

interp_hixd_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$HIX_DOM,
                      xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_hixd_b <- interp2xyz(interp_hixd_b,data.frame=T)
interp_hixd_b$Date <- as.Date(interp_hixd_b$x)

interp_hixp <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$HIX_POM,
                      xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_hixp <- interp2xyz(interp_hixp,data.frame=T)
interp_hixp$Date <- as.Date(interp_hixp$x)

interp_hixp_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$HIX_POM,
                      xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_hixp_b <- interp2xyz(interp_hixp_b,data.frame=T)
interp_hixp_b$Date <- as.Date(interp_hixp_b$x)

interp_bixd <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$BIX_DOM,
                      xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_bixd <- interp2xyz(interp_bixd,data.frame=T)
interp_bixd$Date <- as.Date(interp_bixd$x)

interp_bixd_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$BIX_DOM,
                      xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_bixd_b <- interp2xyz(interp_bixd_b,data.frame=T)
interp_bixd_b$Date <- as.Date(interp_bixd_b$x)

interp_bixp <- interp(x=mydata_s$Date,y=mydata_s$Station,z=mydata_s$BIX_POM,
                      xo = seq(min(mydata_s$Date),max(mydata_s$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_bixp <- interp2xyz(interp_bixp,data.frame=T)
interp_bixp$Date <- as.Date(interp_bixp$x)

interp_bixp_b <- interp(x=mydata_b$Date,y=mydata_b$Station,z=mydata_b$BIX_POM,
                      xo = seq(min(mydata_b$Date),max(mydata_b$Date),by=1),
                      yo = seq(0,180,by = 1),
                      extrap = F, linear = T, duplicate = "strip")
interp_bixp_b <- interp2xyz(interp_bixp_b,data.frame=T)
interp_bixp_b$Date <- as.Date(interp_bixp_b$x)

# Plot heatmaps for HIX and BIX; POM and DOM; S and B
hixds <- ggplot()+
  geom_tile(interp_hixd,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,26))+
  labs(x = "",y="Distance down estuary (km)",fill="HIX")+
  theme_classic(base_size=15)

hixdb <- ggplot()+
  geom_tile(interp_hixd_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,26))+
  labs(x = "",y="Distance down estuary (km)",fill="HIX")+
  theme_classic(base_size=15)

hixps <- ggplot()+
  geom_tile(interp_hixp,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,26))+
  labs(x = "",y="Distance down estuary (km)",fill="HIX")+
  theme_classic(base_size=15)

hixpb <- ggplot()+
  geom_tile(interp_hixp_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,26))+
  labs(x = "",y="Distance down estuary (km)",fill="HIX")+
  theme_classic(base_size=15)

bixds <- ggplot()+
  geom_tile(interp_bixd,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,1.3))+
  labs(x = "",y="Distance down estuary (km)",fill="BIX")+
  theme_classic(base_size=15)

bixdb <- ggplot()+
  geom_tile(interp_bixd_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,1.3))+
  labs(x = "",y="Distance down estuary (km)",fill="BIX")+
  theme_classic(base_size=15)

bixps <- ggplot()+
  geom_tile(interp_bixp,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,1.3))+
  labs(x = "",y="Distance down estuary (km)",fill="BIX")+
  theme_classic(base_size=15)

bixpb <- ggplot()+
  geom_tile(interp_bixp_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,1.3))+
  labs(x = "",y="Distance down estuary (km)",fill="BIX")+
  theme_classic(base_size=15)

fl <- ggarrange(hixds,hixdb,hixps,hixpb,bixds,bixdb,bixps,bixpb,nrow=4,ncol=2,labels=c("A","B","C","D","E","F","G","H"))

ggsave("Fig_Output/heatmaps_Fl.png",fl,width = 8, height = 10)