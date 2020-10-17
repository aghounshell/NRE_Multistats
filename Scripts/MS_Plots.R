# Script to plot box plots for 2015-2016 and calculate long term median (2000-2019)
# A Hounshell, 11 Mar 2020

# Load in libraries
pacman::p_load(tidyverse,PerformanceAnalytics,GGally,dplyr,ggpubr)

# Load in data (Database_DOSat.csv)
my_data <- read.csv(file.choose())
# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]

# Load in historical data: to calculate median (NRWQ_2000to2019)
long_data <- read.csv(file.choose())

# Convert DOC from uM to mg/L
long_data$DOC <- long_data$DOC*12.011/1000

# Convert POC from ug/L to mg/L
long_data$POC <- long_data$POC/1000

# FIlter out S and B
long_b <- long_data %>% filter(Depth == "B")
long_s <- long_data %>% filter(Depth == "S")

sal_b <- median(long_b$YSI_Salinity,na.rm=TRUE)
sal_s <- median(long_s$YSI_Salinity,na.rm=TRUE)

# Separate by season and calculate median
long_winter <- long_data %>% filter(Season=="Winter")
long_spring <- long_data %>% filter(Season=="Spring")
long_summer <- long_data %>% filter(Season=="Summer")
long_fall <- long_data %>% filter(Season=="Fall")

# Calculate median for each season
# Salinity
sal_winter <- median(long_winter$YSI_Salinity,na.rm=TRUE)
sal_spring <- median(long_spring$YSI_Salinity,na.rm=TRUE)
sal_summer <- median(long_summer$YSI_Salinity,na.rm=TRUE)
sal_fall <- median(long_fall$YSI_Salinity,na.rm=TRUE)

# Chla
chla_winter <- median(long_winter$Correct.Chla_IV,na.rm=TRUE)
chla_spring <- median(long_spring$Correct.Chla_IV,na.rm=TRUE)
chla_summer <- median(long_summer$Correct.Chla_IV,na.rm=TRUE)
chla_fall <- median(long_fall$Correct.Chla_IV,na.rm=TRUE)

# DOC
doc_winter <- median(long_winter$DOC,na.rm=TRUE)
doc_spring <- median(long_spring$DOC,na.rm=TRUE)
doc_summer <- median(long_summer$DOC,na.rm=TRUE)
doc_fall <- median(long_fall$DOC,na.rm=TRUE)

# POC
poc_winter <- median(long_winter$POC,na.rm=TRUE)
poc_spring <- median(long_spring$POC,na.rm=TRUE)
poc_summer <- median(long_summer$POC,na.rm=TRUE)
poc_fall <- median(long_fall$POC,na.rm=TRUE)

## Calculate median for the yearly data (2015-2016) for each season
med <- my_data2 %>% select(Season,Sal,Chla,DOC_mg,POC_mg,,HIX_DOM,HIX_POM,Flushing_Time) %>% group_by(Season) %>% 
  summarize_all(funs(median))

median(my_data2$Flushing_Time)

# Plot salinity and chla
my_data2$Season<-factor(my_data2$Season, levels=c("Summer", "Fall", "Winter", "Spring"))


jpeg("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Plots/Figure3.jpg",width=200,height=110,units="mm",res=800)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,2))

boxplot(Sal~Season,data=my_data2,varwidth=TRUE,ylab="Salinity",cex.lab=0.8,cex.axis=0.8,col="white")
segments(0.7,11.03,1.3,11.03,col="#005b96",lwd=2,lty=5) # Summer
segments(1.7,9.39,2.3,9.39,col="#005b96",lwd=2,lty=5) # Fall
segments(2.7,7.53,3.3,7.53,col="#005b96",lwd=2,lty=5) # Winter
segments(3.7,4.73,4.3,4.73,col="#005b96",lwd=2,lty=5) # Spring
text(3.5,19,labels="- - - 2000-2019\nmedian",col="#005b96",cex=0.8)

# Chla plotted w/o outliers
ylab.text=expression(paste("Chla (",mu,"g L"^"-1"*")"))
boxplot(Chla~Season,data=my_data2,varwidth=TRUE,ylab=ylab.text,cex.lab=0.8,cex.axis=0.8,col="white")
segments(0.7,12.27,1.3,12.27,col="#005b96",lwd=2,lty=5) #Summer
segments(1.7,10.40,2.3,10.40,col="#005b96",lwd=2,lty=5) #Fall
segments(2.7,11.86,3.3,11.86,col="#005b96",lwd=2,lty=5) #Winter
segments(3.7,12.23,4.3,12.23,col="#005b96",lwd=2,lty=5) #Spring
text(3.5,125,labels="- - - 2000-2019\nmedian",col="#005b96",cex=0.8)

dev.off()

# Plot stratication index for SI (1100 x 600)
my_data2$Season<-factor(my_data2$Season, levels=c("Summer", "Fall", "Winter", "Spring"))
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,2))
ylab.text=expression(paste("Stratification index (ppt m"^"-1"*")"))
boxplot(Strat_Index~Season,data=my_data2,varwidth=TRUE,ylab=ylab.text,cex.lab=1.5,cex.axis=1.5,col=FALSE)

ylab.text=expression(paste("Flushing time (d"^"-1"*")"))
boxplot(Flushing_Time~Season,data=my_data2,varwidth=TRUE,ylab=ylab.text,cex.lab=1.5,cex.axis=1.5,col=FALSE)

# Calculate median flushing time



# Plot DOM and POM parameter by season as box plots (1100 x 1000)
jpeg("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Plots/Figure4.jpg",width=250,height=250,units="mm",res=800)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
my_data2$Season<-factor(my_data2$Season, levels=c("Summer", "Fall", "Winter", "Spring"))
ylab.text=expression(paste("DOC (mg L"^"-1"*")"))
boxplot(DOC_mg~Season,data=my_data2,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,ylim=c(0,15),col=FALSE)
segments(0.7,6.52,1.3,6.52,col="#005b96",lwd=2,lty=5) #Summer
segments(1.7,7.78,2.3,7.78,col="#005b96",lwd=2,lty=5) #Fall
segments(2.7,6.83,3.3,6.83,col="#005b96",lwd=2,lty=5) #Winter
segments(3.7,7.00,4.3,7.00,col="#005b96",lwd=2,lty=5) #Spring
text(3.4,1,labels="- - - 2000-2019 median",col="#005b96",cex=1.5)

my_data2$Season<-factor(my_data2$Season, levels=c("Summer", "Fall", "Winter", "Spring"))
ylab.text=expression(paste("POC (mg L"^"-1"*")"))
boxplot(POC_mg~Season,data=my_data2,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,col=FALSE)
segments(0.7,1.43,1.3,1.43,col="#005b96",lwd=2,lty=5) #Summer
segments(1.7,1.22,2.3,1.22,col="#005b96",lwd=2,lty=5) #Fall
segments(2.7,1.33,3.3,1.33,col="#005b96",lwd=2,lty=5) #Winter
segments(3.7,1.46,4.3,1.46,col="#005b96",lwd=2,lty=5) #SPring
text(1.7,4.8,labels="- - - 2000-2019 median",col="#005b96",cex=1.5)

boxplot(HIX_DOM~Season,data=my_data2,varwidth=TRUE,ylab="DOM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25),col=FALSE)
abline(h=6,lty=2)
abline(h=16,lty=2)
text(1.7,4,labels="Less humified material",cex=1.5)
text(3.3,23,labels="More humified material",cex=1.5)

boxplot(HIX_POM~Season,data=my_data2,varwidth=TRUE,ylab="POM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25),col=FALSE)
abline(h=6,lty=2)
abline(h=16,lty=2)

dev.off()

# Crossplot of HIX and BIX by season/location
dom_plot <- ggplot()+
  geom_point(my_data2,mapping=aes(BIX_DOM,HIX_DOM,color=Season),size=2.5)+
  xlim(0,1.3)+
  ylim(0,27)+
  scale_color_manual(breaks=c("Winter","Spring","Summer","Fall"),
                     values=c("#1E88E5","#004D40","#FFC107","#D81B60"))+
  ylab("DOM HIX")+
  xlab("DOM BIX")+
  theme_classic(base_size = 15)+
  theme(legend.position=c(0.8,0.8))

pom_plot <- ggplot()+
  geom_point(my_data2,mapping=aes(BIX_POM,HIX_POM,color=Season),size=2.5)+
  xlim(0,1.3)+
  ylim(0,27)+
  scale_color_manual(breaks=c("Winter","Spring","Summer","Fall"),
                     values=c("#1E88E5","#004D40","#FFC107","#D81B60"))+
  xlab("POM BIX")+
  ylab("POM HIX")+
  theme_classic(base_size = 15)+
  theme(legend.position='none')

ggarrange(dom_plot,pom_plot,common.legend=FALSE,ncol=2,nrow=1)

# Separate data by data pool and by surface and bottom
env_s <- my_data2 %>% filter(Depth == "S") %>% select(Temp,Sal,DO_Sat,Turb,Chla)
env_b <- my_data2 %>% filter(Depth == "B") %>% select(Temp,Sal,DO_Sat,Turb,Chla)
dom_s <- my_data2 %>% filter(Depth == "S") %>% select(DOC_mg,DON_mg,DOC_DON,a254_DOM,SUVA_DOM,HIX_DOM,BIX_DOM,B_DOM,T_DOM,A_DOM,C_DOM,N_DOM,M_DOM)
dom_b <- my_data2 %>% filter(Depth == "B") %>% select(DOC_mg,DON_mg,DOC_DON,a254_DOM,SUVA_DOM,HIX_DOM,BIX_DOM,B_DOM,T_DOM,A_DOM,C_DOM,N_DOM,M_DOM)
pom_s <- my_data2 %>% filter(Depth == "S") %>% select(POC_mg,PN_mg,POCtoPN,a254_POM,SUVA_POC,HIX_POM,BIX_POM,B_POM,T_POM,A_POM,C_POM,N_POM,M_POM)
pom_b <- my_data2 %>% filter(Depth == "B") %>% select(POC_mg,PN_mg,POCtoPN,a254_POM,SUVA_POC,HIX_POM,BIX_POM,B_POM,T_POM,A_POM,C_POM,N_POM,M_POM)

env_all <- my_data2 %>% select(Date,Season,Station,Depth,Temp,Sal,DO_Sat,Turb,Chla)
dom_all <- my_data2 %>% select(Date,Season,Station,Depth,DOC_mg,DON_mg,DOC_DON,a254_DOM,SUVA_DOM,HIX_DOM,BIX_DOM,B_DOM,T_DOM,A_DOM,C_DOM,N_DOM,M_DOM) 
pom_all <- my_data2 %>% select(Date,Season,Station,Depth,POC_mg,PN_mg,POCtoPN,a254_POM,SUVA_POC,HIX_POM,BIX_POM,B_POM,T_POM,A_POM,C_POM,N_POM,M_POM)

# Plot Correlation chart - environmental data (1000, 800)
par(mar=c(5.1,4.1,4.1,2.1))
chart.Correlation(env_s, histogram=TRUE, method=c("pearson"))

chart.Correlation(env_b, histogram=TRUE, method=c("pearson"))

# Plot correlation chart - DOM data
chart.Correlation(dom_s, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_b, histogram=TRUE, method=c("pearson"))

res2 <- cor(dom.all,method=c("pearson"))

# Plot correlation chart - POM data
chart.Correlation(pom_s, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_b, histogram=TRUE, method=c("pearson"))

## Plot OM parameters (DOC, HIX_D, BIX_D, POC, HIX_P, HIX_D) down estuary by station
## Plot as combined S and B (no big difference between the two!)
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(3,2))
ylab.text=expression(paste("DOC (mg L"^"-1"*")"))
boxplot(DOC_mg~Station,data=dom_all,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,ylim=c(0,15),col="white")

ylab.text=expression(paste("POC (mg L"^"-1"*")"))
boxplot(POC_mg~Station,data=pom_all,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,col="white")

boxplot(HIX_DOM~Station,data=dom_all,varwidth=TRUE,ylab="DOM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25),col="white")
abline(h=6,lty=2)
abline(h=16,lty=2)
text(3,4,labels="Fresher material",cex=1.5)
text(4,23,labels="More humified material",cex=1.5)

boxplot(HIX_POM~Station,data=pom_all,varwidth=TRUE,ylab="POM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25),col="white")
abline(h=6,lty=2)
abline(h=16,lty=2)

boxplot(BIX_DOM~Station,data=dom_all,varwidth=TRUE,ylab="DOM BIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1.3),col="white")

boxplot(BIX_POM~Station,data=pom_all,varwidth=TRUE,ylab="POM BIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1.3),col="white")

## Plot surface and bottom salinity by station down estuary
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,1))
boxplot(Sal~Depth*Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5,col=c("blue","green"),ylab="Salinity",xlab="Station")
legend("topleft", c("Surface","Bottom"),inset=0.04, fill=c("blue","green"), cex=1.5)

boxplot(YSI_Salinity~Depth*Station,data=long_data3,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5,col=c("blue","green"),ylab="Salinity",xlab="Station")
legend("topleft", c("Surface","Bottom"),inset=0.04, fill=c("blue","green"), cex=1.5)

boxplot(Sal~Depth,data=my_data2,varwidth=TRUE,ylab="Salinity",cex.lab=1.5,cex.axis=1.5,ylim=c(0,30))
segments(0.6,10.9,1.4,10.9,col="#005b96",lwd=2,lty=5)
segments(1.6,5.12,2.4,5.12,col="#005b96",lwd=2,lty=5)
text(1.2,28,labels="- - - 2000-2019 median",col="#005b96",cex=1.5)

