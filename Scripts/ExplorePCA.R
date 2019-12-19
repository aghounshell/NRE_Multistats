# Script to conduct exploratory PCA for Env, DOM, and POM data collected from the NRE
# Following suggestions from S.R. Fegley
# A Hounshell, 09 Oct 2019
# Updated:
#   Database includes DO_Sat and station as a number
#   Re-do existing PCA with DO_Sat
#   Conduct PCA with Station and Depth included as variables
#   Visualize both with 3-axes (if variability < 80% for 2-axes)
# A Hounshell, 19 Dec 2019

# Load in libraries need
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc,ggplot2,tidyverse,vegan3d,
               scatterplot3d,rgl)

# Load in data (Database_DOSat.csv)
my_data <- read.csv(file.choose())
# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]
my_data2$Date <- as.POSIXct(strptime(my_data2$Date, "%m/%d/%Y", tz="EST"))

# Load in surface temperature from KNKT weather station (Cherry Point)
airtemp <- read_csv('C:/Users/ahoun/Dropbox/NRE_MultiStats/NRE_Multistats/Data/AirTemp_KNKT.csv')
airtemp$Date <- as.POSIXct(strptime(airtemp$Date, "%m/%d/%Y", tz="EST"))
airtemp$Temp_C  <- as.numeric(airtemp$Temp_C) 

# Add new categorical variable to describe upper, mid, lower estuary
attach(my_data2)
my_data2$est[Station < 55] <- "Upper"
my_data2$est[Station > 55 & Station < 130] <- "Mid"
my_data2$est[Station > 130] <- "Lower"
detach(my_data2)
my_data2$est <- as.factor(my_data2$est)

# Separate data by data pool: updated to include DO_Sat instead of DO (as mg/L) for Env_all
# Updated selections to reflect new .csv
env_all <- my_data2[,c(1:7,9:11,46)]
dom_all <- my_data2[,c(1:5,15,17:22,24:29,46)]
pom_all <- my_data2[,c(1:5,31,33:38,40:46)]

# Separate into Surface and bottom
env_s <- env_all %>% filter(Depth=='S')
env_b <- env_all %>% filter(Depth=='B')
dom_s <- dom_all %>% filter(Depth=='S')
dom_b <- dom_all %>% filter(Depth=='B')
pom_s <- pom_all %>% filter(Depth=='S')
pom_b <- pom_all %>% filter(Depth=='B')

# Calculate daily avg temp for Air Temp
daily_temp <- airtemp %>% group_by(Date) %>% summarize(avg_c = mean(Temp_C,na.rm=TRUE))

# Plot Env parameters by time

pdf("C:/Users/ahoun/Dropbox/NRE_MultiStats/NRE_Multistats/Plots/Env_Time.pdf", width=12, height=8)

ggplot(daily_temp,aes(Date,avg_c))+
  geom_point()+
  geom_line()+
  ggtitle('Air Temp')+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(Temp~Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(env_s,aes(Date,Temp,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Temp')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(env_b,aes(Date,Temp,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Temp')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(Sal~Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(env_s,aes(Date,Sal,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Sal')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(env_b,aes(Date,Sal,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Sal')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(DO_Sat~Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(env_s,aes(Date,DO_Sat,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface DO Sat')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(env_b,aes(Date,DO_Sat,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom DO Sat')+
  theme_classic(base_size=15)+
  scale_colour_gradientn(colours=rainbow(5))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(Turb~Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(env_s,aes(Date,Turb,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Turb')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(env_b,aes(Date,Turb,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Turb')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(Chla~Station,data=env_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(env_s,aes(Date,Chla,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Chla')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(env_b,aes(Date,Chla,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Chla')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

dev.off()

## DOM Data

pdf("C:/Users/ahoun/Dropbox/NRE_MultiStats/NRE_Multistats/Plots/DOM_Time.pdf", width=12, height=8)

boxplot(DOC_mg~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,DOC_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface DOC')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,DOC_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom DOC')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(DON_mg~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,DON_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface DON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,DON_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom DON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(DOC_DON~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,DOC_DON,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface DOC:DON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,DOC_DON,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom DOC:DON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(a254_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,a254_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface a254')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,a254_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom a254')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(SUVA_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,SUVA_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface SUVA')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,SUVA_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom SUVA')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(HIX_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,HIX_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface HIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,HIX_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom DOC')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(BIX_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,BIX_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface BIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,BIX_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom BIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(B_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,B_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak B')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,B_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak B')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(T_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,T_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak T')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,T_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak T')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(A_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,A_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak A')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,A_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak A')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(C_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,C_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak C')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,C_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak C')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(M_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,M_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak M')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,M_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak M')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(N_DOM~Station,data=dom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(dom_s,aes(Date,N_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak N')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(dom_b,aes(Date,N_DOM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak N')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

dev.off()

## POM Data

pdf("C:/Users/ahoun/Dropbox/NRE_MultiStats/NRE_Multistats/Plots/POM_Time.pdf", width=12, height=8)

boxplot(POC_mg~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,POC_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface POC')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,POC_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom POC')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(PN_mg~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,PN_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface PON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,PN_mg,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom PON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(POCtoPN~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,POCtoPN,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface POC:PON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,POCtoPN,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom POC:PON')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(a254_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,a254_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface a254')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,a254_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom a254')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(SUVA_POC~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,SUVA_POC,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface SUVA')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,SUVA_POC,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom SUVA')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(HIX_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,HIX_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface HIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,HIX_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom HIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(BIX_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,BIX_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface BIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,BIX_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom BIX')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(B_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,B_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak B')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,B_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak B')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(T_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,T_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak T')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,T_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak T')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(A_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,A_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak A')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,A_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak A')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(C_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,C_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak C')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,C_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak C')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(M_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,M_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak M')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,M_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak M')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

boxplot(N_POM~Station,data=pom_all,varwidth=TRUE,cex.axis=1.5,cex.lab=1.5)

ggplot(pom_s,aes(Date,N_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Surface Peak N')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

ggplot(pom_b,aes(Date,N_POM,group=Station))+
  geom_point(aes(colour=Station))+
  geom_line(aes(colour=Station))+
  ggtitle('Bottom Peak N')+
  scale_colour_gradientn(colours=rainbow(5))+
  theme_classic(base_size=15)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle=90))

dev.off()

# Normalize data using 'scale': do not include station or depth as variables
# Separate data by data pool
env_all <- my_data2[,c(6:7,9:11)]
dom_all <- my_data2[,c(15,17:22,24:29)]
pom_all <- my_data2[,c(31,33:38,40:45)]

# Normalize data: subract mean/SD
env_scale <- scale(env_all)
dom_scale <- scale(dom_all)
pom_scale <- scale(pom_all)

# Plot correlation charts for each data matrix
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/Corr_Plots_DOSat.pdf", width=12, height=8)

chart.Correlation(env_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_scale, histogram=TRUE, method=c("pearson"))

dev.off()

## Conduct PCA on each data martix: no transformations; no removal of outliers
env_pca <- rda(env_scale)
summary(env_pca,axes=0)
plot(env_pca)
text(env_pca)
screeplot(env_pca, bstick = TRUE)

dom_pca <- rda(dom_scale)
summary(dom_pca,axes=0)
plot(dom_pca)
text(dom_pca)

pom_pca <- rda(pom_scale)
summary(pom_pca,axes=0)
plot(pom_pca)
text(pom_pca)

# Make screeplot for each of the PCA output
par(mfrow=c(2,2))
screeplot(env_pca,bstick=TRUE)
screeplot(dom_pca,bstick=TRUE)
screeplot(pom_pca,bstick=TRUE)

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
envspe_sc2 <- scores(env_pca, choices=1:3, display="sp", scaling=2)
domspe_sc2 <- scores(dom_pca, choices=1:3, display="sp", scaling=2)
pomspe_sc2 <- scores(pom_pca, choices=1:3, display="sp", scaling=2)

# Try plotting in 3D using ordiplot3d
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/PCA_Norm_3d.pdf", width=12, height=8)

# PCA results plotted in 3D by season
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

env <- ordiplot3d(env_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (43%)",
                   zlab="PC3 (12%)",ylab="PC2 (28%)")
points(env,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(env$xyz.convert(envspe_sc2),rownames(envspe_sc2),cex=0.8,xpd=TRUE)

dom <- ordiplot3d(dom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (66%)",
                  zlab="PC3 (7%)",ylab="PC2 (17%)")
points(dom,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(dom$xyz.convert(domspe_sc2),rownames(domspe_sc2),cex=0.8,xpd=TRUE)

pom <- ordiplot3d(pom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (45%)",
                  zlab="PC3 (10%)",ylab="PC2 (24%)")
points(pom,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(pom$xyz.convert(pomspe_sc2),rownames(pomspe_sc2),cex=0.8,xpd=TRUE)

# PCA results plotted in 3D by depth
my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")

env <- ordiplot3d(env_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (43%)",
                  zlab="PC3 (12%)",ylab="PC2 (28%)")
points(env,"points",pch=sq[my_data2$Depth],col=c("black","black","black","black"),bg=colvec[my_data2$Depth],
       cex=0.7)
text(env$xyz.convert(envspe_sc2),rownames(envspe_sc2),cex=0.8,xpd=TRUE)

dom <- ordiplot3d(dom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (66%)",
                  zlab="PC3 (7%)",ylab="PC2 (17%)")
points(dom,"points",pch=sq[my_data2$Depth],col=c("black","black","black","black"),bg=colvec[my_data2$Depth],
       cex=0.7)
text(dom$xyz.convert(domspe_sc2),rownames(domspe_sc2),cex=0.8,xpd=TRUE)

pom <- ordiplot3d(pom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (45%)",
                  zlab="PC3 (10%)",ylab="PC2 (24%)")
points(pom,"points",pch=sq[my_data2$Depth],col=c("black","black","black","black"),bg=colvec[my_data2$Depth],
       cex=0.7)
text(pom$xyz.convert(pomspe_sc2),rownames(pomspe_sc2),cex=0.8,xpd=TRUE)

# PCA results plotted in 3D as location in the estuary (Upper, Mid, Lower)
my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")

env <- ordiplot3d(env_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (43%)",
                  zlab="PC3 (12%)",ylab="PC2 (28%)")
points(env,"points",pch=sq[my_data2$est],col=c("black","black","black","black"),bg=colvec[my_data2$est],
       cex=0.7)
text(env$xyz.convert(envspe_sc2),rownames(envspe_sc2),cex=0.8,xpd=TRUE)

dom <- ordiplot3d(dom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (66%)",
                  zlab="PC3 (7%)",ylab="PC2 (17%)")
points(dom,"points",pch=sq[my_data2$est],col=c("black","black","black","black"),bg=colvec[my_data2$est],
       cex=0.7)
text(dom$xyz.convert(domspe_sc2),rownames(domspe_sc2),cex=0.8,xpd=TRUE)

pom <- ordiplot3d(pom_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-2,2),xlab="PC1 (45%)",
                  zlab="PC3 (10%)",ylab="PC2 (24%)")
points(pom,"points",pch=sq[my_data2$est],col=c("black","black","black","black"),bg=colvec[my_data2$est],
       cex=0.7)
text(pom$xyz.convert(pomspe_sc2),rownames(pomspe_sc2),cex=0.8,xpd=TRUE)

# Plotting Env data in 2D by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(env_pca,type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env_pca, display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8,
     col = "black")

plot(env_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,3], angle=20, col="black")
text(env_pca,choices=c(1,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

plot(env_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (28% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,2], envspe_sc2[,3], angle=20, col="black")
text(env_pca,choices=c(2,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

# Plotting Env data in 2D by depth
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")

plot(env_pca,type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env_pca, display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8,
     col = "black")

plot(env_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,3], angle=20, col="black")
text(env_pca,choices=c(1,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

plot(env_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (28% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,2], envspe_sc2[,3], angle=20, col="black")
text(env_pca,choices=c(2,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

# Plotting Env data in 2D by location
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")

plot(env_pca,type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env_pca, display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8,
     col = "black")

plot(env_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (43% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,3], angle=20, col="black")

text(env_pca,choices=c(1,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

plot(env_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (28% var. explained)",ylab="PC3 (12% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,2], envspe_sc2[,3], angle=20, col="black")

text(env_pca,choices=c(2,3), display = "species", labels=c("Temp","Sal","DO_Sat","Turb","Chla"), scaling=2, cex = 0.8, col = "black")

# Plotting DOM data in 2D by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (17% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, choices=c(1,2), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                            "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC3 (7% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(1,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (17% var. explained)",ylab="PC3 (6% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,2], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(2,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")
# Plot DOM Data in 2D by depth
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (17% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, choices=c(1,2), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                            "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC3 (7% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(1,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (17% var. explained)",ylab="PC3 (6% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,2], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(2,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")
# Plot DOM data in 2D by location
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (17% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, choices=c(1,2), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                            "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC3 (7% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(1,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")

plot(dom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (17% var. explained)",ylab="PC3 (6% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sc2[,2], domspe_sc2[,3], angle=20, col="black")
text(dom_pca,choices=c(2,3), display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C",
                                                           "M","N"), scaling=2, cex = 0.8, col = "black")
# Plot POM data in 2D space by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, choices=c(1,2), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B",
                                                            "T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(1,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (24% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,2], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(2,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

# Plot POM data in 2D by Depth
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, choices=c(1,2), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B",
                                                            "T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(1,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (24% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Depth],
                     bg=colvec[Depth],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,2], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(2,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

# Plot POM data in 2D by location
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, choices=c(1,2), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B",
                                                            "T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(1,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (24% var. explained)",ylab="PC3 (10% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[est],
                     bg=colvec[est],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sc2[,2], pomspe_sc2[,3], angle=20, col="black")
text(pom_pca,choices=c(2,3), display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

dev.off()

### Include Depth and Station as variables in PCA space
# Normalize data using 'scale': include station and depth as variables
# Separate data by data pool
env_sta <- my_data2[,c(3,5,6:7,9:11)]
dom_sta <- my_data2[,c(3,5,15,17:22,24:29)]
pom_sta <- my_data2[,c(3,5,31,33:38,40:45)]

# Normalize data: subract mean/SD
env_sta_scale <- scale(env_sta)
dom_sta_scale <- scale(dom_sta)
pom_sta_scale <- scale(pom_sta)

## Conduct PCA on each data martix: no transformations; no removal of outliers
env_sta_pca <- rda(env_sta_scale)
summary(env_sta_pca,axes=0)
plot(env_sta_pca)
text(env_sta_pca)

dom_sta_pca <- rda(dom_sta_scale)
summary(dom_sta_pca,axes=0)
plot(dom_sta_pca)
text(dom_sta_pca)

pom_sta_pca <- rda(pom_sta_scale)
summary(pom_sta_pca,axes=0)
plot(pom_sta_pca)
text(pom_sta_pca)

# Make screeplot for each of the PCA output
par(mfrow=c(2,2))
screeplot(env_pca,bstick=TRUE)
screeplot(dom_pca,bstick=TRUE)
screeplot(pom_pca,bstick=TRUE)

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
envspe_sta_sc2 <- scores(env_sta_pca, choices=1:3, display="sp", scaling=2)
domspe_sta_sc2 <- scores(dom_sta_pca, choices=1:3, display="sp", scaling=2)
pomspe_sta_sc2 <- scores(pom_sta_pca, choices=1:3, display="sp", scaling=2)

# Try plotting in 3D using ordiplot3d
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/PCA_Norm_3d_sta.pdf", width=12, height=8)

# PCA results plotted in 3D by season
par(mfrow=c(1,1))
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

env_2 <- ordiplot3d(env_sta_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-3,3),xlab="PC1 (35%)",
                  zlab="PC3 (16%)",ylab="PC2 (28%)")
points(env_2,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(env_2$xyz.convert(envspe_sta_sc2),rownames(envspe_sta_sc2),cex=0.8,xpd=TRUE)

dom_2 <- ordiplot3d(dom_sta_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-3,3),xlab="PC1 (59%)",
                  zlab="PC3 (8%)",ylab="PC2 (15%)")
points(dom_2,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(dom_2$xyz.convert(domspe_sta_sc2),rownames(domspe_sta_sc2),cex=0.8,xpd=TRUE)

pom_2 <- ordiplot3d(pom_sta_pca,display="sites",choices=1:3,scaling=2,xlim=c(-4,4),zlim=c(-3,3),xlab="PC1 (40%)",
                  zlab="PC3 (9%)",ylab="PC2 (22%)")
points(pom_2,"points",pch=sq[my_data2$Season],col=c("black","black","black","black"),bg=colvec[my_data2$Season],
       cex=0.7)
text(pom_2$xyz.convert(pomspe_sta_sc2),rownames(pomspe_sta_sc2),cex=0.8,xpd=TRUE)

# Plot Env in 2D by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(env_sta_pca,type="n",scaling=2,xlab="PC1 (35% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(env_sta_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sta_sc2[,1], envspe_sta_sc2[,2], angle=20, col="black")
text(env_sta_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")

plot(env_sta_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (35% var. explained)",ylab="PC3 (16% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_sta_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sta_sc2[,1], envspe_sta_sc2[,3], angle=20, col="black")
text(env_sta_pca,choices=c(1,3), display = "species", scaling=2, cex = 0.8, col = "black")

plot(env_sta_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (28% var. explained)",ylab="PC3 (16% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_sta_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sta_sc2[,2], envspe_sta_sc2[,3], angle=20, col="black")
text(env_sta_pca,choices=c(2,3), display = "species", scaling=2, cex = 0.8, col = "black")

# Plot DOM data in 2D by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(dom_sta_pca,type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC2 (15% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(dom_sta_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sta_sc2[,1], domspe_sta_sc2[,2], angle=20, col="black")
text(dom_sta_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")

plot(dom_sta_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC3 (8% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_sta_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sta_sc2[,1], domspe_sta_sc2[,3], angle=20, col="black")
text(dom_sta_pca,choices=c(1,3), display = "species", scaling=2, cex = 0.8, col = "black")

plot(dom_sta_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (15% var. explained)",ylab="PC3 (8% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_sta_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_sta_sc2[,2], domspe_sta_sc2[,3], angle=20, col="black")
text(dom_sta_pca,choices=c(2,3), display = "species", scaling=2, cex = 0.8, col = "black")

# Plot POM data in 2D by season
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(pom_sta_pca,type="n",scaling=2,xlab="PC1 (40% var. explained)",ylab="PC2 (22% var. explained)",cex.axis=1.5,
     cex.lab=1.5)
with(my_data2,points(pom_sta_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sta_sc2[,1], pomspe_sta_sc2[,2], angle=20, col="black")
text(pom_sta_pca, display = "species", scaling=2, cex = 0.8,
     col = "black")

plot(pom_sta_pca,choices=c(1,3),type="n",scaling=2,xlab="PC1 (40% var. explained)",ylab="PC3 (9% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_sta_pca,choices=c(1,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sta_sc2[,1], pomspe_sta_sc2[,3], angle=20, col="black")
text(pom_sta_pca,choices=c(1,3), display = "species", scaling=2, cex = 0.8, col = "black")

plot(pom_sta_pca,choices=c(2,3),type="n",scaling=2,xlab="PC2 (22% var. explained)",ylab="PC3 (9% var. explained)",
     cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_sta_pca,choices=c(2,3),display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                     bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                     pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_sta_sc2[,2], pomspe_sta_sc2[,3], angle=20, col="black")
text(pom_sta_pca,choices=c(2,3), display = "species", scaling=2, cex = 0.8, col = "black")

dev.off()

###### Old code for 2D plots ##########
# Not including depth or location in the PCA

# 2D plots
# Season ordered as: winter, spring, summer, fall
pdf("C:/Users/ahoun/Dropbox/NRE_MultiStats/NRE_Multistats/Plots/PCA_Norm.pdf", width=12, height=8)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env_pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (7% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

# PCA results plotted as S vs. B
my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env.pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (7% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

# PCA results plotted as location in the estuary (Upper, Mid, Lower)
my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env_pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
arrows(0, 0, envspe_sc2[,1], envspe_sc2[,2], angle=20, col="black")
text(env_pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

plot(dom_pca,type="n",scaling=2,xlab="PC1 (66% var. explained)",ylab="PC2 (7% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom_pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
arrows(0, 0, domspe_sc2[,1], domspe_sc2[,2], angle=20, col="black")
text(dom_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom_pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

dev.off()

## Explore the impact of outliers on POM data
# Make boxplots of all POM data first
par(mfrow=c(1,1))
boxplot(pom_all,use.cols=TRUE,varwidth=TRUE)

# Take sqrt of POM data then normalize
pom_sqrt <- sqrt(pom_all)
boxplot(pom_sqrt,use.cols=TRUE,varwidth=TRUE)

pom_sq_scale <- scale(pom_sqrt)

# Now conduct PCA on transformed AND normalized data
pom_pca_t <- rda(pom_sq_scale)
summary(pom_pca_t,axes=0)
plot(pom_pca_t)
text(pom_pca_t)

# Plot screeplot for POM data: transformed and un-transformed
par(mfrow=c(1,2))
screeplot(pom_pca,bstick=TRUE)
screeplot(pom_pca_t,bstick=TRUE)

pomspe_sc2_t <- scores(pom_pca_t, display="sp", choices=c(1,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

par(mfrow=c(1,2))
plot(pom_pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,cex.lab=1.5,main='Un-Transformed')
with(my_data2,points(pom_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, pomspe_sc2[,1], pomspe_sc2[,2], angle=20, col="black")
text(pom_pca, display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

plot(pom_pca_t,type="n",scaling=2,xlab="PC1 (47% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5,main='Sqrt')
with(my_data2,points(pom_pca_t,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, pomspe_sc2_t[,1], pomspe_sc2_t[,2], angle=20, col="black")
text(pom_pca_t, display = "species", labels=c("POC","PN","POC:PN","a254","SUVA","HIX","BIX","B","T","A","C","M","N"), scaling=2, cex = 0.8, col = "black")

# Saved rfile as: PCA_3d.RDATA