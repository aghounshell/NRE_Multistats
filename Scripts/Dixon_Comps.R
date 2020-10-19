## Script to plot monthly averaged Chla from Dixon et al., 2014 as compared to the 2015-2016 period
## Include: Mean, Max, Min for each month
## A Hounshell, 12 Mar 2020

# Load in libraries
pacman::p_load(tidyverse,ggpubr)

# Load in data from 2015-2016
data <- read.csv("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Database_DOSat.csv")

# Remove un-complete data rows (any rows that do not have all data associated with them)
data2 <- data[complete.cases(data),]
data2$Date <- as.POSIXct(strptime(data2$Date, "%m/%d/%Y", tz="EST"))
data2$Date <- format(as.POSIXct(data2$Date), "%m")

# Load in data from Dixon et al., 2014
dixon <- read.csv("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Dixon_Data.csv")
dixon$X.1 <- as.POSIXct(strptime(dixon$X.1, "%m/%d/%Y", tz="EST"))
dixon$X.1 <- format(as.POSIXct(dixon$X.1), "%m")

# Calculate mean, min, and max for each month for 2015-2016
mean <- data2 %>% select(Date,Sal,Chla,DOC_mg,SUVA_DOM,HIX_DOM) %>% group_by(Date) %>% summarise_all(funs(mean(.,na.rm=TRUE)))
max <- data2 %>% select(Date,Sal,Chla,DOC_mg,SUVA_DOM,HIX_DOM) %>% group_by(Date) %>% summarise_all(funs(max(.,na.rm=TRUE)))
min <- data2 %>% select(Date,Sal,Chla,DOC_mg,SUVA_DOM,HIX_DOM) %>% group_by(Date) %>% summarise_all(funs(min(.,na.rm=TRUE)))

# Format Dixon et al., Data
dixon_mean <- dixon %>% filter(X.2 == "Mean")
dixon_min <- dixon %>% filter(X.2 == "Min")
dixon_max <- dixon %>% filter(X.2 == "Max")

# Re-plot: Sal, DOC, SUVA, and HIX
sal <- ggplot()+
  geom_errorbar(mapping=aes(x=mean$Date,y=mean$Sal,ymin=min$Sal,ymax=max$Sal,color="2015-2016"),size=1)+
  geom_point(mapping=aes(x=mean$Date,y=mean$Sal,color="2015-2016"),size=3)+
  geom_line(mean,mapping=aes(x=Date,y=Sal,group=1,color="2015-2016"),size=1)+
  geom_errorbar(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$Salinity,ymin=dixon_min$Salinity,ymax=dixon_max$Salinity,color="2010-2011"),size=1)+
  geom_point(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$Salinity,color="2010-2011"),size=3)+
  geom_line(dixon_mean,mapping=aes(x=X.1,y=Salinity,group=1,color="2010-2011"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2010-2011"), values=c('#173f5f','#3caea3'))+
  ylab("Salinity")+
  xlab("Month")+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

doc <- ggplot()+
  geom_errorbar(mapping=aes(x=mean$Date,y=mean$DOC_mg,ymin=min$DOC_mg,ymax=max$DOC_mg,color="2015-2016"),size=1)+
  geom_point(mapping=aes(x=mean$Date,y=mean$DOC_mg,color="2015-2016"),size=3)+
  geom_line(mean,mapping=aes(x=Date,y=DOC_mg,group=1,color="2015-2016"),size=1)+
  geom_errorbar(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$DOC..mg.L.,ymin=dixon_min$DOC..mg.L.,ymax=dixon_max$DOC..mg.L.,color="2010-2011"),size=1)+
  geom_point(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$DOC..mg.L.,color="2010-2011"),size=3)+
  geom_line(dixon_mean,mapping=aes(x=X.1,y=DOC..mg.L.,group=1,color="2010-2011"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2010-2011"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  xlab("Month")+
  ylim(c(0,15))+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

suva <- ggplot()+
  geom_errorbar(mapping=aes(x=mean$Date,y=mean$SUVA_DOM,ymin=min$SUVA_DOM,ymax=max$SUVA_DOM,color="2015-2016"),size=1)+
  geom_point(mapping=aes(x=mean$Date,y=mean$SUVA_DOM,color="2015-2016"),size=3)+
  geom_line(mean,mapping=aes(x=Date,y=SUVA_DOM,group=1,color="2015-2016"),size=1)+
  geom_errorbar(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$SUVA,ymin=dixon_min$SUVA,ymax=dixon_max$SUVA,color="2010-2011"),size=1)+
  geom_point(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$SUVA,color="2010-2011"),size=3)+
  geom_line(dixon_mean,mapping=aes(x=X.1,y=SUVA,group=1,color="2010-2011"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2010-2011"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("SUVA"[254]*" (L mg"^-1*" C m"^-1*")")))+
  xlab("Month")+
  ylim(c(0,8))+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

hix <- ggplot()+
  geom_errorbar(mapping=aes(x=mean$Date,y=mean$HIX_DOM,ymin=min$HIX_DOM,ymax=max$HIX_DOM,color="2015-2016"),size=1)+
  geom_point(mapping=aes(x=mean$Date,y=mean$HIX_DOM,color="2015-2016"),size=3)+
  geom_line(mean,mapping=aes(x=Date,y=HIX_DOM,group=1,color="2015-2016"),size=1)+
  geom_errorbar(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$HIX,ymin=dixon_min$HIX,ymax=dixon_max$HIX,color="2010-2011"),size=1)+
  geom_point(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$HIX,color="2010-2011"),size=3)+
  geom_line(dixon_mean,mapping=aes(x=X.1,y=HIX,group=1,color="2010-2011"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2010-2011"), values=c('#173f5f','#3caea3'))+
  ylab("HIX")+
  xlab("Month")+
  ylim(c(0,25))+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())

all <- ggarrange(sal,doc,suva,hix,common.legend=TRUE,ncol=2,nrow=2)

ggsave("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Plots/Figure7.jpg",all,dpi=800,width=250,height=250,
       units=c("mm"))

# Plot Chla comp between the two study periods
ggplot()+
  geom_errorbar(mapping=aes(x=mean$Date,y=mean$Chla,ymin=min$Chla,ymax=max$Chla,color="2015-2016"),size=1)+
  geom_point(mapping=aes(x=mean$Date,y=mean$Chla,color="2015-2016"),size=3)+
  geom_line(mean,mapping=aes(x=Date,y=Chla,group=1,color="2015-2016"),size=1)+
  geom_errorbar(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$Chla..ug.L.,ymin=dixon_min$Chla..ug.L.,ymax=dixon_max$Chla..ug.L.,color="2010-2011"),size=1)+
  geom_point(mapping=aes(x=dixon_mean$X.1,y=dixon_mean$Chla..ug.L.,color="2010-2011"),size=3)+
  geom_line(dixon_mean,mapping=aes(x=X.1,y=Chla..ug.L.,group=1,color="2010-2011"),size=1)+
  scale_color_manual(breaks=c("2015-2016","2010-2011"), values=c('#173f5f','#3caea3'))+
  ylab(expression(paste("Chla (", mu,"g L"^-1*")")))+
  xlab("Month")+
  theme_classic(base_size = 15)+
  theme(legend.title = element_blank())
