# Script to plot box plots for 2015-2016 and calculate long term median (2000-2019)
# A Hounshell, 11 Mar 2020

# Load in libraries
pacman::p_load(tidyverse,PerformanceAnalytics,GGally,dplyr,ggpubr,ggplot2,akima,lubridate,colorRamps,RColorBrewer)

# Load in data (Database_DOSat.csv)
my_data <- read.csv(file.choose())
# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]
my_data2$Date <- as.POSIXct(strptime(my_data2$Date, "%m/%d/%Y", tz="EST"))

# Load in historical data: to calculate median (NRWQ_2000to2019)
long_data <- read.csv(file.choose())

# Convert DOC from uM to mg/L
long_data$DOC <- long_data$DOC*12.011/1000

# Convert POC from ug/L to mg/L
long_data$POC <- long_data$POC/1000

############ Heatmaps for spatial/temporal visualizations of key parameters? #######################
## In response to paper revisions: 26 Jan 2021

# Meh - not a great visualization....
# Try heatmaps?
# Need to interpret among data first....
# Separate into S and B
mydata_s <- my_data2 %>% 
  filter(Depth == "S") %>% 
  mutate(Date = as.Date(Date))
mydata_b <- my_data2 %>% 
  filter(Depth == "B") %>% 
  mutate(Date=as.Date(Date))

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

# Plot salinity heatmap
sals <- ggplot()+
  geom_tile(interp_sal,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,20))+
  labs(x = "",y="Distance down estuary (km)",fill="Sal")+
  theme_classic(base_size=10)

salb <- ggplot()+
  geom_tile(interp_sal_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,20))+
  labs(x = "",y="Distance down estuary (km)",fill="Sal")+
  theme_classic(base_size=10)

# Plot chla heatmaps
chlas <- ggplot()+
  geom_tile(interp_chla,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,130))+
  labs(x = "",y="Distance down estuary (km)",fill=expression(paste("Chla (",mu,"g L"^-1*")")))+
  theme_classic(base_size=10)

chlab <- ggplot()+
  geom_tile(interp_chla_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,130))+
  labs(x = "",y="Distance down estuary (km)",fill=expression(paste("Chla (",mu,"g L"^-1*")")))+
  theme_classic(base_size=10)

# What if we did a graph of sal, DOC, POC? Then a separate heatmap of FDOM parameters?
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

# DOC heatmaps
docs <- ggplot()+
  geom_tile(interp_doc,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(4,15))+
  labs(x = "",y="Distance down estuary (km)",fill=expression("DOC (mg L"^-1*")"))+
  theme_classic(base_size=10)

docb <- ggplot()+
  geom_tile(interp_doc_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(4,15))+
  labs(x = "",y="Distance down estuary (km)",fill=expression("DOC (mg L"^-1*")"))+
  theme_classic(base_size=10)

# POC heatmaps
pocs <- ggplot()+
  geom_tile(interp_poc,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_s,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,6))+
  labs(x = "",y="Distance down estuary (km)",fill=expression("POC (mg L"^-1*")"))+
  theme_classic(base_size=10)

pocb <- ggplot()+
  geom_tile(interp_poc_b,mapping=aes(x=Date,y=y,fill=z))+
  geom_point(mydata_b,mapping=aes(x=Date,y=Station),color="white",size=0.7)+
  scale_fill_distiller(palette = "YlGnBu",direction = 1,na.value="gray",limits=c(0,6))+
  labs(x = "",y="Distance down estuary (km)",fill=expression("POC (mg L"^-1*")"))+
  theme_classic(base_size=10)

test <- ggarrange(sals,salb,chlas,chlab,docs,docb,pocs,pocb,nrow=4,ncol=2,labels=c("A","B","C","D","E","F","G","H"))

ggsave("Fig_Output/heatmaps_EnvC.png",test,width = 8, height = 10)

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

###### UPDATED EXPLORATORY GRAPHS FOR REVISIONS ######
# Create graphs that capture temporal and spatial variability of various parameters
# Calculate mean and SD for each station for each season (Summer 2015; Fall 2015; Winter 2015/2016; Spring 2015;
# Summer 2016)
mean_s_summer15 <- my_data2 %>% 
  filter(Date >= "2015-07-20" & Date <= "2015-09-01") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Summer15") %>% 
  mutate(Depth = "S")

stdev_s_summer15 <- my_data2 %>% 
  filter(Date >= "2015-07-20" & Date <= "2015-09-01") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Summer15") %>% 
  mutate(Depth = "S")

mean_b_summer15 <- my_data2 %>% 
  filter(Date >= "2015-07-20" & Date <= "2015-09-01") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Summer15") %>% 
  mutate(Depth = "B")

stdev_b_summer15 <- my_data2 %>% 
  filter(Date >= "2015-07-20" & Date <= "2015-09-01") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Summer15") %>% 
  mutate(Depth = "B")

mean_s_fall15 <- my_data2 %>% 
  filter(Date>="2015-09-14" & Date<"2015-12-08") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Fall15") %>% 
  mutate(Depth = "S")

stdev_s_fall15 <- my_data2 %>% 
  filter(Date>="2015-09-14" & Date<"2015-12-08") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Fall15") %>% 
  mutate(Depth = "S")

mean_b_fall15 <- my_data2 %>% 
  filter(Date>="2015-09-14" & Date<"2015-12-08") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Fall15") %>% 
  mutate(Depth = "B")

stdev_b_fall15 <- my_data2 %>% 
  filter(Date>="2015-09-14" & Date<"2015-12-08") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Fall15") %>% 
  mutate(Depth = "B")

mean_s_winter15 <- my_data2 %>% 
  filter(Date>= "2015-12-01" & Date< "2016-03-07") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Winter15") %>% 
  mutate(Depth = "S")

stdev_s_winter15 <- my_data2 %>% 
  filter(Date>= "2015-12-01" & Date< "2016-03-07") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Winter15") %>% 
  mutate(Depth = "S")

mean_b_winter15 <- my_data2 %>% 
  filter(Date>= "2015-12-01" & Date< "2016-03-07") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Winter15") %>% 
  mutate(Depth = "B")

stdev_b_winter15 <- my_data2 %>% 
  filter(Date>= "2015-12-01" & Date< "2016-03-07") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Winter15") %>% 
  mutate(Depth = "B")

mean_s_spring16 <- my_data2 %>% 
  filter(Date>="2016-03-07" & Date<"2016-06-06") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Spring16") %>% 
  mutate(Depth = "S")

stdev_s_spring16 <- my_data2 %>% 
  filter(Date>="2016-03-07" & Date<"2016-06-06") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Spring16") %>% 
  mutate(Depth = "S")

mean_b_spring16 <- my_data2 %>% 
  filter(Date>="2016-03-07" & Date<"2016-06-06") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Spring16") %>% 
  mutate(Depth = "B")

stdev_b_spring16 <- my_data2 %>% 
  filter(Date>="2016-03-07" & Date<"2016-06-06") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Spring16") %>% 
  mutate(Depth = "B")

mean_s_summer16 <- my_data2 %>% 
  filter(Date>="2016-06-06") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Summer16") %>% 
  mutate(Depth = "S")

stdev_s_summer16 <- my_data2 %>% 
  filter(Date>="2016-06-06") %>% 
  filter(Depth == "S") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Summer16") %>% 
  mutate(Depth = "S")

mean_b_summer16 <- my_data2 %>% 
  filter(Date>="2016-06-06") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(mean)) %>% 
  mutate(Season = "Summer16") %>% 
  mutate(Depth = "B")

stdev_b_summer16 <- my_data2 %>% 
  filter(Date>="2016-06-06") %>% 
  filter(Depth == "B") %>% 
  group_by(Station) %>% 
  summarize_all(funs(sd)) %>% 
  mutate(Season = "Summer16") %>% 
  mutate(Depth = "B")

## Use rbind to concatenate into one matrix for surface and bottom, mean and sd
mean_s <- rbind(mean_s_summer15,mean_s_fall15,mean_s_winter15,mean_s_spring16,mean_s_summer16)

stdev_s <- rbind(stdev_s_summer15,stdev_s_fall15,stdev_s_winter15,stdev_s_spring16,stdev_s_summer16)

mean_b <- rbind(mean_b_summer15,mean_b_fall15,mean_b_winter15,mean_b_spring16,mean_b_summer16)

stdev_b <- rbind(stdev_b_summer15,stdev_b_fall15,stdev_b_winter15,stdev_b_spring16,stdev_b_summer16)

## Then plot!
ggplot()+
  geom_line(mean_s,mapping=aes(x=Station,y=Sal,color=Season))+
  geom_point(mean_s,mapping=aes(x=Station,y=Sal,color=Season))+
  geom_errorbar(stdev_s,mapping=aes(x=Station,ymin=mean_s$Sal-Sal,ymax=mean_s$Sal+Sal,color=Season))

ggplot()+
  geom_line(mean_b,mapping=aes(x=Station,y=Sal,color=Season))+
  geom_point(mean_b,mapping=aes(x=Station,y=Sal,color=Season))+
  geom_errorbar(stdev_b,mapping=aes(x=Station,ymin=mean_b$Sal-Sal,ymax=mean_b$Sal+Sal,color=Season))

# Filter out S and B
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

