### Script to plot discharge and average discharge from Ft. Barnwell from July 1, 2015
### to July 31, 2016
### A Hounshell, 11 Mar 2020

# Load in libraries
pacman::p_load(dplyr,ggplot2,tidyverse,lubridate,ggpubr)

# Load daily (averaged) discharge data (in cfs)
q <- read_csv("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Daily_Q.csv")
q$datetime <- as.POSIXct(strptime(q$datetime, "%m/%d/%Y", tz="EST"))

# Convert to cms and scale to area of ungaged watershed (Peierls et al.)
q$`85489_00060_00003` <- q$`85489_00060_00003`*0.02832/0.69

# Load in average Q
q_avg <- read_csv("C:/Users/ahoun/OneDrive/Desktop/NRE_Multistats/Data/Avg_Q.csv")
q_avg$DateTime <- as.POSIXct(strptime(q_avg$DateTime, "%m/%d/%Y", tz="EST"))
q_avg$mean_va <- q_avg$mean_va*0.02832/0.69

ggplot()+
  geom_line(q,mapping=aes(x=datetime,y=q$`85489_00060_00003`,linetype="Q"),size=1)+
  geom_line(q_avg,mapping=aes(x=DateTime,y=mean_va,linetype="Yearly average"),size=1)+
  geom_vline(xintercept = as.POSIXct("2015-07-20"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-08-03"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-08-17"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-08-31"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2015-09-14"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-09-29"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2015-10-12"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-10-29"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2015-11-17"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2015-12-08"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-01-20"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-02-17"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-03-07"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-03-22"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-04-06"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-04-19"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-05-09"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-05-24"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-06-06"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-06-20"),color="#D3D3D3")+
  geom_vline(xintercept = as.POSIXct("2016-07-06"),color="#D3D3D3")+ 
  geom_vline(xintercept = as.POSIXct("2016-07-18"),color="#D3D3D3")+
  xlab("Time")+
  ylab(expression(paste("Discharge (m"^3*" s"^-1*")")))+
  theme_classic(base_size=15)+
  theme(legend.position = c(0.15,0.85),legend.title = element_blank(),
        legend.box.background = element_rect(color="black"), 
        legend.box.margin = margin(0.2,0.2,0.2,0.2))
