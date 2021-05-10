# Script to conduct exploratory PCA for Env, DOM, and POM data collected from the NRE
# Following suggestions from S.R. Fegley
# Updated:
#   Database includes DO_Sat and station as a number
#   Re-do existing PCA with DO_Sat
#   Conduct PCA with Station and Depth included as variables
#   Visualize both with 3-axes (if variability < 80% for 2-axes)
# Following additional suggestions from S.R. Fegley (01/08/202):
#   Separate surface and bottom
#   Include DO Sat (instead of %DO)
#   ONLY include two-axes (third axis did not add valuable information)
#   Include all variables in PCA (but then remove correlated variables)
#   Conduct RDA using selected variables from PCA
# A Hounshell, 01 Jan. 2020

# Updated for MS revisions:
#       Remove fluorescence peak-picking indicators
#       Re-do analyses - did anything substantially change?
# A Hounshell, 02 Feb 2021
# Rfile saved as: MultiStats

# Load in libraries needed
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc,ggplot2,tidyverse,vegan3d,
               scatterplot3d,rgl,car,patchwork)

# Load in data (Database_DOSat.csv)
my_data <- read.csv('C:/Users/ahoun/Desktop/NRE_Multistats/Data/Database_DOSat.csv')
# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]
my_data2$Date <- as.POSIXct(strptime(my_data2$Date, "%m/%d/%Y"))

# Add new categorical variable to describe upper, mid, lower estuary
attach(my_data2)
my_data2$est[Station < 55] <- "Upper"
my_data2$est[Station > 55 & Station < 130] <- "Mid"
my_data2$est[Station > 130] <- "Lower"
detach(my_data2)
my_data2$est <- as.factor(my_data2$est)

# Updated season to reflect summer 2015 and summer 2016
my_data2 <- my_data2 %>% 
        mutate(Season = ifelse(Date < '2015-09-01' & Season == "Summer","Summer15",
                                ifelse(Date > '2016-06-01' & Season == "Summer","Summer16",
                                       Season)))

# Separate data by data pool: updated to include DO_Sat instead of DO (as mg/L) for Env_all
# Updated selections to reflect new .csv
env_all <- my_data2 %>% 
        select("Date","Season","Station","Depth","Depth_num","Temp","Sal","DO_Sat","Turb",
               "Chla","est")
dom_all <- my_data2 %>% 
        select("Date","Season","Station","Depth","Depth_num","DOC_mg","DON_mg",
               "DOC_DON","a254_DOM","SUVA_DOM","HIX_DOM","BIX_DOM","est")
pom_all <- my_data2 %>% 
        select("Date","Season","Station","Depth","Depth_num","POC_mg","PN_mg",
               "POCtoPN","a254_POM","SUVA_POC","HIX_POM","BIX_POM","est")

# Separate into Surface and bottom
env_s <- env_all %>% filter(Depth=='S')
env_b <- env_all %>% filter(Depth=='B')
dom_s <- dom_all %>% filter(Depth=='S')
dom_b <- dom_all %>% filter(Depth=='B')
pom_s <- pom_all %>% filter(Depth=='S')
pom_b <- pom_all %>% filter(Depth=='B')

# Normalize data using 'scale': do not include station or depth as variables
# Separate data by data pool and by surface and bottom
env_s_pca <- env_s[,c(6:10)]
env_b_pca <- env_b[,c(6:10)]
dom_s_pca <- dom_s[,c(6:12)]
dom_b_pca <- dom_b[,c(6:12)]
pom_s_pca <- pom_s[,c(6:12)]
pom_b_pca <- pom_b[,c(6:12)]

# Normalize data: subtract mean/SD
env_s_scale <- scale(env_s_pca)
env_b_scale <- scale(env_b_pca)
dom_s_scale <- scale(dom_s_pca)
dom_b_scale <- scale(dom_b_pca)
pom_s_scale <- scale(pom_s_pca)
pom_b_scale <- scale(pom_b_pca)

# Plot correlation charts for each data matrix
pdf("C:/Users/ahoun/Desktop/NRE_MultiStats/Fig_Output/Corr_Plots_SandB_Revs.pdf", width=12, height=8)

chart.Correlation(env_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(env_b_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_b_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_b_scale, histogram=TRUE, method=c("pearson"))

dev.off()

## Construct correlation plots for publication - use unscaled data for better interpretation
# Calculate correlations for surface environmental data
cor(env_s$Temp,env_s$Sal,method="pearson")
cor(env_s$Temp,env_s$DO_Sat,method="pearson")
cor(env_s$DO_Sat,env_s$Sal,method="pearson")
cor(env_s$Temp,env_s$Turb,method="pearson")
cor(env_s$Sal,env_s$Turb,method="pearson")
cor(env_s$DO_Sat,env_s$Turb,method="pearson")
cor(env_s$Temp,env_s$Chla,method="pearson")
cor(env_s$Sal,env_s$Chla,method="pearson")
cor(env_s$DO_Sat,env_s$Chla,method="pearson")
cor(env_s$Turb,env_s$Chla,method="pearson")

# Surface Environmental
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS4.jpg",width=420,height=400,units="mm",res=800)

p1 <- ggplot(data=env_s,mapping=aes(x=Temp))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab(expression('Temp ('*degree*C*')'))+
        theme_classic(base_size=21)

p2 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=Sal,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=Sal),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=20,label="r = 0.39",size=6)+
        ylab("Sal")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())


p3 <- ggplot(data=env_s,mapping=aes(x=Sal))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        theme_classic(base_size=21)

p4 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=DO_Sat,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=DO_Sat),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=160,label="r = 0.13",size=6)+
        ylab("DO Sat (%)")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))


p5 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Sal,y=DO_Sat,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Sal,y=DO_Sat),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=15,y=150,label="r = 0.34",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))


p6 <- ggplot(data=env_s,mapping=aes(x=DO_Sat))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab("DO Sat (%)")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p7 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=Turb),method=lm,se=FALSE,color="black")+
        ylab("Turb (NTU)")+
        annotate(geom="text",x=25,y=28,label="r = -0.56",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())


p8 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Sal,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Sal,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=30,label="r = -0.57",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())


p9 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=DO_Sat,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=DO_Sat,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=55,y=30,label="r = -0.27",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
                    panel.background=element_rect(fill="grey93"))


p10 <- ggplot(data=env_s,mapping=aes(x=Turb))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab("Turb (NTU)")+
        theme_classic(base_size=21)

p11 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=130,label="r = 0.12",size=6)+
        ylab(expression(paste("Chla (",mu,"g L"^"-1"*")")))+
        xlab(expression('Temp ('*degree*C*')'))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")


p12 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Sal,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Sal,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=130,label="r = 0.29",size=6)+
        xlab("Sal")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())


p13 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=DO_Sat,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=DO_Sat,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=55,y=130,label="r = 0.50",size=6)+
        xlab("DO Sat (%)")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))


p14 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Turb,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Turb,y=Chla),method=lm,se=FALSE,color="black")+
        xlab("Turb (NTU")+
        annotate(geom="text",x=20,y=130,label="r = -0.28",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())


p15 <- ggplot(data=env_s,mapping=aes(x=Chla))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab(expression(paste("Chla (",mu,"g L"^"-1"*")")))+
        theme_classic(base_size=21)

layout <- '
A####
BC###
DEG##
HIJK#
LMNOP'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, design=layout)

dev.off()

# Bottom Environmental
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS5.jpg",width=420,height=400,units="mm",res=800)

p1 <- ggplot(data=env_b,mapping=aes(x=Temp))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab(expression('Temp ('*degree*C*')'))+
        theme_classic(base_size=21)

cor(env_b$Temp,env_b$Sal,method="pearson")

p2 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Temp,y=Sal,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Temp,y=Sal),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=20,label="r = 0.49",size=6)+
        ylab("Sal")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())


p3 <- ggplot(data=env_b,mapping=aes(x=Sal))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        theme_classic(base_size=21)

cor(env_b$Temp,env_b$DO_Sat,method="pearson")

p4 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Temp,y=DO_Sat,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Temp,y=DO_Sat),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=110,label="r = -0.65",size=6)+
        ylab("DO Sat (%)")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(env_b$DO_Sat,env_b$Sal,method="pearson")

p5 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Sal,y=DO_Sat,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Sal,y=DO_Sat),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=15,y=110,label="r = -0.42",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p6 <- ggplot(data=env_b,mapping=aes(x=DO_Sat))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab("DO Sat (%)")+
        theme_classic(base_size=21)

cor(env_b$Temp,env_b$Turb,method="pearson")

p7 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Temp,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Temp,y=Turb),method=lm,se=FALSE,color="black")+
        ylab("Turb (NTU)")+
        annotate(geom="text",x=25,y=28,label="r = -0.46",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(env_b$Sal,env_b$Turb,method="pearson")

p8 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Sal,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Sal,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=30,label="r = -0.66",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(env_b$DO_Sat,env_b$Turb,method="pearson")

p9 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=DO_Sat,y=Turb,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=DO_Sat,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=55,y=30,label="r = 0.36",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p10 <- ggplot(data=env_b,mapping=aes(x=Turb))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab("Turb (NTU)")+
        theme_classic(base_size=21)

cor(env_b$Temp,env_b$Chla,method="pearson")

p11 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Temp,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Temp,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=100,label="r = -0.04",size=6)+
        ylab(expression(paste("Chla (",mu,"g L"^"-1"*")")))+
        xlab(expression('Temp ('*degree*C*')'))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")

cor(env_b$Sal,env_b$Chla,method="pearson")

p12 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Sal,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Sal,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=100,label="r = 0.18",size=6)+
        xlab("Sal")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(env_b$DO_Sat,env_b$Chla,method="pearson")

p13 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=DO_Sat,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=DO_Sat,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=55,y=100,label="r = 0.13",size=6)+
        xlab("DO Sat (%)")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(env_b$Turb,env_b$Chla,method="pearson")

p14 <- ggplot()+
        geom_point(data=env_b,mapping=aes(x=Turb,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_b,mapping=aes(x=Turb,y=Chla),method=lm,se=FALSE,color="black")+
        xlab("Turb (NTU)")+
        annotate(geom="text",x=20,y=100,label="r = -0.19",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())


p15 <- ggplot(data=env_b,mapping=aes(x=Chla))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        xlab(expression(paste("Chla (",mu,"g L"^"-1"*")")))+
        theme_classic(base_size=21)

layout <- '
A####
BC###
DEG##
HIJK#
LMNOP'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, design=layout)

dev.off()

## Surface DOM
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS6.jpg",width=600,height=400,units="mm",res=800)

p1 <- ggplot(data=dom_s,mapping=aes(x=DOC_mg))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab(expression(paste("DOC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)

cor(dom_s$DOC_mg,dom_s$DON_mg,method="pearson")

p2 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=DON_mg,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=DON_mg),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7,y=0.52,label="r = 0.81",size=6)+
        ylab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p3 <- ggplot(data=dom_s,mapping=aes(x=DON_mg))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_mg,dom_s$DOC_DON,method="pearson")

p4 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=DOC_DON,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=DOC_DON),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=13,y=45,label="r = 0.51",size=6)+
        ylab("DOC:DON")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(dom_s$DON_mg,dom_s$DOC_DON,method="pearson")

p5 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DON_mg,y=DOC_DON,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DON_mg,y=DOC_DON),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=45,label="r = -0.08",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p6 <- ggplot(data=dom_s,mapping=aes(x=DOC_DON))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("DOC:DON")+
        theme_classic(base_size=21)

cor(dom_s$DOC_mg,dom_s$a254_DOM,method="pearson")

p7 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7,y=125,label="r = 0.96",size=6)+
        ylab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DON_mg,dom_s$a254_DOM,method="pearson")

p8 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DON_mg,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DON_mg,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.3,y=125,label="r = 0.79",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_DON,dom_s$a254_DOM,method="pearson")

p9 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_DON,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_DON,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=40,y=125,label="r = 0.47",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p10 <- ggplot(data=dom_s,mapping=aes(x=a254_DOM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_mg,dom_s$SUVA_DOM,method="pearson")

p11 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=12.5,y=2.5,label="r = 0.63",size=6)+
        ylab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(dom_s$DON_mg,dom_s$SUVA_DOM,method="pearson")

p12 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DON_mg,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DON_mg,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=2.5,label="r = 0.52",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_DON,dom_s$SUVA_DOM,method="pearson")

p13 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_DON,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_DON,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=40,y=2.5,label="r = 0.47",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(dom_s$a254_DOM,dom_s$SUVA_DOM,method="pearson")

p14 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=a254_DOM,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=a254_DOM,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=110,y=2.5,label="r = 0.82",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p15 <- ggplot(data=dom_s,mapping=aes(x=SUVA_DOM))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab("SUVA")+
        theme_classic(base_size=21)

cor(dom_s$DOC_mg,dom_s$HIX_DOM,method="pearson")

p16 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7,y=21,label="r = 0.80",size=6)+
        ylab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DON_mg,dom_s$HIX_DOM,method="pearson")

p17 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DON_mg,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DON_mg,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.3,y=21,label="r = 0.64",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_DON,dom_s$HIX_DOM,method="pearson")

p18 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_DON,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_DON,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=40,y=21,label="r = 0.43",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$a254_DOM,dom_s$HIX_DOM,method="pearson")

p19 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=a254_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=a254_DOM,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=50,y=21,label="r = 0.88",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$SUVA_DOM,dom_s$HIX_DOM,method="pearson")

p20 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=SUVA_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=SUVA_DOM,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=2.6,y=21,label="r = 0.78",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p21 <- ggplot(data=dom_s,mapping=aes(x=HIX_DOM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_mg,dom_s$BIX_DOM,method="pearson")

p22 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_mg,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_mg,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=13,y=0.8,label="r = -0.54",size=6)+
        ylab("BIX")+
        xlab(expression(paste("DOC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")

cor(dom_s$DON_mg,dom_s$BIX_DOM,method="pearson")

p23 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DON_mg,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DON_mg,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=0.8,label="r = -0.44",size=6)+
        xlab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$DOC_DON,dom_s$BIX_DOM,method="pearson")

p24 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=DOC_DON,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=DOC_DON,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=40,y=0.8,label="r = -0.32",size=6)+
        xlab("DOC:DON")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(dom_s$a254_DOM,dom_s$BIX_DOM,method="pearson")

p25 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=a254_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=a254_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=110,y=0.8,label="r = -0.66",size=6)+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_s$SUVA_DOM,dom_s$BIX_DOM,method="pearson")

p26 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=SUVA_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=SUVA_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3.7,y=0.8,label="r = -0.77",size=6)+
        xlab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(dom_s$HIX_DOM,dom_s$BIX_DOM,method="pearson")

p27 <- ggplot()+
        geom_point(data=dom_s,mapping=aes(x=HIX_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_s,mapping=aes(x=HIX_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=18,y=0.8,label="r = -0.77",size=6)+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p28 <- ggplot(data=dom_s,mapping=aes(x=BIX_DOM))+
        geom_histogram(binwidth = 0.05)+
        ylab("Count")+
        xlab("BIX")+
        theme_classic(base_size=21)

layout <- '
A######
BC#####
DEG####
HIJK###
LMNOP##
QRSUVW#
XYZabcd'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, 
           Q=p16, R=p17, S=p18, U=p19, V=p20, W=p21, X=p22, Y=p23, Z=p24, a=p25, b=p26, c=p27, d=p28, design=layout)

dev.off()

## Bottom DOM
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS7.jpg",width=600,height=400,units="mm",res=800)

p1 <- ggplot(data=dom_b,mapping=aes(x=DOC_mg))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab(expression(paste("DOC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)

cor(dom_b$DOC_mg,dom_b$DON_mg,method="pearson")

p2 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=DON_mg,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=DON_mg),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=6,y=0.52,label="r = 0.77",size=6)+
        ylab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

p3 <- ggplot(data=dom_b,mapping=aes(x=DON_mg))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)

cor(dom_b$DOC_mg,dom_b$DOC_DON,method="pearson")

p4 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=DOC_DON,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=DOC_DON),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=6,y=38,label="r = 0.58",size=6)+
        ylab("DOC:DON")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(dom_b$DON_mg,dom_b$DOC_DON,method="pearson")

p5 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DON_mg,y=DOC_DON,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DON_mg,y=DOC_DON),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=35,label="r = -0.07",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p6 <- ggplot(data=dom_b,mapping=aes(x=DOC_DON))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("DOC:DON")+
        theme_classic(base_size=21)

cor(dom_b$DOC_mg,dom_b$a254_DOM,method="pearson")

p7 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=6,y=100,label="r = 0.94",size=6)+
        ylab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DON_mg,dom_b$a254_DOM,method="pearson")

p8 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DON_mg,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DON_mg,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.3,y=105,label="r = 0.74",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_DON,dom_b$a254_DOM,method="pearson")

p9 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_DON,y=a254_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_DON,y=a254_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=35,y=35,label="r = 0.54",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p10 <- ggplot(data=dom_b,mapping=aes(x=a254_DOM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_mg,dom_b$SUVA_DOM,method="pearson")

p11 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=11,y=2.5,label="r = 0.63",size=6)+
        ylab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DON_mg,dom_b$SUVA_DOM,method="pearson")

p12 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DON_mg,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DON_mg,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=2.5,label="r = 0.49",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_DON,dom_b$SUVA_DOM,method="pearson")

p13 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_DON,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_DON,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=35,y=2.5,label="r = 0.39",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$a254_DOM,dom_b$SUVA_DOM,method="pearson")

p14 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=a254_DOM,y=SUVA_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=a254_DOM,y=SUVA_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=80,y=2.5,label="r = 0.85",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p15 <- ggplot(data=dom_b,mapping=aes(x=SUVA_DOM))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab("SUVA")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_mg,dom_b$HIX_DOM,method="pearson")

p16 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=6,y=18,label="r = 0.79",size=6)+
        ylab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DON_mg,dom_b$HIX_DOM,method="pearson")

p17 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DON_mg,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DON_mg,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=7,label="r = 0.58",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_DON,dom_b$HIX_DOM,method="pearson")

p18 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_DON,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_DON,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=35,y=7,label="r = 0.50",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$a254_DOM,dom_b$HIX_DOM,method="pearson")

p19 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=a254_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=a254_DOM,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=45,y=18,label="r = 0.86",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$SUVA_DOM,dom_b$HIX_DOM,method="pearson")

p20 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=SUVA_DOM,y=HIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=SUVA_DOM,y=HIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=2.6,y=21,label="r = 0.77",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p21 <- ggplot(data=dom_b,mapping=aes(x=HIX_DOM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$DOC_mg,dom_b$BIX_DOM,method="pearson")

p22 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_mg,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_mg,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=0.75,label="r = -0.68",size=6)+
        ylab("BIX")+
        xlab(expression(paste("DOC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")

cor(dom_b$DON_mg,dom_b$BIX_DOM,method="pearson")

p23 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DON_mg,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DON_mg,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.5,y=0.75,label="r = -0.53",size=6)+
        xlab(expression(paste("DON (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(dom_b$DOC_DON,dom_b$BIX_DOM,method="pearson")

p24 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=DOC_DON,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=DOC_DON,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=35,y=0.75,label="r = -0.41",size=6)+
        xlab("DOC:DON")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(dom_b$a254_DOM,dom_b$BIX_DOM,method="pearson")

p25 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=a254_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=a254_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=80,y=0.75,label="r = -0.75",size=6)+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$SUVA_DOM,dom_b$BIX_DOM,method="pearson")

p26 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=SUVA_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=SUVA_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4.5,y=0.75,label="r = -0.70",size=6)+
        xlab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(dom_b$HIX_DOM,dom_b$BIX_DOM,method="pearson")

p27 <- ggplot()+
        geom_point(data=dom_b,mapping=aes(x=HIX_DOM,y=BIX_DOM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=dom_b,mapping=aes(x=HIX_DOM,y=BIX_DOM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=16,y=0.75,label="r = -0.75",size=6)+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p28 <- ggplot(data=dom_b,mapping=aes(x=BIX_DOM))+
        geom_histogram(binwidth = 0.05)+
        ylab("Count")+
        xlab("BIX")+
        theme_classic(base_size=21)

layout <- '
A######
BC#####
DEG####
HIJK###
LMNOP##
QRSUVW#
XYZabcd'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, 
           Q=p16, R=p17, S=p18, U=p19, V=p20, W=p21, X=p22, Y=p23, Z=p24, a=p25, b=p26, c=p27, d=p28, design=layout)

dev.off()

## Surface POM
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS8.jpg",width=600,height=400,units="mm",res=800)

p1 <- ggplot(data=pom_s,mapping=aes(x=POC_mg))+
        geom_histogram(binwidth = 0.5)+
        ylab("Count")+
        xlab(expression(paste("POC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)

cor(pom_s$POC_mg,pom_s$PN_mg,method="pearson")

p2 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=PN_mg,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=PN_mg),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4.5,y=0.1,label="r = 0.91",size=6)+
        ylab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p3 <- ggplot(data=pom_s,mapping=aes(x=PN_mg))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_s$POC_mg,pom_s$POCtoPN,method="pearson")

p4 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=POCtoPN,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=POCtoPN),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=4.5,label="r = 0.02",size=6)+
        ylab("POC:PN")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_s$PN_mg,pom_s$POCtoPN,method="pearson")

p5 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=PN_mg,y=POCtoPN,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=PN_mg,y=POCtoPN),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.6,y=13,label="r = -0.31",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p6 <- ggplot(data=pom_s,mapping=aes(x=POCtoPN))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("POC:PN")+
        theme_classic(base_size=21)

cor(pom_s$POC_mg,pom_s$a254_POM,method="pearson")

p7 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=1,label="r = 0.42",size=6)+
        ylab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_s$PN_mg,pom_s$a254_POM,method="pearson")

p8 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=PN_mg,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=PN_mg,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.2,y=8.5,label="r = 0.31",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_s$POCtoPN,pom_s$a254_POM,method="pearson")

p9 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POCtoPN,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POCtoPN,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=9,label="r = 0.24",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p10 <- ggplot(data=pom_s,mapping=aes(x=a254_POM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("a254")+
        theme_classic(base_size=21)

cor(pom_s$POC_mg,pom_s$SUVA_POC,method="pearson")

p11 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=4.5,label="r = -0.43",size=6)+
        ylab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_s$PN_mg,pom_s$SUVA_POC,method="pearson")

p12 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=PN_mg,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=PN_mg,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.6,y=4.5,label="r = -0.46",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_s$POCtoPN,pom_s$SUVA_POC,method="pearson")

p13 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POCtoPN,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POCtoPN,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=5,y=4.5,label="r = 0.21",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_s$a254_POM,pom_s$SUVA_POC,method="pearson")

p14 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=a254_POM,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=a254_POM,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=2,y=4.5,label="r = 0.51",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p15 <- ggplot(data=pom_s,mapping=aes(x=SUVA_POC))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab("SUVA")+
        theme_classic(base_size=21)

cor(pom_s$POC_mg,pom_s$HIX_POM,method="pearson")

p16 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=21,label="r = -0.62",size=6)+
        ylab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_s$PN_mg,pom_s$HIX_POM,method="pearson")

p17 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=PN_mg,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=PN_mg,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.6,y=21,label="r = -0.64",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_s$POCtoPN,pom_s$HIX_POM,method="pearson")

p18 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POCtoPN,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POCtoPN,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=12,y=21,label="r = 0.16",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_s$a254_POM,pom_s$HIX_POM,method="pearson")

p19 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=a254_POM,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=a254_POM,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7.5,y=21,label="r = -0.16",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_s$SUVA_POC,pom_s$HIX_POM,method="pearson")

p20 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=SUVA_POC,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=SUVA_POC,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4.5,y=3,label="r = 0.46",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p21 <- ggplot(data=pom_s,mapping=aes(x=HIX_POM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("HIX")+
        theme_classic(base_size=21)

cor(pom_s$POC_mg,pom_s$BIX_POM,method="pearson")

p22 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POC_mg,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POC_mg,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=0.25,label="r = 0.32",size=6)+
        ylab("BIX")+
        xlab(expression(paste("POC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")

cor(pom_s$PN_mg,pom_s$BIX_POM,method="pearson")

p23 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=PN_mg,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=PN_mg,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.2,y=1.2,label="r = 0.35",size=6)+
        xlab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_s$POCtoPN,pom_s$BIX_POM,method="pearson")

p24 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=POCtoPN,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=POCtoPN,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=11,y=1.2,label="r = -0.17",size=6)+
        xlab("POC:PN")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_s$a254_POM,pom_s$BIX_POM,method="pearson")

p25 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=a254_POM,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=a254_POM,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7.5,y=1.2,label="r = 0.05",size=6)+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_s$SUVA_POC,pom_s$BIX_POM,method="pearson")

p26 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=SUVA_POC,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=SUVA_POC,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4.5,y=1.2,label="r = -0.25",size=6)+
        xlab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_s$HIX_POM,pom_s$BIX_POM,method="pearson")

p27 <- ggplot()+
        geom_point(data=pom_s,mapping=aes(x=HIX_POM,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_s,mapping=aes(x=HIX_POM,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=15,y=1.2,label="r = -0.43",size=6)+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

p28 <- ggplot(data=pom_s,mapping=aes(x=BIX_POM))+
        geom_histogram(binwidth = 0.05)+
        ylab("Count")+
        xlab("BIX")+
        theme_classic(base_size=21)

layout <- '
A######
BC#####
DEG####
HIJK###
LMNOP##
QRSUVW#
XYZabcd'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, 
           Q=p16, R=p17, S=p18, U=p19, V=p20, W=p21, X=p22, Y=p23, Z=p24, a=p25, b=p26, c=p27, d=p28, design=layout)

dev.off()

## Bottom POM
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/FigureS9.jpg",width=600,height=400,units="mm",res=800)

p1 <- ggplot(data=pom_b,mapping=aes(x=POC_mg))+
        geom_histogram(binwidth = 0.5)+
        ylab("Count")+
        xlab(expression(paste("POC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)

cor(pom_b$POC_mg,pom_b$PN_mg,method="pearson")

p2 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=PN_mg,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=PN_mg),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=0.1,label="r = 0.89",size=6)+
        ylab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p3 <- ggplot(data=pom_b,mapping=aes(x=PN_mg))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_b$POC_mg,pom_b$POCtoPN,method="pearson")

p4 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=POCtoPN,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=POCtoPN),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=4.5,label="r = 0.00",size=6)+
        ylab("POC:PN")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_b$PN_mg,pom_b$POCtoPN,method="pearson")

p5 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=PN_mg,y=POCtoPN,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=PN_mg,y=POCtoPN),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.35,y=15,label="r = -0.37",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

p6 <- ggplot(data=pom_b,mapping=aes(x=POCtoPN))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("POC:PN")+
        theme_classic(base_size=21)

cor(pom_b$POC_mg,pom_b$a254_POM,method="pearson")

p7 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=1,label="r = 0.28",size=6)+
        ylab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_b$PN_mg,pom_b$a254_POM,method="pearson")

p8 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=PN_mg,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=PN_mg,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.4,y=8.5,label="r = 0.06",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_b$POCtoPN,pom_b$a254_POM,method="pearson")

p9 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POCtoPN,y=a254_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POCtoPN,y=a254_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=6,y=9,label="r = 0.41",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p10 <- ggplot(data=pom_b,mapping=aes(x=a254_POM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("a254")+
        theme_classic(base_size=21)

cor(pom_b$POC_mg,pom_b$SUVA_POC,method="pearson")

p11 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=4.5,label="r = -0.41",size=6)+
        ylab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_b$PN_mg,pom_b$SUVA_POC,method="pearson")

p12 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=PN_mg,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=PN_mg,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.35,y=4.5,label="r = -0.52",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_b$POCtoPN,pom_b$SUVA_POC,method="pearson")

p13 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POCtoPN,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POCtoPN,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=17,y=4.5,label="r = 0.31",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_b$a254_POM,pom_b$SUVA_POC,method="pearson")

p14 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=a254_POM,y=SUVA_POC,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=a254_POM,y=SUVA_POC),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=2,y=4.5,label="r = 0.68",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p15 <- ggplot(data=pom_b,mapping=aes(x=SUVA_POC))+
        geom_histogram(binwidth = 0.1)+
        ylab("Count")+
        xlab("SUVA")+
        theme_classic(base_size=21)

cor(pom_b$POC_mg,pom_b$HIX_POM,method="pearson")

p16 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=21,label="r = -0.38",size=6)+
        ylab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank())

cor(pom_b$PN_mg,pom_b$HIX_POM,method="pearson")

p17 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=PN_mg,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=PN_mg,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.35,y=22,label="r = -0.38",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_b$POCtoPN,pom_b$HIX_POM,method="pearson")

p18 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POCtoPN,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POCtoPN,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=17,y=21,label="r = 0.17",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_b$a254_POM,pom_b$HIX_POM,method="pearson")

p19 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=a254_POM,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=a254_POM,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7.5,y=21,label="r = -0.04",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

cor(pom_b$SUVA_POC,pom_b$HIX_POM,method="pearson")

p20 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=SUVA_POC,y=HIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=SUVA_POC,y=HIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=3,label="r = 0.27",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())

p21 <- ggplot(data=pom_b,mapping=aes(x=HIX_POM))+
        geom_histogram(binwidth = 1)+
        ylab("Count")+
        xlab("HIX")+
        theme_classic(base_size=21)

cor(pom_b$POC_mg,pom_b$BIX_POM,method="pearson")

p22 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POC_mg,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POC_mg,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=3,y=1.2,label="r = 0.19",size=6)+
        ylab("BIX")+
        xlab(expression(paste("POC (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")

cor(pom_b$PN_mg,pom_b$BIX_POM,method="pearson")

p23 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=PN_mg,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=PN_mg,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=0.38,y=1.2,label="r = 0.21",size=6)+
        xlab(expression(paste("PN (mg L"^"-1"*")")))+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank(),
              plot.background=element_rect(fill = "grey93"),
              panel.background=element_rect(fill="grey93"))

cor(pom_b$POCtoPN,pom_b$BIX_POM,method="pearson")

p24 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=POCtoPN,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=POCtoPN,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=15,y=1.2,label="r = -0.11",size=6)+
        xlab("POC:PN")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_b$a254_POM,pom_b$BIX_POM,method="pearson")

p25 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=a254_POM,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=a254_POM,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=7.5,y=1.2,label="r = -0.26",size=6)+
        xlab("a254")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_b$SUVA_POC,pom_b$BIX_POM,method="pearson")

p26 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=SUVA_POC,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=SUVA_POC,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=4,y=1.2,label="r = -0.38",size=6)+
        xlab("SUVA")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

cor(pom_b$HIX_POM,pom_b$BIX_POM,method="pearson")

p27 <- ggplot()+
        geom_point(data=pom_b,mapping=aes(x=HIX_POM,y=BIX_POM,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=pom_b,mapping=aes(x=HIX_POM,y=BIX_POM),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=17,y=1.2,label="r = -0.35",size=6)+
        xlab("HIX")+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())

p28 <- ggplot(data=pom_b,mapping=aes(x=BIX_POM))+
        geom_histogram(binwidth = 0.05)+
        ylab("Count")+
        xlab("BIX")+
        theme_classic(base_size=21)

layout <- '
A######
BC#####
DEG####
HIJK###
LMNOP##
QRSUVW#
XYZabcd'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, 
           Q=p16, R=p17, S=p18, U=p19, V=p20, W=p21, X=p22, Y=p23, Z=p24, a=p25, b=p26, c=p27, d=p28, design=layout)

dev.off()

#############################################################################
## Conduct PCA on each scaled data matrix: no transformations; no removal of outliers
env_s_pca <- rda(env_s_scale)
summary(env_s_pca,axes=0)
plot(env_s_pca)
text(env_s_pca)
screeplot(env_s_pca, bstick = TRUE)

env_b_pca <- rda(env_b_scale)
summary(env_b_pca,axes=0)
plot(env_b_pca)
text(env_b_pca)
screeplot(env_b_pca, bstick = TRUE)

dom_s_pca <- rda(dom_s_scale)
summary(dom_s_pca,axes=0)
plot(dom_s_pca)
text(dom_s_pca)
screeplot(dom_s_pca, bstick = TRUE)

dom_b_pca <- rda(dom_b_scale)
summary(dom_b_pca,axes=0)
plot(dom_b_pca)
text(dom_b_pca)
screeplot(dom_b_pca, bstick = TRUE)

pom_s_pca <- rda(pom_s_scale)
summary(pom_s_pca,axes=0)
plot(pom_s_pca)
text(pom_s_pca)
screeplot(pom_s_pca, bstick=TRUE)

pom_b_pca <- rda(pom_b_scale)
summary(pom_b_pca,axes=0)
plot(pom_b_pca)
text(pom_b_pca)
screeplot(pom_b_pca, bstick=TRUE)

# Make screeplot for each of the PCA output
pdf("C:/Users/ahoun/Desktop/NRE_MultiStats/Fig_Output/Scree_Plots_SandB_All.pdf", width=12, height=8)

par(mfrow=c(2,3))
screeplot(env_s_pca,bstick=TRUE)
screeplot(dom_s_pca,bstick=TRUE)
screeplot(pom_s_pca,bstick=TRUE)
screeplot(env_b_pca,bstick=TRUE)
screeplot(dom_b_pca,bstick=TRUE)
screeplot(pom_b_pca,bstick=TRUE)

dev.off()

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
envspe_s_sc2 <- scores(env_s_pca, choices=1:3, display="sp", scaling=2)
domspe_s_sc2 <- scores(dom_s_pca, choices=1:3, display="sp", scaling=2)
pomspe_s_sc2 <- scores(pom_s_pca, choices=1:3, display="sp", scaling=2)

envspe_b_sc2 <- scores(env_b_pca, choices=1:3, display="sp", scaling=2)
domspe_b_sc2 <- scores(dom_b_pca, choices=1:3, display="sp", scaling=2)
pomspe_b_sc2 <- scores(pom_b_pca, choices=1:3, display="sp", scaling=2)

# Plot in 2D (PC1 and PC2)

pdf("C:/Users/ahoun/Desktop/NRE_MultiStats/Fig_Output/PCA_SandB_Revs.pdf", width=12, height=8)

# Env_S
# Season ordered as: winter, spring, summer, fall
env_s$Season<-factor(env_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(env_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(env_s,levels(Season))
sq<-c(24,21,22,23,25)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_s_pca,type="n",scaling=2,xlab="PC1 (48% var. explained)",ylab="PC2 (23% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_S")
with(env_s,points(env_s_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_s_sc2[,1], envspe_s_sc2[,2], angle=20, col="black")
text(env_s_pca, display = "species", labels=c("Temp","Sal","%DO","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_s
dom_s$Season<-factor(dom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_s_pca,type="n",scaling=2,xlab="PC1 (65% var. explained)",ylab="PC2 (18% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_S")
with(dom_s,points(dom_s_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_s_sc2[,1], domspe_s_sc2[,2], angle=20, col="black")
text(dom_s_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA",
                                              "HIX","BIX"), scaling=2, cex = 0.8, col = "black")

# pom_s
pom_s$Season<-factor(pom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(pom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(pom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(pom_s_pca,type="n",scaling=2,xlab="PC1 (48% var. explained)",ylab="PC2 (21% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_S")
with(pom_s,points(pom_s_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_s,legend("bottomright",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_s_sc2[,1], pomspe_s_sc2[,2], angle=20, col="black")
text(pom_s_pca, display = "species", labels=c("POC","PON","POC:PON","a254","SUVA",
                                              "HIX","BIX"), scaling=2, cex = 0.8, col = "black")


## Then plot bottom results
# Env_B
# Season ordered as: winter, spring, summer, fall
env_b$Season<-factor(env_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(env_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(env_b,levels(Season))
sq<-c(24,21,22,23,25)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_b_pca,type="n",scaling=2,xlab="PC1 (51% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_B")
with(env_b,points(env_b_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_b,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_b_sc2[,1], envspe_b_sc2[,2], angle=20, col="black")
text(env_b_pca, display = "species", labels=c("Temp","Sal","%DO","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_b
dom_b$Season<-factor(dom_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_b_pca,type="n",scaling=2,xlab="PC1 (64% var. explained)",ylab="PC2 (18% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_B")
with(dom_b,points(dom_b_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_b_sc2[,1], domspe_b_sc2[,2], angle=20, col="black")
text(dom_b_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX"), scaling=2, cex = 0.8, col = "black")

# pom_b
pom_b$Season<-factor(pom_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(pom_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(pom_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(pom_b_pca,type="n",scaling=2,xlab="PC1 (42% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_B")
with(pom_b,points(pom_b_pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_b_sc2[,1], pomspe_b_sc2[,2], angle=20, col="black")
text(pom_b_pca, display = "species", labels=c("POC","PON","POC:PON","a254","SUVA","HIX","BIX"), scaling=2, cex = 0.8, col = "black")

dev.off()

## Remove correlated variables from PCA: need to think about removing variables from PCA
# Using a combination of: 1. PCA shows variables are collinear, 2. Data shows variables are correlated >75-80%,
# 3. Use chemical knowledge that the two variables are related (aka: a perfectly, imperfect science)
# Go through each data set individually:

# ENV_S: 1. PCA shows %DO and Chla are closely collinear; 2. Have an r2 about 50%; 3. %DO is partly controlled by
# Chla - in the surface I would expect %DO to be relatively constant but Chla to fluctuate, especially seasonally
# REMOVE: %DO from subsequent analyses - very collinear via PCA, so isn't adding additional information to the
# analysis
# KEEP: Temp, Sal, Turb, Chla
env_s_scale_2 <- env_s_scale[,-c(3)]

# DOM_S: 1. PCA shows lots of collinear variables; 2. Several variables have high r2; 3. Lots of the variables
# should be similar chemically
# REMOVE: DON, a254, HIX
# KEEP: DOC, DOC:DON, SUVA, BIX
dom_s_scale_2 <- dom_s_scale[,c(1,3,5,7)]
chart.Correlation(dom_s_scale_2, histogram=TRUE, method=c("pearson"))

# POM_S: 1. PCA shows lots of collinear variables; 2. Several variables have high r2; 3. Lots of variables are
# chemically similar
# REMOVE: PN
# KEEP: POC, PN:POC, a254, SUVA, HIX, BIX
pom_s_scale_2 <- pom_s_scale[,-c(2)]
chart.Correlation(pom_s_scale_2, histogram=TRUE, method=c("pearson"))

# ENV_B: 1. PCA shows nothing collinear; 2. No variables with high correlations
# KEEP: ALL
env_b_scale_2 <- env_b_scale

# DOM_B: 1. Lots of collinear variables on PCA; 2. Several variables with high r2; 3. Lots of variables chemically
# similar
# REMOVE: a254, SUVA, (highly correlated w/ DOC in PCA space),HIX (highly correlated w/ DOC in PCA space)
# KEEP: DOC, DON, DOC:DON, BIX
dom_b_scale_2 <- dom_b_scale[,c(1,2,3,7)]
chart.Correlation(dom_b_scale_2, histogram=TRUE, method=c("pearson"))

# POM_B: 1. Lots of collinear variables via PCA; 2. Several high r2; 3. Lots of chemically similar variables
# REMOVE: PN
# KEEP: POC, POC:PN, a254, SUVA, HIX, BIX
pom_b_scale_2 <- pom_b_scale[,-c(2)]
chart.Correlation(pom_b_scale_2, histogram=TRUE, method=c("pearson"))

## Conduct PCA on updated data matrices (i.e., with collinear variables removed)
## Conduct PCA on each scaled data martix: no transformations; no removal of outliers
env_s_pca_2 <- rda(env_s_scale_2)
summary(env_s_pca_2,axes=0)
plot(env_s_pca_2)
text(env_s_pca_2)
screeplot(env_s_pca_2, bstick = TRUE)

env_b_pca_2 <- rda(env_b_scale_2)
summary(env_b_pca_2,axes=0)
plot(env_b_pca_2)
text(env_b_pca_2)
screeplot(env_b_pca_2, bstick = TRUE)

dom_s_pca_2 <- rda(dom_s_scale_2)
summary(dom_s_pca_2,axes=0)
plot(dom_s_pca_2)
text(dom_s_pca_2)
screeplot(dom_s_pca_2, bstick = TRUE)

dom_b_pca_2 <- rda(dom_b_scale_2)
summary(dom_b_pca_2,axes=0)
plot(dom_b_pca_2)
text(dom_b_pca_2)
screeplot(dom_b_pca_2, bstick = TRUE)

pom_s_pca_2 <- rda(pom_s_scale_2)
summary(pom_s_pca_2,axes=0)
plot(pom_s_pca_2)
text(pom_s_pca_2)
screeplot(pom_s_pca_2, bstick=TRUE)

pom_b_pca_2 <- rda(pom_b_scale_2)
summary(pom_b_pca_2,axes=0)
plot(pom_b_pca_2)
text(pom_b_pca_2)
screeplot(pom_b_pca_2, bstick=TRUE)

# Make screeplot for each of the PCA output
pdf("C:/Users/ahoun/Desktop/NRE_MultiStats/Fig_Output/Scree_Plots_SandB_2.pdf", width=12, height=8)

par(mfrow=c(2,3))
screeplot(env_s_pca_2,bstick=TRUE)
screeplot(dom_s_pca_2,bstick=TRUE)
screeplot(pom_s_pca_2,bstick=TRUE)
screeplot(env_b_pca_2,bstick=TRUE)
screeplot(dom_b_pca_2,bstick=TRUE)
screeplot(pom_b_pca_2,bstick=TRUE)

dev.off()

# Construct PCA biplot that will be divided by season and will include objects and variables 
# The graph will be in scaling 2
# Extract species scores for scaling 2
envspe_s_sc2 <- scores(env_s_pca_2, choices=1:3, display="sp", scaling=2)
domspe_s_sc2 <- scores(dom_s_pca_2, choices=1:3, display="sp", scaling=2)
pomspe_s_sc2 <- scores(pom_s_pca_2, choices=1:3, display="sp", scaling=2)

envspe_b_sc2 <- scores(env_b_pca_2, choices=1:3, display="sp", scaling=2)
domspe_b_sc2 <- scores(dom_b_pca_2, choices=1:3, display="sp", scaling=2)
pomspe_b_sc2 <- scores(pom_b_pca_2, choices=1:3, display="sp", scaling=2)

# Plot in 2D (PC1 and PC2)

# pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/PCA_SandB_2.pdf", width=12, height=8)


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
        region <- match.arg(region, c("figure", "plot", "device"))
        pos <- match.arg(pos, c("topleft", "top", "topright", 
                                "left", "center", "right", 
                                "bottomleft", "bottom", "bottomright"))
        if(region %in% c("figure", "device")) {
                ds <- dev.size("in")
                # xy coordinates of device corners in user coordinates
                x <- grconvertX(c(0, ds[1]), from="in", to="user")
                y <- grconvertY(c(0, ds[2]), from="in", to="user")
                # fragment of the device we use to plot
                if(region == "figure") {
                        # account for the fragment of the device that 
                        # the figure is using
                        fig <- par("fig")
                        dx <- (x[2] - x[1])
                        dy <- (y[2] - y[1])
                        x <- x[1] + dx * fig[1:2]
                        y <- y[1] + dy * fig[3:4]
                } 
        }
        # much simpler if in plotting region
        if(region == "plot") {
                u <- par("usr")
                x <- u[1:2]
                y <- u[3:4]
        }
        sw <- strwidth(text, cex=cex) * 60/100
        sh <- strheight(text, cex=cex) * 60/100
        x1 <- switch(pos,
                     topleft     =x[1] + sw, 
                     left        =x[1] + sw,
                     bottomleft  =x[1] + sw,
                     top         =(x[1] + x[2])/2,
                     center      =(x[1] + x[2])/2,
                     bottom      =(x[1] + x[2])/2,
                     topright    =x[2] - sw,
                     right       =x[2] - sw,
                     bottomright =x[2] - sw)
        y1 <- switch(pos,
                     topleft     =y[2] - sh,
                     top         =y[2] - sh,
                     topright    =y[2] - sh,
                     left        =(y[1] + y[2])/2,
                     center      =(y[1] + y[2])/2,
                     right       =(y[1] + y[2])/2,
                     bottomleft  =y[1] + sh,
                     bottom      =y[1] + sh,
                     bottomright =y[1] + sh)
        old.par <- par(xpd=NA)
        on.exit(par(old.par))
        text(x1, y1, text, cex=cex, ...)
        return(invisible(c(x,y)))
}

## Re-arranged for publication: 6 panel graph with each PCA result
# Env_S
# Season ordered as: winter, spring, summer, fall
jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure7_r2.jpg",width=300,height=350,units="mm",res=800)

env_s$Season<-factor(env_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(env_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(env_s,levels(Season))
sq<-c(24,21,22,23,25)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(3,2))

plot(env_s_pca_2,type="n",scaling=2,xlab="PC1 (54% var. explained)",ylab="PC2 (23% var. explained)",cex.axis=2,
     cex.lab=2,xlim=c(-3,3))
with(env_s,points(env_s_pca_2,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
arrows(0, 0, envspe_s_sc2[,1], envspe_s_sc2[,2], angle=20, col="black")
text(env_s_pca_2, display = "species", labels=c("","","",""), scaling=2, cex = 1, col = "black")
text(-2.55,1.1,labels="Temp",cex=2,col="black")
text(-2.55,-0.15,labels="Sal",cex=2,col="black")
text(2.85,-0.35,labels="Turb",cex=2,col="black")
text(-1.8,-2.3,labels="Chla",cex=2,col="black")
fig_label("A.",region="figure",pos="topleft",cex=2.5)

# Env_B
# Season ordered as: winter, spring, summer, fall
env_b$Season<-factor(env_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(env_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(env_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(env_b_pca_2,type="n",scaling=2,xlab="PC1 (51% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=2,
     cex.lab=2,xlim=c(-3,3))
with(env_b,points(env_b_pca_2,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
arrows(0, 0, envspe_b_sc2[,1], envspe_b_sc2[,2], angle=20, col="black")
text(env_b_pca_2, display = "species", labels=c("","","","",""), scaling=2, cex = 1, col = "black")
text(-2.55,1.0,labels="Temp",cex=2,col="black")
text(-2.5,-0.65,labels="Sal",cex=2,col="black")
text(2.55,-1.1,labels="%DO",cex=2,col="black")
text(2.55,0.8,labels="Turb",cex=2,col="black")
text(-0.3,-2.6,labels="Chla",cex=2,col="black")
fig_label("B.",region="figure",pos="topleft",cex=2.5)

# dom_s
dom_s$Season<-factor(dom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_s_pca_2,type="n",scaling=2,xlab="PC1 (65% var. explained)",ylab="PC2 (20% var. explained)",cex.axis=2,
     cex.lab=2,xlim=c(-3,3))
with(dom_s,points(dom_s_pca_2,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
with(dom_s,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=2))
arrows(0, 0, domspe_s_sc2[,1], domspe_s_sc2[,2], angle=20, col="black")
text(dom_s_pca_2, display = "species", labels=c("","","",""), scaling=2, cex = 1, 
     col = "black")
text(2.7,0.5,labels="DOC",cex=2,col="black")
text(1.95,2.2,labels="DOC:DON",cex=2,col="black")
text(2.8,-1.0,labels="SUVA",cex=2,col="black")
text(-2.6,1.1,labels="BIX",cex=2,col="black")
fig_label("C.",region="figure",pos="topleft",cex=2.5)

# dom_b
dom_b$Season<-factor(dom_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_b_pca_2,type="n",scaling=2,xlab="PC1 (63% var. explained)",ylab="PC2 (27% var. explained)",cex.axis=2,
     cex.lab=2,xlim=c(-3,3))
with(dom_b,points(dom_b_pca_2,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
arrows(0, 0, domspe_b_sc2[,1], domspe_b_sc2[,2], angle=20, col="black")
text(dom_b_pca_2, display = "species", labels=c("","","",""), scaling=2, cex = 1, 
     col = "black")
text(-2.9,-0.2,labels="DOC",cex=2,col="black")
text(-2.5,1.7,labels="DON",cex=2,col="black")
text(-2.2,-2.3,labels="DOC:DON",cex=2,col="black")
text(2.6,0.1,labels="BIX",cex=2,col="black")
fig_label("D.",region="figure",pos="topleft",cex=2.5)

# pom_s
pom_s$Season<-factor(pom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(pom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(pom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(pom_s_pca_2,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (27% var. explained)",cex.axis=2,
     cex.lab=2,xlim=c(-3,3))
with(pom_s,points(pom_s_pca_2,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
arrows(0, 0, pomspe_s_sc2[,1], pomspe_s_sc2[,2], angle=20, col="black")
text(pom_s_pca_2, display = "species", labels=c("","","","","",""), scaling=2, 
     cex = 1, col = "black")
text(2.0,-1.2,labels="POC",cex=2,col="black")
text(-0.7,-1.5,labels="POC:PN",cex=2,col="black")
text(0.2,-2.5,labels="a254",cex=2,col="black")
text(-1.9,-1.5,labels="SUVA",cex=2,col="black")
text(-2.4,0.3,labels="HIX",cex=2,col="black")
text(1.9,0.0,labels="BIX",cex=2,col="black")
fig_label("E.",region="figure",pos="topleft",cex=2.5)

# pom_b
pom_b$Season<-factor(pom_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(pom_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(pom_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(pom_b_pca_2,type="n",scaling=2,xlab="PC1 (38% var. explained)",ylab="PC2 (26% var. explained)",cex.axis=2,
     cex.lab=2)
with(pom_b,points(pom_b_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=2))
arrows(0, 0, pomspe_b_sc2[,1], pomspe_b_sc2[,2], angle=20, col="black")
text(pom_b_pca_2, display = "species", labels=c("","","","","",""), scaling=2, cex = 1, col = "black")
text(1.3,-1.9,labels="POC",cex=2,col="black")
text(-1.8,-0.9,labels="POC:PN",cex=2,col="black")
text(-1.7,-1.9,labels="a254",cex=2,col="black")
text(-2.5,-0.2,labels="SUVA",cex=2,col="black")
text(-1.5,1.4,labels="HIX",cex=2,col="black")
text(1.8,-0.5,labels="BIX",cex=2,col="black")
fig_label("F.",region="figure",pos="topleft",cex=2.5)

dev.off()

###################################### RDA ######################################
# Now conduct RDA on each DOM or POM data matrix using pared down variables
# DOM_S
dom_s_scale_2 <- as.data.frame(dom_s_scale_2)
env_s_scale_2 <- as.data.frame(env_s_scale_2)
dom_s_rda <- rda(dom_s_scale_2~.,env_s_scale_2,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a_all <- RsquareAdj(dom_s_rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(dom_s_rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: selected - Sal, Temp, Turb
dom_s_rda_forsel <- forward.sel(dom_s_scale_2,env_s_scale_2,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,
                                adjR2thresh=R2a_all)

# Check backward elimination
dom_s_rda_back <- ordistep(dom_s_rda,direction="backward",permutations=how(nperm=999))

# DOM RDA model: Sal, Temp, Turb
dom_s_rda_final <- rda(dom_s_scale_2 ~ env_s_scale_2$Sal + env_s_scale_2$Temp + env_s_scale_2$Turb,scale=FALSE)

# Plot
dom_s_rspe_sc2 <- scores(dom_s_rda_final, display="sp", choices=c(1,2), scaling=2)
dom_s_rbp.sc2 <- scores(dom_s_rda_final,display="bp",choices=c(1,2),scaling=2)

# POM_S
pom_s_rda <- rda(pom_s_scale_2~.,env_s_scale_2,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a_all <- RsquareAdj(pom_s_rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(pom_s_rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: selected - Chla, Turb, Sal
pom_s_rda_forsel <- forward.sel(pom_s_scale_2,env_s_scale_2,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,
                                adjR2thresh=R2a_all)

# Check backward elimination
pom_s_rda_back <- ordistep(pom_s_rda,direction="backward",permutations=how(nperm=999))

# POM RDA model: Chla, Turb, Sal
pom_s_rda_final <- rda(pom_s_scale_2 ~ env_s_scale_2$Chla + env_s_scale_2$Turb + env_s_scale_2$Sal,scale=FALSE)

# Plot
pom_s_rspe_sc2 <- scores(pom_s_rda_final, display="sp", choices=c(1,2), scaling=2)
pom_s_rbp.sc2 <- scores(pom_s_rda_final,display="bp",choices=c(1,2),scaling=2)

# Convert to data.frame
dom_b_scale_2 <- data.frame(dom_b_scale_2)
pom_b_scale_2 <- data.frame(pom_b_scale_2)
env_b_scale_2 <- data.frame(env_b_scale_2)
# DOM_B
dom_b_rda <- rda(dom_b_scale_2~.,env_b_scale_2,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a_all <- RsquareAdj(dom_b_rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(dom_b_rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: selected - Sal, Temp, DO_Sat, Chla
dom_b_rda_forsel <- forward.sel(dom_b_scale_2,env_b_scale_2,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,
                                adjR2thresh=R2a_all)

# Check backward elimination
dom_b_rda_back <- ordistep(dom_b_rda,direction="backward",permutations=how(nperm=999))

# DOM RDA model: Sal, Temp, DO_Sat, Chla
dom_b_rda_final <- rda(dom_b_scale_2 ~ env_b_scale_2$Sal + env_b_scale_2$Temp + env_b_scale_2$DO_Sat + 
                               env_b_scale_2$Chla, scale=FALSE)

# Plot
dom_b_rspe_sc2 <- scores(dom_b_rda_final, display="sp", choices=c(1,2), scaling=2)
dom_b_rbp.sc2 <- scores(dom_b_rda_final,display="bp",choices=c(1,2),scaling=2)

# POM_B
pom_b_rda <- rda(pom_b_scale_2~.,env_b_scale_2,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a_all <- RsquareAdj(pom_b_rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(pom_b_rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: selected - Chla, Turb, Sal
pom_b_rda_forsel <- forward.sel(pom_b_scale_2,env_b_scale_2,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,
                                adjR2thresh=R2a_all)

# Check backward elimination
pom_b_rda_back <- ordistep(pom_b_rda,direction="backward",permutations=how(nperm=999))

# POM RDA model: Chla, Turb, Sal
## CHANGED: Turb, Chla, Sal
pom_b_rda_final <- rda(pom_b_scale_2 ~ env_b_scale_2$Turb + env_b_scale_2$Chla + env_b_scale_2$Sal,scale=FALSE)

# Plot
pom_b_rspe_sc2 <- scores(pom_b_rda_final, display="sp", choices=c(1,2), scaling=2)
pom_b_rbp.sc2 <- scores(pom_b_rda_final,display="bp",choices=c(1,2),scaling=2)

## Plot RDA results
#pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/RDA_SandB.pdf", width=12, height=8)

jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure8.jpg",width=250,height=250,units="mm",res=800)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

dom_s$Season<-factor(dom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (93% fitted,46% total var.)",
     ylab="RDA2 (6% fitted, 3% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-2.5,5))
with(dom_s,points(dom_s_rda_final,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_s,legend("topleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_s_rspe_sc2[,1], dom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_s_rda_final,display = "species", labels=c("","","SUVA","BIX"), scaling=2, cex = 1.5, 
     col = "black")
text(1.8,0.3,labels="DOC",cex=1.5,col="black")
rect(0.5,0.5,3.1,1.1,col="white",border=NA,density=70)
text(1.6,0.8,labels="DOC:DON",cex=1.5,col="black")
#text(2.6,0.1,labels="SUVA",cex=1.5,col="black")
#text(-2.6,0.4,labels="BIX",cex=1.5,col="black")
text(dom_s_rda_final,display="bp",labels=c("Sal","Temp","Turb"),scaling=2,cex=1.5,col="black")
fig_label("A.",region="figure",pos="topleft",cex=2)

plot(pom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (63% fitted, 28% total var.)",
     ylab="RDA2 (35% fitted, 16% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-3,1.5),xlim=c(-2,3))
with(dom_s,points(pom_s_rda_final,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_s_rspe_sc2[,1], pom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_s_rda_final,display = "species", labels=c("POC","","","SUVA","HIX",""), scaling=2, 
     cex = 1.5, col = "black")
#text(2.55,-0.2,labels="POC",cex=1.5,col="black")
rect(-1.2,-0.85,0,-0.55,col="white",border=NA,density=70)
text(-0.6,-0.7,labels="POC:PN",cex=1.5,col="black")
text(0.7,-1.85,labels="a254",cex=1.5,col="black")
#text(-1.65,-1.3,labels="SUVA",cex=1.5,col="black")
#text(-1.8,-0.3,labels="HIX",cex=1.5,col="black")
text(1.1,0.3,labels="BIX",cex=1.5,col="black")
text(pom_s_rda_final,display="bp",labels=c("Chla","Turb","Sal"),scaling=2,cex=1.5,col="black")
fig_label("B.",region="figure",pos="topleft",cex=2)

## Bottom
env_b$Season<-factor(env_b$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(env_b,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(env_b,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_b_rda_final,scaling=2,display="sites",xlab="RDA1 (90% fitted, 46% total var.)",
     ylab="RDA2 (10% fitted, 5% total var.)",cex.axis=1.5,cex.lab=1.5,xlim=c(-4,3.5))
with(env_b,points(dom_b_rda_final,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
#with(env_b,legend("topleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
#                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_b_rspe_sc2[,1], dom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_b_rda_final,display = "species", labels=c("","DON","","BIX"), scaling=2, cex = 1.5, 
     col = "black")
text(-2.1,0.1,labels="DOC",cex=1.5,col="black")
#text(-2.65,-0.7,labels="DON",cex=1.5,col="black")
text(-1.8,-1.0,labels="DOC:DON",cex=1.5,col="black")
#text(2.6,-0.2,labels="BIX",cex=1.5,col="black")
text(dom_b_rda_final,display="bp",labels=c("Sal","Temp","DO","Chla"),scaling=2,cex=1.5,col="black")
fig_label("C.",region="figure",pos="topleft",cex=2)

# POM Bottom
plot(pom_b_rda_final,scaling=2,display="sites",xlab="RDA1 (62% fitted, 26% total var.)",
     ylab="RDA2 (36% fitted, 15% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,3),xlim=c(-2.5,3))
with(env_b,points(pom_b_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_b_rspe_sc2[,1], pom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_b_rda_final,display = "species", labels=c("","","a254","SUVA","",""), scaling=2, 
     cex = 1.5, col = "black")
text(1.2,1.8,labels="POC",cex=1.5,col="black")
text(-1.5,0.2,labels="POC:PN",cex=1.5,col="black")
#text(-0.6,1.9,labels="a254",cex=1.5,col="black")
#text(-1.8,1.05,labels="SUVA",cex=1.5,col="black")
rect(-1.1,-0.35,-0.6,0,col="white",border=NA,density=70)
text(-0.9,-0.2,labels="HIX",cex=1.5,col="black")
text(pom_b_rda_final,display="bp",labels=c("Turb","Chla","Sal"),scaling=2,cex=1.5,col="black")
rect(0.6,-0.15,1.0,0.15,col="white",border=NA,density=70)
text(0.8,0,labels="BIX",cex=1.5,col="black")
fig_label("D.",region="figure",pos="topleft",cex=2)

dev.off()

