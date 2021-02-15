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
        theme_classic(base_size=21)

p2 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=Sal,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=Sal),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=20,label="r = 0.39",size=6)+
        ylab(expression('Temp ('*degree*C*')'))+
        theme_classic(base_size=21)+
        theme(legend.position = "none")


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
        theme(legend.position = "none",plot.background=element_rect(fill = "grey93"),
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
        theme_classic(base_size=21)

p11 <- ggplot()+
        geom_point(data=env_s,mapping=aes(x=Temp,y=Chla,color=Season,shape=Season,fill=Season),size=3)+
        scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_fill_manual(values=c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB"))+
        scale_shape_manual(values=c(24,21,22,23,25))+
        geom_smooth(data=env_s,mapping=aes(x=Temp,y=Turb),method=lm,se=FALSE,color="black")+
        annotate(geom="text",x=10,y=130,label="r = 0.12",size=6)+
        ylab(expression(paste("Surface Chla (",mu,"g L"^"-1"*")")))+
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
        ylab("Turb (NTU")+
        annotate(geom="text",x=20,y=130,label="r = -0.28",size=6)+
        theme_classic(base_size=21)+
        theme(legend.position = "none",axis.title.y=element_blank())


p15 <- ggplot(data=env_s,mapping=aes(x=Chla))+
        geom_histogram(binwidth=1)+
        ylab("Count")+
        theme_classic(base_size=21)

layout <- '
A####
BC###
DEG##
HIJK#
LMNOP'

wrap_plots(A=p1, B=p2, C=p3, D=p4, E=p5, G=p6, H=p7, I=p8, J=p9, K=p10, L=p11, M=p12, N=p13, O=p14, P=p15, design=layout)

dev.off()

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


## Re-arranged for publication: 6 panel graph with each PCA result
# Env_S
# Season ordered as: winter, spring, summer, fall

jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure5_Rev.jpg",width=300,height=350,units="mm",res=800)

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
with(dom_s,legend("topleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=2))
arrows(0, 0, domspe_s_sc2[,1], domspe_s_sc2[,2], angle=20, col="black")
text(dom_s_pca_2, display = "species", labels=c("DOC","DOC:DON","SUVA","BIX"), scaling=2, cex = 1, 
     col = "black")
#text(2.6,-0.7,labels="DOC",cex=2,col="black")
#text(1.95,1.15,labels="DOC:DON",cex=2,col="black")
#text(2.6,0.4,labels="a254",cex=2,col="black")
#text(-1.8,-1.6,labels="BIX",cex=2,col="black")
#text(1.3,-2.4,labels="T",cex=2,col="black")

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
text(dom_b_pca_2, display = "species", labels=c("DOC","DON","DOC:DON","BIX"), scaling=2, cex = 1, 
     col = "black")
#text(2.5,-0.7,labels="DOC",cex=2,col="black")
#text(1.7,1.4,labels="DOC:DON",cex=2,col="black")
#text(2.4,0.6,labels="SUVA",cex=2,col="black")
#text(-2.55,-0.4,labels="BIX",cex=2,col="black")
#text(1.1,-2.5,labels="T",cex=2,col="black")

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
text(pom_s_pca_2, display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX"), scaling=2, 
     cex = 1, col = "black")
#text(2.6,-0.3,labels="POC",cex=2,col="black")
#text(-0.2,-1.4,labels="POC:PN",cex=2,col="black")
#text(0.9,-2.3,labels="a254",cex=2,col="black")
#text(-1.7,-1.9,labels="SUVA",cex=2,col="black")
#text(-2.3,-0.3,labels="HIX",cex=2,col="black")
#text(1.9,0.4,labels="BIX",cex=2,col="black")
#text(1.9,-0.8,labels="T",cex=2,col="black")

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
text(pom_b_pca_2, display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX"), scaling=2, cex = 1, col = "black")
#text(2.2,-0.2,labels="POC",cex=2,col="black")
#text(0,1.5,labels="POC:PN",cex=2,col="black")
#text(0.9,2.2,labels="a254",cex=2,col="black")
#text(-0.6,2.1,labels="SUVA",cex=2,col="black")
#text(-1.4,0.8,labels="HIX",cex=2,col="black")
#text(1.2,-1.1,labels="BIX",cex=2,col="black")
#text(2.2,0.3,labels="T",cex=2,col="black")
#text(2.2,0.7,labels="N",cex=2,col="black")

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

jpeg("C:/Users/ahoun/Desktop/NRE_Multistats/Fig_Output/Figure6_Rev.jpg",width=250,height=250,units="mm",res=800)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

dom_s$Season<-factor(dom_s$Season, levels=c("Summer15","Fall","Winter", "Spring", "Summer16"))
with(dom_s,levels(Season))
colvec<-c("#D81B60","#FFC107","#1E88E5","#004D40","#F7C1BB")
with(dom_s,levels(Season))
sq<-c(24,21,22,23,25)

plot(dom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (93% fitted,46% total var.)",
     ylab="RDA2 (6% fitted, 3% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-2.5,5))
with(env_s,points(dom_s_rda_final,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
#with(env_s,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
#                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_s_rspe_sc2[,1], dom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_s_rda_final,display = "species", labels=c("DOC","DOC:DON","SUVA","BIX"), scaling=2, cex = 1.5, 
     col = "black")
#text(1.6,0.6,labels="DOC",cex=1.5,col="black")
#text(2.4,-0.55,labels="DOC:DON",cex=1.5,col="black")
#text(2.6,0.1,labels="a254",cex=1.5,col="black")
#text(-2.6,0.4,labels="BIX",cex=1.5,col="black")
#text(0.4,1.0,labels="T",cex=1.5,col="white")
text(dom_s_rda_final,display="bp",labels=c("Sal","Temp","Turb"),scaling=2,cex=1.5,col="black")

plot(pom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (63% fitted, 28% total var.)",
     ylab="RDA2 (35% fitted, 16% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-3,1.5),xlim=c(-2,3))
with(env_s,points(pom_s_rda_final,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_s_rspe_sc2[,1], pom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_s_rda_final,display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX"), scaling=2, 
     cex = 1.5, col = "black")
#text(2.55,-0.2,labels="POC",cex=1.5,col="black")
#text(-0.4,-0.8,labels="POC:PN",cex=1.5,col="black")
#text(0.8,-1.85,labels="a254",cex=1.5,col="black")
#text(-1.65,-1.3,labels="SUVA",cex=1.5,col="black")
#text(-1.8,-0.3,labels="HIX",cex=1.5,col="black")
#text(1.1,0.4,labels="BIX",cex=1.5,col="black")
#text(1.3,-0.5,labels="T",cex=1.5,col="black")
text(pom_s_rda_final,display="bp",labels=c("Chla","Turb","Sal"),scaling=2,cex=1.5,col="black")

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
with(env_b,legend("topleft",legend=levels(Season),bty="n",col=c("black","black","black","black","black"),
                  pch=c(24,21,22,23,25),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_b_rspe_sc2[,1], dom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_b_rda_final,display = "species", labels=c("DOC","DON","DOC:DON","BIX"), scaling=2, cex = 1.5, 
     col = "black")
#text(-2.6,0.35,labels="DOC",cex=1.5,col="black")
#text(-2.65,-0.7,labels="DOC:DON",cex=1.5,col="black")
#text(-2.8,-0.2,labels="SUVA",cex=1.5,col="black")
#text(2.6,-0.2,labels="BIX",cex=1.5,col="black")
#text(-0.7,0.9,labels="T",cex=1.5,col="black")
text(dom_b_rda_final,display="bp",labels=c("Sal","Temp","DO","Chla"),scaling=2,cex=1.5,col="black")

# POM Bottom
plot(pom_b_rda_final,scaling=2,display="sites",xlab="RDA1 (62% fitted, 26% total var.)",
     ylab="RDA2 (36% fitted, 15% total var.)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,3),xlim=c(-2.5,3))
with(env_b,points(pom_b_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_b_rspe_sc2[,1], pom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_b_rda_final,display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX"), scaling=2, 
     cex = 1.5, col = "black")
#text(2.1,1,labels="POC",cex=1.5,col="black")
#text(-0.9,0.75,labels="POC:PN",cex=1.5,col="black")
#text(-0.6,1.9,labels="a254",cex=1.5,col="black")
#text(-1.8,1.05,labels="SUVA",cex=1.5,col="black")
#text(-0.8,0.3,labels="HIX",cex=1.5,col="white")
#text(0.8,-0.5,labels="BIX",cex=1.5,col="black")
#text(-0.7,0.9,labels="T",cex=1.5,col="black")
#text(0.9,1.1,labels="N",cex=1.5,col="black")
text(pom_b_rda_final,display="bp",labels=c("Turb","Chla","Sal"),scaling=2,cex=1.5,col="black")

dev.off()

