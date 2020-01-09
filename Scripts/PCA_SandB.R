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

# Save rfile as: PCA_SandB

# Load in libraries need
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc,ggplot2,tidyverse,vegan3d,
               scatterplot3d,rgl,car)

# Load in data (Database_DOSat.csv)
my_data <- read.csv(file.choose())
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

# Normalize data using 'scale': do not include station or depth as variables
# Separate data by data pool and by surface and bottom
env_s_pca <- env_s[,c(6:10)]
env_b_pca <- env_b[,c(6:10)]
dom_s_pca <- dom_s[,c(6:18)]
dom_b_pca <- dom_b[,c(6:18)]
pom_s_pca <- pom_s[,c(6:18)]
pom_b_pca <- pom_b[,c(6:18)]

# Normalize data: subract mean/SD
env_s_scale <- scale(env_s_pca)
env_b_scale <- scale(env_b_pca)
dom_s_scale <- scale(dom_s_pca)
dom_b_scale <- scale(dom_b_pca)
pom_s_scale <- scale(pom_s_pca)
pom_b_scale <- scale(pom_b_pca)

# Plot correlation charts for each data matrix
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/Corr_Plots_SandB.pdf", width=12, height=8)

chart.Correlation(env_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(env_b_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(dom_b_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_s_scale, histogram=TRUE, method=c("pearson"))
chart.Correlation(pom_b_scale, histogram=TRUE, method=c("pearson"))

dev.off()

## Conduct PCA on each scaled data martix: no transformations; no removal of outliers
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
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/Scree_Plots_SandB.pdf", width=12, height=8)

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

pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/PCA_SandB.pdf", width=12, height=8)

# Env_S
# Season ordered as: winter, spring, summer, fall
env_s$Season<-factor(env_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_s,levels(Season))
sq<-c(21,22,23,24)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_s_pca,type="n",scaling=2,xlab="PC1 (48% var. explained)",ylab="PC2 (23% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_S")
with(env_s,points(env_s_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_s_sc2[,1], envspe_s_sc2[,2], angle=20, col="black")
text(env_s_pca, display = "species", labels=c("Temp","Sal","%DO","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_s
dom_s$Season<-factor(dom_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(dom_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(dom_s,levels(Season))
sq<-c(21,22,23,24)

plot(dom_s_pca,type="n",scaling=2,xlab="PC1 (65% var. explained)",ylab="PC2 (18% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_S")
with(dom_s,points(dom_s_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_s_sc2[,1], domspe_s_sc2[,2], angle=20, col="black")
text(dom_s_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C","M",
                                              "N"), scaling=2, cex = 0.8, col = "black")

# pom_s
pom_s$Season<-factor(pom_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(pom_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(pom_s,levels(Season))
sq<-c(21,22,23,24)

plot(pom_s_pca,type="n",scaling=2,xlab="PC1 (48% var. explained)",ylab="PC2 (21% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_S")
with(pom_s,points(pom_s_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_s,legend("bottomright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_s_sc2[,1], pomspe_s_sc2[,2], angle=20, col="black")
text(pom_s_pca, display = "species", labels=c("POC","PON","POC:PON","a254","SUVA","HIX","BIX","B","T","A","C","M",
                                              "N"), scaling=2, cex = 0.8, col = "black")


## Then plot bottom results
# Env_B
# Season ordered as: winter, spring, summer, fall
env_b$Season<-factor(env_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_b,levels(Season))
sq<-c(21,22,23,24)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_b_pca,type="n",scaling=2,xlab="PC1 (51% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_B")
with(env_b,points(env_b_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_b,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_b_sc2[,1], envspe_b_sc2[,2], angle=20, col="black")
text(env_b_pca, display = "species", labels=c("Temp","Sal","%DO","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_b
dom_b$Season<-factor(dom_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(dom_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(dom_b,levels(Season))
sq<-c(21,22,23,24)

plot(dom_b_pca,type="n",scaling=2,xlab="PC1 (64% var. explained)",ylab="PC2 (18% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_B")
with(dom_b,points(dom_b_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_b_sc2[,1], domspe_b_sc2[,2], angle=20, col="black")
text(dom_b_pca, display = "species", labels=c("DOC","DON","DOC:DON","a254","SUVA","HIX","BIX","B","T","A","C","M",
                                              "N"), scaling=2, cex = 0.8, col = "black")

# pom_b
pom_b$Season<-factor(pom_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(pom_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(pom_b,levels(Season))
sq<-c(21,22,23,24)

plot(pom_b_pca,type="n",scaling=2,xlab="PC1 (42% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_B")
with(pom_b,points(pom_b_pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_b_sc2[,1], pomspe_b_sc2[,2], angle=20, col="black")
text(pom_b_pca, display = "species", labels=c("POC","PON","POC:PON","a254","SUVA","HIX","BIX","B","T","A","C","M",
                                              "N"), scaling=2, cex = 0.8, col = "black")

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
# REMOVE: DON, SUVA, HIX, A, C, B, M, N
# KEEP: DOC, DOC:DON, a254, BIX, T
dom_s_scale_2 <- dom_s_scale[,c(1,3:4,7,9)]

# POM_S: 1. PCA shows lots of collinear variables; 2. Several variables have high r2; 3. Lots of variables are
# chemically similar
# REMOVE: PN, A, C, M, B, N
# KEEP: POC, PN:POC, a254, SUVA, HIX, BIX, T
pom_s_scale_2 <- pom_s_scale[,-c(2,8,10:13)]

# ENV_B: 1. PCA shows nothing collinear; 2. No variables with high correlations
# KEEP: ALL
env_b_scale_2 <- env_b_scale

# DOM_B: 1. Lots of collinear variables on PCA; 2. Several variables with high r2; 3. Lots of variables chemically
# similar
# REMOVE: DON, A, C, M, N, HIX, a254, B
# KEEP: DOC, DOC:DON, SUVA, BIX, T
dom_b_scale_2 <- dom_b_scale[,c(1,3,5,7,9)]

# POM_B: 1. Lots of collinear variables via PCA; 2. Several high r2; 3. Lots of chemically similar variables
# REMOVE: PN, A, C, M, B
# KEEP: POC, POC:PN, a254, SUVA, HIX, BIX, T, N
pom_b_scale_2 <- pom_b_scale[,-c(2,8,10:12)]

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
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/Scree_Plots_SandB_2.pdf", width=12, height=8)

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

pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/PCA_SandB_2.pdf", width=12, height=8)

# Env_S
# Season ordered as: winter, spring, summer, fall
env_s$Season<-factor(env_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_s,levels(Season))
sq<-c(21,22,23,24)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_s_pca_2,type="n",scaling=2,xlab="PC1 (54% var. explained)",ylab="PC2 (23% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_S")
with(env_s,points(env_s_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_s_sc2[,1], envspe_s_sc2[,2], angle=20, col="black")
text(env_s_pca_2, display = "species", labels=c("Temp","Sal","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_s
dom_s$Season<-factor(dom_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(dom_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(dom_s,levels(Season))
sq<-c(21,22,23,24)

plot(dom_s_pca_2,type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC2 (22% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_S")
with(dom_s,points(dom_s_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_s,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_s_sc2[,1], domspe_s_sc2[,2], angle=20, col="black")
text(dom_s_pca_2, display = "species", labels=c("DOC","DOC:DON","a254","BIX","T"), scaling=2, cex = 0.8, 
     col = "black")

# pom_s
pom_s$Season<-factor(pom_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(pom_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(pom_s,levels(Season))
sq<-c(21,22,23,24)

plot(pom_s_pca_2,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_S")
with(pom_s,points(pom_s_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_s,legend("bottomright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_s_sc2[,1], pomspe_s_sc2[,2], angle=20, col="black")
text(pom_s_pca_2, display = "species", labels=c("POC","POC:PON","a254","SUVA","HIX","BIX","T"), scaling=2, 
     cex = 0.8, col = "black")


## Then plot bottom results
# Env_B
# Season ordered as: winter, spring, summer, fall
env_b$Season<-factor(env_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_b,levels(Season))
sq<-c(21,22,23,24)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env_b_pca_2,type="n",scaling=2,xlab="PC1 (51% var. explained)",ylab="PC2 (24% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="ENV_B")
with(env_b,points(env_b_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_b,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe_b_sc2[,1], envspe_b_sc2[,2], angle=20, col="black")
text(env_b_pca_2, display = "species", labels=c("Temp","Sal","%DO","Turb","Chla"), scaling=2, cex = 0.8, col = "black")
# USE EVENTUALLY
# text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
# text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
# text(2.3,-1.8,labels="DO",cex=1.5,col="black")
# text(2.2,0.6,labels="Turb",cex=1.5,col="black")
# text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

# dom_b
dom_b$Season<-factor(dom_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(dom_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(dom_b,levels(Season))
sq<-c(21,22,23,24)

plot(dom_b_pca_2,type="n",scaling=2,xlab="PC1 (60% var. explained)",ylab="PC2 (17% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="DOM_B")
with(dom_b,points(dom_b_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(dom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, domspe_b_sc2[,1], domspe_b_sc2[,2], angle=20, col="black")
text(dom_b_pca_2, display = "species", labels=c("DOC","DOC:DON","SUVA","BIX","T"), scaling=2, cex = 0.8, 
     col = "black")

# pom_b
pom_b$Season<-factor(pom_b$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(pom_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(pom_b,levels(Season))
sq<-c(21,22,23,24)

plot(pom_b_pca_2,type="n",scaling=2,xlab="PC1 (38% var. explained)",ylab="PC2 (28% var. explained)",cex.axis=1.5,
     cex.lab=1.5,main="POM_B")
with(pom_b,points(pom_b_pca_2,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(pom_b,legend("topright",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, pomspe_b_sc2[,1], pomspe_b_sc2[,2], angle=20, col="black")
text(pom_b_pca_2, display = "species", labels=c("POC","POC:PON","a254","SUVA","HIX","BIX","T",
                                              "N"), scaling=2, cex = 0.8, col = "black")

dev.off()

###################################### RDA ######################################
# Now conduct RDA on each DOM or POM data matrix using pared down variables
# DOM_S
dom_s_rda <- rda(dom_s_scale_2~.,env_s_scale_2,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a_all <- RsquareAdj(dom_s_rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(dom_s_rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: selected - Sal, Temp, Turb
dom_s_rda_forsel <- forward.sel(dom_s_scale_2,env_s_scale_2,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,
                                adjR2thresh=R2a_all)

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

# DOM RDA model: Sal, Temp, Turb
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

# POM RDA model: Chla, Turb, Sal
pom_b_rda_final <- rda(pom_b_scale_2 ~ env_b_scale_2$Chla + env_b_scale_2$Turb + env_b_scale_2$Sal,scale=FALSE)

# Plot
pom_b_rspe_sc2 <- scores(pom_b_rda_final, display="sp", choices=c(1,2), scaling=2)
pom_b_rbp.sc2 <- scores(pom_b_rda_final,display="bp",choices=c(1,2),scaling=2)

## Plot RDA results
pdf("C:/Users/ahoun/OneDrive/Desktop/NRE_MultiStats/Plots/RDA_SandB.pdf", width=12, height=8)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

env_s$Season<-factor(env_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_s,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_s,levels(Season))
sq<-c(21,22,23,24)

plot(dom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (86% fitted,35% total var.)",
     ylab="RDA2 (9% fitted, 4% total var.)",main="DOM_S",cex.axis=1.5,cex.lab=1.5)
with(env_s,points(dom_s_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_s,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_s_rspe_sc2[,1], dom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_s_rda_final,display = "species", labels=c("DOC","DOC:DON","a254","BIX","T"), scaling=2, cex = 1.5, 
     col = "black")
# For later
# text(-2.7,0.1,labels="DOC",cex=1.5,col="black")
# text(-2.5,-0.4,labels="SUVA",cex=1.5,col="black")
# text(-0.24,-1,labels="B",cex=1.5,col="white")
# text(-0.7,-1.2,labels="T",cex=1.5,col="white")
text(dom_s_rda_final,display="bp",labels=c("Sal","Temp","Turb"),scaling=2,cex=1.5,col="black")
# text(0.2,-0.6,labels="Chla",cex=1.5,col="white")

plot(pom_s_rda_final,scaling=2,display="sites",xlab="RDA1 (65% fitted, 28% total var.)",
     ylab="RDA2 (33% fitted, 14% total var.)",main="POM_S",cex.axis=1.5,cex.lab=1.5,xlim=c(-2,3),ylim=c(-2,3))
with(env_s,points(pom_s_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_s_rspe_sc2[,1], pom_s_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_s_rda_final,display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX","T"), scaling=2, 
     cex = 1.5, col = "black")
text(pom_s_rda_final,display="bp",labels=c("Chla","Turb","Sal"),scaling=2,cex=1.5,col="black")

## Bottom
env_b$Season<-factor(env_s$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(env_b,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(env_b,levels(Season))
sq<-c(21,22,23,24)

plot(dom_b_rda_final,scaling=2,display="sites",xlab="RDA1 (89% fitted, 43% total var.)",
     ylab="RDA2 (9% fitted, 4% total var.)",main="DOM_B",cex.axis=1.5,cex.lab=1.5)
with(env_b,points(dom_b_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
with(env_b,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),
                  pch=c(21,22,23,24),pt.bg=colvec,cex=1.5))
arrows(0,0,dom_b_rspe_sc2[,1], dom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom_b_rda_final,display = "species", labels=c("DOC","DOC:DON","SUVA","BIX","T"), scaling=2, cex = 1.5, 
     col = "black")
# For later
# text(-2.7,0.1,labels="DOC",cex=1.5,col="black")
# text(-2.5,-0.4,labels="SUVA",cex=1.5,col="black")
# text(-0.24,-1,labels="B",cex=1.5,col="white")
# text(-0.7,-1.2,labels="T",cex=1.5,col="white")
text(dom_b_rda_final,display="bp",labels=c("Sal","Temp","DO_Sat","Chla"),scaling=2,cex=1.5,col="black")
# text(0.2,-0.6,labels="Chla",cex=1.5,col="white")

plot(pom_b_rda_final,scaling=2,display="sites",xlab="RDA1 (52% fitted, 20% total var.)",
     ylab="RDA2 (44% fitted, 17% total var.)",main="POM_B",cex.axis=1.5,cex.lab=1.5,xlim=c(-2,3),ylim=c(-2,3))
with(env_b,points(pom_b_rda_final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],
                  bg=colvec[Season],cex=1.5))
arrows(0,0,pom_b_rspe_sc2[,1], pom_b_rspe_sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom_b_rda_final,display = "species", labels=c("POC","POC:PN","a254","SUVA","HIX","BIX","T","N"), scaling=2, 
     cex = 1.5, col = "black")
text(pom_b_rda_final,display="bp",labels=c("Chla","Turb","Sal"),scaling=2,cex=1.5,col="black")

dev.off()

