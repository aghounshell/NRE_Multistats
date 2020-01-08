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

## Remove correlated variables from PCA: using VIF's
# Env_s
env_s_scale <- data.frame(env_s_scale)

temp_m <- lm(Temp~.,data=env_s_scale)
vif(temp_m)
# Total vif = 1.47

sal_m <- lm(Sal~.,data=env_s_scale)
vif(sal_m)

# DOM_s
dom_s_scale <- data.frame(dom_s_scale)

doc_m <- lm(DOC_mg~.,data=dom_s_scale)
vif(doc_m)
# Total vif = 250

don_m <- lm(DON_mg~.,data=dom_s_scale)
vif(don_m)
# Total vif = 70.92

# Remove DOC, A_DOM, C_DOM, M_DOM, a254_DOM from DOM_S (highest VIF; then re-model)
dom_s_scale_2 <- dom_s_scale[,-c(1,4,10:12)]

don_m <- lm(DON_mg~.,data=dom_s_scale_2)
vif(don_m)

cn_m <- lm(DOC_DON~.,data=dom_s_scale_2)
vif(cn_m)

## Now use dom_s_scale_2 for subsequent analyses??

## Remove variables from POM_S
pom_s_scale <- data.frame(pom_s_scale)

poc_m <- lm(POC_mg~.,data=pom_s_scale)
vif(poc_m)

pn_m <- lm(PN_mg~.,data=pom_s_scale)
vif(pn_m)

# Remove N_POM
pom_s_scale_2 <- pom_s_scale[,-c(9,10,12,13)]

poc_m <- lm(POC_mg~.,data=pom_s_scale_2)
vif(poc_m)

pn_m <- lm(PN_mg~.,data=pom_s_scale_2)
vif(pn_m)
