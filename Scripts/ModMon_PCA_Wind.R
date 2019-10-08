# Script to plot PCA for averaged time points and for all time points
# NRE Monitoring project
# A Hounshell, 28 Feb 2019
# Add to GitHub Repo: 08 Oct 2019
# Updated to include PCA using averaged data that excludes the wind data

# Load in Averaged data: Database_Avg.csv
avg_data <- read.csv(file.choose())
# Remove incomplete cases
avg_data2 <- avg_data[complete.cases(avg_data),]

# Load in libraries
library(vegan)
library(corrplot)
library(PerformanceAnalytics)
library(Hmisc)

## Look at correlations with the averaged dataset
# Start with Env data: Temp, Sal, DO, Turb, Chla, Flushing_Time, Strat_index, Avg_wind, Max_Wind, Wind_Dir
res <- cor(as.matrix(avg_data2[,c(3:7,10:13)]),method=c("pearson"))
corrplot(res,type="upper",order="hclust",tl.col="black",tl.srt=45)
chart.Correlation(avg_data2[,c(3:7,10:13)],histogram=TRUE,pch=19,method=c("pearson"))
env <- avg_data2[,c(4,6:7,10:13)]

# Look at DOM data: 
res2 <- cor(as.matrix(avg_data2[,c(15,17,18:22,24:29)]),method=c("pearson"))
corrplot(res2,type="upper",order="hclust",tl.col="black",tl.srt=45)
chart.Correlation(avg_data2[,c(15,17,18:22,24:29)],histogram=TRUE,pch=19,method=c("pearson"))
dom <- avg_data2[,c(15,17,20,22,25)]

# Look at POM data:
res3 <- cor(as.matrix(avg_data2[,c(31,33:38,40:45)]),method=c("pearson"))
corrplot(res3,type="upper",order="hclust",tl.col="black",tl.srt=45)
chart.Correlation(avg_data2[,c(31,33:38,40:45)],histogram=TRUE,pch=19)
pom <- avg_data2[,c(31,34,36:38,41:42)]

# Look at correlations between all parameters
res4 <- rcorr(as.matrix(avg_data2[,c(4,6:7,10:13,15,17,20,22,25,31,34,36:38,41:42)]),type=c("pearson"))
corrplot(res4$r, type="upper", tl.col = "black", tl.srt=45, p.mat = res4$P, sig.level = 0.05, insig = "blank")
res4$P

# Conduct PCA on all parameters (i.e., don't separate by designation?)
pcadat <- avg_data2[,c(4,6:7,10:13,15,17,20,22,25,31,34,36:38,41:42)]
dat.pca <- rda(pcadat,scale=TRUE)
# Look at summary of results
summary(dat.pca,axes=0)
# Plot using vegan biplots
par(mfrow = c(1, 2))
biplot(dat.pca, scaling = 1, main = "PCA - scaling 1")
biplot(dat.pca, main = "PCA - scaling 2") # Default scaling 2

# Construct PCA biplot that
# The graph will be in scaling 2
# Extract species scores for scaling 2
spe.sc2 <- scores(dat.pca, display="sp", choices=c(1,2))

avg_data2$Season<-factor(avg_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(avg_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

par(mfrow=c(1,2))
plot(dat.pca,type="n",scaling=2,xlab="PC1 (42.9% var. explained)",ylab="PC2 (16.8% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(avg_data2,points(dat.pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.1))
with(avg_data2,legend("topright",legend=levels(Season),bty="n",col="black",pch=c(21,22,23,24),pt.bg=colvec,cex=1.1))

plot (spe.sc2[,1],spe.sc2[,2],pch=20,cex=2,xlab="PC1(42.9% var. explained)",ylab="PC2 (16.8% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1),xlim=c(-1,1.8))
abline(0,0,lty=3)
abline(v=0,lty=3)
text(dat.pca, display="species", labels = c("","","","","","","","","","","","","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(-0.95,0.35,labels="Sal",cex=1.2,col="black")
text(0.85,-0.45,labels="Turb",cex=1.2,col="black")
text(-0.35,0.7,labels="Chla",cex=1.2,col="black")
text(0.15,-0.35,labels="Avg_Wind",cex=1.2,col="black")
text(1.4,-0.13,labels="Max_Wind",cex=1.2,col="black")
text(0.2,-0.85,labels="Wind_Dir",cex=1.2,col="black")
text(1.1,0.02,labels="Q",cex=1.2,col="black")
text(1.03,0.6,labels="DOC",cex=1.2,col="black")
text(0.6,0.85,labels="DON",cex=1.2,col="black")
text(1.2,0.15,labels="SUVA_D",cex=1.2,col="black")
text(-0.78,-0.03,labels="BIX_D",cex=1.2,col="black")
text(0.01,1,labels="T_DOM",cex=1.2,col="black")
text(-0.5,0.2,labels="POC",cex=1.2,col="black")
text(0.19,0.12,labels="POC:PN",cex=1.2,col="black")
text(1.4,-0.33,labels="SUVA_P",cex=1.2,col="black")
text(0.12,0.35,labels="HIX_P",cex=1.2,col="black")
text(-0.55,-0.5,labels="BIX_P",cex=1.2,col="black")
text(-0.05,-0.17,labels="T_POM",cex=1.2,col="black")
text(1.15,0.35,labels="A_POM",cex=1.2,col="black")

## Conduct PCA on averaged variables w/o wind data
nowind.pca <- rda(pcadat[,c(-4,-5,-6)],scale=TRUE)
# Look at summary of results
summary(nowind.pca,axes=0)
# Plot using vegan biplots
par(mfrow = c(1, 2))
biplot(nowind.pca, scaling = 1, main = "PCA - scaling 1")
biplot(nowind.pca, main = "PCA - scaling 2") # Default scaling 2

# Construct PCA biplot that
# The graph will be in scaling 2
# Extract species scores for scaling 2
spe.sc2 <- scores(nowind.pca, display="sp", choices=c(1,2))

avg_data2$Season<-factor(avg_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(avg_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

par(mfrow=c(1,2))
plot(nowind.pca,type="n",scaling=2,xlab="PC1 (45.1% var. explained)",ylab="PC2 (17.2% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(avg_data2,points(nowind.pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.1))
with(avg_data2,legend("topright",legend=levels(Season),bty="n",col="black",pch=c(21,22,23,24),pt.bg=colvec,cex=1.1))

plot (spe.sc2[,1],spe.sc2[,2],pch=20,cex=2,xlab="PC1(45.1% var. explained)",ylab="PC2 (17.2% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1),xlim=c(-1,1.8))
abline(0,0,lty=3)
abline(v=0,lty=3)
text(nowind.pca, display="species", labels = c("Sal","Turb","Chla","Q","DOC","DON","SUVA_D","BIX_D","T_DOM","POC","POC:PN","SUVA_P","HIX_P","BIX_P","T_POM","A_POM"), scaling=2, cex = 0.8, col = "black")

# Make wind direction a categorical variable
attach(avg_data2)
avg_data2$card[Wind_Dir > 0 & Wind_Dir < 22.5] <- "N"
avg_data2$card[Wind_Dir > 22.5 & Wind_Dir < 67.5] <- "NE"
avg_data2$card[Wind_Dir > 67.5 & Wind_Dir < 112.5] <- "E"
avg_data2$card[Wind_Dir > 112.5 & Wind_Dir < 157.5] <- "SE"
avg_data2$card[Wind_Dir > 157.5 & Wind_Dir < 202.5] <- "S"
avg_data2$card[Wind_Dir > 202.5 & Wind_Dir < 247.5] <- "SW"
avg_data2$card[Wind_Dir > 247.5 & Wind_Dir < 292.5] <- "W"
avg_data2$card[Wind_Dir > 292.5 & Wind_Dir < 337.5] <- "NW"
avg_data2$card[Wind_Dir > 337.5] <- "N"
detach(avg_data2)

# Concuct PCA w/o wind direction (then plot samples by wind direction)
pcawind <- avg_data2[,c(4,6:7,10:11,13,15,17,20,22,25,31,34,36:38,41:42)]
wind.pca <- rda(pcawind,scale=TRUE)
# Look at summary of results
summary(wind.pca,axes=0)
# Plot using vegan biplots
par(mfrow = c(1, 2))
biplot(wind.pca, scaling = 1, main = "PCA - scaling 1")
biplot(wind.pca, main = "PCA - scaling 2") # Default scaling 2

# Construct PCA biplot that
# The graph will be in scaling 2
# Extract species scores for scaling 2
spe.sc2 <- scores(wind.pca, display="sp", choices=c(1,2))

avg_data2$card<-factor(avg_data2$card, levels=c("N", "NE","S","SW","NW"))
with(avg_data2,levels(card))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60","#9C5113")
with(avg_data2,levels(card))
sq<-c(21,22,23,24,25)

par(mfrow=c(1,2))
plot(wind.pca,type="n",scaling=2,xlab="PC1 (45% var. explained)",ylab="PC2 (15.7% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(avg_data2,points(wind.pca,display="sites",col=c("black","black","black","black","black"),scaling=2,pch=sq[card],bg=colvec[card],cex=1.1))
with(avg_data2,legend("topright",legend=levels(card),bty="n",col=c("black","black","black","black","black"),pch=c(21,22,23,24,25),pt.bg=colvec,cex=1.1))

plot (spe.sc2[,1],spe.sc2[,2],pch=20,cex=2,xlab="PC1(45% var. explained)",ylab="PC2 (15.7% var. explained)",cex.axis=1.5,cex.lab=1.5,xlim=c(-1,1.9))
abline(0,0,lty=3)
abline(v=0,lty=3)
text(wind.pca, display="species", labels = c("","","","","","","","","","","","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(-0.9,0.38,labels="Sal",cex=1.2,col="black")
text(1.1,-0.4,labels="Turb",cex=1.2,col="black")
text(-0.35,0.5,labels="Chla",cex=1.2,col="black")
text(0.15,-0.11,labels="Avg_Wind",cex=1.2,col="black")
text(1.45,-0.20,labels="Max_Wind",cex=1.2,col="black")
text(1.1,0.1,labels="Q",cex=1.2,col="black")
text(1.07,0.6,labels="DOC",cex=1.2,col="black")
text(0.65,0.83,labels="DON",cex=1.2,col="black")
text(1.4,0,labels="SUVA_D",cex=1.2,col="black")
text(-0.78,0.05,labels="BIX_D",cex=1.2,col="black")
text(-0.45,0.93,labels="T_DOM",cex=1.2,col="black")
text(-0.37,0.25,labels="POC",cex=1.2,col="black")
text(0.15,0,labels="POC:PN",cex=1.2,col="black")
text(1.43,-0.32,labels="SUVA_P",cex=1.2,col="black")
text(0.42,0.43,labels="HIX_P",cex=1.2,col="black")
text(-0.55,-0.32,labels="BIX_P",cex=1.2,col="black")
text(0.3,-0.3,labels="T_POM",cex=1.2,col="black")
text(0.9,0.33,labels="A_POM",cex=1.2,col="black")

# Saved as 'ModMon_PCA_2'