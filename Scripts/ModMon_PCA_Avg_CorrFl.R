# Script to plot PCA for averaged time points and for all time points
# NRE Monitoring project
# Includes corrected fluorescence data
# A Hounshell, 04 Mar 2019
# Added to Github 06 Sep 2019

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
env <- avg_data2[,c(4,6:7,10:12)]

# Look at DOM data: 
res2 <- cor(as.matrix(avg_data2[,c(15,17,18:22,24:29)]),method=c("pearson"))
corrplot(res2,type="upper",order="hclust",tl.col="black",tl.srt=45)
chart.Correlation(avg_data2[,c(15,17,18:22,24:29)],histogram=TRUE,pch=19,method=c("pearson"))
dom <- avg_data2[,c(15,17:18,20,25)]

# Look at POM data:
res3 <- cor(as.matrix(avg_data2[,c(31,33:38,40:45)]),method=c("pearson"))
corrplot(res3,type="upper",order="hclust",tl.col="black",tl.srt=45)
chart.Correlation(avg_data2[,c(31,33:38,40:45)],histogram=TRUE,pch=19,method=c("pearson"))
pom <- avg_data2[,c(31,34,36:38,41:42)]

# Look at correlations between all parameters
res4 <- rcorr(as.matrix(avg_data2[,c(4,6:7,10:12,15,17:18,20,25,31,34,36:38,41:42)]),type=c("spearman"))
corrplot(res4$r, type="upper", tl.col = "black", tl.srt=45, p.mat = res4$P, sig.level = 0.05, insig = "blank")

# Conduct PCA on all parameters (i.e., don't separate by designation?)
pcadat <- avg_data2[,c(4,6:7,10:12,15,17:18,20,25,31,34,36:38,41:42)]
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
colvec<-c("#396AB1","#3E9651","#DA7C30","#CC2529")

par(mfrow=c(1,2))
plot(dat.pca,type="n",scaling=2,xlab="PC1 (41.1% var. explained)",ylab="PC2 (17.4% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-2,1.5),xlim=c(-1.5,3))
with(avg_data2,points(dat.pca,display="sites",col=colvec[Season],scaling=2,pch=22,bg=colvec[Season],cex=1.1))
with(avg_data2,legend("bottomrigh",legend=levels(Season),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))

plot (spe.sc2[,1],spe.sc2[,2],pch=20,cex=2,xlab="PC1(41.1% var. explained)",ylab="PC2 (17.4% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1),xlim=c(-1,1.5))
abline(0,0,lty=3)
abline(v=0,lty=3)
text(dat.pca, display="species", labels = c("","","","","","","","","","","","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(-0.9,0.35,labels="Temp",cex=1.2,col="black")
text(-0.95,-0.3,labels="Sal",cex=1.2,col="black")
text(0.7,0.3,labels="Turb",cex=1.2,col="black")
text(-0.4,-0.8,labels="Chla",cex=1.2,col="black")
text(1.05,-0.15,labels="Avg_Wind",cex=1.2,col="black")
text(0.2,0.85,labels="Wind_Dir",cex=1.2,col="black")
text(0.6,-0.6,labels="DOC",cex=1.2,col="black")
text(0.5,0.1,labels="DOC:DON",cex=1.2,col="black")
text(1.25,0,labels="SUVA_D",cex=1.2,col="black")
text(-0.8,-0.7,labels="B_DOM",cex=1.2,col="black")
text(0.2,-1,labels="N_DOM",cex=1.2,col="black")
text(-0.33,-0.4,labels="POC",cex=1.2,col="black")
text(0.23,-0.15,labels="POC:PN",cex=1.2,col="black")
text(1.25,0.28,labels="SUVA_P",cex=1.2,col="black")
text(0.2,-0.4,labels="HIX_P",cex=1.2,col="black")
text(-0.55,0.65,labels="BIX_P",cex=1.2,col="black")
text(-0.1,-0.05,labels="T_POM",cex=1.2,col="black")
text(1,-0.35,labels="A_POM",cex=1.2,col="black")

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
pcawind <- avg_data2[,c(3:4,6:7,10,15,18,20,24,29,31,34,36:38,41:42)]
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
colvec<-c("#396AB1","#3E9651","#DA7C30","#CC2529","#59656D")

par(mfrow=c(1,2))
plot(wind.pca,type="n",scaling=2,xlab="PC1 (43.8% var. explained)",ylab="PC2 (16.0% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-2,1.5),xlim=c(-1.5,3))
with(avg_data2,points(wind.pca,display="sites",col=colvec[card],scaling=2,pch=22,bg=colvec[card],cex=1.1))
with(avg_data2,legend("bottomrigh",legend=levels(card),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))

plot (spe.sc2[,1],spe.sc2[,2],pch=20,cex=2,xlab="PC1(43.8% var. explained)",ylab="PC2 (16.0% var. explained)",cex.axis=1.5,cex.lab=1.5,ylim=c(-1.5,1),xlim=c(-1,1.5))
abline(0,0,lty=3)
abline(v=0,lty=3)
text(wind.pca, display="species", labels = c("","","","","","","","","","","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(-0.9,0.35,labels="Temp",cex=1.2,col="black")
text(-0.95,-0.4,labels="Sal",cex=1.2,col="black")
text(0.67,0.3,labels="Turb",cex=1.2,col="black")
text(-0.3,-0.8,labels="Chla",cex=1.2,col="black")
text(0.3,-0.05,labels="Avg_Wind",cex=1.2,col="black")
text(1,-0.35,labels="DOC",cex=1.2,col="black")
text(1.2,-0.07,labels="DOC:DON",cex=1.2,col="black")
text(1.27,0.08,labels="SUVA_D",cex=1.2,col="black")
text(-0.75,-0.75,labels="B_DOM",cex=1.2,col="black")
text(0.4,-0.85,labels="N_DOM",cex=1.2,col="black")
text(-0.35,-0.4,labels="POC",cex=1.2,col="black")
text(0.4,-0.25,labels="POC:PN",cex=1.2,col="black")
text(1.25,0.28,labels="SUVA_P",cex=1.2,col="black")
text(0.23,-0.52,labels="HIX_P",cex=1.2,col="black")
text(-0.55,0.7,labels="BIX_P",cex=1.2,col="black")
text(0.15,0.2,labels="T_POM",cex=1.2,col="black")
text(1,-0.6,labels="A_POM",cex=1.2,col="black")

# Saved as Wind_PCA