# Script to conduct PCA, RDA, and CoIA on Env, DOM, and POM data sets
# Using corrected fluorescence data
# A Hounshell, 04 Mar 2019
# Added to GitHub Repo 06 Sep 2019

# Load in data (Database.csv)
my_data <- read.csv(file.choose())
# Remove un-complete data rows (any rows that do not have all data associated with them)
my_data2 <- my_data[complete.cases(my_data),]

# Add new categorical variable to describe upper, mid, lower estuary
attach(my_data2)
my_data2$est[Station < 55] <- "Upper"
my_data2$est[Station > 55 & Station < 130] <- "Mid"
my_data2$est[Station > 130] <- "Lower"
detach(my_data2)
my_data2$est <- as.factor(my_data2$est)

# Load in libraries need
pacman::p_load(vegan,adespatial,ade4,PerformanceAnalytics,corrplot,Hmisc)

# Plot DOM and POM parameter by season as box plots
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
my_data$Season<-factor(my_data$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
ylab.text=expression(paste("DOC (mg L"^"-1"*")"))
boxplot(DOC_mg~Season,data=my_data,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,ylim=c(0,15))

my_data$Season<-factor(my_data$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
ylab.text=expression(paste("POC (mg L"^"-1"*")"))
boxplot(POC_mg~Season,data=my_data,varwidth=TRUE,ylab=ylab.text,cex.axis=1.5,cex.lab=1.5,outline=FALSE,ylim=c(0,4))

boxplot(HIX_DOM~Season,data=my_data,varwidth=TRUE,ylab="DOM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25))
abline(h=6,lty=2)
abline(h=16,lty=2)
text(1.7,4,labels="Fresher material",cex=1.5)
text(1.7,23,labels="More humified material",cex=1.5)

boxplot(HIX_POM~Season,data=my_data,varwidth=TRUE,ylab="POM HIX",cex.axis=1.5,cex.lab=1.5,ylim=c(0,25))
abline(h=6,lty=2)
abline(h=16,lty=2)

# Separate data by data pool
env.all <- my_data2[,c(5:9)]
dom.all <- my_data2[,c(13,15:20,22:27)]
pom.all <- my_data2[,c(29,31:36,38:43)]

# Plot Correlation chart - environmental data
par(mar=c(5.1,4.1,4.1,2.1))
chart.Correlation(env.all, histogram=TRUE, method=c("pearson"))

# Plot correlation chart - DOM data
chart.Correlation(dom.all, histogram=TRUE, method=c("pearson"))
res2 <- cor(dom.all,method=c("pearson"))

# Plot correlation chart - POM data
chart.Correlation(pom.all, histogram=TRUE, method=c("pearson"))

# Separate each data matrix: Env, DOM, POM
env <- my_data2[,c(5:9)]
dom <- my_data2[,c(13,16,18,20,22:23)]
pom <- my_data2[,c(29,32:36,39:40)]

# Take the square root of the DOM and POM data - to 'correct' for the right-skewdness
env.sqrt <- sqrt(env)
dom.sqrt <- sqrt(dom)
pom.sqrt <- sqrt(pom)

# Check correlation plots
chart.Correlation(env.sqrt,histogram=TRUE,method=c("pearson"))
chart.Correlation(dom.sqrt,histogram=TRUE,method=c("pearson"))
chart.Correlation(pom.sqrt,histogram=TRUE,method=c("pearson"))

# Combine Env, DOM, and POM data into a single matrix
all_data <- cbind(env,dom.sqrt,pom.sqrt)

# Plot correlation matrix for all data
res3 <- rcorr(as.matrix(all_data),type=c("pearson"))
corrplot(res3$r, type="upper", tl.col = "black", tl.srt=45, p.mat = res3$P, sig.level = 0.05, insig = "blank")
res3$P

## Conduct PCA on env, dom.sqrt, and pom.sqrt data
env.pca <- rda(env,scale=TRUE)
summary(env.pca,axes=0)
plot(env.pca)
text(env.pca)
dom.pca <- rda(dom.sqrt,scale=TRUE)
summary(dom.pca,axes=0)
plot(dom.pca)
text(dom.pca)
pom.pca <- rda(pom.sqrt,scale=TRUE)
summary(pom.pca,axes=0)
plot(pom.pca)
text(pom.pca)
# Construct PCA biplot that will be divided by season and will include objects and variables as well as a circle of equilibrium contribution
# The graph will be in scaling 2
# Extract species scores for scaling 2
envspe.sc2 <- scores(env.pca,display="sp",choices=c(1,2))
domspe.sc2 <- scores(dom.pca, display="sp", choices=c(1,2))
pomspe.sc2 <- scores(pom.pca, display="sp", choices=c(1,2))
# Season ordered as: winter, spring, summer, fall
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env.pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env.pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),pch=c(21,22,23,24),pt.bg=colvec,cex=1.3))
arrows(0, 0, envspe.sc2[,1], envspe.sc2[,2], angle=20, col="black")
text(env.pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")

plot(dom.pca,type="n",scaling=2,xlab="PC1 (48% var. explained)",ylab="PC2 (29% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom.pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, domspe.sc2[,1], domspe.sc2[,2], angle=20, col="black")
text(dom.pca, display = "species", labels=c("","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.5,-0.4,labels="DOC",cex=1.5,col="black")
text(1.2,1.15,labels="DOC:DON",cex=1.5,col="black")
text(2.4,0.6,labels="SUVA",cex=1.5,col="black")
text(-2.1,-1.35,labels="BIX",cex=1.5,col="black")
text(1.1,-2.65,labels="B",cex=1.5,col="black")
text(1.7,-2.35,labels="T",cex=1.5,col="black")

plot(pom.pca,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (30% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom.pca,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0, 0, pomspe.sc2[,1], pomspe.sc2[,2], angle=20, col="black")
text(pom.pca, display = "species", labels=c("","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.3,0.95,labels="POC",cex=1.5,col="black")
text(0.5,-1.5,labels="POC:PN",cex=1.5,col="black")
text(1.9,-2.1,labels="a254",cex=1.5,col="black")
text(-0.65,-2.3,labels="SUVA",cex=1.5,col="black")
text(-1.85,-1.4,labels="HIX",cex=1.5,col="black")
text(1.15,1.55,labels="BIX",cex=1.5,col="black")
text(2.6,0.3,labels="T",cex=1.5,col="black")
text(2.2,-1.6,labels="A",cex=1.5,col="black")

# PCA results plotted as S vs. B
my_data2$Depth<-factor(my_data2$Depth, levels=c("S", "B"))
with(my_data2,levels(Depth))
colvec<-c("#CC2529","#396AB1")
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env.pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env.pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
with(my_data2,legend("bottomleft",legend=levels(Depth),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
arrows(0, 0, envspe.sc2[,1], envspe.sc2[,2], angle=20, col="black")
text(env.pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")
plot(dom.pca,type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC2 (19% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom.pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
arrows(0, 0, domspe.sc2[,1], domspe.sc2[,2], angle=20, col="black")
text(dom.pca, display = "species", labels=c("","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.5,-0.4,labels="DOC",cex=1.5,col="black")
text(1.2,1.05,labels="DOC:DON",cex=1.5,col="black")
text(2.4,0.6,labels="SUVA",cex=1.5,col="black")
text(-2.1,-1.35,labels="BIX",cex=1.5,col="black")
text(1.1,-2.65,labels="B",cex=1.5,col="black")
text(1.7,-2.35,labels="T",cex=1.5,col="black")
plot(pom.pca,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (30% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom.pca,display="sites",col=colvec[Depth],scaling=2,pch=22,bg=colvec[Depth],cex=1.1))
arrows(0, 0, pomspe.sc2[,1], pomspe.sc2[,2], angle=20, col="black")
text(pom.pca, display = "species", labels=c("","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.5,1,labels="POC",cex=1.5,col="black")
text(0.65,-1.5,labels="POC:PN",cex=1.5,col="black")
text(1.9,-2.1,labels="a254",cex=1.5,col="black")
text(-0.7,-2.3,labels="SUVA",cex=1.5,col="black")
text(-2.05,-1.4,labels="HIX",cex=1.5,col="black")
text(1.3,1.5,labels="BIX",cex=1.5,col="black")
text(2.7,0.3,labels="T",cex=1.5,col="black")
text(2.2,-1.6,labels="A",cex=1.5,col="black")

# PCA results plotted as location in the estuary (Upper, Mid, Lower)
my_data2$est<-factor(my_data2$est, levels=c("Upper", "Mid","Lower"))
with(my_data2,levels(est))
colvec<-c("#396AB1","#3E9651","#DA7C30")
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
plot(env.pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env.pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
with(my_data2,legend("bottomleft",legend=levels(est),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
arrows(0, 0, envspe.sc2[,1], envspe.sc2[,2], angle=20, col="black")
text(env.pca, display = "species", labels=c("","","","",""), scaling=2, cex = 0.8, col = "black")
text(-2.2,0.65,labels="Temp",cex=1.5,col="black")
text(-2.35,-0.7,labels="Sal",cex=1.5,col="black")
text(2.3,-1.8,labels="DO",cex=1.5,col="black")
text(2.2,0.6,labels="Turb",cex=1.5,col="black")
text(0.2,-2.7,labels="Chla",cex=1.5,col="black")
plot(dom.pca,type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC2 (19% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom.pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
arrows(0, 0, domspe.sc2[,1], domspe.sc2[,2], angle=20, col="black")
text(dom.pca, display = "species", labels=c("","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.5,-0.4,labels="DOC",cex=1.5,col="black")
text(1.2,1.05,labels="DOC:DON",cex=1.5,col="black")
text(2.4,0.6,labels="SUVA",cex=1.5,col="black")
text(-2.1,-1.35,labels="BIX",cex=1.5,col="black")
text(1.1,-2.65,labels="B",cex=1.5,col="black")
text(1.7,-2.35,labels="T",cex=1.5,col="black")
plot(pom.pca,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (30% var. explained)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom.pca,display="sites",col=colvec[est],scaling=2,pch=22,bg=colvec[est],cex=1.1))
arrows(0, 0, pomspe.sc2[,1], pomspe.sc2[,2], angle=20, col="black")
text(pom.pca, display = "species", labels=c("","","","","","","",""), scaling=2, cex = 0.8, col = "black")
text(2.5,1,labels="POC",cex=1.5,col="black")
text(0.65,-1.5,labels="POC:PN",cex=1.5,col="black")
text(1.9,-2.1,labels="a254",cex=1.5,col="black")
text(-0.7,-2.3,labels="SUVA",cex=1.5,col="black")
text(-2.05,-1.4,labels="HIX",cex=1.5,col="black")
text(1.3,1.5,labels="BIX",cex=1.5,col="black")
text(2.7,0.3,labels="T",cex=1.5,col="black")
text(2.2,-1.6,labels="A",cex=1.5,col="black")

## Compare PCA results to NMDS results (for appendix/SI) - show that the assumptions made by PCA are valid for the data
# Start with environmental data
env.scale <- scale(env)
env.mds <- metaMDS(env.scale,distance="euclidean",k=3,autotransform=FALSE,noshare=FALSE,trymax=40,plot=TRUE)
env.mds

# Use Procrustes analysis to determine how similar the PCA and NMDS results are
env.pro <- procrustes(env.pca,env.mds,symmetric=TRUE)
summary(env.pro)
env.protest <- protest(env.pca,env.mds,symmetric=TRUE)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
colvec<-c("#396AB1","#3E9651","#DA7C30","#CC2529")
plot(env.pca,type="n",scaling=2,xlab="PC1 (49% var. explained)",ylab="PC2 (25% var. explained)",cex.axis=1.5,cex.lab=1.5,xlim=c(-1,2),ylim=c(-2,1))
with(my_data2,points(env.pca,display="sites",col=colvec[Season],scaling=2,pch=22,bg=colvec[Season],cex=1.1))
with(my_data2,legend("bottomright",legend=levels(Season),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
plot(env.mds,cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(env.mds,col=colvec[Season],pch=22,bg=colvec[Season],cex=1.1))
text(-3,-4.5,labels="Stress = 0.06",cex=1.5,col="black")
plot(env.pro,cex.axis=1.5,cex.lab=1.5,main="")
text(0.1,-0.1,labels="Corr. = 0.92",cex=1.5,col="black")
text(0.1,-0.13,labels="p = 0.001",cex=1.5,col="black")

# Now DOM data
dom.scale <- scale(dom.sqrt)
dom.mds <- metaMDS(dom.scale,distance="euclidean",k=3,autotransform=FALSE,noshare=FALSE,trymax=40,plot=TRUE)
dom.mds
# Use Procrustes analysis to determine how similar the PCA and NMDS results are
dom.pro <- procrustes(dom.pca,dom.mds,symmetric=TRUE)
summary(dom.pro)
dom.protest <- protest(dom.pca,dom.mds,symmetric=TRUE)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
colvec<-c("#396AB1","#3E9651","#DA7C30","#CC2529")
plot(dom.pca,type="n",scaling=2,xlab="PC1 (59% var. explained)",ylab="PC2 (19% var. explained)",cex.axis=1.5,cex.lab=1.5,xlim=c(-1,2),ylim=c(-2,1))
with(my_data2,points(dom.pca,display="sites",col=colvec[Season],scaling=2,pch=22,bg=colvec[Season],cex=1.1))
with(my_data2,legend("bottomright",legend=levels(Season),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
plot(dom.mds,cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom.mds,col=colvec[Season],pch=22,bg=colvec[Season],cex=1.1))
text(-2.3,5.8,labels="Stress = 0.05",cex=1.5,col="black")
plot(dom.pro,cex.axis=1.5,cex.lab=1.5,main="")
text(-0.08,-0.1,labels="Corr. = 0.92",cex=1.5,col="black")
text(-0.08,-0.13,labels="p = 0.001",cex=1.5,col="black")

# PCA and NMDS for POM data
pom.scale <- scale(pom.sqrt)
pom.mds <- metaMDS(pom.scale,distance="euclidean",k=2,autotransform=FALSE,noshare=FALSE,trymax=40,plot=TRUE)
pom.mds
# Use Procrustes analysis to determine how similar the PCA and NMDS results are
pom.pro <- procrustes(pom.pca,pom.mds,symmetric=TRUE)
summary(pom.pro)
pom.protest <- protest(pom.pca,pom.mds,symmetric=TRUE)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))
my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
colvec<-c("#396AB1","#3E9651","#DA7C30","#CC2529")
plot(pom.pca,type="n",scaling=2,xlab="PC1 (39% var. explained)",ylab="PC2 (30% var. explained)",cex.axis=1.5,cex.lab=1.5,xlim=c(-1,2),ylim=c(-2,1))
with(my_data2,points(pom.pca,display="sites",col=colvec[Season],scaling=2,pch=22,bg=colvec[Season],cex=1.1))
with(my_data2,legend("bottomright",legend=levels(Season),bty="n",col=colvec,pch=22,pt.bg=colvec,cex=1.1))
plot(pom.mds,cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(pom.mds,col=colvec[Season],pch=22,bg=colvec[Season],cex=1.1))
text(9,-5,labels="Stress = 0.14",cex=1.5,col="black")
plot(pom.pro,cex.axis=1.5,cex.lab=1.5,main="")
text(0.25,-0.14,labels="Corr. = 0.97",cex=1.5,col="black")
text(0.25,-0.18,labels="p = 0.001",cex=1.5,col="black")

## Then conduct RDA using the Env data as the explanatory matrix
dom.scale <- data.frame(dom.scale)
pom.scale <- data.frame(pom.scale)
env.scale <- data.frame(env.scale)
dom.rda <- rda(dom.scale~.,env.scale,scale=FALSE)
# Global adjusted R^2 (0.41)
(R2a.all <- RsquareAdj(dom.rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(dom.rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
# Scale both the Env parameters and the DOM parameters: Sal, Temp, Turb, Chla
dom.rda.forsel <- forward.sel(dom.scale,env.scale,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,adjR2thresh=R2a.all)

# Use ordistep for forward selection
mod0 <- rda(dom.scale ~ 1, data=env.scale,scale=FALSE)
step.both <- ordistep(mod0, scope = formula(dom.rda),direction = "both",permutations = how(nperm = 499))
step.both$anova
RsquareAdj(step.both)
anova(step.both, by = "axis", permutations = how(nperm = 999))

# Then backward elimination using ordistep
step.back <- ordistep(dom.rda,direction="backward",permutations = how(nperm = 499))
RsquareAdj(step.back)
anova(step.back, by = "axis", permutations = how(nperm = 999))

# DOM RDA model: Temp, Sal, Turb, Chla
dom.rda.final <- rda(dom.scale ~ env.scale$Sal + env.scale$Temp + env.scale$Turb + env.scale$Chla,scale=FALSE)

# Plot
domrspe.sc2 <- scores(dom.rda.final, display="sp", choices=c(1,2), scaling=2)
domrbp.sc2 <- scores(dom.rda.final,display="bp",choices=c(1,2),scaling=2)

# Conduct RDA on POM data using Env data as the explanatory matrix
pom.rda <- rda(pom.scale~.,env.scale,scale=FALSE)
# Global adjusted R^2 (0.47)
(R2a.all <- RsquareAdj(pom.rda)$adj.r.squared)
# Test of all canoical axes from full rda
anova(pom.rda,by="axis", permutation=how(nperm=999))

# Use forward selection to select the model with the best number of explanatory variables
pom.rda.forsel <- forward.sel(pom.scale,env.scale,Xscale=FALSE,Yscale=FALSE,Ycenter=FALSE,adjR2thresh=R2a.all)

# Use ordistep for forward selection
mod0 <- rda(pom.scale ~ 1, data=env.scale,scale=FALSE)
step.both <- ordistep(mod0, scope = formula(pom.rda),direction = "both",permutations = how(nperm = 499))
step.both$anova
RsquareAdj(step.both)
anova(step.both, by = "axis", permutations = how(nperm = 999))

# Then backward elimination using ordistep
step.back <- ordistep(pom.rda,direction="backward",permutations = how(nperm = 499))
RsquareAdj(step.back)
anova(step.back, by = "axis", permutations = how(nperm = 999))

# Conduct and plot final POM rda
pom.rda.final <- rda(pom.scale ~ env.scale$Chla + env.scale$Turb + env.scale$Sal + env.scale$Temp, scale=FALSE)
summary(pom.rda.final,axes=0)
R2a.all <- RsquareAdj(pom.rda.final)$adj.r.squared
# Plot
pomrspe.sc2 <- scores(pom.rda.final, display="sp", choices=c(1,2), scaling=2)
pomrbp.sc2 <- scores(pom.rda.final,display="bp",choices=c(1,2),scaling=2)

# Plot POM and DOM RDA results on one graph
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,2))

my_data2$Season<-factor(my_data2$Season, levels=c("Winter", "Spring", "Summer", "Fall"))
with(my_data2,levels(Season))
colvec<-c("#1E88E5","#004D40","#FFC107","#D81B60")
with(my_data2,levels(Season))
sq<-c(21,22,23,24)

plot(dom.rda.final,scaling=2,display="sites",xlab="RDA1 (82% fitted,33% total var.)",ylab="RDA2 (13% fitted, 5% total var.)",cex.axis=1.5,cex.lab=1.5)
with(my_data2,points(dom.rda.final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
with(my_data2,legend("bottomleft",legend=levels(Season),bty="n",col=c("black","black","black","black"),pch=c(21,22,23,24),pt.bg=colvec,cex=1.5))
arrows(0,0,domrspe.sc2[,1], domrspe.sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(dom.rda.final,display = "species", labels=c("","DOC:DON","","BIX","",""), scaling=2, cex = 1.5, col = "black")
text(-2.7,0.1,labels="DOC",cex=1.5,col="black")
text(-2.5,-0.4,labels="SUVA",cex=1.5,col="black")
text(-0.24,-1,labels="B",cex=1.5,col="white")
text(-0.7,-1.2,labels="T",cex=1.5,col="white")
text(dom.rda.final,display="bp",labels=c("Sal","Temp","Turb",""),scaling=2,cex=1.5,col="black")
text(0.2,-0.6,labels="Chla",cex=1.5,col="white")

plot(pom.rda.final,scaling=2,display="sites",xlab="RDA1 (57% fitted, 27% total var.)",ylab="RDA2 (41% fitted, 19% total var.)",cex.axis=1.5,cex.lab=1.5,xlim=c(-2,3),ylim=c(-2,3))
with(my_data2,points(pom.rda.final,display="sites",col=c("black","black","black","black"),scaling=2,pch=sq[Season],bg=colvec[Season],cex=1.5))
arrows(0,0,pomrspe.sc2[,1], pomrspe.sc2[,2],lty="dashed",col="black",adj=0.5,length=0)
text(pom.rda.final,display = "species", labels=c("POC","POC:PN","a254","","HIX","BIX","T","A"), scaling=2, cex = 1.5, col = "black")
text(-0.8,2.2,labels="SUVA",cex=1.5,col="black")
text(pom.rda.final,display="bp",labels=c("Chla","","Sal","Temp"),scaling=2,cex=1.5,col="black")
text(0.2,2,labels="Turb",cex=1.5,col="black")

## Now conduct CoIA on DOMvEnv; POMvEnv; and DOMvPOM
# Data sets will be scaled during PCA prior to CoIA
# PCA on all matrices using ade4 functions
dudi.env <- dudi.pca(as.matrix(env[,c(1:2,4:5)]),scale=TRUE,center=TRUE,scannf=FALSE)
dudi.dom <- dudi.pca(dom.sqrt,scale = TRUE,center=TRUE,scannf=FALSE)
dudi.pom <- dudi.pca(pom.sqrt,scale = TRUE,center=TRUE,scannf=FALSE)

# Conduct CoIA on DOMvEnv
# Cumulated relative variation of eigenvalues
cumsum(dudi.env$eig / sum(dudi.env$eig))
# Cumulated relative variation of eigenvalues
cumsum(dudi.dom$eig / sum(dudi.dom$eig))
# Are the row weights equal in the 2 analyses?
all.equal(dudi.env$lw, dudi.dom$lw)

coia.env.dom <- coinertia(dudi.env, dudi.dom,scannf = FALSE,nf = 2)
summary(coia.env.dom)
# Relative variation on first eigenvalue
coia.env.dom$eig[1] / sum(coia.env.dom$eig)
# Permutation test
randtest(coia.env.dom, nrepet = 999)

# Plot results
plot(coia.env.dom)

par(mfrow=c(2,2))
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.dom$aX[,1],coia.env.dom$aX[,2],col="black",adj=0.5,length=0)
text(coia.env.dom$aX,labels=c("Ax1","Ax2"),cex=1.5,col="black")
text(-0.1,-0.95,labels="Unconstrained axes (Env)",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.dom$aY[,1],coia.env.dom$aY[,2],col="black",adj=0.5,length=0)
text(coia.env.dom$aY,labels=c("Ax1","Ax2"),cex=1.5,col="black")
text(-0.1,0.95,labels="Unconstrained axes (DOM)",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.dom$l1[,1], coia.env.dom$l1[,2],col="black",adj=0.5,length=0)
text(coia.env.dom$l1, labels=c("DOC","","SUVA","BIX","B","T"),cex = 1.5, col = "black")
text(0.2,-0.4,labels="DOC:DON",cex=1.5,col="black")
arrows(0,0,coia.env.dom$c1[,1], coia.env.dom$c1[,2],lty="dashed",col="black",adj=0.5,length=0)
text(coia.env.dom$c1,labels=c("Temp","Sal","Turb","Chla"),cex=1.5,col="black")
text(0.7,0.95,labels="Loadings",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(coia.env.dom$lX[,1],coia.env.dom$lX[,2],cex.axis=1.3,xlab="Factor 1",ylab="Factor 2",cex.lab=1.5,xlim=c(-4,4),ylim=c(-4,7))
points(coia.env.dom$lY[,1],coia.env.dom$lY[,2],pch=19)
arrows(coia.env.dom$lX[,1],coia.env.dom$lX[,2],coia.env.dom$lY[,1],coia.env.dom$lY[,2],length=0.1)
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
text(-3,6.5,labels="RV = 0.43",cex=1.3,col="black")
text(-3,5.5,labels="p = 0.001",cex=1.3,col="black")

## Conduct CoIA on Env vs. POM data
# Cumulated relative variation of eigenvalues
cumsum(dudi.env$eig / sum(dudi.env$eig))
# Cumulated relative variation of eigenvalues
cumsum(dudi.pom$eig / sum(dudi.pom$eig))
# Are the row weights equal in the 2 analyses?
all.equal(dudi.env$lw, dudi.pom$lw)

coia.env.pom <- coinertia(dudi.env, dudi.pom,scannf = FALSE,nf = 2)
summary(coia.env.pom)
# Relative variation on first eigenvalue
coia.env.pom$eig[1] / sum(coia.env.pom$eig)
# Permutation test
randtest(coia.env.pom, nrepet = 999)

# Plot results
plot(coia.env.pom)

par(mfrow=c(2,2))
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.pom$aX[,1],coia.env.pom$aX[,2],col="black",adj=0.5,length=0)
text(coia.env.pom$aX,labels=c("",""),cex=1.5,col="black")
text(0.95,-0.25,labels="Ax1",cex=1.5,col="black")
text(0.3,0.95,labels="Ax2",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
text(-0.3,-0.95,labels="Unconstrained axes (Env)",cex=1.5,col="black")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.pom$aY[,1],coia.env.pom$aY[,2],col="black",adj=0.5,length=0)
text(coia.env.pom$aY,labels=c("",""),cex=1.5,col="black")
text(-0.45,-0.95,labels="Ax1",cex=1.5,col="black")
text(0.95,-0.4,labels="Ax2",cex=1.5,col="black")
text(-0.3,0.95,labels="Unconstrained axes (POM)",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.env.pom$l1[,1], coia.env.pom$l1[,2],col="black",adj=0.5,length=0)
text(coia.env.pom$l1, labels=c("","","","","","","",""),cex = 1.5, col = "black")
text(-0.5,-0.5,labels="POC",cex=1.5,col="black")
text(0.3,-0.2,labels="POC:PN",cex=1.5,col="black")
text(0.5,-0.5,labels="a254",cex=1.5,col="black")
text(0.85,-0.1,labels="SUVA",cex=1.5,col="black")
text(0.4,0.35,labels="HIX",cex=1.5,col="black")
text(-0.4,-0.05,labels="BIX",cex=1.5,col="black")
text(-0.1,-0.5,labels="T",cex=1.5,col="black")
text(0.15,-0.6,labels="A",cex=1.5,col="black")
arrows(0,0,coia.env.pom$c1[,1], coia.env.pom$c1[,2],lty="dashed",col="black",adj=0.5,length=0)
text(coia.env.pom$c1,labels=c("","","",""),cex=1.5,col="black")
text(-0.4,0.4,labels="Temp",cex=1.5,col="black")
text(-0.68,0.19,labels="Sal",cex=1.5,col="black")
text(0.8,-0.3,labels="Turb",cex=1.5,col="black")
text(-0.5,-0.95,labels="Chla",cex=1.5,col="black")
text(-0.7,0.95,labels="Loadings",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(coia.env.pom$lX[,1],coia.env.pom$lX[,2],xlim=c(-5,5),ylim=c(-12,5),cex.axis=1.3,xlab="Factor 1",ylab="Factor 2",cex.lab=1.5)
points(coia.env.pom$lY[,1],coia.env.pom$lY[,2],pch=19)
arrows(coia.env.pom$lX[,1],coia.env.pom$lX[,2],coia.env.pom$lY[,1],coia.env.pom$lY[,2],length=0.1)
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
text(-3.5,4.5,labels="RV = 0.53",cex=1.3,col="black")
text(-3.5,3,labels="p = 0.001",cex=1.3,col="black")

## Then conduct CoIA on DOM vs. POM data - how well are these two pools aligned?
# Cumulated relative variation of eigenvalues
cumsum(dudi.dom$eig / sum(dudi.dom$eig))
# Cumulated relative variation of eigenvalues
cumsum(dudi.pom$eig / sum(dudi.pom$eig))
# Are the row weights equal in the 2 analyses?
all.equal(dudi.dom$lw, dudi.pom$lw)

coia.dom.pom <- coinertia(dudi.dom, dudi.pom,scannf = FALSE,nf = 2)
summary(coia.dom.pom)
# Relative variation on first eigenvalue
coia.dom.pom$eig[1] / sum(coia.dom.pom$eig)
# Permutation test
randtest(coia.dom.pom, nrepet = 999)

# Plot results
plot(coia.dom.pom)

par(mfrow=c(2,2))
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.dom.pom$aX[,1],coia.dom.pom$aX[,2],col="black",adj=0.5,length=0)
text(coia.dom.pom$aX,labels=c("Ax1","Ax2"),cex=1.5,col="black")
text(-0.1,0.9,labels="Unconstrained axes (DOM)",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.dom.pom$aY[,1],coia.dom.pom$aY[,2],col="black",adj=0.5,length=0)
text(coia.dom.pom$aY,labels=c("Ax1",""),cex=1.5,col="black")
text(-0.9,-0.1,labels="Ax2",cex=1.5,col="black")
text(-0.1,0.9,labels="Unconstrained axes (POM)",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")

plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.dom.pom$c1[,1], coia.dom.pom$c1[,2],col="black",adj=0.5,length=0)
text(coia.dom.pom$c1,labels=c("","","","","",""),cex=1.5,col="black")
text(-0.5,0.2,labels="DOC",cex=1.5,col="black")
text(-0.3,0.9,labels="DOC:DON",cex=1.5,col="black")
text(-0.6,0,labels="SUVA",cex=1.5,col="black")
text(0.65,0.55,labels="BIX",cex=1.5,col="black")
text(0,-0.2,labels="B",cex=1.5,col="black")
text(-0.2,-0.35,labels="T",cex=1.5,col="black")
#text(-0.25,-0.3,labels="T",cex=1.5,col="black")
#text(0.5,-0.95,labels="DOM Loadings",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
plot(0,0,ylab="",xlab="",cex.axis=1.3)
arrows(0,0,coia.dom.pom$l1[,1], coia.dom.pom$l1[,2],col="black",adj=0.5,length=0)
text(coia.dom.pom$l1, labels=c("","","","","","","",""),cex = 1.5, col = "black")
text(0.3,0.6,labels="POC",cex=1.5,col="black")
text(-0.11,-0.17,labels="POC:PN",cex=1.5,col="black")
text(-0.6,0.4,labels="a254",cex=1.5,col="black")
text(-0.85,-0.15,labels="SUVA",cex=1.5,col="black")
text(-0.2,-0.45,labels="HIX",cex=1.5,col="black")
text(0.4,0.15,labels="BIX",cex=1.5,col="black")
text(-0.05,0.55,labels="T",cex=1.5,col="black")
text(-0.35,0.55,labels="A",cex=1.5,col="black")
#text(-0.5,0.95,labels="POM Loadings",cex=1.5,col="black")
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")

par(mfrow=c(1,1))
plot(coia.dom.pom$lX[,1],coia.dom.pom$lX[,2],xlim=c(-6,5),ylim=c(-5,10),cex.axis=1.3,xlab="Factor 1",ylab="Factor 2",cex.lab=1.5)
points(coia.dom.pom$lY[,1],coia.dom.pom$lY[,2],pch=19)
arrows(coia.dom.pom$lX[,1],coia.dom.pom$lX[,2],coia.dom.pom$lY[,1],coia.dom.pom$lY[,2],length=0.1)
abline(h=0,col="black",lty="dotted")
abline(v=0,col="black",lty="dotted")
text(-5,9.5,labels="RV = 0.20",cex=1.3,col="black")
text(-5,8,labels="p = 0.001",cex=1.3,col="black")

# Saved as MultiStats_corr