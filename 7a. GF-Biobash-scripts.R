#sdm: a reproducible and extensible R package for species distribution modelling: Tutorial
#### R
#Clean up!
rm(list=ls())

library("devtools")
library("sdm")
library("dismo")
library("dplyr")
library("tidyr")
library("mapview")
library("raster")
library("rgdal")
library("gradientForest")

#Species for tutorial 
Meehania montis-koyae Ohwi

#Load the 7 points occurrence from sampling sites using dismo
sp <- paste0("C:/Users/Bashir/Desktop/Meehania/Population results/Gradient Forest/7 Pop Coordinates in decimals.csv")
sp

montis <- read.csv(sp)
# or do: montis <- read.table(sp,  header=TRUE,  sep=",")
# first rows
head(montis)


#-------
#Convert the spg dataframe to a special point dataframe
class(montis)

coordinates(montis) <- c("lon","lat")

#Check the number of rows to confirm that the occurrence has reduced
nrow(montis)

#To confirm the new special point dataframe
class(montis)

plot(montis)
#---------
##Download the bioclim data
#Valid resolutions are 0.5, 2.5, 5, and 10 (minutes of a degree). 

?raster
raster::getD


bio <- raster::getData("worldclim", var = "bio", res=2.5)
bio #information about the downloaded bioclim

# lets also get the elevational data associated with the climate data
bioalt <- raster::getData("worldclim", var = "alt", res=2.5)
bioalt
plot(bioalt)
points(montis)

names(bio) #check the names of the climate data which should be bio1-19

#Plot one of them e.g. World temperature
plot(bio[[1]])

#To add species occurrence data on the temperature
points(montis)

#If you have knowledge about the species, then you can remove or subset the data to cover a particular region and fit your data to an area
#We will use a rectangle using Extent in the raster package

e <- drawExtent()

#Drag and select the area, then use crop fxn to crop it out
montis <- crop(montis, e)

plot(montis)
#To check the cropped part using another colour (red)
points(montis, col="red")

#To do the same for the raster bio
bioc <- crop(bio, e)

#To view the subset bio raster
plot(bioc[[1]])

#-------

#To extract the workflow the species distribution model, we need to test the climatic variables for multicollinearity (Pearson correlation, vif or hybrid approach (Pearson and vif))

#We can use usdm package functions vifstep, vifcor
library("usdm")
#vifstep = rule of thumb, if the vif is <10, then there is multicollinearity
#vifcor = checks the coefficient if <0.9, then you'll decide which one to exclude and can exclude the one with high vif

#To check the vif of the raster bioc
vif(bioc)

#To introduce a dataframe as the input and extract the raster (bioc) using the species occurrence points
pointexbio <- raster::extract(bioc, montis)
head(pointexbio)
nrow(pointexbio)

plot(pointexbio)

#Write the extracted point-bioclimate variables to a csv file.
write.csv(x = pointexbio, file = "pointexbio.csv", row.names = TRUE)

#Extract point elevation data form the global data
pointexbioalt <- raster::extract(bioalt, montis)
head(pointexbioalt)
nrow(pointexbioalt)

#Write the extracted point-bioclimate variables to a csv file.
write.csv(x = pointexbioalt, file = "pointexbioalt.csv", row.names = TRUE)

#To run the vifstep and vifcor
v <- vifstep(pointexbio)
#Run v to check the list of remaining variables
v

library("stats")
cor(x, y = NULL, use = "everything",
    method = c("pearson", "kendall", "spearman"))
corpoint <- cor(pointexbio, y = NULL, use = "everything",
    method = c("pearson"))

corpoint
plot(corpoint)

symnum(corpoint)

library("corrplot")

corrplot(corpoint, method = 'number')

#To remove the unwanted variables
bioc <- exclude(bioc, v)
bioc

#### install.packages(c("psych","vegan"), dependencies=TRUE) ####
install.packages(c("psych","vegan"), dependencies=TRUE)
# Load packages
# -------------
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RD

env <- read.csv("C:/Users/ENDALE/Desktop/Gradient Forest/envcorrcoef0.8.csv")
str(env) # Look at the structure of the data frame

pairs.panels(env[,1:20], scale=T)
pairs.panels(env[,2:21], scale=T)

v <- vifstep(env)
v
#Variables      VIF bio16 1.355769, bio17 1.790235, alt 1.543053
install.packages("GGally")
library("GGally")

ggcorr(env,
       nbreaks = 6,
       label = TRUE,
       label_size = 3,
       color = "grey50")

#------
#Convert the .ped file to acceptable format for GF using LEA package
install.packages("gradientForest", repos="http://R-Forge.R-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LEA")

library("BiocManager")
library(LEA)
library(readr)
library(gradientForest)

ped2lfmm(input.file = "out.ped", output.file = "outGFinput.ad.lfmm")

alle<- read.lfmm("outGFinput.ad.lfmm")

######Plot gradientForest result##########

env <- read.csv("envcorrcoef0.8.csv", header = T)
View(env)

env.new <- data.frame(env$bio3, env$bio5, env$bio7, env$bio8, env$bio9, env$bio17, env$bio18, env$alt)

ped2lfmm(input.file = "snp2_nomiss_noLD_ad.ped", output.file = "snp2_nomiss.nold.ad.lfmm")
#ped2lfmm(input.file = "snp2_nomiss_or.ped", output.file = "snp2_nomiss_or.lfmm")
#ped2lfmm(input.file = "snp2_nomiss_ad.ped", output.file = "snp2_nomiss_ad.lfmm")
alle<- read.lfmm("snp2_nomiss.nold.ad.lfmm")
#alle_or<- read.lfmm("snp2_nomiss_or.lfmm")
#alle_ad<- read.lfmm("snp2_nomiss_ad.lfmm")

nSites <- dim(alle)[1]
lev <- floor(log2(nSites * 0.368/2))
lev

gf<- gradientForest(data = cbind(env.new, alle), predictor.vars = colnames(env.new), 
                    response.vars = colnames(alle), ntree = 500, transform = NULL, 
                    maxLevel = lev, corr.threshold = 0.5, compact = T, nbin = 201)
gf <- gradientForest(data = cbind(env.new, alle), 
                     predictor.vars = colnames(env.new), response.vars = colnames(alle),
                     ntree = 500, transform = NULL, compact = T,
                     nbin = 201, maxLevel = lev, corr.threshold = 0.5)

gf

pdf("GFMeehaniaCorCoef0.8.pdf")
plot(gf, plot.type = "O")
most_important <- names(importance(gf))[1:9]

plot(gf, plot.type = "C", imp.vars = most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                              0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
plot(gf, plot.type = "C", imp.vars = most_important, show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), 
                                                     mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
plot(gf, plot.type = "P", show.names = F, horizontal = F,
     cex.axis = 1, cex.labels = 0.7, line = 2.5)
dev.off()

######Plot RDA result##########

#env.vip<- subset(env,select = most_important)
rda <- rda(alle ~ ., data=env.new, scale=T)
signif.all <- anova.cca(rda, parallel=getOption("mc.cores"))
rda #######proportion=27.06%, p=0.001
summary(eigenvals(rda, model = "constrained")) ####Proportion Explained RDA1=25.80%    RDA2=17.03%    RDA3=14.20% 
group<- c(rep("QL",15),rep("MS",25),rep("DQ",20))
env.new<- data.frame(env.new,group)
eco <- env.new$group
bg <- c("#E07F35","#5F9A8A","#66A51E") # 3 nice colors for our ecotypes
pdf("RDA_all_alle_new.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the species
text(rda, scaling=3, display="bp", col="#0071BC", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(rda, type="n", scaling=3, choices=c(1,3))
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(1,3))
text(rda, scaling=3, display="bp", col="#0071BC", cex=1, choices=c(1,3))
legend("topright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

#######Each enviromental variables######

bio02= subset(env,select = c(bio02))
bio04= subset(env,select = c(bio04))
bio10= subset(env,select = c(bio10))
bio11= subset(env,select = c(bio11))
bio12= subset(env,select = c(bio12))
bio13= subset(env,select = c(bio13))
bio15= subset(env,select = c(bio15))



bio02.rda <- rda(alle ~ ., data=bio02, scale=T)
bio02.rda#######proportion=5.14%
signif.bio02 <- anova.cca(bio02.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio04.rda <- rda(alle ~ ., data=bio04, scale=T)
bio04.rda#####proportion=6.50%
signif.bio04 <- anova.cca(bio04.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio10.rda <- rda(alle ~ ., data=bio10, scale=T)
bio10.rda#######proportion=3.97%
signif.bio10 <- anova.cca(bio10.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio11.rda <- rda(alle ~ ., data=bio11, scale=T)
bio11.rda#######proportion=4.02%
signif.bio11 <- anova.cca(bio11.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio12.rda <- rda(alle ~ ., data=bio12, scale=T)
bio12.rda#######proportion=3.94%
signif.bio12 <- anova.cca(bio12.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio13.rda <- rda(alle ~ ., data=bio13, scale=T)
bio13.rda#######proportion=4.48%
signif.bio13 <- anova.cca(bio13.rda, parallel=getOption("mc.cores")) #####p=0.001 ***
bio15.rda <- rda(alle ~ ., data=bio15, scale=T)
bio15.rda#######proportion=6.27%
signif.bio15 <- anova.cca(bio15.rda, parallel=getOption("mc.cores")) #####p=0.011 ***