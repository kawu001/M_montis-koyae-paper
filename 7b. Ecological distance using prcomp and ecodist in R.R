#############################################################################
#How to calculate ecodist using prcomp for pca and ecodist for distance matrix
##############################################################################

#1. The ranked bioclims from gradient forest including altitude will be used as the env data
#### R
#Clean up!
rm(list=ls())
library("dplyr")
library("stats")
library("tidyverse")
library("ecodist")
library(stats)

#Read in the data (bioclims per population)
bio.env <- read.csv(file = "GF5envbio_0.8.csv", header = TRUE)
bio.env

#Transpose the data from population as column, and row as bioclims to column (bioclims) and rows (populations)
t.env <- t(bio.env)
t.env
write.csv(t.env, file = "transformed_bioclim.csv")

#[,1]    [,2]     [,3]    [,4]    [,5]    [,6]    [,7]   
#Population "JP.CS" "JP.XYS" "CH.FJ" "CH.JL" "CH.SY" "CH.WY" "CH.XE"
#I deleted the header, leaving the population as rows
#Read in the data (transposed bioclims)
env <- read.table(file = "transformed_bioclim.txt", header = TRUE)
env

#Explore Env data
dim(env)
str(env)
head(env)
summary(env)
pairs(env[2:8], main = "Bivariate Plots of the Environmental Data")

pca.env <- prcomp(env[2:8], scale. = TRUE, center = TRUE)
pca.env$rotation
names(pca.env)
summary(pca.env)
str(pca.env)

attributes(pca.env)

pca.env$rotation

#Get the eigenvalues of principal components
#Adopt the one with eigenvalue greater than 1
eig <- (pca.env$sdev)^2
eig
# 6.843617e+00 1.401097e-01 1.589881e-02 3.740380e-04 3.860153e-33 this shows that pcs1-3.

#Get variances of principal components
variance <- eig*100/sum(eig)
variance

#Cummulative variances
cumvar <- cumsum(variance)
cumvar

#Combine to dataframe
eig.var.cum <- data.frame(eig=eig, variance=variance, cumvariance=cumvar)
eig.var.cum

#barplot of principal components eigenvalues
barplot(eig.var.cum[1:5,1], names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"))

#barplot of principal components variances
barplot(eig.var.cum[1:5,2], names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"))

#The PCs were copied manually, since I was unable to transform the data from prcomp dataframe to normal
#Load the copied data and select the first 3 with cummulative variance summing up to 999.99466

env.pcs <- read.table(file = "3 PCs from prcomp (transformed).txt", header = TRUE)
env.pcs
str(env.pcs)
summary(env.pcs)

ncol(env.pcs) #3 columns

nrow(env.pcs) #7 rows

library(stats)
library(ecodist)

#Calculate the distance using rows 2-5, since row 1 is alphabets
eco.dist <- distance(env.pcs, method = "euclidean", sprange=NULL, spweight=NULL, icov)
eco.dist

eco.dist <- data.frame(as.matrix(eco.dist))
eco.dist

#write.csv(eco.dist, file = "ecodist_dissimilarity.csv")
write.csv(eco.dist, file = "ecodist_dissimilarity_transformed.csv")