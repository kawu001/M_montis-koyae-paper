
#=============================================================================================================#
#
# Script created by Angelica Cuevas
# 
# Study published in Molecular Ecology :
# Cite as : Cuevas et al., 2020. Intraspecific genomic variation and local adaptation in a young hybrid species. Molecular Ecology
#
#
# This script: run all analysis in .... Produce statistical analysis, figures and supplementary material
#	
# Usage notes: run line by line
#=============================================================================================================#

#===========================================================================================================================================================================####
#-------------------------------------- Genomic variation for House, Spanish and Italian sparrows --------------------------------------------------------------------------####

rm(list=ls())
working.directory <- setwd("C:/Users/Bashir/Desktop/Meehania/Population results/IBE by MMRR/5 GF Seperate bioclims")
working.directory
source('multiplot_function.R', chdir = F)
source('multiplot_function.R', chdir = F)

library(adegenet)
library(ggfortify)
library(dplyr)
library(ggplot2)
library("lattice")


#=============================================================================================================================================================================####
#----------------------- Local adaptation to climate. Isolation by environment (IBE) and Isolation by distance (IBD) in the Italian sparrow ----------------------------------####

#######------------------------------------------------------------------------------------------#######
#######--------------------------  Environmental variables across Italy   -----------------------#######
#######--------------------- Maps - Loading values for the climatic variables -------------------#######
#######------------------------------------------------------------------------------------------#######
# setwd("/Users/angelimc/Documents/Angelica/RAD_seq.Analysis/RAD_2019/Inland.Italians/Environment")
#### packages ####

# CLIMATE
# Altitude = Elevation
# BIO7 = Temperature Annual Range
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO17 = Precipitation of Driest Quarter

# Read Geodistance and Ecological distance data
library(tseries)

# Read in all the bioclims seperately
#Alt
Alt <- read.table("./5 GF Seperate bioclims/BIO.Elevation.txt", header = TRUE)
#TEMP.A
bio7 <-read.table("./5 GF Seperate bioclims/BIO7.Temperature Annual Range.txt", header = T)
#TEMP.W
bio8 <-read.table("./5 GF Seperate bioclims/BIO8.Mean Temperature of Wettest Quarter.txt", header = T)
#TEMP.D
bio9 <-read.table("./5 GF Seperate bioclims/BIO9.Mean Temperature of Driest Quarter.txt", header = T)
#PREP.D
bio17 <-read.table("./5 GF Seperate bioclims/BIO17.Precipitation of Driest Quarter.txt", header = T)

geo.lat <- read.table("C:/Users/Bashir/Desktop/Meehania/Population results/IBE by MMRR/7.Pop.Coordinates.in.decimals.txt", header = T)
geo.lat
class(geo.lat)

######### Create distance matrices ########
####	Standarizing Distance Matrices
#functions for standardizing the values so coefficients comparable as beta-weights
library("sp")

standardize = function(x){tmp=(x-mean(x,na.rm=T))/sd(x,na.rm=T);diag(tmp)=0;return(tmp)}

prep.standardize <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x=(x-mean(x))/sqrt(var(x)) 
x}

prep <- function(x){x=data.matrix(x)
x=x[lower.tri(x, diag = FALSE)]
x}

####	GEOGRAPHIC DISTANCE - GREAT CIRCLE DISTANCE MATRIX 
tmp = geo.lat[,c('Population','Longitude','Latitude')]
tmp
pop.loc = as.matrix(tmp[,c('Longitude','Latitude')])
pop.loc
rownames(pop.loc) = tmp[,'Population']
rownames(pop.loc)
geo = spDists(pop.loc,longlat=TRUE) ; colnames(geo) = rownames(geo) = rownames(pop.loc)
geo

GEO = as.matrix(standardize(geo))
GEO

# Altitude
tmp = Alt[c('Population','alt')]
tmp
ma.alt = as.matrix(tmp[,'alt'])
rownames(ma.alt) = tmp[,'Population']
order = rownames(pop.loc)
ma.alt.raw = ma.alt[order,]
ALT = as.matrix(standardize(as.matrix(dist(ma.alt.raw))))
ALT
class(ALT)
write.table(ALT, "Matrix_Altitude.txt", sep = "\t")

# BIO7 = Temperature Annual Range TEMP.A
tmp = bio7[,c('Population','bio7')]
temp.a = as.matrix(tmp[,'bio7'])
rownames(temp.a) = tmp[,'Population']
order= rownames(pop.loc)
temp.a.raw=temp.a[order,]
TEMP.A = as.matrix(standardize(as.matrix(dist(temp.a.raw))))
TEMP.A
write.table(TEMP.A, "Matrix_Temp.A.txt", sep = "\t")

# BIO8 = Mean Temperature of Wettest Quarter #TEMP.W
tmp = bio8[,c('Population','bio8')]
temp.w = as.matrix(tmp[,'bio8'])
rownames(temp.w) = tmp[,'Population']
order= rownames(pop.loc)
temp.w.raw=temp.w[order,]
TEMP.W = as.matrix(standardize(as.matrix(dist(temp.w.raw))))
TEMP.W
write.table(TEMP.W, "Matrix_Temp.W.txt", sep = "\t")

# BIO9 = Mean Temperature of Driest Quarter #TEMP.D
tmp = bio9[,c('Population','bio9')]
temp.D = as.matrix(tmp[,'bio9'])
rownames(temp.D) = tmp[,'Population']
order= rownames(pop.loc)
temp.D.raw=temp.D[order,]
TEMP.D = as.matrix(standardize(as.matrix(dist(temp.D.raw))))
TEMP.D
write.table(TEMP.D, "Matrix_Temp.D.txt", sep = "\t")

# BIO17 = Precipitation of Driest Quarter #PREP.D
tmp = bio17[,c('Population','bio17')]
prep.D = as.matrix(tmp[,'bio17'])
rownames(prep.D) = tmp[,'Population']
order= rownames(pop.loc)
prep.D.raw=prep.D[order,]
PREP.D = as.matrix(standardize(as.matrix(dist(prep.D.raw))))
PREP.D
write.table(PREP.D, "Matrix_Prep.D.txt", sep = "\t")

#### reading FST the matrix
fst.m<-as.matrix(read.table("C:/Users/Bashir/Desktop/Meehania/Population results/IBE by MMRR/Univariate, Multivariate, CA/genetic distance.txt"))
class(fst.m)
head(fst.m)
fst.m[upper.tri(fst.m)] = fst.m[lower.tri(fst.m)]
head(fst.m)
str(fst.m)


#######-------------------------------------------------------------------------------------------#######
#######-------------------------------  IBD and IBE   --------------------------------------------#######
#######-- Univariate, Multivariate Regression with radomization (MMRR) and Commonality Analyses --#######
#######-------------------------------------------------------------------------------------------#######
source("C:/Users/Bashir/Desktop/Meehania/Population results/IBE by MMRR/MMRR_function.R", chdir = F)
source("C:/Users/Bashir/Desktop/Meehania/Population results/IBE by MMRR/multiplot_function.R", chdir = F)

#install.packages("yhat")
library(yhat)

# defining predictors
#predictors=list(GEO=GEO, A.TEMP=A.TEMP, A.PREC=A.PREC, TEMP.S=TEMP.S, PREC.S=PREC.S, ALT=ALT)
# predictors=list(GEO=GEO, A.TEMP=A.TEMP, A.PREC=A.PREC, TEMP.S=TEMP.S, PREC.S=PREC.S, ALT=ALT, BEAK.H=BEAK.H, BEAK.L=BEAK.H)
# predictors=list(GEO=GEO, A.PREC=A.PREC, TEMP.S=TEMP.S, PREC.S=PREC.S, BEAK.H=BEAK.H, BEAK.L=BEAK.L)
predictors=list(GEO=GEO, ALT=ALT, TEMP.A=TEMP.A, TEMP.W=TEMP.W, TEMP.D=TEMP.D, PREP.D=PREP.D)

####------------------                   Table 1.                  ------------------####
###	Table 1 - UNIVARIATE TESTS with randomization ####
colnames = c('predictor','r2','coefficient','tstatistic','pvalue')#
uni.table = data.frame(matrix(ncol=length(colnames),nrow=length(predictors)))#
colnames(uni.table) = colnames#
i = 1#
for(i in 1:length(predictors)){#
  cat(i,'\n')#
  tmp.fst.m = fst.m[colnames(predictors[[i]]),colnames(predictors[[i]])]#
  fit = MMRR(tmp.fst.m,predictors[i],nperm=1000)#
  uni.table[i,'predictor'] = names(predictors)[i]#
  uni.table[i,'r2'] = sprintf("%.3f", round(fit$r.squared,3))#
  uni.table[i,'coefficient'] = sprintf("%.3f", round(fit$coefficients[2],3))#
  uni.table[i,'tstatistic'] = sprintf("%.3f", round(fit$tstatistic[2],3))#
  uni.table[i,'pvalue'] = sprintf("%.3f", round(fit$Fpvalue,3))#
}#

colnames(uni.table) = c('','r2','B','t','p')#
uni.table[uni.table[,'p'] < 0.001,'p'] = '0.001'#
uni.table = uni.table[rev(order(uni.table$r2)),]#
uni.table

# write.table(uni.table,file = "table.uni_italianMF.0.02.txt",quote=F,col.names=T,row.names=F,sep='\t')#
write.table(uni.table, file = "table.uni_Meehania.txt", quote = F, col.names = T, row.names = F, sep = '\t')

#### Plotting of univariant regressions ###### 

# FST vs Temperature Seasonality_Distance
TemperatureS<- read.table("FST.v.s.TEMP.S.txt", header = T) 
ts.fst<-ggplot(TemperatureS, aes(TEMP.S,Mean.FST)) + geom_point(size = 3) +
  stat_smooth(method = "lm", col = "black", size = 0.3)+
  labs(y="Mean Pairwise FST",x="Pairwise Temperature Seasonality Distance", size=0.10)+
  theme(axis.title=element_text(size=12))
ts.fst

#regression
linear_regres<-lm(Mean.FST~TEMP.S, TemperatureS)
summary(linear_regres) 
# Multiple R-squared:  0.163,	Adjusted R-squared:  0.1309 
# F-statistic: 5.065 on 1 and 26 DF,  p-value: 0.0331*

# Mantel test is significant as well, see below

####------------------                   Table 2.                  ------------------####
####	Table 2 - Multivariate Matrix Regression with Randomization - MMRR ####

multi.table = c()
resBootstrap.list = list()

# MMRR to get regression coefficients and significance values for multivariate model
# data for MMRR
fit = MMRR(Y = fst.m, X = predictors, nperm=1000)
fit

write.table(fit, file = "MMRR_Meehania.txt", sep = '\t')

#write.table(fit, file = "MMRR_Meehania.txt", sep = '\t')

# 
####	Distance matrices 
# 

FST = prep.standardize(fst.m) 
C.GEO = prep(GEO)
C.ALT = prep(ALT)
C.TEMP.A = prep(TEMP.A)
C.TEMP.W = prep(TEMP.W)
C.TEMP.D = prep(TEMP.D)
C.PREP.D = prep(PREP.D)


# Do Commonality Analysis to get Unique, Common, and Total contribution of each predictor variable
data=data.frame(FST = FST, 
                C.GEO,
                C.ALT,
                C.TEMP.A,
                C.TEMP.W,
                C.TEMP.D,
                C.PREP.D)

ca_i = regr(lm(FST ~ C.GEO + 
                 C.ALT + C.TEMP.A + C.TEMP.W + C.TEMP.D + C.PREP.D, data=data))

ca_i

#add in the common
colnames = c('model','predictor','coefficient','tstatistic','pvalue','Unique','Common','Total')
table = data.frame(matrix(ncol=length(colnames),nrow=length(predictors)))
colnames(table) = colnames
model = paste0('Fst ~ ',paste(names(predictors),collapse=' + '))
table[,'model'] = c(model,paste0('r2 = ',round(fit$r.squared,2)),'','','','')
#table[,'model'] = c(model,paste0('r2 = ',round(fit$r.squared,2)),'')
table[,'predictor'] = names(predictors)
table[,'coefficient'] = round(fit$coefficients[-1],3)
table[,'tstatistic'] = round(fit$tstatistic[-1],2)
table[,'pvalue'] = round(fit$tpvalue[-1],2)
table[,'Unique'] = round(ca_i$Commonality_Data$CCTotalbyVar[,'Unique'],3)
table[,'Common'] = round(ca_i$Commonality_Data$CCTotalbyVar[,'Common'],2)
table[,'Total']  = round(ca_i$Commonality_Data$CCTotalbyVar[,'Total'],2)
table = apply(table,2,as.character)
table

#multi.table = rbind(multi.table,table,rep('',ncol(table)))
multi.table = rbind(multi.table,table,rep('',ncol(table)))
write.table(multi.table, file = "MMRR_Multitable_Meehania_model1.txt",quote=F,col.names=T,row.names=F,sep='\t')

####------------------                   Figure. S3.                  ------------------####
# Commonality analysis with confidence intervals ####
#bootstrap procedure modified from code in supplementary of Prunier et al. (2015), which was based on methods in Peterman et al. (2014)
#	Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. Molecular Ecology, 24, 263–283.
#	Peterman WE, Connette GM, Semlitsch RD, Eggert LS (2014) Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology, 23, 2402–2413.

nperm=1000
n.predictors = 6
ncombos = ((2^n.predictors)-1)
boot=matrix(data = 0, nrow = nperm, ncol = ncombos)
resBootstrap=data.frame(n=rep(0,ncombos),o=rep(0,ncombos),l=rep(0,ncombos),u=rep(0,ncombos),p=rep(0,ncombos))
n=ncol(fst.m)
sn=0.9*n

for (i in 1:nperm){
  rarray= sort(sample(n,sn,replace=F))
  mmFst = fst.m[rarray,rarray][lower.tri(fst.m[rarray,rarray],diag=F)]
  mmGEO= prep(GEO[rarray,rarray])
  mmALT = prep(ALT[rarray,rarray])
  mmTEMP.A = prep(TEMP.A[rarray,rarray])
  mmTEMP.W =prep(TEMP.W[rarray,rarray])
  mmTEMP.D =prep(TEMP.D[rarray,rarray])
  mmPREP.D =prep(PREP.D[rarray,rarray])
  
  comm=regr(lm(mmFst ~ mmGEO + mmALT + mmTEMP.A + mmTEMP.W + mmTEMP.D + mmPREP.D))
  boot[i,]=comm$Commonality_Data$CC[c(1:ncombos),1]
  print(i)
}

for (f in 1:ncombos){
  q=quantile(boot[,f], c(.025,.975))
  resBootstrap[f,1]=f
  resBootstrap[f,3]=q[1]
  resBootstrap[f,4]=q[2]
}

resBootstrap[,2]=comm$Commonality_Data$CC[c(1:ncombos),1]
resBootstrap[,5]=comm$Commonality_Data$CC[c(1:ncombos),2]
rownames(resBootstrap) = rownames(ca_i$Commonality_Data$CC)[1:ncombos]
resBootstrap[,c('o','l','u')] = resBootstrap[,c('o','l','u')]
ca_i$Commonality_Data$CC[c(1:ncombos),1]
new.rownames = gsub('\\s|[?!Unique$|Common$|to$|and$]','',rownames(resBootstrap))
rownames(resBootstrap) = new.rownames

tmp = resBootstrap
tmp = round(tmp,3)
tmp$n = rev(tmp$n)

par(mar=c(4,2,4,30))
xlim = c(min(tmp[,3]), max(tmp[,4]))
plot(tmp[,1],xlim=xlim,font=5,lab=c(ncombos, 7, 1),xaxt="n",yaxt="n",cex.lab=1,xlab='',ylab='',cex=1.5) 
arrows(tmp[,3],tmp[,1],tmp[,4],tmp[,1],code=3,length=0.05,angle=90,lwd=1)
abline(v=0,lty=2)
axis(1,cex.axis=0.9)
mtext('Confidence Intervals of Correlation Coefficient',1,line=2.5,cex=0.9)
#yaxis
x=cbind(1:63,rev(rownames(tmp)))
axis(4,at=x[,1],labels=x[,2],cex.axis=0.75,las=2)
mtext('Predictor Sets',2,cex=0.9,line=.5)
#upper xaxis
at = seq(0,round(max(tmp[,4]),2),by=round(max(tmp[,4]),2)/5)
labels= round(at/sum(tmp[,'o']),2)
axis(3,at=at,labels=labels,cex.axis=0.75)
mtext('% Total',3,line=2.5,cex=1)
#coefficients
mtext('Coefficient',4,at = 65,las=2,line=17,adj=.2)
axis(4,line=17,at=tmp[,1],labels=tmp[,'o'],tick=F,cex.axis=0.8,las=2,hadj=1)
mtext('% Total',4,at = 65,las=2,line=24,adj=.5)
axis(4,line=24,at=tmp[,1],labels=tmp[,'p'],tick=F,cex.axis=0.75,las=2,hadj=1)

mtext('Total',4,at = -.5,las=2,line=14,adj=.5, cex=0.9)
mtext(sum(tmp[,'o']),4,at = -.5,las=2,line=17,adj=.5, cex=0.8)
mtext(100,4,at = -.5,las=2,line=24.5,adj=.5, cex=0.8)
par(mfrow=c(1,1))
dev.off()
#######-------------------------------------------------------------------------------------------#######
#######-------------------------------  Mantel test Beak/TempS v.s. Fst   -------------------------------#######
#######-------------------------------------------------------------------------------------------#######
# install.packages("ade4") #=> this required objects of class *.dist to run the mantel test
library(ade4)            #   does not accept only the matrix I created for Geographic and Genomic (fst) Distance
# install.packages("vegan") # -> for mantel test
library(vegan)

# FST matrix
fst.m

# Beak matrices
BEAK.H
BEAK.L

# TEMP.S matrix
TEMP.S

mantel_test_BH<-mantel(fst.m, BEAK.H, method="pearson", permutations=999)
mantel_test_BH # Mantel statistic r: 0.1908 - Significance: 0.266 

mantel_test_BL<-mantel(fst.m, BEAK.L, method="pearson", permutations=999)
mantel_test_BL # Mantel statistic r: 0.1723 - Significance: 0.208 

mantel_test_TS<-mantel(fst.m, TEMP.S, method="pearson", permutations=999)
mantel_test_TS # Mantel statistic r: 0.4038 - Significance: 0.029 *



######
#
#
#
######----------------------------------------------------------------------------------------------------------------------------------
# Cite as : Cuevas et al., 2020. Intraspecific genomic variation and local adaptation in a young hybrid species. Molecular Ecology
######----------------------------------------------------------------------------------------------------------------------------------
#
#
#
#======================================================================================================================================================================================######
#############################################################################################################################################################################################
#======================================================================================================================================================================================######