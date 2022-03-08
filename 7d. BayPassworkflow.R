#!/usr/bin/Rscript
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona. If not, see <http://www.gnu.org/licenses/>.

require(corrplot)
require(ape)
require(mvtnorm)
require(geigen)


## Source the baypass R functions

# Full path to "baypass_utils.R" provided with baypass.
# Type: str
# Example: "~/Software/Science/baypass_2.1/utils/baypass_utils.R"
source("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/utils/baypass_utils.R")

####################################
#1. Core Model
####################################
## upload estimate of omega
omega=as.matrix(read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/1.Core Model/Meehania_mat_omega.out"))
pop.names = c("JP-CS","JP-XYS","CH-FJ","CH-JL","CH-SY","CH-WY","CH-XE")
#
dimnames(omega)=list(pop.names,pop.names)

## Compute and visualize the correlation matrix using SVD decomposition
cor.mat=cov2cor(omega)
pdf(file=paste0("./1.Core Model/","omega_corr.pdf"),width = 16, height = 16)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
dev.off()

## Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
pdf(file=paste0("./1.Core Model/", "Hier_clust_tree.pdf"),width = 16, height = 8)
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

## Estimates of the XtX differentiation measures
anacore.snp.res=read.table("./1.Core Model/Meehania_summary_pi_xtx.out",h=T)
pdf(file=paste0("./1.Core Model/", "XtX_diff.pdf"),width = 16, height = 8)
plot(anacore.snp.res$M_XtX)
dev.off()

# Get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution
pi.beta.coef=read.table("./1.Core Model/Meehania_summary_beta_params.out",h=T)$Mean

# Upload the original data to obtain total allele count
current.data<-geno2YN("7populations.BAYPASS")

# Create the POD
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=22561,
                           sample.size=current.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")

file.rename("G.btapods", paste("./1.Core Model/", "G.btapods", sep=""))


#######################################################
# 1b. Sanity Check: Compare POD and original data estimates
#######################################################
# Get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table("./1.Core Model/Meehaniapod_mat_omega.out"))
pdf(file=paste0("./1.Core Model/", "Omega_estimate_from_POD.pdf"),width = 16, height = 8)
plot(pod.omega,omega) ; abline(a=0,b=1)
dev.off()
fmd.dist(pod.omega,omega)

# Get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table("./1.Core Model/Meehaniapod_summary_beta_params.out",h=T)$Mean
pdf(file=paste0("./1.Core Model/", "POD_estimates.pdf"),width = 16, height = 8)
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)
dev.off()

#
#######################################################
# XtX calibration
#######################################################
# Get the pod XtX
pod.xtx=read.table("./1.Core Model/Meehaniapod_summary_pi_xtx.out",h=T)$M_XtX

# Compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)

# Add the thresh to the actual XtX plot
snp.info=read.csv("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/snp.info.csv")
snp.info
head(snp.info)
str(snp.info)

chr.id=as.factor(snp.info$CHROM)
mycolors=rep(rainbow(2),20)

pdf(file=paste0("./1.Core Model/", "XtX_POD_diff.pdf"),width = 16, height = 8)
plot(anacore.snp.res$M_XtX,col=mycolors[chr.id])
abline(h=pod.thresh,lty=2)
dev.off()

anacore.snp.res.1 <- cbind(snp.info,anacore.snp.res)

anacore.snp.res.1

write.table(anacore.snp.res.1,"anacore1.snp.restults",row.names=T,sep=" ",quote=F)
anacore.snp.res.row <- read.table("anacore1.snp.restults",sep = " ",header = T)

anacore.snp.res.row

outliers <- anacore.snp.res.1[,c("ID")][which(anacore.snp.res.1$M_XtX>pod.thresh)]
outliers

anacore.snp.res.1[,"ID"][which(anacore.snp.res.1$M_XtX>pod.thresh)]
xtx.outliers <- anacore.snp.res.1[,"ID"][which(anacore.snp.res.1[,"M_XtX"]>pod.thresh)]
write.table(xtx.outliers,file=paste0("./1.Core Model/","Meehania.xtx.ouliers.snp.txt"),sep = " ",row.names = F,col.names = F,quote = F)

############################################
#2 Analysis under the IS covariate mode (MCMC is run under the core model)
############################################
covis.snp.res=read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/2.IS covariate mode/Meehaniacovis_summary_betai_reg.out",h=T)
graphics.off()

pdf(file=paste0("./2.IS covariate mode/", "BFs_layout.pdf"),width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covis.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

###################################################
#2b. Omega matrix Importance Sampling estimates of the Bayes Factor
###################################################
covis2.snp.res=read.table("./2.IS covariate mode/Meehaniacovis2_summary_betai_reg.out", h=T)
graphics.off()

pdf(file=paste0("./2.IS covariate mode/", "BFs_layout_pass2.pdf"),width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covis2.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis2.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis2.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

#################################################
#3  Analysis under the MCMC covariate mode (MCMC is run under the STD model)
#################################################
covmcmc.snp.res=read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/3.MCMC covariate mode/Meehaniacovmcmc_summary_betai.out", h=T)
covmcmc.snp.xtx=read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/3.MCMC covariate mode/Meehaniacovmcmc_summary_pi_xtx.out",h=T)$M_XtX
graphics.off()

pdf(file=paste0("./3.MCMC covariate mode/", "BFs_layout.pdf"), width = 16, height = 8)
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res$eBPmc,xlab="SNP",ylab="eBPmc")
plot(covmcmc.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covmcmc.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
dev.off()

#################################################
#4 Analysis under the AUX covariate mode: MCMC is run under the AUX model
#################################################
covaux.snp.res.raw=read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/4.AUX model/Meehaniacovaux_summary_betai.out", h=T)
covaux.snp.xtx.raw=read.table("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/4.AUX model/Meehaniacovaux_summary_pi_xtx.out", h=T)$M_XtX
snp.info=read.csv("C:/Users/Bashir/Desktop/Meehania/Population results/Outlier using BayPass2/snp.info.csv")
snp.info
str(snp.info)
head(snp.info)

covaux.snp.res=cbind(snp.info,covaux.snp.res.raw)
covaux.snp.xtx=cbind(snp.info,covaux.snp.xtx.raw)

write.table(covaux.snp.res,"./4.AUX model/covaux.snp.res.results",row.names=T,sep=" ",quote=F)
write.table(covaux.snp.xtx,"./4.AUX model/covaux.snp.xtx.results",row.names=T,sep=" ",quote=F)
mycolors=rep(rainbow(2),20)


pdf(file=paste0("./4.AUX model/", "BFs_layout.pdf"),width = 16, height = 16)
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)",col=mycolors[chr.id])
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"),col=mycolors[chr.id])
plot(covaux.snp.xtx$covaux.snp.xtx.raw,xlab="SNP",ylab="XtX corrected for SMS",col=mycolors[chr.id])
dev.off()

covaus.outliers=covaux.snp.res[,"ID"][which(covaux.snp.res[,"BF.dB."]>20)]
write.table(covaus.outliers,file=paste0("./4.AUX model/","Meehaina.GEA.covaus.ouliers.BF20.snp.txt"),sep = " ",row.names = F,col.names = F,quote = F)
