#sdm: a reproducible and extensible R package for species distribution modelling: Tutorial
#### R
#Clean up!
rm(list=ls())

library("ggplot2")
library("devtools")
library("sdm")
library("dismo")
library("dplyr")
library("tidyr")
library("mapview")
library("raster")
library("rgdal")


#Load in combined locality data
sp <- read.csv("meehania_montis_f_thin1.csv", header=T)
sp

class(sp)

dim(sp)

#check the number of rows to see if filtering worked
nrow(sp)


#To keep the coordinate column name will change from sp to spg
spg <- sp %>% select(lon,lat)
head(spg)

#To make the data presence only, we need to add 1 to species column
spg$species <- 1

#To check if it has added a new column with species and assigned 1
head(spg)

#-------
#Convert the spg dataframe to a special point dataframe
class(spg)

coordinates(spg) <- c("lon","lat")

#This gives this error message "Error in `coordinates<-`(`*tmp*`, value = c("lon", "lat")) : coordinates are not allowed to contain missing values"

#To remove the data with missing coordinates, we'll use "drop" tidyr package
spg <- spg %>% drop_na()

#Check the number of rows to confirm that the occurrence has reduced
nrow(spg)

#Now convert the spg dataframe to a special point dataframe
coordinates(spg) <- c("lon","lat")

#To confirm the new special point dataframe
class(spg)

#---------
##Download the bioclim data
#Valid resolutions are 0.5, 2.5, 5, and 10 (minutes of a degree). 

raster::getD
bio <- raster::getData("worldclim", var = "bio", res=2.5)
bio #information about the downloaded bioclim

names(bio) #check the names of the climate data which should be bio1-19

#Plot one of them e.g. World temperature
plot(bio[[1]])

#To add species occurrence data on the temperature
points(spg)

#If you have knowledge about the species, then you can remove or subset the data to cover a particular region and fit your data to an area
#We will use a rectangle using Extent in the raster package

e <- drawExtent()

#Drag and select the area, then use crop fxn to crop it out
spg <- crop(spg, e)

#To check the cropped part using another colour (red)
points(spg, col="red")

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
ex <- raster::extract(bioc, spg)
head(ex)

#To run the vifstep and vifcor
v <- vifstep(ex)
#Run v to check the list of remaining variables
v
#To remove the unwanted variables
bioc <- exclude(bioc, v)
bioc
#----------
#Now to explore functions in sdm package
library("sdm")
#d <- sdmData(species~., spg, predictors = bioc) #report of the dataset as  presence only

#report of the dataset as presence with pseudo-absence data created using bg = list(method="gRandom", n=1000)
d <- sdmData(species~., spg, predictors = bioc, bg = list(method="gRandom", n=10000))
d

#run the sdm function to fit the model, and you can use getmethodNames() to see available methods
getmethodNames()
#n=3 means 3 replicates for each methods, and specify cores for high performance computing
m <- sdm(species~., d, methods = c('maxent'),
         replication=c('sub','boot'),test..p=25,n=10)

#m <- sdm(species~., d, methods = c('glm','brt','rf','fda'), replication=c('sub','boot'),test..p=30,n=3, parallelSettings=list(ncore=4,method='parallel')) Just to add threads and other settings
#sdm with F1 and click on index to check out the functions in the pacakage

m

#To extract a specific output from any of the models (m@models$species$brt$`8`@object)
#m from the previous step, select the type, and the species (if more than 1)and model type, model number and object
#To view m in the graphical user interface
gui(m)

#To generate the prediction using the predict function and you can specify different format like tiff, gif and so on
p1 <- predict(m, bioc, filename='pr.img') #This generate 20 outputs
p1
names(p1)

#plot different predicted probability for 4 different methods (glm, brt, rf, and fda)
plot(p1[[c(1,7,13,20)]])

#While ensemble generates the 20 output combined
#en1 <- ensemble(m, bio, file name = 'en.img', setting = list(method='weighted',stat='tss',opt=2))
#To generate ensemble for some of the models and not all the 24
#en1 <- ensemble(m, bio, file name = 'en.img', setting = list(id=c(1,3,4,6,12,22),method='weighted',stat='tss'))

#To differentiate the models, as there are some inconsistencies, we'll use the ensemble function. 
#Since, we have generated the predicted prob., then we will insert p1.
en1 <- ensemble(m, p1, filename = 'en.img', setting = list(method='weighted',stat='AUC',opt=2))

#Plot the raster output for the potential distribution for the current time.
plot(en1)

#############
#To generate for the future prediction, we will download future climate CMIP5 for different GCM models and for the target year 2070, resolutions, and variables
#Also, there is different RCP scenarios
#############

biof <- raster::getData('CMIP5', var='bio', res=2.5, rcp=85, model='CN', year=70)
biof
plot(biof[[1]])
names(biof)

#There is need to change the names of the biof to conform with the previous current time.
names(biof) <- names(bio)
#To confirm the changed names
names(biof)

#To do the same for the raster of future biockim variables
bioff <- crop(biof, e)
plot(bioff[[1]])

#To remove the unwanted variables
bio_future <- exclude(bioff, v)
bio_future

#To generate the predicted output for the future climate
en2 <- ensemble(m, bio_future, filename = 'enf.img', setting = list(method='weighted',stat='AUC',opt=2))

#Plot the output for the future
plot(en2)

#############
#To generate for the future prediction, we will download future climate CMIP5 for different GCM models and for the target year 2070, resolutions, and variables
#Also, there is different RCP scenarios
#############
#load the LGM climatic scenario (about 22,000 years ago).
##load the bioclim data
envL <- list.files("D:/MeehaniaLGM_SDM/LGM_2-5m", pattern='tif', full.names=TRUE)
envL
names(envL) #check the names of the climate data which should be bio1-19
bioL <- stack(envL)

plot(bioL[[10]])
names(bioL)

#To do the same for the raster of future biockim variables
bioLGM <- crop(bioL, e)
plot(bioLGM[[3]])

#To generate the predicted output for the future climate
enLGM <- ensemble(m, bioLGM, filename = 'enLGM.img', setting = list(method='weighted',stat='AUC',opt=2))

#Plot the output for the future
plot(enLGM)

#############
#To generate for the future prediction, we will download future climate CMIP5 for different GCM models and for the target year 2070, resolutions, and variables
#Also, there is different RCP scenarios
#############
#load the Mid Holocene climatic scenario (about 6,000 years ago).
##load the bioclim data
envMH <- list.files("D:/MeehaniaLGM_SDM/MIDHL_2-5m", pattern='tif', full.names=TRUE)
envMH
names(envMH) #check the names of the climate data which should be bio
bioMHL <- stack(envMH)

names(bioMHL)
plot(bioMHL)

#To do the same for the raster of future biockim variables
bioMDHL <- crop(bioMHL, e)
plot(bioMDHL[[3]])

#To generate the predicted output for the future climate
enMHL <- ensemble(m, bioMDHL, filename = 'enMDHL.img', setting = list(method='weighted',stat='AUC',opt=2))

#Plot the output for the future
plot(enMHL)
#-------------

#To observe the changes from present to future climate prediction, we will use different colours
#Specify the colours
cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred','blue'))

plot(en1, col=cl(200)) #plot will be restricted to using only the specified colours
plot(en2, col=cl(200))
plot(enLGM, col=cl(200))
plot(enMHL, col=cl(200))

#To interatively explore the colours and merge present and future predictions together
library("mapview")
#Stack the two layers together
stack(en1,en2,enLGM)
mapview(stack(en1,en2,enLGM), col.regions=cl(200))

#You can add the special data points to the mapview
mapview(stack(en1,en2,enLGM), col.regions=cl(200)) + spg #You will get error as the crs for spg is not defined

#To define the crs for the spg
proj4string(spg)
projection(en1)
#Since the projection for en1 is same as that of spg
proj4string(spg) <- projection(en1)

mapview(en1,en2,enLGM, col.regions=cl(200)) + spg #Another error, so you will have to view the en1 one after the other
mapview(en1, col.regions=cl(200)) + spg
mapview(en2, col.regions=cl(200)) + spg
mapview(enLGM, col.regions=cl(200)) + spg

#--------
#Now to quantify the changes in both present and future predictions
ch1 <- en2 - en1
cl2 <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch1,col=cl2(200))

#Now to quantify the changes in both LGM and present predictions
ch2 <- enLGM - en1
plot(ch2,col=cl2(200))
#--------
#We may calculate the changes based on presence and absence data
#To do this, we need to add presence/absence format to the en predictions using thresholds
#To do this we will use the evaluates function in sdm package
d #the data format with presence-background

#The probability of occurrence can be extracted from the ensemble (en)
df <- as.data.frame(d)
head(df)

#Then we need to add the coordinates to be extracted from object (d) and add to df
df <- as.data.frame(species=df$species,coordinates(d))
head(df)

evaluates(df$species,extract(en1,df[]))
df$species
xy <- as.matrix(df[,c('lon','lat')])
head(xy)

#Probability of occurrence
p <- raster::extract(en1,xy)
head(p)
nrow(df)
length(p)

#Comparing the probability results with the ensemble output
ev <- evaluates(df$species,p)

#To see the statistics such nas AUC, and others
ev@statistics

ev@threshold_based
#After picking a threshold, then less say, I picked threshold 3.
th <- ev@threshold_based$threshold[2]

#create a raster from ensemble 1 (en1) and 2, then use the threshold to compare
pa1 <- raster(en1)
pa1[] <- ifelse(en1[] >= th, 1, 0)
plot(pa1)

pa2 <- raster(en2)
pa2[] <- ifelse(en2[] >= th, 1, 0)
plot(pa2)

#Compare pa1 and pa2 by subtraction
chp <- pa2 - pa1
plot(chp) #I will get 3 possible values 1 = presence, 0 = no change either in presence/absence and -1 = absence

#In future, if we have 1 = colonization, and -1 in future = extinction
#Let quantify using colours red = extinction, gray = no change and blue = colonization
plot(chp,col=c('red','gray','blue'))

#---------
#for interpretation, response curve for all the 24 models with confidence interval
rcurve(m)
rcurve(m,id=1:6)
rcurve(m,id=7:12)

#Percentage of importance of the variables with the confidence interval and it can be plotted
getVarImp(m,id=1)
plot(getVarImp(m))
plot(getVarImp(m,method='glm'))
plot(getVarImp(m,method='rf'))

#Mapping the suitability over the environmental space
niche(bio,en1,n=c('bio9','bio15'),col=cl(200))

# measure variable importance
dir.create("D:/MeehaniaLGM_SDM/variable_importance")
for (i in 1:20) {
  lm <- getVarImp(m,id={i})
  capture.output(lm,file = paste0("D:/MeehaniaLGM_SDM/variable_importance/var_imp_",i-1,".txt"))
}


writeRaster(en1, filename='D:/MeehaniaLGM_SDM/asc/en.asc', overwrite=TRUE)
writeRaster(en2, filename='D:/MeehaniaLGM_SDM/asc/en2.asc', overwrite=TRUE)
writeRaster(enLGM, filename='D:/MeehaniaLGM_SDM/asc/enLGM.asc', overwrite=TRUE)
plot(midholocene)
writeRaster(midholocene, filename='enMD_holocene.asc', overwrite=TRUE)