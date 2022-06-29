## O.L. Pescott
## 29/06/2022
## Make French neighbourhoods for Isle-de-France Frescalo application
#rm(list=ls())
#.pardefault <- par()
#par(.pardefault)
library(sf)
library(terra)
library(raster)
## Helper function
rbind.fill <- function(x) {
  nam <- sapply(x, names)
  unam <- unique(unlist(nam))
  len <- sapply(x, length)
  out <- vector("list", length(len))
  for (i in seq_along(len)) {
    out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
  }
  setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
} ## End

#############################################################
## Part 1 -- extract info for estimating commune similarities
#############################################################
# Read in shapefile from Jeanne
idfShp <- st_read(dsn = "data/COMMUNES-IDF/COMMUNE-IDF.shp")
crs(idfShp) # EPSG:2154; https://epsg.io/2154
plot(idfShp[1])

# European landcover
clc <- rast(x = "W:/PYWELL_SHARED/Pywell Projects/BRC/Oli Pescott/___SCRATCH/maps/CLC2018_EuropeLandCover/CLC2018/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")
plot(clc, legend = FALSE)

## Extract relevant part of Euro landcover (possibly not necessary, could just skip to extract below, but maybe this speeds up extract as it is not working with
# the entire European raster)
idfBuf <- st_buffer(st_union(idfShp[1]), 1000)
plot(idfBuf)
# create bounding box coords as well for visualisations later
#idfBuf4326 <- st_transform(st_buffer(st_union(idfShp[1]), 1000), crs = 4326)
#idfBB <- st_bbox(idfBuf4326)
#save(idfBB, file = "outputs/idfBB.rdata")
# carry on with land cover crop
#clc2Crop <- crop(x = clc, y = st_transform(idfBuf, crs = "epsg:3035"))
#plot(clc2Crop)
#writeRaster(clc2Crop, filename=file.path("outputs/", "100mCLC2018Crop_IDF.tif"), overwrite=TRUE)
#clc2Crop <- rast(x = file.path("outputs/", "100mCLC2018Crop.tif"))
#clcCrop2154 <- project(clc2Crop, "epsg:2154")
#writeRaster(clcCrop2154, filename=file.path("outputs/", "100mCLC2018Crop_2154.tif"), overwrite=TRUE)
clcCrop2154 <- rast(x = file.path("outputs/", "100mCLC2018Crop_2154.tif"))
plot(clcCrop2154, legend = FALSE)
plot(idfBuf, bg = "transparent", add = TRUE)

# raster::extract needs raster to be RasterLayer not SpatRaster, so read with raster::raster(), not terra::rast()
#clcCrop2154_rl <- raster(x = file.path("outputs/", "100mCLC2018Crop_2154.tif"))
#e <- extract(clcCrop2154_rl, idfShp) # slow
#class.counts <- lapply(e, table)
#save(class.counts, file = "outputs/clc2018idfCounts.Rdata")
load(file = "outputs/clc2018idfCounts.Rdata")

p.counts <- rbind.fill(class.counts)
p.counts[is.na(p.counts)] <- 0

head(as.data.frame(idfShp))
idfInf <- data.frame(idfShp[,c(1:3)])
idfInfC <- cbind(idfInf, p.counts); head(idfInfC)
idfLcM <- as.matrix(idfInfC[,c(5:29)]); head(idfLcM)  # just keep land cover info for matrix operations
## normalise relative to total commune area (otherwise commune similarity might just be based on size rather than composition of land covers)
idfLcMN <- t(round(apply(idfLcM, 1, function(x) x/sum(x)), digits = 4))

## for the moment just stick with land cover and distance as the weights, although could add in geology, soil etc. later on
## join everything back together
idfShp2 <- st_drop_geometry(idfShp)
allDat <- cbind(idfShp2[,c(1:3,5:6)], idfLcMN); head(allDat)
#write.csv(allDat, file = "outputs/finalIDF_communeNormalisedLandcoverMatrix.csv")

#################################################################
## Part 2 -- Create lists of weighted neighbours for each commune
#################################################################
library(rgdal)
library(spdep)
library(plyr)
library(useful) # for corner views of large matrices
##
#allDat <- read.csv(file = "outputs/finalIDF_communeNormalisedLandcoverMatrix.csv", header = T, stringsAsFactors = F)
#head(allDat)
index <- allDat[,c(1:5)] # ids
lcNum <- allDat[,c(6:30)] # normalised land cover data
str(lcNum) # all numeric
## Calculate 100 nearest neighbours spatially -- this could be any number, but go with 100 initially
## re-add centroids, as Paris ones are missing
idfShpXY <- cbind(idfShp, st_coordinates(st_centroid(idfShp)))
numbers <- idfShpXY[,c("X","Y")]
numbers <- as(numbers, Class = "Spatial")
#coordinates(numbers) = ~X+Y
#proj4string(numbers) = CRS("+init=epsg:2154")

numbers.df <- SpatialPointsDataFrame(numbers, data.frame(id=1:length(numbers), nom = index$NOM_COMM))
# spdep functions 'knearneigh' -- get 100 nearest neighbours
numbers.knn <- spdep::knearneigh(numbers.df, k = 100) # calculate k nearest neighbours
# For each hectad, from closest 100, select the 50 with the closest land cover similarity
#numbers.knn$nn # matrix, rows hectads, columns neighbours
# str(numbers.knn$nn)
#int [1:1300, 1:100] 6 3 21 3 8 1 10 9 8 7 ...
data_matrix <- numbers.knn$nn
## add target hectad to first column (otherwise hectad itself is not considered in Frescalo -- see Hill 2012), and delete 201st
data_matrix <- cbind(1:1300, data_matrix); data_matrix <- data_matrix[,c(1:100)]

## need to iterate by row, and by cell within row, grabbing relevant rows of hectad env data, calculating similarity
## and writing to a matrix of the same dimensions
topleft(data_matrix)
dim(data_matrix)
out <- outR <- outD <- matrix(nrow = 1300, ncol = 100)
for (j in 1:1300){
  for (i in 1:100){
    c <- data_matrix[j,i] # lookup hectad reference for each neighbour of j (j = row = hectad in data_matrix)
    Mat <- as.matrix(lcNum[c(j,c),]) # grab relevant rows of env data from hectad table
    sim <- Mat / sqrt(rowSums(Mat * Mat))
    sim <- sim %*% t(sim)
    out[j,i] <- sim[2,1]
    outD[j,i] <- i
  }
}
#save(out, file = "outputs/envCosineSimMatrix_idfCommunes.rdata")
#load(file = "outputs/envCosineSimMatrix_idfCommunes.rdata")

## Now get ranks of environmental similarities for each 50 hectad neighbourhood
for (j in 1:1300){
  outR[j,] <- rank(-out[j,], ties.method = "first") # minus negates rank so ordering is reversed (larger similarity = highest rank)
}
topleft(outR)
## and the original hectad IDs of those <= 50
outTop50_id <- outTop50_envRank <- outTop50_distRank <- matrix(nrow = 1300, ncol = 50)
for (j in 1:1300){
  outTop50_id[j,] <- subset(data_matrix[j,], outR[j,]<51) # get original hectad IDs, in distance order from target, for env similarity ranks <= 50
  outTop50_envRank[j,] <- subset(outR[j,], outR[j,]<51) # also get actual environmental similarity ranks in distance order, for env similarity ranks <= 50
  outTop50_distRank[j,] <- subset(outD[j,], outR[j,]<51) # also get actual distance similarity ranks in distance order, for env similarity ranks <= 50
}
topleft(outTop50_id)
topleft(outTop50_envRank)
topleft(outTop50_distRank)

## Distance-based weights -- but this should only be for the 50 environmentally similar hectads
# the actual numbers of data_matrix are the orignal hectad IDs based on distance
outTop50 <- as.data.frame(t(outTop50_id))
names(outTop50) <- rep(1:1300,1)
topleft(outTop50)
outTop50 <- reshape(outTop50, direction = "long", varying = list(1:1300))
outTop50$time <- as.factor(outTop50$time) # change to factor so can use by() in split-apply-combine style to interate function across target neighbourhoods
head(outTop50)

outTop50E <- as.data.frame(t(outTop50_envRank))
names(outTop50E) <- rep(1:1300,1)
topleft(outTop50E)
outTop50E <- reshape(outTop50E, direction = "long", varying = list(1:1300))
outTop50E$time <- as.factor(outTop50E$time) # change to factor so can use by() in sp
head(outTop50E)

outTop50D <- as.data.frame(t(outTop50_distRank))
names(outTop50D) <- rep(1:1300,1)
topleft(outTop50D)
outTop50D <- reshape(outTop50D, direction = "long", varying = list(1:1300))
outTop50D$time <- as.factor(outTop50D$time) # change to factor so can use by() in sp
head(outTop50D)

outTop50$rankDist <- outTop50D[,"1"]
outTop50$rankEnv <- outTop50E[,"1"]
outTop50 <- outTop50[,c(1,2,4,5)] # drop ID which is auto-generated by reshape
names(outTop50) <- c("target", "neighID", "rankDist", "rankEnv" )
head(outTop50)


########################################################################################################################
## Now calculate actual weights from ranks as per Hill (2012)
########################################################################################################################
dlist <- by(outTop50, outTop50$target, function(x) ((1-((x[,"rankDist"]-1)^2)/150^2)^4) * ((1-((x[,"rankEnv"]-1)^2)/100^2)^4) ) # change formula to match Hill (2012)
outTop50$weights <- unlist(dlist)
head(outTop50)
longformOrder <- do.call(rbind, by(outTop50, outTop50$target, function(x) x[order(x$weights, decreasing = T),]))
head(longformOrder)
## few quick checks:
head(outTop50D[outTop50D$target==5,])
head(longformOrder[longformOrder$target==5,])
(max(longformOrder[longformOrder$target==2,]$rankDist))
max(longformOrder$rankDist)

####################################################################################
## Finally, add in commune IDs, as makes live easier, and Frescalo uses site names
####################################################################################
labels <- numbers.df@data

longform2 <- merge(longformOrder, labels, by.x = 'target', by.y = 'id', all.x = TRUE)
names(longform2)[6] <- c('target.let') 

longform3 <- merge(longform2, labels, by.x = 'neighID', by.y = 'id', all.x = TRUE)
longform3 <- longform3[order(longform3$target, longform3$weights),]
names(longform3)[7] <- c('neighbour.let') 

final_weights <- longform3[,c('target.let','neighbour.let','weights')]
final_weights <- do.call(rbind, by(final_weights, final_weights$target.let, function(x) x[order(x$weights, decreasing = T),]))
head(final_weights)
nrow(final_weights)
row.names(final_weights) <- 1:nrow(final_weights)
###
#write.csv(final_weights, file = "outputs/idfCommuneWeights_2022-06-29_v0.csv")
###

### Create files for Frescalo visualisation
idfNew_v0 <- final_weights
#save(idfNew_v0, file = "outputs/idfCommuneWeights_2022-06-29_v0.rdata")
load(file = "outputs/idfCommuneWeights_2022-06-29_v0.rdata")

# add in x/y for visualisation purposes here
idfGrLookup <- data.frame(letGr = idfShpXY$NOM_COMM, x = idfShpXY$X, y = idfShpXY$Y, region = "IDF")
idfGrLookup <- idfGrLookup[order(idfGrLookup$letGr),]; head(idfGrLookup)
#save(idfGrLookup, file = "outputs/idfCentroidXY.rdata")

