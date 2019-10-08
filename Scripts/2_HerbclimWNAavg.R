###################################################################################
## PROGRAM FUNCTIONS: EXTRACT CLIMATIC VARIABLES FOR POINT COORDINATES!
###################################################################################
# EXTRACT CLIMATIC VARIABLES FOR POINT COORDINATES:
# Extract ClimateWNA variables from point coordinates and convert to bioclim variables
# Climate records are averaged from a set time period 1981-2010
# ClimateWNA run initially with the time series function from 1901 - 2012 for all records
# ClimateWNA run with DEM 90m res obtained from HydroSHEDS
# (http://hydrosheds.cr.usgs.gov/index.php)

# PART 1: Extract elevation values for records
# PART 2: Calculate 30-year climate normals and coefficients of variation
# PART 3: Average climate within 62 km of focal populations, clipped to NHD waterways
# PART 4: Average climate within 62 km of all occurences, clipped to NHD waterways

######## INITIALIZATION ########################### 

path.root <- "~/Google Drive/CardAdapt" 
# path.root <- "~/Documents/CardAdapt" # Amy's computer
setwd(path.root)
source("Data/Climate/Rcode/1_Setup_directories.R")
source("Rcode/functions.R")

########################################################
## Source Climate Data 
## (Matthew Bayly: April 7, 2014)
## (Last updated by Chris Muir: December 30, 2016)
########################################################
## Use the time series option in ClimateWNA to get annual
## records for the total range of the data (ie. 1901 - 2012)

library(dismo)
library(KernSmooth)
library(maptools)
library(ggplot2)
library(rgdal)
library(rgeos)
library(gam)
library(randomForest)
library(gbm)
library(ncdf4)
library(ncdf4.helpers)

# PART 1: Extract elevation values for records

	# Focal populations: 16 populations used in growth chamber study
	foc.pres <- read.csv(file = paste(path.root, "/Data/Populations.csv", sep = ""))
	
	# Comprehensive set of occurences from Angert ENM paper
	all.pres <- read.csv("~/Google Drive/cardinalisENM/data files/all.records.aug.31.csv")
	
	# Select field-based occurences and herbarium records since 2000
	all.pres <- subset(all.pres[, c("ID", "DATASET", "Latitude", "Longitude")], 
		all.pres$PRESABS == 1 & all.pres$YEAR >= 2000)
	
	# Select correct DEM file
	foc.pres$demFile <- character(nrow(foc.pres))
	foc.pres$demFile[which(foc.pres$Lat > 30 & foc.pres$Lat < 35 & foc.pres$Lon > -120 & 
		foc.pres$Lon < -115)] <- "n30w120_dem_bil"
	foc.pres$demFile[which(foc.pres$Lat > 30 & foc.pres$Lat < 35 & foc.pres$Lon > -125 & 
		foc.pres$Lon < -120)] <- "n30w125_dem_bil"
	foc.pres$demFile[which(foc.pres$Lat > 35 & foc.pres$Lat < 40 & foc.pres$Lon > -120 & 
		foc.pres$Lon < -115)] <- "n35w120_dem_bil"
	foc.pres$demFile[which(foc.pres$Lat > 35 & foc.pres$Lat < 40 & foc.pres$Lon > -125 & 
		foc.pres$Lon < -120)] <- "n35w125_dem_bil"
	foc.pres$demFile[which(foc.pres$Lat > 40 & foc.pres$Lat < 45 & foc.pres$Lon > -120 & 
		foc.pres$Lon < -115)] <- "n40w120_dem_bil"
	foc.pres$demFile[which(foc.pres$Lat > 40 & foc.pres$Lat < 45 & foc.pres$Lon > -125 & 
		foc.pres$Lon < -120)] <- "n40w125_dem_bil"

	all.pres$demFile <- character(nrow(all.pres))
	all.pres$demFile[which(all.pres$Lat > 30 & all.pres$Lat < 35 & all.pres$Lon > -120 & 
		all.pres$Lon < -115)] <- "n30w120_dem_bil"
	all.pres$demFile[which(all.pres$Lat > 30 & all.pres$Lat < 35 & all.pres$Lon > -125 & 
		all.pres$Lon < -120)] <- "n30w125_dem_bil"
	all.pres$demFile[which(all.pres$Lat > 35 & all.pres$Lat < 40 & all.pres$Lon > -120 & 
		all.pres$Lon < -115)] <- "n35w120_dem_bil"
	all.pres$demFile[which(all.pres$Lat > 35 & all.pres$Lat < 40 & all.pres$Lon > -125 & 
		all.pres$Lon < -120)] <- "n35w125_dem_bil"
	all.pres$demFile[which(all.pres$Lat > 40 & all.pres$Lat < 45 & all.pres$Lon > -120 & 
		all.pres$Lon < -115)] <- "n40w120_dem_bil"
	all.pres$demFile[which(all.pres$Lat > 40 & all.pres$Lat < 45 & all.pres$Lon > -125 & 
		all.pres$Lon < -120)] <- "n40w125_dem_bil"

	coordinates(foc.pres) <- ~Lon+Lat
	projection(foc.pres) <- CRS('+proj=longlat');

	coordinates(all.pres) <- ~Longitude+Latitude
	projection(all.pres) <- CRS('+proj=longlat');

	foc.pres$Elev <- numeric(nrow(foc.pres))
	for (i in unique(foc.pres$demFile))
	{
		# load 90m DEM
		dem <- raster(paste("DEM/", i, "/", gsub("_bil", ".bil", i), sep = ""))
		
		# extract coordinate values from raster stack
		foc.pres$Elev[foc.pres$demFile == i] <- 
			extract(dem, foc.pres[foc.pres$demFile == i, ]) 
	}

	out1 <- data.frame(ID1 = foc.pres$Site, ID2 = foc.pres$ID, lat = foc.pres$Lat,
		long = foc.pres$Lon, elev = foc.pres$Elev)
	# For some reason, this has to be copied and saved in Excel after export to get ClimateWNA to work properly
	# See 'Data/Climate/ClimateWNA/Instructions for using ClimateWNA.txt' for directions
    # write.csv(out1, file = "climateWNA_input1.csv", row.names = FALSE, quote = F)

	all.pres$Elev <- numeric(nrow(all.pres))
	for (i in unique(all.pres$demFile))
	{
		# load 90m DEM
		dem <- raster(paste("DEM/", i, "/", gsub("_bil", ".bil", i), sep = "")) 
		
		# extract coordinate values from raster stack
		all.pres$Elev[all.pres$demFile == i] <- 
			extract(dem, all.pres[all.pres$demFile == i, ]) 
	}

	out2 <- data.frame(ID1 = all.pres$ID, ID2 = all.pres$DATASET, 
										 lat = all.pres$Latitude, long = all.pres$Longitude, 
										 elev = all.pres$Elev)
	# For some reason, this has to be copied and saved in Excel after export to get ClimateWNA to work properly
    # See 'Data/Climate/ClimateWNA/Instructions for using ClimateWNA.txt' for directions
    # write.csv(out2, file = "climateWNA_input2.csv", row.names = FALSE, quote = F)

#########################################################
# Export & run in climateWNA
	
# PART 2: Load time.series records & calculate biovars
# run climateWNA! external to this script.


	
	setwd()
	timeSer1 <- read.csv(file = 'climateWNA/climateWNA_input1_1981-2010MSYT.csv')
	timeSer2 <- read.csv(file = 'climateWNA/climateWNA_input2_1981-2010MSYT.csv')
  ClimateWNA_input1_1981-2010MSYT.csv
	timeSer1 <- ClimateWNA_3sites
	# Focal Populations
	
	tmax1 <- as.matrix(timeSer1[, paste("Tmax", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	tmin1 <- as.matrix(timeSer1[, paste("Tmin", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	prec1 <- as.matrix(timeSer1[, paste("PPT", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	bio1 <- biovars(prec1, tmin1, tmax1)

	# Add biovars to data.frame
	
	timeSer1 <- cbind(timeSer1, bio1)
	
	# Replace values of -9999 with NA
	
	timeSer1[timeSer1 == -9999] <- NA

	# Average by site for 30 year normal (1981 - 2010)

	yearmin <- 1981
	yearmax <- 2010
	
	timeSer1 <- subset(timeSer1, timeSer1$Year >= yearmin & timeSer1$Year <= yearmax)

	# Change all temperatures to Kelvin
	tempColumns <- c(grep("T[max|min|ave]", colnames(timeSer1)), which(colnames(timeSer1) 
		%in% c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11")))
	timeSer1[, colnames(timeSer1)[tempColumns]] <- timeSer1[, 
		colnames(timeSer1)[tempColumns]] + 273.16

	# 30-year mean
	foc.avg30 <- as.data.frame(apply(timeSer1[, 4:ncol(timeSer1)], 2, function(X) tapply(X,
		timeSer1$ID2, mean, na.rm = TRUE)))
	foc.avg30$PopID <- rownames(foc.avg30)
	
	# 30-year variance
	foc.var30 <- as.data.frame(apply(timeSer1[, 4:ncol(timeSer1)], 2, function(X) tapply(X, 
		timeSer1$ID2, var, na.rm = TRUE)))
	foc.var30[, c("Lat", "Long", "Elev")] <- foc.avg30[, c("Lat", 
		"Long", "Elev")]
	foc.var30$PopID <- rownames(foc.var30)

	# 30-year coefficient of variation
	foc.cva30 <- as.data.frame(apply(timeSer1[, 4:ncol(timeSer1)], 2, function(X) tapply(X, 
		timeSer1$ID2, function(X) sd(X, na.rm = TRUE) / mean(X, na.rm = TRUE))))
	foc.cva30[, c("Lat", "Long", "Elev")] <- foc.avg30[, c("Lat", 
		"Long", "Elev")]
	foc.cva30$PopID <- rownames(foc.cva30)

	# Adjust column names
	
	x <- which(colnames(foc.avg30) %in% c("Lat", "Long", "Elev", "PopID")) 
	colnames(foc.avg30)[-x] <- paste(colnames(foc.avg30)[-x], "_avg", sep = "")
	colnames(foc.var30)[-x] <- paste(colnames(foc.var30)[-x], "_var", sep = "")
	colnames(foc.cva30)[-x] <- paste(colnames(foc.cva30)[-x], "_cva", sep = "")

# Historical California Basin Characterization Model Downscaled Hydrology 

	# Function to get hydrological data from set of occurences
	getHydroVar <- function(hvar, occ, path = "~/Google Drive/BCM_monthly_climate_layers/")
	{
	
	# Get current working directory
	wd <- getwd()
	
	# Set working directory for hydrological data
	setwd(paste(path, hvar, sep = "/"))
	
	# Years 1981-2010, all months
	yrs <- as.character(1981:2010)
	mns <- as.character(1:12)
	mns[1:9] <- paste("0", mns[1:9], sep = "")
	
	# Occurences must be class SpatialPoints* with an ID column
	dat <- expand.grid(occ$ID, mns, yrs)
	colnames(dat) <- c("ID", "month", "year")
	dat[hvar] <- numeric(nrow(dat))
	
	# Extract data for each month and year combination for each population
	for (i in 1:(length(yrs) * length(mns)))
	{
		y <- yrs[ceiling(i / 12)]
		m <- mns[ifelse(i %% 12 != 0, i %% 12, 12)]
		r <- readRDS(sprintf("CA_BCM_HST_Monthly_%s_%s_%s.Rdata", hvar, y, m))
		occ <- spTransform(occ, CRS(r[[2]]$proj4string))
		dat[dat$year == y & dat$month == m, hvar] <- extract(r[[1]], occ)
		# Update status
		cat(hvar, "\t", m, ",", y, "(", 
				round(i / (length(yrs) * length(mns)) * 100, 1), "% )\n")
	}
	
	# Annual mean, variance, and coefficient of variation
	annDat <- as.data.frame(tapply(dat[, hvar], list(dat$ID, dat$year), sum))
	annDat[, paste(hvar, "_avg", sep = "")] <- apply(annDat[, yrs], 1, mean, 
																									 na.rm = T)
	annDat[, paste(hvar, "_var", sep = "")] <- apply(annDat[, yrs], 1, var, 
																									 na.rm = T)
	annDat[, paste(hvar, "_cva", sep = "")] <- apply(annDat[, yrs], 1, 
																									 function(X) sd(X, na.rm = TRUE) / mean(X, na.rm = TRUE))
	
	# Return
	col2ret <- paste(hvar, c("_avg", "_var", "_cva"), sep = "")
	return(annDat[, col2ret])
	
	# Reset working directory
	setwd(wd)
	
}

	# Focal populations - takes several minutes to run
	# foc.aet <- getHydroVar(hvar = "aet", occ = foc.pres)
	# foc.cwd <- getHydroVar(hvar = "cwd", occ = foc.pres)
	# foc.pet <- getHydroVar(hvar = "pet", occ = foc.pres)
	# foc.rch <- getHydroVar(hvar = "rch", occ = foc.pres)
	# foc.run <- getHydroVar(hvar = "run", occ = foc.pres)

	# hst.foc <- cbind(foc.aet, foc.cwd, foc.pet, foc.rch, foc.run)
	# write.csv(hst.foc, file = 'HST/hst.foc.csv', row.names = TRUE, quote = F)
	hst.foc <- read.csv('HST/hst.foc.csv', row.names = 1)

	if (all(row.names(foc.avg30) == row.names(hst.foc)))
	{
		foc.avg30 <- cbind(foc.avg30, hst.foc[, grep("_avg", colnames(hst.foc))])
		foc.var30 <- cbind(foc.var30, hst.foc[, grep("_var", colnames(hst.foc))])
		foc.cva30 <- cbind(foc.cva30, hst.foc[, grep("_cva", colnames(hst.foc))])
	}

	save(timeSer1, file = paste(path.obj, "/FocalOcc_timeSer.RData", sep = ""))
	save(foc.avg30, file = paste(path.obj, "/FocalOcc_avg1981-2010.RData", sep = ""))
	save(foc.var30, file = paste(path.obj, "/FocalOcc_var1981-2010.RData", sep = ""))
	save(foc.cva30, file = paste(path.obj, "/FocalOcc_cva1981-2010.RData", sep = ""))

	# All Populations
	
	tmax2 <- as.matrix(timeSer2[, paste("Tmax", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	tmin2 <- as.matrix(timeSer2[, paste("Tmin", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	prec2 <- as.matrix(timeSer2[, paste("PPT", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	bio2 <- biovars(prec2, tmin2, tmax2)

	# Add biovars to data.frame
	
	timeSer2 <- cbind(timeSer2, bio2)

	# Replace values of -9999 with NA
	
	timeSer2[timeSer2 == -9999] <- NA

	# Average by site for 30 year normal (1981 - 2010)

	yearmin <- 1981
	yearmax <- 2010

	timeSer2 <- subset(timeSer2, timeSer2$Year >= yearmin & timeSer2$Year <= yearmax)

	# Change all temperatures to Kelvin
	tempColumns <- c(grep("T[max|min|ave]", colnames(timeSer2)), which(colnames(timeSer2) 
		%in% c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11")))
	timeSer2[, colnames(timeSer2)[tempColumns]] <- timeSer2[, 
		colnames(timeSer2)[tempColumns]] + 273.16

	# 30-year mean
	all.avg30 <- as.data.frame(apply(timeSer2[, 4:ncol(timeSer2)], 2, function(X) tapply(X, 
		timeSer2$ID1, mean, na.rm = TRUE)))
	all.avg30$ID <- rownames(all.avg30)

	# 30-year variance
	all.var30 <- as.data.frame(apply(timeSer2[, 4:ncol(timeSer2)], 2, function(X) tapply(X, 
		timeSer2$ID1, var, na.rm = TRUE)))
	all.var30[, c("Latitude", "Longitude", "Elevation")] <- all.avg30[, c("Latitude", 
		"Longitude", "Elevation")]
	all.var30$ID <- rownames(all.var30)

	# 30-year coefficient of variation
	all.cva30 <- as.data.frame(apply(timeSer2[, 4:ncol(timeSer2)], 2, function(X) tapply(X, 
		timeSer2$ID1, function(X) sd(X, na.rm = TRUE) / mean(X, na.rm = TRUE))))
	all.cva30[, c("Latitude", "Longitude", "Elevation")] <- all.avg30[, c("Latitude", 
		"Longitude", "Elevation")]
	all.cva30$ID <- rownames(all.cva30)

	# Adjust column names
	
	x <- which(colnames(all.avg30) %in% c("Latitude", "Longitude", "Elevation", "ID")) 
	colnames(all.avg30)[-x] <- paste(colnames(all.avg30)[-x], "_avg", sep = "")
	colnames(all.var30)[-x] <- paste(colnames(all.var30)[-x], "_var", sep = "")
	colnames(all.cva30)[-x] <- paste(colnames(all.cva30)[-x], "_cva", sep = "")

# Historical California Basin Characterization Model Downscaled Hydrology 

	# All populations - takes several minutes to run
	# Remove population IDs (in ID factor) that are not used
	# all.pres1 <- all.pres
	# all.pres1$ID <- as.character(all.pres1$ID)
	# all.aet <- getHydroVar(hvar = "aet", occ = all.pres1)
	# all.cwd <- getHydroVar(hvar = "cwd", occ = all.pres1)
	# all.pet <- getHydroVar(hvar = "pet", occ = all.pres1)
	# all.rch <- getHydroVar(hvar = "rch", occ = all.pres1)
	# all.run <- getHydroVar(hvar = "run", occ = all.pres1)

	# hst.all <- cbind(all.aet, all.cwd, all.pet, all.rch, all.run)
	# write.csv(hst.all, file = 'HST/hst.all.csv', row.names = TRUE, quote = F)
	hst.all <- read.csv('HST/hst.all.csv', row.names = 1)

    if (all(row.names(all.avg30) %in% row.names(hst.all)))
    {
        x <- match(row.names(all.avg30), row.names(hst.all))    
        all.avg30 <- cbind(all.avg30, hst.all[x, grep("_avg", colnames(hst.all))])
        all.var30 <- cbind(all.var30, hst.all[x, grep("_var", colnames(hst.all))])
        all.cva30 <- cbind(all.cva30, hst.all[x, grep("_cva", colnames(hst.all))])
    }

    save(timeSer2, file = paste(path.obj, "/AllOcc_timeSer.RData", sep = ""))
	save(all.avg30, file = paste(path.obj, "/AllOcc_avg1981-2010.RData", sep = ""))
	save(all.var30, file = paste(path.obj, "/AllOcc_var1981-2010.RData", sep = ""))
	save(all.cva30, file = paste(path.obj, "/AllOcc_cva1981-2010.RData", sep = ""))

####################################################

# PART 3: Average climate within 62 km of focal populations, clipped to NHD waterways

	bufferSize <- 6.2e4 # Buffer size in meters
	# Combine separate DEM files into one (do not repeat)
	# dems <- c("n30w120_dem", "n30w125_dem", "n35w120_dem", "n35w125_dem", "n40w120_dem", 
		# "n40w125_dem")
	# dems <- paste("DEM/", dems, "_bil/", dems, ".bil", sep = "")
	# dem1 <- raster(dems[1])
	# for (i in 2:length(dems))
	# {
		# dem2 <- raster(dems[i])
		# dem1 <- mosaic(dem1, dem2, fun = mean)
	# }

	# writeRaster(dem1, filename = "DEM/dem_Lat30-45_Lon115-125.grd")
	
	# Start from here if you need to choose a different buffer than 10km
	# Otherwise do not rerun as this takes awhile
	dem <- raster("DEM/dem_Lat30-45_Lon115-125.grd")
  foc.pres <- read.csv(file = "~/Google Drive/CardAdapt/Data/Populations.csv")
  coordinates(foc.pres) <- ~Lon+Lat
	projection(foc.pres) <- CRS(projection(dem));

	# Raster of distance from nearest focal population
	# distFoc <- distanceFromPoints(dem, foc.pres)
	
	# Make all points greater than 62 km NA
	# distFoc[distFoc > bufferSize] <- NA
	# writeRaster(distFoc, filename = paste(path.obj, "/distFoc_62km.grd", sep = ""), overwrite = T)
	distFoc <- raster(paste(path.obj, "/distFoc_62km.grd", sep = ""))

	# Clip to streams

	# Get coordinates, elevation, and distance to nearest focal population from cells
	# within 62 km of focal populations
	# cells <- which(!is.na(values(distFoc)))
	# cellCoords <- xyFromCell(distFoc, cells, spatial = T)
	# saveRDS(cellCoords,  file = paste(path.obj, "/cellCoords_62km.rds", sep = ""))
	# cellCoords <- readRDS(file = paste(path.obj, "/cellCoords_62km.rds", sep = ""))
	
	# cellElevation <- extract(dem, cellCoords) # Takes ~6 min on MacBook Air
	# saveRDS(cellElevation,  file = paste(path.obj, "/cellElevation_62km.rds", sep = ""))
	# cellElevation <- readRDS(file = paste(path.obj, "/cellElevation_62km.rds", sep = ""))

	# for each cell, find closest population
	# pointDist <- pointDistance(cellCoords, foc.pres, lonlat = T) # takes ~10 min on MacBook Air
	# saveRDS(pointDist,  file = paste(path.obj, "/pointDist_62km.rds", sep = ""))
	# pointDist <- readRDS(file = paste(path.obj, "/pointDist_62km.rds", sep = ""))
	# x <- which(pointDist > bufferSize)
	# pointDist[x] <- NA
	
	# NEW VERSION:
	# extentByPop <- list()
	# for (i in 1:16) extentByPop[[i]] <- extent(cellCoords[which(!is.na(pointDist[, i]))]) #Takes a minute
	# saveRDS(extentByPop, file = paste(path.obj, "/extentByPop_62km.rds", sep = ""))
	extentByPop <- readRDS(file = paste(path.obj, "/extentByPop_62km.rds", sep = ""))

	# GET HUC_08 (needs to already be downloaded)
	# to download goto Terminal: 	
	# $> wget ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_08_usa.zip
  HUC8 <- rgdal::readOGR("NHD/huc_08_USA", layer = "huc_08_USA", verbose = FALSE)
	HUC8@proj4string <- sp::CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=WGS84")

	# Subset of HUC8 subbasins in CA and OR (HUC2 = {17, 18})
	subHUC8 <- HUC8[HUC8@data$HUC4 %in% c(1710, 1802:1810), ]
	
	# Code to roughly find out which subbasins need to be downloaded
	# This doesn't quite work because there is some discrepency between HUC8 codes in 
	# "NHD/huc_08_USA" and that in the WBD available from the National Map. This might be
	# resolved once I can unzip "WBD_National.zip"
	uniqueHUC8 <- character()
	for (i in 1:length(extentByPop))
	{
		r1 <- crop(distFoc, extentByPop[[i]])
		d1 <- distanceFromPoints(r1, foc.pres[i, ])
		d1[d1 > bufferSize] <- NA
		
		# pole O
		centroidCell <- which.min(values(d1))
	
		b1 <- boundaries(d1)
		b1coords <- xyFromCell(d1, which(values(b1) == 1))
		
		# Polar coordinates (r, phi) of boundary
		x <- b1coords[, "x"] - xyFromCell(d1, centroidCell)[, "x"]
		y <- b1coords[, "y"] - xyFromCell(d1, centroidCell)[, "y"]
		r <- sqrt(x ^ 2 + y ^ 2)
		phi <- atan2(y, x)
	
		# Order boundary by polar angle (= azimuth)
		b1coords <- b1coords[order(phi), ]
		b1coords <- rbind(b1coords, b1coords[1, ])

		# Make spatial polygon
		p1 <- SpatialPolygons(list(Polygons(list(Polygon(b1coords)), "p1")), 
			proj4string = CRS(projection(d1)))
		p1 <- spTransform(p1, CRS(projection(HUC8)))
	
		# Find out which HUC8 subbasins intersect with polygon around focal population
		whichHUC8 <- logical(length(subHUC8))
		for (j in 1:length(subHUC8))
		{
			try(whichHUC8[j] <- !is.null(intersect(extent(subHUC8[j, ]), extent(p1))))
		}
	
		uniqueHUC8 <- c(uniqueHUC8, subHUC8@data$HUC[whichHUC8])
	}
	
	alreadyHave <- list.dirs("NHD/SUBREGIONS")
	alreadyHave <- substr(alreadyHave[grepl("[0-9]{8}", alreadyHave)], 20, 27)
	sort(unique(uniqueHUC8[!(uniqueHUC8 %in% alreadyHave)]))

	# to download subregion data, goto http://viewer.nationalmap.gov/viewer/nhd.html?p=nhd
	# and download shape files of needed subregions based on HUC4/HUC8 codes from area.list
	
	# Graph to make sure that all focal populations + buffer are encompassed by NHD data
	# I have downloaded.
	
	# Get list of extents for available NHD shapefiles
	subregs <- list.dirs("NHD/SUBREGIONS", recursive = F)
	subregs <- subregs[grepl("[0-9]{8}", subregs)]
	HUC8 <- substr(subregs, 20, 27)

	extentByHUC8 <- list()
	for (i in 1:length(subregs))
	{
		wbd <- readShapeLines(paste(path.dat, "/", subregs[i], "/WBDHU8", sep = ""),
			proj4string = CRS(projection(distFoc)))
		extentByHUC8[[i]] <- extent(wbd[which(wbd@data$HUC8 == HUC8[i]), ])
	}	
	
	#  Which subregions needed for a given focal population + buffer
	pdf("NHD/Subregions_Pop.pdf", 5, 5)
	whichSubregs <- list()
	for (i in 1:length(extentByPop))
	{
		whichHUC8 <- logical(length(extentByHUC8))
		for (j in 1:length(extentByHUC8))
		{
			try(whichHUC8[j] <- !is.null(intersect(extentByHUC8[[j]], extentByPop[[i]])))
		}
		
		whichSubregs[[i]] <- subregs[whichHUC8]
		codes <- HUC8[whichHUC8]
		
		# Bind subregions and plot with focal population + buffer
		wbd <- readShapePoly(paste(path.dat, "/", whichSubregs[[i]][1], "/WBDHU8", sep = ""),
			proj4string = CRS(projection(distFoc)), IDvar = "HUC8")
		wbd <- wbd[wbd@data$HUC8 == codes[1], ]
		
		if (length(whichSubregs[[i]]) > 1)
		{
			for (k in 2:length(whichSubregs[[i]]))
			{
				x <- readShapePoly(paste(path.dat, "/", whichSubregs[[i]][k], "/WBDHU8", 
					sep = ""), proj4string = CRS(projection(distFoc)), IDvar = "HUC8")
				x <- x[x@data$HUC8 == codes[k], ]
				wbd <- rbind(wbd, x)
			}
		}
		
		# Make focal population i + buffer into polygon
		r1 <- crop(distFoc, extentByPop[[i]])
		d1 <- distanceFromPoints(r1, foc.pres[i, ])
		d1[d1 > bufferSize] <- NA
		
		# pole O
		centroidCell <- which.min(values(d1))
	
		b1 <- boundaries(d1)
		b1coords <- xyFromCell(d1, which(values(b1) == 1))
	
		# Polar coordinates (r, phi) of boundary
		x <- b1coords[, "x"] - xyFromCell(d1, centroidCell)[, "x"]
		y <- b1coords[, "y"] - xyFromCell(d1, centroidCell)[, "y"]
		r <- sqrt(x ^ 2 + y ^ 2)
		phi <- atan2(y, x)
	
		# Order boundary by polar angle (= azimuth)
		b1coords <- b1coords[order(phi), ]
		b1coords <- rbind(b1coords, b1coords[1, ])
		
		# Make spatial polygon
		p1 <- SpatialPolygons(list(Polygons(list(Polygon(b1coords)), "p1")), 
			proj4string = CRS(projection(r1)))

		plot(wbd, col = rgb(0, 0, 1, 0.5), main = foc.pres@data$Site[i])
		plot(p1, add = T, col = rgb(1, 0, 0, 0.5))
	}
	dev.off()

	# Now, import Flowline Shapefile, bind, rasterize, and save
	for (i in 1:length(foc.pres))
	{
		# Update and estimated time remaing
		t1 <- Sys.time()
		cat(as.character(foc.pres@data[i, "Site"]), "\t")
		
		# Make focal population i + buffer into polygon
		r1 <- crop(distFoc, extentByPop[[i]])
		d1 <- distanceFromPoints(r1, foc.pres[i, ])
		d1[d1 > bufferSize] <- NA
		
		# pole O
		centroidCell <- which.min(values(d1))
		
		b1 <- boundaries(d1)
		b1coords <- xyFromCell(d1, which(values(b1) == 1))
		
		# Polar coordinates (r, phi) of boundary
		x <- b1coords[, "x"] - xyFromCell(d1, centroidCell)[, "x"]
		y <- b1coords[, "y"] - xyFromCell(d1, centroidCell)[, "y"]
		r <- sqrt(x ^ 2 + y ^ 2)
		phi <- atan2(y, x)
		
		# Order boundary by polar angle (= azimuth)
		b1coords <- b1coords[order(phi), ]
		b1coords <- rbind(b1coords, b1coords[1, ])
		
		# Make spatial polygon
		p1 <- SpatialPolygons(list(Polygons(list(Polygon(b1coords)), "p1")), 
			proj4string = CRS(projection(r1)))

		# Bind subregion flowlines
		flowLine <- readShapeLines(paste(path.dat, "/", whichSubregs[[i]][1], "/NHDFlowline", 
			sep = ""), proj4string = CRS(projection(distFoc)))#, IDvar = "PERMANENT_")
		flowLine <- spChFIDs(flowLine, paste("pol", LETTERS[1], "_", 1:length(flowLine), 
			sep = ""))
		if (length(whichSubregs[[i]]) > 1)
		{
			for (k in 2:length(whichSubregs[[i]]))
			{
				x <- readShapeLines(paste(path.dat, "/", whichSubregs[[i]][k], "/NHDFlowline", 
					sep = ""), proj4string = CRS(projection(distFoc)))
				x <- spChFIDs(x, paste("pol", LETTERS[k], "_", 1:length(x), sep = ""))
				flowLine <- rbind(flowLine, x)
			}
		}
		
		# Rasterize and save
		tmp1 <- crop(distFoc, extentByPop[[i]])
		tmp1 <- distanceFromPoints(tmp1, foc.pres[1, ])
		tmp1[tmp1 > bufferSize] <- NA
		values(tmp1) <- !is.na(values(tmp1))
		tmp2 <- crop(dem, extentByPop[[i]])
		tmp3 <- overlay(tmp1, tmp2, fun = "*")
		values(tmp3) <- ifelse(values(tmp3) == 0, NA, values(tmp3))
		#plot(tmp3)
		rFlowline <- rasterize(flowLine, tmp3)
		writeRaster(rFlowline, filename = paste("NHD/Flowline_rasters_62km_buffer/", 
																						foc.pres@data$Site[i], sep = ""), overwrite = T)
		
		# Update and estimated time remaing
		t2 <- Sys.time()
		z <- difftime(t2, t1, units = "hours")[1] * (length(foc.pres) - i)
		cat("Time remaining:", format(z), "\n")
	}

	# Crop Flowline rasters and extract elevation and other information

	n <- 16000 # Number of population x 1e3 points per species
	cellCoords <- data.frame(long = numeric(n), lat = numeric(n))
	cellElevation <- numeric(n)
	cellDist <- numeric(n)
	closestPop <- character(n)

	pdf("DEM_Clipped2Stream.pdf", 5, 5)
	set.seed(746401)
	for (i in 1:length(foc.pres))
	{
		tmp1 <- crop(distFoc, extentByPop[[i]])
		tmp1 <- tmp5 <- distanceFromPoints(tmp1, foc.pres[i, ])
		tmp1[tmp1 > bufferSize] <- NA
		values(tmp1) <- !is.na(values(tmp1))
		tmp2 <- crop(dem, extentByPop[[i]])
		tmp3 <- overlay(tmp1, tmp2, fun = "*")
		values(tmp3) <- ifelse(values(tmp3) == 0, NA, values(tmp3))
		# plot(tmp3)
		rFlowline <- raster(paste("NHD/Flowline_rasters_62km_buffer/", foc.pres@data$Site[i], sep = ""))
		rFlowCrop <- crop(rFlowline, extent(tmp3))
		# plot(rFlowCrop)
		values(rFlowCrop) <- !is.na(values(rFlowCrop))
		tmp4 <- overlay(tmp3, rFlowCrop, fun = "*")
		values(tmp4) <- ifelse(values(tmp4) == 0, NA, values(tmp4))
		
		# Plot
		plot(tmp4, main = foc.pres@data$Site[i], 
				 col = colorRampPalette(c("tomato", "steelblue"))(100))

		# Get coordinates, elevation, and distance to nearest focal population from cells
		# within 100 km of focal populations, clipped to streams
		cells <- which(!is.na(values(tmp4)))

		# sample 1000 points within each population
		cells <- cells[sample(length(cells), 1e3)]
		cellCoordsSp <- xyFromCell(tmp4, cells, spatial = T)
		cellCoords[(1e3 * (i - 1) + 1):(i * 1e3), ] <- as.data.frame(cellCoordsSp)
		cellElevation[(1e3 * (i - 1) + 1):(i * 1e3)] <- extract(tmp4, cellCoordsSp)
		cellDist[(1e3 * (i - 1) + 1):(i * 1e3)] <- extract(tmp5, cellCoordsSp)
		closestPop[(1e3 * (i - 1) + 1):(i * 1e3)] <- rep(as.character(foc.pres@data$Site[i]), 
			length(cells))
		
	}
	dev.off()
	
	out3 <- data.frame(ID1 = closestPop, ID2 = cellDist, lat = cellCoords$lat, 
		long = cellCoords$long, elev = cellElevation)
	# For some reason, this has to be copied and saved in Excel after export to get 
	# ClimateWNA to work properly
  # See 'Data/Climate/ClimateWNA/Instructions for using ClimateWNA.txt' for directions
  # write.csv(out3, file = "ClimateWNA/climateWNA_input3_62km.csv", row.names = FALSE, quote = F)

#########################################################
# Export & run in climateWNA
	
	timeSer3 <- read.csv(file = 'climateWNA/climateWNA_input3_62km_1981-2010MSYT.csv')

	# Change ID2 (NOTE: this would be better done before exporting for climateWNA)
	timeSer3$ID2 <- rep(1:16000, 30)
	
	# calculate biovars and export
	
	tmax3 <- as.matrix(timeSer3[, paste("Tmax", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	tmin3 <- as.matrix(timeSer3[, paste("Tmin", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	prec3 <- as.matrix(timeSer3[, paste("PPT", paste(ifelse(1:12 < 10, "0", ""), 
		as.character(1:12), sep = ""), sep = "")], ncol = 12)

	bio3 <- biovars(prec3, tmin3, tmax3)

	# Add biovars to data.frame
	
	timeSer3 <- cbind(timeSer3, bio3)
	
	# Replace values of -9999 with NA
	
	timeSer3[timeSer3 == -9999] <- NA

	# Average by site for 30 year normal (1981 - 2010)

	yearmin <- 1981
	yearmax <- 2010

	# Change all temperatures to Kelvin
	tempColumns <- c(grep("T[max|min|ave]", colnames(timeSer3)), which(colnames(timeSer3) 
		%in% c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11")))
	timeSer3[, colnames(timeSer3)[tempColumns]] <- timeSer3[, 
		colnames(timeSer3)[tempColumns]] + 273.16
	timeSer3 <- subset(timeSer3, timeSer3$Year >= yearmin & timeSer3$Year <= yearmax)
	
	tmp <- read.csv(paste(path.root, "/Data/Populations.csv", sep = ""))
	timeSer3$PopID <- tmp$ID[match(timeSer3$ID1, tmp$Site)]

	# SA = spatial average
	# 30-year mean
	focSA.avg30 <- as.data.frame(apply(timeSer3[, 4:(ncol(timeSer3) - 1)], 2, function(X) tapply(X,
		timeSer3$ID2, mean, na.rm = TRUE)))
	
	# 30-year variance
	focSA.var30 <- as.data.frame(apply(timeSer3[, 4:(ncol(timeSer3) - 1)], 2, function(X) tapply(X, 
		timeSer3$ID2, var, na.rm = TRUE)))
	focSA.var30[, c("Latitude", "Longitude", "Elevation")] <- focSA.avg30[, c("Latitude", 
		"Longitude", "Elevation")]

	# 30-year coefficient of variation
	focSA.cva30 <- as.data.frame(apply(timeSer3[, 4:(ncol(timeSer3) - 1)], 2, function(X) tapply(X, 
		timeSer3$ID2, function(X) sd(X, na.rm = TRUE) / mean(X, na.rm = TRUE))))

for (i in 4:ncol(focSA.cva30)) {

	# Make Inf's NA
		focSA.cva30[is.infinite(focSA.cva30[, i]), i] <- NA
		
	# Replace NA's with 0's for cases where all values are 0
	if (any(is.na(focSA.cva30[, i]), na.rm = T) & any(focSA.avg30[, i] == 0, na.rm = T)) {
			focSA.cva30[which(focSA.avg30[, i] == 0), i] <- 0

		}

	}
	focSA.cva30[, c("Latitude", "Longitude", "Elevation")] <- focSA.avg30[, c("Latitude", 
		"Longitude", "Elevation")]

	focSA.avg30$PopID <- timeSer3$PopID[1:16000]
	# focSA.var30$PopID <- timeSer3$PopID[1:16000]
	focSA.cva30$PopID <- timeSer3$PopID[1:16000]

	# Import predicted suitability to weight points for calculating spatially-integrated climate
	# some biovars need to be converted from K to C, divided by 10, and log-transformed
	# bio2, bio3, bio4, bio10, bio11, bio12, bio14, bio15
	# log(x + 0.5) for bio3, bio10, bio12, bio14

  setwd("~/Google Drive/cardinalisENM/R objects/")

  # Load fitted models
  for (i in 1:10) {
    mod.lr = get(load(paste("LR.mod2.",i,".pseudo11.Rda", sep="")))
    assign(paste("LR.mod2.",i, sep=""), mod.lr)
    mod.gam = get(load(paste("GAM.mod4.",i,".pseudo11.Rda", sep="")))
    assign(paste("GAM.mod4.",i, sep=""), mod.gam)
    mod.rf = get(load(paste("RF.mod1.",i,".pseudo11.Rda", sep="")))  
    assign(paste("RF.mod1.",i, sep=""), mod.rf)
    mod.brt = get(load(paste("BRT.mod4.",i,".pseudo11.Rda", sep="")))
    assign(paste("BRT.mod4.",i, sep=""), mod.brt)  
    mod.max = get(load(paste("MAX.mod1.",i,".pseudo11.Rda", sep="")))
    assign(paste("MAX.mod1.",i, sep=""), mod.max)  
    }

  ext.accs.lr = get(load("LR.mod2.extaccs.pseudo11.Rda"))
  ext.accs.gam = get(load("GAM.mod4.extaccs.pseudo11.Rda"))
  ext.accs.rf = get(load("RF.mod1.extaccs.pseudo11.Rda"))
  ext.accs.brt = get(load("BRT.mod4.extaccs.pseudo11.Rda"))
  ext.accs.max = get(load("MAX.mod1.extaccs.pseudo11.Rda"))

  lr.cuts = ext.accs.lr[ext.accs.lr$thresh=="SensSpec", "threshold"]
  gam.cuts = ext.accs.gam[ext.accs.gam$thresh=="SensSpec", "threshold"]
  rf.cuts = ext.accs.rf[ext.accs.rf$thresh=="SensSpec", "threshold"]
  brt.cuts = ext.accs.brt[ext.accs.brt$thresh=="SensSpec", "threshold"]
  max.cuts = ext.accs.max[ext.accs.max$thresh=="SensSpec", "threshold"]

  cuts <- cbind(lr.cuts, gam.cuts, rf.cuts, brt.cuts, max.cuts)

  # Bioclimatic variables (30-year normal 1981 - 2010) for points are in focSA.avg30
  pred <- focSA.avg30[, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]

  # Convert bio10 and bio1 from K to C
  pred$bio10 <- pred$bio10 - 273.16
  pred$bio11 <- pred$bio11 - 273.16

  # log-transform variables, consistent with ENM models
  pred$bio3 <- log(pred$bio3 + 0.5)
  pred$bio10 <- log(pred$bio10 + 0.5)
  pred$bio12 <- log(pred$bio12 + 0.5)
  pred$bio14 <- log(pred$bio14 + 0.5)

  # data.frame for suitability scores
  ENMprob <- data.frame(matrix(NA, ncol = 50, nrow = nrow(pred)))
  colnames(ENMprob) <- c(paste("LRprob", 1:10, sep = ""), 
    paste("GAMprob", 1:10, sep = ""), paste("RFprob", 1:10, sep = ""), 
    paste("BRTprob", 1:10, sep = ""), paste("MAXprob", 1:10, sep = ""))

  for (i in 1:10) {
    ## LR predictions
    mod = get(paste("LR.mod2.", i, sep=""))
    ENMprob[, i] = predict(mod, pred, type = "response")

    ## GAM predictions
    mod = get(paste("GAM.mod4.", i, sep = ""))
    ENMprob[, 10 + i] = predict(mod, pred, type = "response")

    ## RF predictions
    mod = get(paste("RF.mod1.", i, sep = ""))
    ENMprob[, 20 + i] = predict(mod, pred, type = "prob", fun = predict, index = 2, overwrite = T)[, 2]
 
    ## BRT predictions
    mod = get(paste("BRT.mod4.", i, sep = ""))
    ENMprob[, 30 + i] = predict(mod, pred, n.trees = mod$gbm.call$best.trees, type = "response") 
 
    ## MAX predictions
    mod = get(paste("MAX.mod1.", i, sep = ""))
    ENMprob[, 40 + i] = predict(mod, pred, overwrite = T) 
    
  }   

  # Calculate ensemble average of probabilities
  ENMprob$ensAvgProb <- apply(ENMprob, 1, mean)

	# Save occurence probabilities
  save(ENMprob, file = paste(path.obj, "/ENMprob_62km.RData", sep = ""))

  # Calculate spatial average climate (weighted by suitability)
	# Do for all climate columns
	cols <- which(!colnames(focSA.avg30) %in% c("Latitude", "Longitude", "Elevation", "PopID"))
	focSA.avg <- focSA.cva <- data.frame(matrix(rep(NA, 16 * ncol(focSA.avg30)), nrow = 16))
	colnames(focSA.avg) <- colnames(focSA.avg30)
	colnames(focSA.cva) <- colnames(focSA.cva30)
	focSA.avg$PopID <- rownames(focSA.avg) <- levels(focSA.avg30$PopID)
	focSA.cva$PopID <- rownames(focSA.cva) <- levels(focSA.cva30$PopID)
	for (i in cols)
	{
		for (j in levels(focSA.avg30$PopID)) 
  	{
			focSA.avg[j, i] <- weighted.mean(focSA.avg30[focSA.avg30$PopID == j, i], 
				ENMprob$ensAvgProb[focSA.avg30$PopID == j])
			focSA.cva[j, i] <- weighted.mean(focSA.cva30[focSA.cva30$PopID == j, i], 
				ENMprob$ensAvgProb[focSA.cva30$PopID == j])
		}
	}
	
	# Adjust column names and export
	# Big data.frame with all 16000 points
	x <- which(colnames(focSA.avg30) %in% c("Latitude", "Longitude", "Elevation", "PopID")) 
	colnames(focSA.avg30)[-x] <- paste(colnames(focSA.avg30)[-x], "_avg", sep = "")
	colnames(focSA.var30)[-x] <- paste(colnames(focSA.var30)[-x], "_var", sep = "")
	colnames(focSA.cva30)[-x] <- paste(colnames(focSA.cva30)[-x], "_cva", sep = "")

	save(timeSer3, file = paste(path.obj, "/FocalOcc1000_timeSer.RData", sep = ""))
	save(focSA.avg30, file = paste(path.obj, "/FocalOcc1000_avg1981-2010.RData", sep = ""))
	save(focSA.var30, file = paste(path.obj, "/FocalOcc1000_var1981-2010.RData", sep = ""))
	save(focSA.cva30, file = paste(path.obj, "/FocalOcc1000_cva1981-2010.RData", sep = ""))

	# Small data.frame with all 16 points from focal populations
	x <- which(colnames(focSA.avg) %in% c("Latitude", "Longitude", "Elevation", "PopID")) 
	colnames(focSA.avg)[-x] <- paste(colnames(focSA.avg)[-x], "_avg", sep = "")
	colnames(focSA.cva)[-x] <- paste(colnames(focSA.cva)[-x], "_cva", sep = "")	

	save(focSA.avg, file = paste(path.obj, "/FocalOccSA_avg1981-2010_62km.RData", sep = ""))
	save(focSA.cva, file = paste(path.obj, "/FocalOccSA_cva1981-2010_62km.RData", sep = ""))
	
	####################################################
	
	# PART 4: Average climate within 62 km of focal populations, clipped to NHD waterways
	
	bufferSize <- 6.2e4 # Buffer size in meters
	# Combine separate DEM files into one (do not repeat)
	# dems <- c("n30w120_dem", "n30w125_dem", "n35w120_dem", "n35w125_dem", "n40w120_dem", 
	# "n40w125_dem")
	# dems <- paste("DEM/", dems, "_bil/", dems, ".bil", sep = "")
	# dem1 <- raster(dems[1])
	# for (i in 2:length(dems))
	# {
	# dem2 <- raster(dems[i])
	# dem1 <- mosaic(dem1, dem2, fun = mean)
	# }
	
	# writeRaster(dem1, filename = "DEM/dem_Lat30-45_Lon115-125.grd")
	
	# Start from here if you need to choose a different buffer than 10km
	# Otherwise do not rerun as this takes awhile
	dem <- raster("DEM/dem_Lat30-45_Lon115-125.grd")

	# Comprehensive set of occurences from Angert ENM paper
	all.pres <- read.csv("~/Google Drive/cardinalisENM/data files/all.records.aug.31.csv")
	# all.pres <- read.csv("~/Documents/all.records.aug.31.csv") # Amy's computer
	
	# Select field-based occurences and herbarium records since 2000
	all.pres <- subset(all.pres[, c("ID", "DATASET", "Latitude", "Longitude")], 
	                   all.pres$PRESABS == 1 & all.pres$YEAR >= 2000)
	coordinates(all.pres) <- ~Longitude+Latitude
	projection(all.pres) <- CRS(projection(dem));
	
	# Crop 62-km circle around every occurence
	# circleCrops <- mclapply(1:nrow(all.pres), function(row) {
	#  makeCircleCrop(all.pres@coords[row, ], 
	#                 r = dem, bufferSize = bufferSize)
	# }, mc.cores = 3)
	# circleCrops <- lapply(circleCrops, cropNA) # remove extra extent
  # saveRDS(circleCrops, paste0(path.obj, "/circleCrops.rds"))
	circleCrops <- readRDS(paste0(path.obj, "/circleCrops.rds"))
	
	# Graph to make sure that all focal populations + buffer are encompassed by NHD data
	# I have downloaded.
	
	# Get list of extents for available NHD shapefiles
	subregs <- list.dirs("NHD/SUBREGIONS", recursive = F)
	subregs <- subregs[grepl("[0-9]{8}", subregs)]
	HUC8 <- substr(subregs, 20, 27)
	anyDuplicated(HUC8) # should be 0
	
	extentByPop <- lapply(circleCrops, extent)
	extentByHUC8 <- list()
	for (i in 1:length(subregs))
	{
	  wbd <- readShapeLines(paste(path.dat, "/", subregs[i], "/WBDHU8", sep = ""),
	                        proj4string = CRS(projection(dem)))
	  extentByHUC8[[i]] <- extent(wbd[which(wbd@data$HUC8 == HUC8[i]), ])
	}	
	
	#  Which subregions needed for a given focal population + buffer
	pdf("NHD/Subregions_Pop.pdf", 5, 5)
	whichSubregs <- list()
	needsMoreNHD <- logical(length(extentByPop)) # get rid of once all are download
	for (i in 1:length(extentByPop))
	  {
	  whichHUC8 <- logical(length(extentByHUC8))
	  for (j in 1:length(extentByHUC8))
	  {
	    try(whichHUC8[j] <- !is.null(intersect(extentByHUC8[[j]], extentByPop[[i]])))
	  }
	  
	  whichSubregs[[i]] <- subregs[whichHUC8]
	  codes <- HUC8[whichHUC8]
	  
	  # Bind subregions and plot with focal population + buffer
	  wbd <- readShapePoly(paste(path.dat, "/", whichSubregs[[i]][1], "/WBDHU8", sep = ""),
	                       proj4string = CRS(projection(dem)), IDvar = "HUC8")
	  wbd <- wbd[wbd@data$HUC8 == codes[1], ]
	  
	  if (length(whichSubregs[[i]]) > 1)
	  {
	    for (k in 2:length(whichSubregs[[i]]))
	    {
	      x <- readShapePoly(paste(path.dat, "/", whichSubregs[[i]][k], "/WBDHU8", 
	                               sep = ""), proj4string = CRS(projection(dem)), IDvar = "HUC8")
	      x <- x[x@data$HUC8 == codes[k], ]
	      # wbd <- rbind(wbd, x[names(wbd)])
	      wbd <- rbind(wbd["HUC8"], x["HUC8"])
	    }
	  }
	  
	  # Dissolve borders
	  wbd.dissolve <- gUnaryUnion(wbd)

	  # Make focal population i + buffer into polygon
	  r1 <- circleCrops[[i]]
	  d1 <- distanceFromPoints(r1, all.pres[i, ])
	  d1[d1 > bufferSize] <- NA
	  d1[is.na(r1)] <- NA

	  # pole O
	  centroidCell <- which.min(values(d1))
	  
	  b1 <- boundaries(d1, type = "outer")
	  b1coords <- xyFromCell(d1, which(values(b1) == 1))
	  
	  # Polar coordinates (r, phi) of boundary
	  x <- b1coords[, "x"] - xyFromCell(d1, centroidCell)[, "x"]
	  y <- b1coords[, "y"] - xyFromCell(d1, centroidCell)[, "y"]
	  r <- sqrt(x ^ 2 + y ^ 2)
	  phi <- atan2(y, x)
	  
	  # Order boundary by polar angle (= azimuth)
	  b1coords <- b1coords[order(phi), ]
	  b1coords <- rbind(b1coords, b1coords[1, ])
	  
	  # Make spatial polygon
	  p1 <- SpatialPolygons(list(Polygons(list(Polygon(b1coords)), "p1")), 
	                        proj4string = CRS(projection(r1)))
	  
	  plot(wbd.dissolve, col = rgb(0, 0, 1, 0.5), main = all.pres@data$Site[i])
	  centroids <- getSpPPolygonsLabptSlots(wbd)
	  text(centroids, labels = wbd$HUC8)
	  plot(p1, add = T, col = rgb(1, 0, 0, 0.5))
	  axis(1); axis(2, las = 1)
	  needsMoreNHD[i] <- !gWithin(p1, wbd.dissolve)
	}
	dev.off()

	# saveRDS(whichSubregs, paste0(path.obj, "/whichSubregs.rds"))
	whichSubregs <- readRDS(paste0(path.obj, "/whichSubregs.rds"))
	
	# Now, import Flowline Shapefile, bind, rasterize, and save (this takes DAYS!)
	for (i in 1:length(extentByPop)) {
	  # Update and estimated time remaing
	  t1 <- Sys.time()
	  cat(i, format(Sys.time(), "%X"), "\n")
	  
	  # Make focal population i + buffer into polygon
	  r1 <- circleCrops[[i]]
	  d1 <- distanceFromPoints(r1, all.pres[i, ])
	  d1[d1 > bufferSize] <- NA
	  d1[is.na(r1)] <- NA
	  
	  # pole O
	  centroidCell <- which.min(values(d1))
	  
	  b1 <- boundaries(d1)
	  b1coords <- xyFromCell(d1, which(values(b1) == 1))
	  
	  # Polar coordinates (r, phi) of boundary
	  x <- b1coords[, "x"] - xyFromCell(d1, centroidCell)[, "x"]
	  y <- b1coords[, "y"] - xyFromCell(d1, centroidCell)[, "y"]
	  r <- sqrt(x ^ 2 + y ^ 2)
	  phi <- atan2(y, x)
	  
	  # Order boundary by polar angle (= azimuth)
	  b1coords <- b1coords[order(phi), ]
	  b1coords <- rbind(b1coords, b1coords[1, ])
	  
	  # Make spatial polygon
	  p1 <- SpatialPolygons(list(Polygons(list(Polygon(b1coords)), "p1")), 
	                        proj4string = CRS(projection(r1)))
	  
	  # Bind subregion flowlines
	  flowLine <- readShapeLines(paste(path.dat, "/", whichSubregs[[i]][1], "/NHDFlowline", 
	                                   sep = ""), proj4string = CRS(projection(dem)))#, IDvar = "PERMANENT_")
	  flowLine <- spChFIDs(flowLine, paste("pol", LETTERS[1], "_", 1:length(flowLine), 
	                                       sep = ""))
	  if (length(whichSubregs[[i]]) > 1)
	  {
	    for (k in 2:length(whichSubregs[[i]]))
	    {
	      x <- readShapeLines(paste(path.dat, "/", whichSubregs[[i]][k], "/NHDFlowline", 
	                                sep = ""), proj4string = CRS(projection(dem)))
	      x <- spChFIDs(x, paste("pol", LETTERS[k], "_", 1:length(x), sep = ""))
	      n <- intersect(names(x), names(flowLine)) # only bind column names that are present in all spatial objects
	      flowLine <- rbind(flowLine[n], x[n])
	    }
	  }
	  
	  # Rasterize and save
	  tmp1 <- r1
	  values(tmp1) <- !is.na(values(tmp1))
	  tmp <- crop(flowLine, extent(r1))
	  # system.time(rFlowline <- rasterize(flowLine, r1)) # this takes a long time, but need for allpres[230] for some reason
	  system.time(rFlowline <- rasterize(tmp, r1)) # 2X as fast
    tmp2 <- overlay(tmp1, rFlowline, fun = function(x, y) x * y)
    values(tmp2) <- ifelse(values(tmp2) == 0, NA, values(tmp2))
	  # writeRaster(tmp2, filename = paste0("NHD/Flowline_rasters_62km_buffer/allpres", i), 
	  #             overwrite = T)
    writeRaster(tmp2, filename = paste0("~/Google Drive/AmyComputer/Flowline_rasters_62km_buffer/allpres", i), 
                overwrite = T)

	  # Update and estimated time remaing
	  t2 <- Sys.time()
	  z <- difftime(t2, t1, units = "hours")[1] * (length(all.pres) - i)
	  cat("Time remaining:", format(z), "\n")
	}
	
	x <- list.files(paste0("~/Google Drive/AmyComputer/Flowline_rasters_62km_buffer/"))
	x <- gsub("[[:alpha:]]", "", x)
	x <- gsub("[.]", "", x)
	x <- as.numeric(x)
	all(table(x) == 2) # should be true
	all(1:358 %in% x) # should be true

	# Crop Flowline rasters and extract elevation and other information
	n <- length(all.pres) * 1e3
	cellCoords <- data.frame(long = numeric(n), lat = numeric(n))
	cellElevation <- numeric(n)
	closestPop <- character(n)
	
	pdf("DEM_Clipped2Stream_allpres_62km.pdf", 5, 5)
	set.seed(469727)
	for (i in 1:length(all.pres)) {
	  # Import saved flowline raster
	  # r  <- raster(paste0("NHD/Flowline_rasters_62km_buffer/allpres", i))
	  r <- raster(paste0("~/Google Drive/AmyComputer/Flowline_rasters_62km_buffer/allpres", i))
	  values(r) <- ifelse(values(r) == 0, NA, values(r))
	  values(r) <- ifelse(is.na(values(r)), NA, 1)
	  
	  r1 <- circleCrops[[i]]
	  
	  r2 <- overlay(r, r1, fun = function(x, y) x * y)
	  
	  # Plot
	  plot(r2, main = i, col = colorRampPalette(c("tomato", "steelblue"))(100))
	  
	  # Get coordinates, elevation, and distance to nearest focal population from cells
	  # within 100 km of focal populations, clipped to streams
	  cells <- which(!is.na(values(r2)))
	  
	  # sample 1000 points within each population
	  cells <- cells[sample(length(cells), 1e3)]
	  cellCoordsSp <- xyFromCell(r2, cells, spatial = T)
	  cellCoords[(1e3 * (i - 1) + 1):(i * 1e3), ] <- as.data.frame(cellCoordsSp)
	  cellElevation[(1e3 * (i - 1) + 1):(i * 1e3)] <- extract(r2, cellCoordsSp)
	  closestPop[(1e3 * (i - 1) + 1):(i * 1e3)] <- rep(i, 1e3)
	  
	}
	
	
	out4 <- data.frame(ID1 = closestPop, ID2 = rep(1:1000, length(all.pres)), 
	                   lat = cellCoords$lat, long = cellCoords$long, elev = cellElevation)
	# For some reason, this has to be copied and saved in Excel after export to get 
	# ClimateWNA to work properly
	# See 'Data/Climate/ClimateWNA/Instructions for using ClimateWNA.txt' for directions
	write.csv(out4, file = "ClimateWNA/climateWNA_input4_62km.csv", row.names = FALSE, quote = F)

	# Break up into smaller parts for ClimateWNA; recombine later
	for (i in 1:36) {
	  start <- 1e4 * (i - 1) + 1
	  end <- 1e4 * i
	  if (end > nrow(out4)) end <- nrow(out4)
	  write.csv(out4[start:end, ], row.names = F,
	            file = sprintf("ClimateWNA/climateWNA_input4_62km_part%s.csv", i))
	}
	
#########################################################
# Export & run in climateWNA

	# Recombine (takes a few hours)
	# n <- length(all.pres) * 1e3
	# tmp <- read.csv(file = 'climateWNA/climateWNA_input4_62km_part1_1981-2010MSYT.csv')
	# timeSer4 <- data.frame(matrix(NA, ncol = ncol(tmp), nrow = length(all.pres) * 3e4))
	# colnames(timeSer4) <- colnames(tmp)
	# start <- seq(1, n * 30, 3e5)
	# end <- seq(3e5, n * 30, 3e5)
	# end[36] <- length(all.pres) * 3e4
	# for (i in 1:36) {
	#  cat(i, "\n")
	#  tmp <- read.csv(file = sprintf('climateWNA/climateWNA_input4_62km_part%s_1981-2010MSYT.csv', i))
	#  timeSer4[start[i]:end[i], ] <- tmp
	# }
	
	# saveRDS(timeSer4, file = paste0(path.obj, "/timeSer4.rds"))
  # timeSer4 <- readRDS(file = paste0(path.obj, "/timeSer4.rds"))
  
	# Change ID2 (NOTE: this would be better done before exporting for climateWNA)
	# for (i in 1:35) {
	#  timeSer4$ID2[(3e5 * (i - 1) + 1):(3e5 * i)] <- rep((1e4 * (i - 1) + 1):(1e4 * i), 30)
	# }
	#	i <- 36
	# timeSer4$ID2[(3e5 * (i - 1) + 1):(358 * 3e4)] <- rep((1e4 * (i - 1) + 1):(1e4 * i - 2e3), 30)
	
	# calculate biovars and export
	
	# tmax4 <- as.matrix(timeSer4[, paste("Tmax", paste(ifelse(1:12 < 10, "0", ""), 
	#                                                  as.character(1:12), sep = ""), sep = "")], ncol = 12)
	
	# tmin4 <- as.matrix(timeSer4[, paste("Tmin", paste(ifelse(1:12 < 10, "0", ""), 
	#                                                  as.character(1:12), sep = ""), sep = "")], ncol = 12)
	
	# prec4 <- as.matrix(timeSer4[, paste("PPT", paste(ifelse(1:12 < 10, "0", ""), 
	#                                                 as.character(1:12), sep = ""), sep = "")], ncol = 12)
	
	# bio4 <- biovars(prec4, tmin4, tmax4)
	
	# Add biovars to data.frame
	# rm(list = c("tmax4", "tmin4", "prec4"))
	# timeSer4 <- cbind(timeSer4, bio4)
  # rm("bio4")

	# Replace values of -9999 with NA
	# timeSer4[timeSer4 == -9999] <- NA
	# saveRDS(timeSer4, file = paste0(path.obj, "/timeSer4.rds"))
	timeSer4 <- readRDS(file = paste0(path.obj, "/timeSer4.rds"))
	
	# Average by site for 30 year normal (1981 - 2010)
	# Change all temperatures to Kelvin
	tempColumns <- c(grep("T[max|min|ave]", colnames(timeSer4)), 
	                 which(colnames(timeSer4) %in% c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11")))
	timeSer4[, colnames(timeSer4)[tempColumns]] <- 
	  timeSer4[, colnames(timeSer4)[tempColumns]] + 273.16

	# SA = spatial average
	# 30-year mean
	n <- length(all.pres) * 1e3
	start <- seq(1, n, 1e4)
	end <- seq(1e4, n, 1e4)
  end[36] <- n
  
	for (i in 1:length(start)) {
  	x <- timeSer4$ID2 %in% start[i]:end[i]
	  allSA.avg30 <- as.data.frame(apply(timeSer4[x, 4:ncol(timeSer4)], 2, 
	                                     function(X) tapply(X, timeSer4$ID2[x], mean, na.rm = TRUE)))
	  saveRDS(allSA.avg30, file = paste0(path.obj, sprintf("/allSA.avg30.%s.rds", i)))
	}
	
	# Recombine, erase small files at end, then save large file
	allSA.avg30 <- readRDS(file = paste0(path.obj, "/allSA.avg30.1.rds"))
	allSA.avg30 <- data.frame(matrix(NA, nrow = n, ncol = ncol(allSA.avg30)))
	for (i in 1:length(start)) {
	  tmp <- readRDS(file = paste0(path.obj, sprintf("/allSA.avg30.%s.rds", i)))
	  allSA.avg30[start[i]:end[i], ] <- tmp
	  file.remove(file = paste0(path.obj, sprintf("/allSA.avg30.%s.rds", i)))
	  rm("tmp")
	}
	
	colnames(allSA.avg30) <- colnames(timeSer4)[4:ncol(timeSer4)]
	saveRDS(allSA.avg30, file = paste0(path.obj, "/allSA.avg30.rds"))
	allSA.avg30 <- readRDS(file = paste0(path.obj, "/allSA.avg30.rds"))
	
	# Check that 100 random entries give correct answer
	out <- logical(100)
	for (i in 1:100) {
	  col <- sample(colnames(allSA.avg30), 1)
	  row <- sample(1:nrow(allSA.avg30), 1)
  	out[i] <- mean(subset(timeSer4[, col], timeSer4$ID2 == i), na.rm = T) == allSA.avg30[i, col]
	}
	all(out) # should be true
	
	# 30-year coefficient of variation
	cva <- function(x, ...) sd(x, ...) / mean(x, ...)
	for (i in 1:length(start)) {
	  x <- timeSer4$ID2 %in% start[i]:end[i]
	  allSA.cva30 <- as.data.frame(apply(timeSer4[x, 4:ncol(timeSer4)], 2, 
	                                     function(X) tapply(X, timeSer4$ID2[x], cva, na.rm = TRUE)))
	  saveRDS(allSA.cva30, file = paste0(path.obj, sprintf("/allSA.cva30.%s.rds", i)))
	}
	
	# Recombine, erase small files at end, then save large file:
	allSA.cva30 <- readRDS(file = paste0(path.obj, "/allSA.cva30.1.rds"))
	allSA.cva30 <- data.frame(matrix(NA, nrow = n, ncol = ncol(allSA.cva30)))
	for (i in 1:length(start)) {
	  tmp <- readRDS(file = paste0(path.obj, sprintf("/allSA.cva30.%s.rds", i)))
	  allSA.cva30[start[i]:end[i], ] <- tmp
	  file.remove(file = paste0(path.obj, sprintf("/allSA.cva30.%s.rds", i)))
	  rm("tmp")
	}
	colnames(allSA.cva30) <- colnames(timeSer4)[4:ncol(timeSer4)]

	for (i in 4:ncol(allSA.cva30)) {
	  
	  # Make Inf's NA
	  allSA.cva30[is.infinite(allSA.cva30[, i]), i] <- NA
	  
	  # Replace NA's with 0's for cases where all values are 0
	  if (any(is.na(allSA.cva30[, i]), na.rm = T) & any(allSA.avg30[, i] == 0, na.rm = T)) {
	    allSA.cva30[which(allSA.avg30[, i] == 0), i] <- 0
	    
	  }
	  
	}
	allSA.cva30[, c("Latitude", "Longitude", "Elevation")] <- 
	  allSA.avg30[, c("Latitude", "Longitude", "Elevation")]
	
	# Check that 100 random entries give correct answer
	out <- logical(100)
	for (i in 1:100) {
	  col <- sample(colnames(allSA.cva30), 1)
	  row <- sample(1:nrow(allSA.cva30), 1)
	  names(out)[i] <- sprintf("%s,%s", col, row)
	  out[i] <- cva(subset(timeSer4[, col], timeSer4$ID2 == i), na.rm = T) == allSA.cva30[i, col]
	}
	all(out, na.rm = T) # should be true (some might be NA, long and lat will be FALSE)
	
	saveRDS(allSA.cva30, file = paste0(path.obj, "/allSA.cva30.rds"))
	allSA.cva30 <- readRDS(file = paste0(path.obj, "/allSA.cva30.rds"))
	
	# Import predicted suitability to weight points for calculating spatially-integrated climate
	# some biovars need to be converted from K to C, divided by 10, and log-transformed
	# bio2, bio3, bio4, bio10, bio11, bio12, bio14, bio15
	# log(x + 0.5) for bio3, bio10, bio12, bio14
	
	setwd("~/Google Drive/cardinalisENM/R objects/")
	
	# Load fitted models
	for (i in 1:10) {
	  mod.lr = get(load(paste("LR.mod2.",i,".pseudo11.Rda", sep="")))
	  assign(paste("LR.mod2.",i, sep=""), mod.lr)
	  mod.gam = get(load(paste("GAM.mod4.",i,".pseudo11.Rda", sep="")))
	  assign(paste("GAM.mod4.",i, sep=""), mod.gam)
	  mod.rf = get(load(paste("RF.mod1.",i,".pseudo11.Rda", sep="")))  
	  assign(paste("RF.mod1.",i, sep=""), mod.rf)
	  mod.brt = get(load(paste("BRT.mod4.",i,".pseudo11.Rda", sep="")))
	  assign(paste("BRT.mod4.",i, sep=""), mod.brt)  
	  mod.max = get(load(paste("MAX.mod1.",i,".pseudo11.Rda", sep="")))
	  assign(paste("MAX.mod1.",i, sep=""), mod.max)  
	}
	
	ext.accs.lr = get(load("LR.mod2.extaccs.pseudo11.Rda"))
	ext.accs.gam = get(load("GAM.mod4.extaccs.pseudo11.Rda"))
	ext.accs.rf = get(load("RF.mod1.extaccs.pseudo11.Rda"))
	ext.accs.brt = get(load("BRT.mod4.extaccs.pseudo11.Rda"))
	ext.accs.max = get(load("MAX.mod1.extaccs.pseudo11.Rda"))
	
	lr.cuts = ext.accs.lr[ext.accs.lr$thresh=="SensSpec", "threshold"]
	gam.cuts = ext.accs.gam[ext.accs.gam$thresh=="SensSpec", "threshold"]
	rf.cuts = ext.accs.rf[ext.accs.rf$thresh=="SensSpec", "threshold"]
	brt.cuts = ext.accs.brt[ext.accs.brt$thresh=="SensSpec", "threshold"]
	max.cuts = ext.accs.max[ext.accs.max$thresh=="SensSpec", "threshold"]
	
	cuts <- cbind(lr.cuts, gam.cuts, rf.cuts, brt.cuts, max.cuts)
	
	# Bioclimatic variables (30-year normal 1981 - 2010) for points are in allSA.avg30
	pred <- allSA.avg30[, c("bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]
	
	# Convert bio10 and bio1 from K to C
	pred$bio10 <- pred$bio10 - 273.16
	pred$bio11 <- pred$bio11 - 273.16
	
	# log-transform variables, consistent with ENM models
	pred$bio3 <- log(pred$bio3 + 0.5)
	pred$bio10 <- log(pred$bio10 + 0.5)
	pred$bio12 <- log(pred$bio12 + 0.5)
	pred$bio14 <- log(pred$bio14 + 0.5)
	
	# data.frame for suitability scores
	ENMprob <- data.frame(matrix(NA, ncol = 50, nrow = nrow(pred)))
	colnames(ENMprob) <- c(paste("LRprob", 1:10, sep = ""), 
	                       paste("GAMprob", 1:10, sep = ""), paste("RFprob", 1:10, sep = ""), 
	                       paste("BRTprob", 1:10, sep = ""), paste("MAXprob", 1:10, sep = ""))
	
	for (i in 1:10) {
	  ## LR predictions
	  mod = get(paste("LR.mod2.", i, sep=""))
	  ENMprob[, i] = predict(mod, pred, type = "response")
	  
	  ## GAM predictions
	  mod = get(paste("GAM.mod4.", i, sep = ""))
	  ENMprob[, 10 + i] = predict(mod, pred, type = "response")
	  
	  ## RF predictions
	  mod = get(paste("RF.mod1.", i, sep = ""))
	  ENMprob[, 20 + i] = predict(mod, pred, type = "prob", fun = predict, index = 2, overwrite = T)[, 2]
	  
	  ## BRT predictions
	  mod = get(paste("BRT.mod4.", i, sep = ""))
	  ENMprob[, 30 + i] = predict(mod, pred, n.trees = mod$gbm.call$best.trees, type = "response") 
	  
	  ## MAX predictions
	  mod = get(paste("MAX.mod1.", i, sep = ""))
	  ENMprob[, 40 + i] = predict(mod, pred, overwrite = T) 
	  
	}   
	
	# Calculate ensemble average of probabilities
	ENMprob$ensAvgProb <- apply(ENMprob, 1, mean)
	
	# Save occurence probabilities
	save(ENMprob, file = paste(path.obj, "/ENMprob_62km_all.RData", sep = ""))
	
	# Calculate spatial average climate (weighted by suitability)
	# Do for all climate columns
	cols <- which(!colnames(allSA.avg30) %in% c("Latitude", "Longitude", "Elevation", "PopID"))
	allSA.avg <- allSA.cva <- data.frame(matrix(rep(NA, 358 * ncol(allSA.avg30)), nrow = 358))
	colnames(allSA.avg) <- colnames(allSA.avg30)
	colnames(allSA.cva) <- colnames(allSA.cva30)
	allSA.avg$PopID <- rownames(allSA.avg) <- as.factor(as.character(all.pres$ID))
	allSA.cva$PopID <- rownames(allSA.cva) <- as.factor(as.character(all.pres$ID))
	allSA.avg30$PopID <- allSA.cva30$PopID <- as.factor(rep(as.character(all.pres$ID), each = 1e3))
	for (i in cols)
	{
	  for (j in levels(allSA.avg30$PopID)) 
	  {
	    allSA.avg[j, i] <- weighted.mean(allSA.avg30[allSA.avg30$PopID == j, i], 
	                                     ENMprob$ensAvgProb[allSA.avg30$PopID == j])
	    allSA.cva[j, i] <- weighted.mean(allSA.cva30[allSA.cva30$PopID == j, i], 
	                                     ENMprob$ensAvgProb[allSA.cva30$PopID == j])
	  }
	}
	
	# Adjust column names and export
	# Big data.frame with all 358000 points
	x <- which(colnames(allSA.avg30) %in% c("Latitude", "Longitude", "Elevation", "PopID")) 
	colnames(allSA.avg30)[-x] <- paste(colnames(allSA.avg30)[-x], "_avg", sep = "")
	colnames(allSA.cva30)[-x] <- paste(colnames(allSA.cva30)[-x], "_cva", sep = "")
	
	save(timeSer4, file = paste(path.obj, "/AllOcc1000_timeSer.RData", sep = ""))
	save(allSA.avg30, file = paste(path.obj, "/AllOcc1000_avg1981-2010.RData", sep = ""))
	save(allSA.cva30, file = paste(path.obj, "/AllOcc1000_cva1981-2010.RData", sep = ""))
	
	# Small data.frame with all 358 points from all populations
	x <- which(colnames(allSA.avg) %in% c("Latitude", "Longitude", "Elevation", "PopID")) 
	colnames(allSA.avg)[-x] <- paste(colnames(allSA.avg)[-x], "_avg", sep = "")
	colnames(allSA.cva)[-x] <- paste(colnames(allSA.cva)[-x], "_cva", sep = "")	
	
	save(allSA.avg, file = paste(path.obj, "/AllOccSA_avg1981-2010_62km.RData", sep = ""))
	save(allSA.cva, file = paste(path.obj, "/AllOccSA_cva1981-2010_62km.RData", sep = ""))
	
