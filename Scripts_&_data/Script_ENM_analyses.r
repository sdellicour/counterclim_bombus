# NOTES:
	# - for the sampling of pseudo-absences, some improvements should be performed as compared to the previous study: (i) we should use a ratio 1:1, (ii) and each replicate should be based
	#	on an independent sampling of pseudo-absences within the background (= set of raster cells where there is at least one other record for another Bombus species); and (iii) finally, 
	#	for the species for which there is not enough background cells remaining to reach a 1:1 ratio, we should sample the remaining pseudo-absences in the rest of the study area

library(ade4)
library(ape)
library(blockCV)
library(diagram)
library(dismo)
library(gbm)
library(geosphere)
library(ggplot2)
library(HDInterval)
library(lubridate)
library(ncdf4)
library(ncf)
library(picante)
library(pgirmess)
library(phytools)
library(RColorBrewer)
library(raster)
# library(rgdal)
# library(rgeos)
library(seqinr)
library(sp)
library(vioplot)

directory = "Bombus_obs_111224"; savingPlots = FALSE
timePeriods = c("2000_2019"); periods = list(c(2000,2019))
data = readRDS(paste0(directory,"/",directory,".rds"))
species = gsub(" ","_",unique(data$TAXON))
species = data.frame(species[order(species)])
extentOfStudyArea = extent(-12,29,36,72)

envVariableNames = c("temperature_winter","temperature_spring","temperature_summer","temperature_inFall",
					 "precipitation_winter","tprecipitation_spring","precipitation_summer","precipitation_inFall",
					 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_inFall",
					 "primary_forest_areas","primary_non-forest_areas","secondary_forest_areas","secondary_non-forest_areas",
					 "croplands_all_categories","managed_pasture_and_rangeland","human_pop_density_log10")
models_isimip3a = c("gswp3-w5e5","20crv3","20crv3-era5","20crv3-w5e5")
models_isimip3a_names = c("GSWP3-W5E5", "20CRv3", "20CRv3-ERA5", "20CRv3-W5E5")

# 1. Preparation of land cover and climatic environmental rasters

	# 1.1. Preparation of the European shapefile that will be used as a mask

if (!file.exists("Continents_shapefile/Simplified_European_contour.shp"))
	{
		coastlines = shapefile("Continents_shapefile/NaturalEarth_10m_continents.shp")
		countries = shapefile("Continents_shapefile/NaturalEarth_10m_countries.shp")
		countries_EU = subset(crop(countries,extentOfStudyArea), CONTINENT=="Europe")
		europe1 = crop(coastlines, countries_EU); polygons = list(); c = 0
		for (i in 1:length(europe1@polygons))
			{
				for (j in 1:length(europe1@polygons[[i]]@Polygons))
					{
						if (europe1@polygons[[i]]@Polygons[[j]]@area > 0.5)
							{
								c = c+1; polygons[[c]] = europe1@polygons[[i]]@Polygons[[j]]
							}
					}
			}
		pols = Polygons(polygons, 1); pols_list = list(); pols_list[[1]] = pols
		europe2 = SpatialPolygons(pols_list); europe3 = gSimplify(europe2, 0.1)
		europe2@proj4string=europe1@proj4string; europe3@proj4string=europe1@proj4string
		europe3_spdf = SpatialPolygonsDataFrame(europe3, data.frame(IS=1:length(europe3@polygons)))
		writeOGR(europe3_spdf, dsn="Continents_shapefile", layer="Simplified_European_contour", driver="ESRI Shapefile")
	}	else		{
		europe3 = shapefile("Continents_shapefile/Simplified_European_contour.shp")
	}

	# 1.2. Preparing the different environmental rasters for the ENM analyses

envVariableNames = c("temperature_winter","temperature_spring","temperature_summer","temperature_inFall",
					 "precipitation_winter","tprecipitation_spring","precipitation_summer","precipitation_inFall",
					 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_inFall",
					 "primary_forest_areas","primary_non-forest_areas","secondary_forest_areas","secondary_non-forest_areas",
					 "croplands_all_categories","managed_pasture_and_rangeland","human_pop_density_log10")
envVariables_list = list(); i = 1
for (i in 1:length(models_isimip3a))
	{
		temperature = brick(paste0("Environmental_rasters/ISIMIP3a/tas_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
		precipitation = brick(paste0("Environmental_rasters/ISIMIP3a/pr_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
		relative_humidity = brick(paste0("Environmental_rasters/ISIMIP3a/hurs_day_obsclim_historical_",models_isimip3a[i],"_2000_2019_ymonmean.nc"))
		# land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2015_timmean.nc4"))
		land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2019_timmean.nc"))
		population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_2000_2019_timmean.nc4"))		
		temperature_winter = mean(temperature[[12]],temperature[[1]],temperature[[2]])-273.15 # conversion to Celcius degrees
		temperature_spring = mean(temperature[[3]],temperature[[4]],temperature[[5]])-273.15 # conversion to Celcius degrees
		temperature_summer = mean(temperature[[6]],temperature[[7]],temperature[[8]])-273.15 # conversion to Celcius degrees
		temperature_inFall = mean(temperature[[9]],temperature[[10]],temperature[[11]])-273.15 # conversion to Celcius degrees
		precipitation_winter = mean(precipitation[[12]],precipitation[[1]],precipitation[[2]])*60*60*24 # conversion to kg/m2/day
		precipitation_spring = mean(precipitation[[3]],precipitation[[4]],precipitation[[5]])*60*60*24 # conversion to kg/m2/day
		precipitation_summer = mean(precipitation[[6]],precipitation[[7]],precipitation[[8]])*60*60*24 # conversion to kg/m2/day
		precipitation_inFall = mean(precipitation[[9]],precipitation[[10]],precipitation[[11]])*60*60*24 # conversion to kg/m2/day
		relative_humidity_winter = mean(relative_humidity[[12]],relative_humidity[[1]],relative_humidity[[2]])
		relative_humidity_spring = mean(relative_humidity[[3]],relative_humidity[[4]],relative_humidity[[5]])
		relative_humidity_summer = mean(relative_humidity[[6]],relative_humidity[[7]],relative_humidity[[8]])
		relative_humidity_inFall = mean(relative_humidity[[9]],relative_humidity[[10]],relative_humidity[[11]])
		landCoverVariableNames = as.character(read.csv("Environmental_rasters/Luse.csv")[1:12,2])
		landCoverVariableIDs = as.character(read.csv("Environmental_rasters/Luse.csv")[1:12,1])
		land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
		for (j in 1:12)
			{
				land_covers1[[j]] = brick(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_2000_2019_timmean.nc"), varname=landCoverVariableIDs[j])
			}
		variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
		variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
				  		   "potentially forested secondary land","potentially non-forested secondary land")
		for (j in 1:length(variable_names))
			{
				names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[j])
				if (length(indices) == 0) indices = which(grepl(variable_names[j],names))
				if (variable_names[j] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
				land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[j]; # print(indices)
				if (length(indices) > 1)
					{
						for (k in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[k]]][]
					}
				land_covers2[[j]] = land_cover[[1]]; land_covers3[[j]] = raster::aggregate(land_cover[[1]],2)
			}
		envVariables = list(); population_log10 = population; population_log10[] = log10(population_log10[]+1) 
		envVariables[[1]] = temperature_winter; envVariables[[2]] = temperature_spring
		envVariables[[3]] = temperature_summer; envVariables[[4]] = temperature_inFall
		envVariables[[5]] = precipitation_winter; envVariables[[6]] = precipitation_spring
		envVariables[[7]] = precipitation_summer; envVariables[[8]] = precipitation_inFall
		envVariables[[9]] = relative_humidity_winter; envVariables[[10]] = relative_humidity_spring
		envVariables[[11]] = relative_humidity_summer; envVariables[[12]] = relative_humidity_inFall
		envVariables[[13]] = land_covers2[[4]] # primary forest areas
		envVariables[[14]] = land_covers2[[5]] # primary non-forest areas
		envVariables[[15]] = land_covers2[[6]] # secondary forest areas
		envVariables[[16]] = land_covers2[[7]] # secondary non-forest areas
		envVariables[[17]] = land_covers2[[1]] # croplands (all catergories)
		envVariables[[18]] = land_covers2[[2]] # managed pasture + rangeland
		envVariables[[19]] = population_log10 # human population (log-transformed)
		for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], extentOfStudyArea, snap="out")
		for (j in 1:12) envVariables[[j]][which(is.na(envVariables[[19]][]))] = NA 
		for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], europe3, snap="out")
		for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], europe3)
		envVariables_list[[1]] = envVariables
		if ((i == 1)&(savingPlots))
			{
				pdf("All_the_figures_&_SI/Rasters_GSWP3_t0_NEW.pdf", width=8, height=5.8); par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
				r = envVariables[[1]]; r[!is.na(r[])] = 1; contour1 = rasterToPolygons(r>0, dissolve=T)
				r = envVariables[[13]]; r[!is.na(r[])] = 1; contour2 = rasterToPolygons(r>0, dissolve=T)	
				j = 1; colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[1]][],envVariables[[2]][],envVariables[[3]][],envVariables[[4]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 2; colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[1]][],envVariables[[2]][],envVariables[[3]][],envVariables[[4]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 3; colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[1]][],envVariables[[2]][],envVariables[[3]][],envVariables[[4]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 4; colourScale = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[1]][],envVariables[[2]][],envVariables[[3]][],envVariables[[4]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Air temperature", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, °C)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 13; colourScale = c("#E5E5E5",colorRampPalette(c("gray97","chartreuse4"),bias=1)(121)[21:121]); contour = contour2
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Primary forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 14; colourScale = c("#E5E5E5",colorRampPalette(c("gray97","darkseagreen4"),bias=1)(121)[21:121]); contour = contour2
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Primary non-forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 5; colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101]; contour = contour1
				vS = c(envVariables[[5]][],envVariables[[6]][],envVariables[[6]][],envVariables[[7]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 6; colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101]; contour = contour1
				vS = c(envVariables[[5]][],envVariables[[6]][],envVariables[[6]][],envVariables[[7]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 7; colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101]; contour = contour1
				vS = c(envVariables[[5]][],envVariables[[6]][],envVariables[[6]][],envVariables[[7]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 8; colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101]; contour = contour1
				vS = c(envVariables[[5]][],envVariables[[6]][],envVariables[[6]][],envVariables[[7]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Precipitation", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, kg/m2/day)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 18; colourScale = c("#E5E5E5",colorRampPalette(c("gray97","burlywood3"),bias=1)(121)[21:121]); contour = contour2
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Pastures and", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("rangeland", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 15; colourScale = c("#E5E5E5",colorRampPalette(c("gray97","olivedrab3"),bias=1)(121)[21:121]); contour = contour2
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Secondary forested", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("areas", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 9; colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[9]][],envVariables[[10]][],envVariables[[11]][],envVariables[[12]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(winter, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 10; colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[9]][],envVariables[[10]][],envVariables[[11]][],envVariables[[12]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(spring, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 11; colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[9]][],envVariables[[10]][],envVariables[[11]][],envVariables[[12]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(summer, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 12; colourScale = colorRampPalette(brewer.pal(9,"PuBu"))(151)[1:101]; contour = contour1
				vS = c(envVariables[[9]][],envVariables[[10]][],envVariables[[11]][],envVariables[[12]][]); vS = vS[which(!is.na(vS))]
				index1 = (((min(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1; index2 = (((max(envVariables[[j]][],na.rm=T)-min(vS))/(max(vS)-min(vS)))*100)+1
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScale[index1:index2]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T))))
				plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Relative humidity", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(fall, %)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 17; colourScale = c("#E5E5E5",colorRampPalette(c("gray97","navajowhite4"),bias=1)(121)[21:121]); contour = contour2
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Croplands", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("(all categories)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				j = 19; colourScale = colorRampPalette(brewer.pal(9,"BuPu"))(191)[21:121]; contour = contour1
				vS = envVariables[[j]][]; rast = raster(as.matrix(c(min(vS,na.rm=T),max(vS,na.rm=T)))); cols = colourScale
				plot(europe3, lwd=0.1, border=NA, col=NA); plot(envVariables[[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Human population", side=3, line=-1.1, at=5, cex=0.37, col="gray30"); mtext("density (log10/km2)", side=3, line=-1.6, at=5, cex=0.37, col="gray30")
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
					 axis.args=list(cex.axis=0.55, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.30,0)), alpha=1, side=3)
				dev.off()
			}
	}

# 2. Loading and ploting all the Bombus records by continent

nullRaster = envVariables_list[[1]][[1]]; nullRaster[!is.na(nullRaster[])] = 1; names(nullRaster) = "nullRaster"
for (i in 1:dim(species)[1])
	{
		if (!file.exists(paste0(directory,"/",species[i,1],".csv")))
			{
				tab = data[which(data$TAXON==gsub("_"," ",species[i,1])),]
				tab = tab[,c("YEAR_2","LONGITUDE","LATITUDE")]
				colnames(tab) = c("year","longitude","latitude")
				tab = tab[which(!is.na(tab[,"year"])),]
				write.csv(tab, paste0(directory,"/",species[i,1],".csv"), row.names=F, quote=F)
			}
	}
observations_list = list(); t = 1
for (t in 1:length(periods))
	{
		species = gsub(" ","_",unique(data$TAXON))
		species = species[which(species!="Bombus_cullumanus")] # discarded because mainly in Middle-Eas/Asia
		species = species[which(species!="Bombus_haematurus")] # discarded because mainly in Middle-Eas/Asia
		species = species[which(species!="Bombus_inexspectatus")] # discarded because associated with a distribution
		species = species[which(species!="Bombus_laesus")] # discarded because mainly in Middle-Eas/Asia
		species = species[which(species!="Bombus_magnus")] # discarded because insular (in Corsica)
		species = species[which(species!="Bombus_schrencki")] # discarded because insular (in Corsica)
		species = species[which(species!="Bombus_xanthopus")] # discarded because insular (in Corsica)
		species = data.frame(species[order(species)]); indices = c(); c = 0 # too restricted compared to its host
		observations = list(); minYear = periods[[t]][1]; maxYear = periods[[t]][2]
		for (i in 1:dim(species)[1])
			{
				tab1 = read.csv(paste0(directory,"/",species[i,1],".csv"), header=T)
				tab2 = tab1[which((tab1[,"year"]>=minYear)&(tab1[,"year"]<=maxYear)),c("longitude","latitude")]
				tab3 = tab2[which(!is.na(raster::extract(nullRaster,tab2))),]; coordinates = rep(NA, dim(tab3)[1])
				for (j in 1:dim(tab3)[1]) coordinates[j] = paste0(tab3[j,"longitude"],"_",tab3[j,"latitude"])
				if (length(unique(coordinates)) >= 30)
					{
						c = c+1; indices = c(indices, i); observations[[c]] = tab3
					}
			}
		species = data.frame(species[indices,]); colnames(species) = "species"; observations_list[[t]] = observations
		if (savingPlots == TRUE)
			{
				pdf(paste0("All_the_figures_&_SI/Bombus_data_1on2_NEW.pdf"), width=8, height=((5.8/3)*4))
				par(mfrow=c(4,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.4, col="gray30")
				for (i in 1:24)
					{
						# plot(nullRaster, col=NA, axes=F, ann=F, box=F, legend=F)
						plot(europe3, lwd=0.8, border="gray50", col="gray90")
						col1 = rgb(13,110,111,220,maxColorValue=255); col2 = rgb(203,185,68,220,maxColorValue=255)
						points(observations[[i]], col=col1, pch=3, cex=0.3, lwd=0.3)
						mtext(gsub("Bombus_","B. ",species[i,1]), side=3, line=-2, at=0, cex=0.50, col="gray30")
					}
				dev.off()
				pdf(paste0("All_the_figures_&_SI/Bombus_data_2on2_NEW.pdf"), width=8, height=((5.8/3)*4))
				par(mfrow=c(4,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.4, col="gray30")
				for (i in 25:47)
					{
						# plot(nullRaster, col=NA, axes=F, ann=F, box=F, legend=F)
						plot(europe3, lwd=0.8, border="gray50", col="gray90")
						col1 = rgb(13,110,111,220,maxColorValue=255); col2 = rgb(203,185,68,220,maxColorValue=255)
						points(observations[[i]], col=col1, pch=3, cex=0.3, lwd=0.3)
						mtext(gsub("Bombus_","B. ",species[i,1]), side=3, line=-2, at=0, cex=0.50, col="gray30")
					}
				dev.off()
			}
	}

