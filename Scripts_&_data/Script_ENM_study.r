library(ade4)
library(ape)
library(beeswarm)
library(blockCV)
library(colorspace)
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
library(phytools)
library(RColorBrewer)
library(raster)
library(rgdal)
library(rgeos)
library(seqinr)
library(sp)
library(vioplot)

# 1. Preparation of the different land cover and climatic environmental rasters
# 2. Loading and ploting all the occurrence records for each Bombus species in Europe
# 3. Boosted regression tree analyses with a spatial cross-validation procedure
# 4. Computation of the prevalence-pseudoabsence-calibrated Sørensen index (SI_pcc)
# 5. Computation and analysis of the relative influence of each environmental factor
# 6. BRT projections based on historical and counterfactual climate simulations
# 7. Computation and mapping of the ESI and SRI metrics for the different time periods
# 8. Computation and mapping of the evolution of the Boyce index (BI) through time

directory = "Bombus_obs_111224"; savingPlots = FALSE
timePeriods = c("2000_2019"); periods = list(c(2000,2019))
data = readRDS(paste0(directory,"/",directory,".rds"))
species = gsub(" ","_",unique(data$TAXON))
species = data.frame(species[order(species)])
extentOfStudyArea = extent(-12,29,36,72)

envVariableNames = c("temperature_winter","temperature_spring","temperature_summer","temperature_inFall",
					 "precipitation_winter","precipitation_spring","precipitation_summer","precipitation_inFall",
					 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_inFall",
					 "primary_forest_areas","primary_non-forest_areas","secondary_forest_areas","secondary_non-forest_areas",
					 "croplands_all_categories","managed_pasture_and_rangeland","human_pop_density_log10")
models_isimip3a = c("gswp3-w5e5","20crv3","20crv3-era5","20crv3-w5e5")
models_isimip3a_names = c("GSWP3-W5E5", "20CRv3", "20CRv3-ERA5", "20CRv3-W5E5")
scenarios = c("obsclim","counterclim"); selectedCV = "SCV2"
pastPeriods = c("1901-1919","1920-1939","1940-1959","1960-1979","1980-1999","2000-2019")

# 1. Preparation of the different land cover and climatic environmental rasters

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
	}	else	{
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
		population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_2000_2019_timmean.nc4"), varname="total-population")		
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
		envVariables[[13]] = land_covers3[[4]] # primary forest areas
		envVariables[[14]] = land_covers3[[5]] # primary non-forest areas
		envVariables[[15]] = land_covers3[[6]] # secondary forest areas
		envVariables[[16]] = land_covers3[[7]] # secondary non-forest areas
		envVariables[[17]] = land_covers3[[1]] # croplands (all catergories)
		envVariables[[18]] = land_covers3[[2]] # managed pasture + rangeland
		envVariables[[19]] = population_log10 # human population (log-transformed)
		for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], extentOfStudyArea, snap="out")
		for (j in 1:12) envVariables[[j]][which(is.na(envVariables[[19]][]))] = NA 
		for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], europe3, snap="out")
		for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], europe3)
		envVariables_list[[i]] = envVariables
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

# 2. Loading and ploting all the occurrence records for each Bombus species in Europe

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
		species = species[which(species!="Bombus_cullumanus")] # discarded because mainly in Middle-East/Asia
		species = species[which(species!="Bombus_haematurus")] # discarded because mainly in Middle-East/Asia
		species = species[which(species!="Bombus_inexspectatus")] # discarded because associated with a distribution
		species = species[which(species!="Bombus_laesus")] # discarded because mainly in Middle-East/Asia
		species = species[which(species!="Bombus_magnus")] # discarded because associated with a sampling bias
		species = species[which(species!="Bombus_schrencki")] # discarded because with an eastern European distribution
		species = species[which(species!="Bombus_xanthopus")] # discarded because insular (in Corsica)
		species = data.frame(species[order(species)]); indices = c(); c = 0
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

# 3. Boosted regression tree analyses with a spatial cross-validation procedure

nullRaster = envVariables_list[[1]][[1]]; nullRaster[!is.na(nullRaster[])] = 1
names(nullRaster) = "nullRaster"; allObservationsOnTheContinent = c()
for (i in 1:dim(species)[1])
	{
		allObservationsOnTheContinent = rbind(allObservationsOnTheContinent, observations_list[[1]][[i]])
	}
backgroundCells = unique(raster::extract(nullRaster, allObservationsOnTheContinent, cellnumbers=T))
background = nullRaster; background[!(1:length(background[]))%in%backgroundCells] = NA
newAnalyses = FALSE; spatialCrossValidation1 = FALSE; spatialCrossValidation2 = TRUE
nberOfReplicates = 10; occurrence_data_summary = matrix(nrow=dim(species)[1], ncol=2)
row.names(occurrence_data_summary) = species[,1]; colnames(occurrence_data_summary) = c("n","n_filtered")
if (!file.exists("BRT_data_frames.rds"))
	{
		all_data = list(); rasters_stack = stack(envVariables_list[[1]]); names(rasters_stack) = envVariableNames
		for (i in 1:dim(species)[1])
			{
				observations = observations_list[[1]][[i]]; data_to_discard = c()
				for (j in 1:length(rasters_stack@layers)) # to discard the observations that do not fall within the rasters
					{
						data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[j]],observations))))
					}
				data_to_discard = unique(data_to_discard); data_to_discard = data_to_discard[order(data_to_discard)]
				if (length(data_to_discard) > 0) observations = observations[which(!c(1:dim(observations)[1])%in%data_to_discard),]
				presenceCells = unique(raster::extract(nullRaster, observations, cellnumbers=T))
				targetSpeciesBackground = background; targetSpeciesBackground[(1:length(targetSpeciesBackground[]))%in%presenceCells] = NA
				studyAreaMinusBackground = nullRaster; studyAreaMinusBackground[which(targetSpeciesBackground[]==1)] = NA
				cellIDs = unique(cellFromXY(nullRaster, observations)); presences = xyFromCell(nullRaster, cellIDs)
				occurrence_data_summary[i,c("n","n_filtered")] = cbind(dim(observations)[1], length(cellIDs))
				n_PAs_1 = length(cellIDs); n_background = length(which(!is.na(targetSpeciesBackground[]))); datas = list()
				for (j in 1:nberOfReplicates)
					{
						if (n_background > n_PAs_1)
							{
								pseudoAbsences = xyFromCell(targetSpeciesBackground, sample(which(!is.na(values(targetSpeciesBackground))), n_PAs_1, replace=F))
							}	else	{
								n_PAs_2 = n_PAs_1-n_background
								pseudoAbsences1 = xyFromCell(targetSpeciesBackground, sample(which(!is.na(values(targetSpeciesBackground))), n_background, replace=F))
								pseudoAbsences2 = xyFromCell(studyAreaMinusBackground, sample(which(!is.na(values(studyAreaMinusBackground))), n_PAs_2 , replace=F))
								pseudoAbsences = rbind(pseudoAbsences1, pseudoAbsences2); colnames(pseudoAbsences) = colnames(observations)
							}
						data1 = cbind(rep(1,dim(presences)[1]),presences); colnames(data1) = c("response","x","y")
						data2 = cbind(rep(0,dim(pseudoAbsences)[1]),pseudoAbsences); colnames(data2) = c("response","x","y")
						data = rbind(data1, data2); data = cbind(data, raster::extract(rasters_stack, data[,2:3]))
						colnames(data) = c("response","x","y",envVariableNames); datas[[j]] = data
					}
				all_data[[i]] = datas
			}
		saveRDS(all_data, "BRT_data_frames.rds")
	}
if (!file.exists(paste0("Occurrence_data.csv")))
	{
		write.csv(occurrence_data_summary, "Occurrence_data.csv", quote=F)
	}
all_data = readRDS("BRT_data_frames.rds")
savingCorrelogram = FALSE; blockSizeExplorationByBlockCV = FALSE
if (savingCorrelogram == TRUE)
	{
		for (i in 1:dim(species)[1])
			{
				if (!file.exists(paste0("Correlogram_graphics/",species[i,1],".pdf")))
					{
						datas = all_data[[i]]; correlograms = list()
						for (j in 1:nberOfReplicates)
							{
								correlograms[[j]] = ncf::correlog(datas[[j]][,"x"], datas[[j]][,"y"], datas[[j]][,"response"], na.rm=T, increment=100, resamp=0, latlon=T)
							}
						pdf(paste0("Correlogram_graphics/",species[i,1],".pdf"), width=4.5, height=3); par(mar=c(2.7,2.8,1.2,1.2))
						plot(correlograms[[1]]$mean.of.class[-1], correlograms$correlations[[1]][-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.7,0.5), xlim=c(0,2500))
						for (j in 1:nberOfReplicates)
							{
								lines(correlograms[[j]]$mean.of.class[-1], correlograms[[j]]$correlation[-1], lwd=0.1, col="gray30")
							}
						abline(h=0, lwd=0.5, col="red", lty=2)
						axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.020, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,3000,500), labels=c("0","500","1000","1500","2000","2500","3000"))
						axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.020, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.8,0.6,0.2), labels=c("","-0.6","-0.4","-0.2","0","0.2","0.4",""))
						title(xlab="Distance (km)", cex.lab=0.7, mgp=c(0.9,0,0), col.lab="gray30")
						title(ylab="Correlation", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
						dev.off()
					}
			}
	}
if (blockSizeExplorationByBlockCV == TRUE)
	{
		buffer = list(); ranges1 = matrix(nrow=dim(species)[1], ncol=1); row.names(ranges1) = species[,1]
		for (i in 2:length(envVariables_list[[1]])) buffer[[i-1]] = terra::rast(envVariables_list[[1]][[i]])
		for (i in 2:length(envVariables_list[[1]])) names(buffer[[i-1]]) = envVariableNames[i]
		test = cv_spatial_autocor(r=terra::rast(buffer), # a SpatRaster object or path to files
								  num_sample=10000, plot=T) # num_sample: number of cells to be used
		for (i in 1:dim(species)[1])
			{
				ranges2 = rep(NA, nberOfReplicates) 
				for (j in 1:nberOfReplicates)
					{
						pa_data = sf::st_as_sf(data.frame(all_data[[i]][[j]]), coords=c("x","y"), crs=4326)
						trycatch = tryCatch(
							{
								ranges2[j] = cv_spatial_autocor(x=pa_data, column="response", plot=F)$range
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
					}
				cat(i,", ",gsub("_","",species[i,]),": ",length(which(is.na(ranges2)))," \"NA\"","\n",sep="")
				ranges1[i,1] = round(mean(ranges2, na.rm=T))
			}
		ranges = ranges1
	}
samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
				selection = which(dists[[(i-1)]] > minDist)
				if (length(selection) == 0)
					{
						stop("Restarts the function with a smaller minimum distance")
					}
				selection_list[[i]] = selection
				test = table(unlist(selection_list))
				indices_minDist = as.numeric(names(which(test==i)))
				indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}
if (newAnalyses == TRUE) { for (h in 1:length(models_isimip3a)) { for (i in 1:dim(species)[1]) {
		rasters_stack = stack(envVariables_list[[h]])
		names(rasters_stack) = envVariableNames
		theRanges = c(500,500)*1000 # distance in meters
		gbm.x = gsub("\\.","-",names(rasters_stack))
		gbm.y = "response"
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.005 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(all_data[[i]][[1]])[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 10 # initial number of trees
		step.size = 5 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE
		brt_model_ccvs = list() # classic cross-validations (CCVs)
		brt_model_scv1 = list() # spatial cross-validations 1 (SCV1)
		brt_model_scv2 = list() # spatial cross-validations 2 (SCV2)
		if (spatialCrossValidation2 == TRUE)	
			{
				if (spatialCrossValidation1 == TRUE)
					{
						AUCs = matrix(nrow=nberOfReplicates, ncol=3); colnames(AUCs) = c("CCV_AUC","SCV1_AUC","SCV2_AUC")
					}	else	{
						AUCs = matrix(nrow=nberOfReplicates, ncol=2); colnames(AUCs) = c("CCV_AUC","SCV2_AUC")
					}
			}	else	{
				if (spatialCrossValidation1 == TRUE)
					{
						AUCs = matrix(nrow=nberOfReplicates, ncol=2); colnames(AUCs) = c("CCV_AUC","SCV1_AUC")
					}	else	{
						AUCs = matrix(nrow=nberOfReplicates, ncol=1); colnames(AUCs) = c("CCV_AUC")
					}
			}
		for (j in 1:nberOfReplicates)
			{
				data = as.data.frame(all_data[[i]][[j]]); n.trees = 10; learning.rate = 0.005; step.size = 5; fold.vector = NULL; worked = FALSE
				# BRT with classic (standard) cross-validation (CCV):
				pdf(file=paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_CCV_",j,".pdf"))
				# while (worked == FALSE)
					# {
						# trycatch = tryCatch(
							# {
								brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							# },	error = function(cond) {
							# },	warning = function(cond) {
							# },	finally = {
							# })
						# if (length(brt_model_ccvs) == j) worked = TRUE
					# }
				dev.off()
				AUCs[j,1] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
				if (spatialCrossValidation1 == TRUE)
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the folds generation of Dhingra, Artois et al. (2016, eLife):
						folds_with_similar_sizes = FALSE; c = 0
						while (folds_with_similar_sizes == FALSE) # while loop to select a partition where the x folds gather at least
							{									  # proportion = (1/(x+1)) of the total number of presence points
								data_presence = data[which(data[,1]==1),]; c = c+1; # print(c)
								fivePoints = samplingPtsMinDist(data_presence[,c("x","y")], minDist=200, nberOfPoints=n.folds)
								fold.vector = foldSelection(data[,c("x","y")], selectedPoints=data_presence[fivePoints,c("x","y")])
								fold.vector_presences = fold.vector[which(data[,1]==1)]
								counts = hist(fold.vector_presences, plot=F)$counts
								props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
								if (min(props) > (1/(n.folds*2))) folds_with_similar_sizes = TRUE
							}
						if (showingFoldsPlot == TRUE)
							{
								par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
								cols = c("olivedrab3","tan3","steelblue3","orange1","tomato2","mediumseagreen")[fold.vector]
								plot(backgrounds[[1]], col="gray90", useRaster=T, colNA=NA, box=F, axes=F, legend=F)
								pchs = c(16,3)[data[,1]+1]; cexs = c(0.25,0.5)[data[,1]+1]
								points(data[,c("x","y")], col=cols, pch=pchs, cex=cexs, lwd=0.7)
							}
						n.trees = 10; learning.rate = 0.001; step.size = 2; worked = FALSE
						pdf(file=paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_SCV1_",j,".pdf"))
						# while (worked == FALSE)
							# {
								# trycatch = tryCatch(
									# {
										brt_model_scv1[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
											var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
											verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
									# },	error = function(cond) {
									# },	warning = function(cond) {
									# },	finally = {
									# })
								# if (length(brt_model_scv1) == j) worked = TRUE
							# }
						dev.off()
						AUCs[j,"SCV1_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean # mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)		
					}
				if (spatialCrossValidation2 == TRUE)
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
						spdf = SpatialPointsDataFrame(data[c("x","y")], data[,c(1,4:dim(data)[2])], proj4string=crs(nullRaster)); worked = FALSE
						# while (worked == FALSE)
							# {
								# trycatch = tryCatch(
									# {
										myblocks = NULL
										# myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRaster, k=n.folds, theRange=theRanges[1], selection="random")
										myblocks = cv_spatial(spdf, column="response", k=n.folds, size=theRanges[1], selection="random")
									# },	error = function(cond) {
									# },	finally = {
									# })
								# if (!is.null(myblocks)) worked = TRUE
							# }
						n.trees = 10; learning.rate = 0.005; step.size = 5; fold.vector = myblocks$folds_ids # "myblocks$foldID" not correct!
						pdf(file=paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_SCV2_",j,".pdf"))
						brt_model_scv2[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.off()
						AUCs[j,"SCV2_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean # mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
					}
			}
		saveRDS(brt_model_ccvs, paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_",nberOfReplicates,"_CCV.rds"))
		if (spatialCrossValidation1 == TRUE)	 saveRDS(brt_model_scv1, paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_",nberOfReplicates,"_SCV1.rds"))
		if (spatialCrossValidation2 == TRUE) saveRDS(brt_model_scv2, paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_",nberOfReplicates,"_SCV2.rds"))
		write.csv(AUCs, paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_AUCs.csv"), row.names=F, quote=F)
	}}}
if (!file.exists(paste0("All_AUC_values.csv")))
	{
		AUC_values = matrix(nrow=dim(species)[1], ncol=8); row.names(AUC_values) = species[,"species"]; colNames = c()
		for (h in 1:length(models_isimip3a))
			{
				colNames = c(colNames, paste0("CCV_",models_isimip3a_names[h]), paste0("SCV2_",models_isimip3a_names[h]))
				for (i in 1:dim(species)[1])
					{
						tab = read.csv(paste0("BRT_projection_files/BRT_models/",species[i,"species"],"_",models_isimip3a_names[h],"_AUCs.csv"), head=T)
						for (j in 1:dim(tab)[2])
							{
								meanV = round(mean(tab[,j]),3); sdV = round(sd(tab[,j]),3)
								if (nchar(meanV) == 4) meanV = paste0(meanV,"0")
								if (nchar(meanV) == 3) meanV = paste0(meanV,"00")
								if (nchar(meanV) == 1) meanV = paste0(meanV,".000")
								if (nchar(sdV) == 4) sdV = paste0(sdV,"0")
								if (nchar(sdV) == 3) sdV = paste0(sdV,"00")
								if (nchar(sdV) == 1) sdV = paste0(sdV,".000")
								AUC_values[i,((h-1)*2)+j] = paste0(meanV," (",sdV,")")
							}
					}
			}
		colnames(AUC_values) = colNames; write.csv(AUC_values, "All_AUC_values.csv", quote=F)
	}
AUC_values = read.csv("All_AUC_values.csv", head=T)

# 4. Computation of the prevalence-pseudoabsence-calibrated Sørensen index (SI_pcc)

	# - computation performed according to the formulas of Leroi et al. (2018, J. Biogeography)
	# - optimisation of the threshold with a 0.01 step increment according to Li & Guo (2013, Ecography)

if (!file.exists(paste0("All_SIppc_values.csv")))
	{
		tab = matrix(nrow=dim(species)[1], ncol=8); row.names(tab) = species[,1]; tabs_list1 = list(); colNames = c()
		for (h in 1:length(models_isimip3a))
			{
				colNames = c(colNames, paste0("SIppc_",models_isimip3a_names[h]), paste0("OptTres_",models_isimip3a_names[h]))
				rasters_stack = stack(envVariables_list[[h]]); tabs_list2 = list(); background_cells = dim(backgroundCells)[1]
				for (i in 1:dim(species)[1])
					{
						brt_model_scv2 = readRDS(paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_10_SCV2.rds"))
						sorensen_ppcs = rep(NA, length(brt_model_scv2)); thresholds = rep(NA, length(brt_model_scv2)); tabs_list3 = list()
						for (j in 1:length(brt_model_scv2))
							{
								tmp = matrix(nrow=101, ncol=2); tmp[,1] = seq(0,1,0.01)
								df = brt_model_scv2[[j]]$gbm.call$dataframe
								responses = df$response; data = df[,4:dim(df)[2]]
								n.trees = brt_model_scv2[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
								prediction = predict.gbm(brt_model_scv2[[j]], data, n.trees, type, single.tree)		
								N = backgroundCells
								P = sum(responses==1); A = sum(responses==0)
								prev = P/(P+A) # proportion of recorded sites where the species is present
								x = (P/A)*((1-prev)/prev)
								sorensen_ppc = 0
								for (threshold in seq(0,1,0.01))
									{
										TP = length(which((responses==1)&(prediction>=threshold))) # true positives
										FN = length(which((responses==1)&(prediction<threshold))) # false negatives
										FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
										sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
										tmp[which(tmp[,1]==threshold),2] = sorensen_ppc_tmp
										if (sorensen_ppc < sorensen_ppc_tmp)
											{
												sorensen_ppc = sorensen_ppc_tmp
												optimised_threshold = threshold
											}
									}
								tabs_list3[[j]] = tmp
								sorensen_ppcs[j] = sorensen_ppc
								thresholds[j] = optimised_threshold
							}
						tabs_list2[[i]] = tabs_list3
						medianV = round(median(sorensen_ppcs),2)
						minV = round(min(sorensen_ppcs),2)
						maxV = round(max(sorensen_ppcs),2)
						if (nchar(medianV) == 3) medianV = paste0(medianV,"0")
						if (nchar(medianV) == 1) medianV = paste0(medianV,".00")
						if (nchar(minV) == 3) minV = paste0(minV,"0")
						if (nchar(minV) == 1) minV = paste0(minV,".00")
						if (nchar(maxV) == 3) maxV = paste0(maxV,"0")
						if (nchar(maxV) == 1) maxV = paste0(maxV,".00")
						tab[i,((h-1)*2)+1] = paste0(medianV," [",minV,"-",maxV,"]")
						medianV = round(median(thresholds),2)
						minV = round(min(thresholds),2)
						maxV = round(max(thresholds),2)
						if (nchar(medianV) == 3) medianV = paste0(medianV,"0")
						if (nchar(medianV) == 1) medianV = paste0(medianV,".00")
						if (nchar(minV) == 3) minV = paste0(minV,"0")
						if (nchar(minV) == 1) minV = paste0(minV,".00")
						if (nchar(maxV) == 3) maxV = paste0(maxV,"0")
						if (nchar(maxV) == 1) maxV = paste0(maxV,".00")
						tab[i,((h-1)*2)+2] = paste0(medianV," [",minV,"-",maxV,"]")
					}
				tabs_list1[[h]] = tabs_list2
			}
		colnames(tab) = colNames; write.csv(tab, "All_SIppc_values.csv", quote=F)		
		pdf(paste0("All_the_figures_&_SI/All_SI_ppc_curves_NEW.pdf"), width=10, height=12)
		par(mfrow=c(8,6), oma=c(0,0,0,0), mar=c(2.5,2.5,0.5,0.5), lwd=0.4, col="gray30")
		for (i in 1:dim(species)[1]) # only for the first ISIMIP3a model (GSWP3-W5R5)
			{
				plot(tabs_list1[[1]][[i]][[1]], col=NA, ann=F, axes=F, xlim=c(0,1), ylim=c(0,1))
				for (j in 1:length(tabs_list1[[1]][[i]])) lines(tabs_list1[[2]][[i]][[j]], lwd=0.3, col="gray50", lty=1)
				axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.07,0), at=c(0,0.2,0.4,0.6,0.8,1), label=c("0","0.2","0.4","0.6","0.8","1"))
				axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.30,0), at=c(0,0.2,0.4,0.6,0.8,1), label=c("0","0.2","0.4","0.6","0.8","1"))
				if (i %in% c(1,7,13,19,25,31,37,43)) title(ylab=expression("SI"["ppc"]), cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
				if (i %in% c(43,44,45,46,47)) title(xlab="Threshold", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
				box(lwd=0.2, col="gray30"); SIppc = unlist(strsplit(tab[i,1]," "))[1]; threshold = unlist(strsplit(tab[i,2]," "))[1]
				mtext(paste0(gsub("Bombus_","B. ",species[i,1])), side=3, line=-6.9, at=0.03, cex=0.55, col="gray30", adj=0)
				mtext(paste0("SIppc = ",SIppc," (",threshold,")"), side=3, line=-7.8, at=0.03, cex=0.50, col="gray30", adj=0)
			}
		dev.off()
	}

# 5. Computation and analysis of the relative influence of each environmental factor

for (h in 1:length(models_isimip3a))
	{
		fileName = paste0("BRT_projection_files/RI_estimates/RI_",models_isimip3a_names[h],".csv")
		if (!file.exists(paste0(fileName)))
			{
				relativeInfluences = matrix(0, nrow=dim(species)[1], ncol=length(envVariables_list[[h]]))
				row.names(relativeInfluences) = species[,"species"]
				for (i in 1:dim(species)[1])
					{
						brt_models = readRDS(paste0("BRT_projection_files/BRT_models/",species[i,1],"_",models_isimip3a_names[h],"_10_SCV2.rds"))				
						for (j in 1:length(brt_models))
							{
								tab = summary(brt_models[[j]]); row.names(tab) = gsub("`","",row.names(tab))
								for (k in 1:length(envVariables_list[[h]]))
									{
										relativeInfluences[i,k] = as.numeric(relativeInfluences[i,k]) + tab[envVariableNames[k],"rel.inf"]
									}
							}
						if (i == 1) colnames(relativeInfluences) = envVariableNames
						vS = round(as.numeric(relativeInfluences[i,])/length(brt_models),2)
						for (j in 1:length(vS))
							{
								if (nchar(vS[j]) == 3) vS[j] = paste0(vS[j],"0")
								if (!grepl("\\.",vS[j])) vS[j] = paste0(vS[j],".00")
							}
						relativeInfluences[i,] = vS
					}
				write.table(relativeInfluences, fileName, quote=F, sep=",")
			}
	}

# 6. BRT projections based on historical and counterfactual climate simulations

envVariables_list_1 = list()
for (g in 1:length(models_isimip3a))
	{
		envVariables_list_2 = list()
		for (h in 1:length(scenarios))
			{
				envVariables_list_3 = list()
				for (i in 1:length(pastPeriods))
					{
						temperature = brick(paste0("Environmental_rasters/ISIMIP3a/tas_day_",scenarios[h],"_historical_",models_isimip3a[g],"_",gsub("-","_",pastPeriods[i]),"_ymonmean.nc"))
						precipitation = brick(paste0("Environmental_rasters/ISIMIP3a/pr_day_",scenarios[h],"_historical_",models_isimip3a[g],"_",gsub("-","_",pastPeriods[i]),"_ymonmean.nc"))
						relative_humidity = brick(paste0("Environmental_rasters/ISIMIP3a/hurs_day_",scenarios[h],"_historical_",models_isimip3a[g],"_",gsub("-","_",pastPeriods[i]),"_ymonmean.nc"))
						if (pastPeriods[i] != "2000-2019")
							{
								land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_",gsub("-","_",pastPeriods[i]),"_timmean.nc4"))
							}	else	{
								land_cover = nc_open(paste0("Environmental_rasters/ISIMIP3a/landcover_annual_",gsub("-","_",pastPeriods[i]),"_timmean.nc"))
							}
						population = brick(paste0("Environmental_rasters/ISIMIP3a/population_histsoc_0p5deg_annual_",gsub("-","_",pastPeriods[i]),"_timmean.nc4"), varname="total-population")		
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
						envVariables[[13]] = land_covers3[[4]] # primary forest areas
						envVariables[[14]] = land_covers3[[5]] # primary non-forest areas
						envVariables[[15]] = land_covers3[[6]] # secondary forest areas
						envVariables[[16]] = land_covers3[[7]] # secondary non-forest areas
						envVariables[[17]] = land_covers3[[1]] # croplands (all catergories)
						envVariables[[18]] = land_covers3[[2]] # managed pasture + rangeland
						envVariables[[19]] = population_log10 # human population (log-transformed)
						for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], extentOfStudyArea, snap="out")
						for (j in 1:12) envVariables[[j]][which(is.na(envVariables[[19]][]))] = NA 
						for (j in 1:length(envVariables)) envVariables[[j]] = crop(envVariables[[j]], europe3, snap="out")
						for (j in 1:length(envVariables)) envVariables[[j]] = mask(envVariables[[j]], europe3)
						envVariables_list_3[[i]] = envVariables
					}
				envVariables_list_2[[h]] = envVariables_list_3
			}
		envVariables_list_1[[g]] = envVariables_list_2
	}
r = envVariables[[1]]; r[!is.na(r[])] = 1
contour = rasterToPolygons(r>0, dissolve=T)

projections_list_1 = list()
for (g in 1:length(models_isimip3a))
	{
		projections_list_2 = list()
		for (h in 1:length(scenarios))
			{
				projections_list_3 = list()
				for (i in 1:length(pastPeriods))
					{
						projections_list_4 = list(); errorInRasterName = FALSE
						rasters_stack = stack(envVariables_list_1[[g]][[h]][[i]]); names(rasters_stack) = envVariableNames
						first_brt_model = readRDS(paste0("BRT_projection_files/BRT_models/",species[1,1],"_",models_isimip3a_names[g],"_10_SCV2.rds"))[[1]]
						if ("tprecipitation_spring" %in% colnames(first_brt_model$data$x.order)) errorInRasterName = TRUE
						for (j in 1:length(names(rasters_stack)))
							{
								if (names(rasters_stack)[j] == "primary_non.forest_areas") names(rasters_stack)[j] = "primary_non-forest_areas"
								if (names(rasters_stack)[j] == "secondary_non.forest_areas") names(rasters_stack)[j] = "secondary_non-forest_areas"
								if (errorInRasterName)
									{
										if (names(rasters_stack)[j] == "precipitation_spring") names(rasters_stack)[j] = "tprecipitation_spring"
									}
							}
						for (j in 1:dim(species)[1])
							{
								replicates = list()
								brt_models = readRDS(paste0("BRT_projection_files/BRT_models/",species[j,1],"_",models_isimip3a_names[g],"_10_SCV2.rds"))
								for (k in 1:length(brt_models))
									{
										df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df)))
										for (l in 1:dim(df)[2])
											{
												if (colnames(df)[l] == "primary_non.forest_areas") colnames(df)[l] = "primary_non-forest_areas"
												if (colnames(df)[l] == "secondary_non.forest_areas") colnames(df)[l] = "secondary_non-forest_areas"
											}
										newdata = df[not_NA,]; n.trees = brt_models[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
										projection = predict.gbm(brt_models[[k]], newdata, n.trees, type, single.tree)
										rast = rasters_stack[[1]]; rast[!is.na(rast[])] = projection; replicates[[k]] = rast
									}
								rasts = stack(replicates); projections_list_4[[j]] = mean(rasts)
							}
						projections_list_3[[i]] = projections_list_4
					}
				projections_list_2[[h]] = projections_list_3
			}
		projections_list_1[[g]] = projections_list_2
	}

colourScales = list()
colourScales[[1]] = c(rep("#F2F4F4",10),rev(hcl.colors(100,palette="terrain2")[1:90]))
colourScales[[2]] = c(rep("#F2F4F4",10),rev(hcl.colors(100,palette="terrain2")[1:90]))
colourScales[[3]] = colorRampPalette(brewer.pal(11,"RdBu"))(100)
for (g in 1:length(models_isimip3a))
	{
		for (j in 1:dim(species)[1])
			{
				pdf(paste0("Species_estimations/",species[j,1],"_",models_isimip3a_names[g],".pdf"), width=8, height=5.8)
				par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); vS1 = c(); vS2 = c()
				for (h in 1:length(scenarios))
					{
						for (i in 1:length(pastPeriods)) vS1 = c(vS1, projections_list_1[[g]][[h]][[i]][[j]][])
					}
				for (i in 1:length(pastPeriods))
					{
						vS2 = c(vS2, projections_list_1[[g]][[1]][[i]][[j]][]-projections_list_1[[g]][[2]][[i]][[j]][])
					}
				vS1 = vS1[which(!is.na(vS1))]; vS2 = vS2[which(!is.na(vS2))]
				minVS2 = min(vS2); maxVS2 = max(vS2)
				if (abs(minVS2) < abs(maxVS2)) minVS2 = -maxVS2
				if (abs(maxVS2) < abs(minVS2)) maxVS2 = -minVS2
				for (h in 1:length(scenarios))
					{
						for (i in 1:length(pastPeriods))
							{				
								index1 = (((min(projections_list_1[[g]][[h]][[i]][[j]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
								index2 = (((max(projections_list_1[[g]][[h]][[i]][[j]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
								plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[h]][index1:index2]; rast = raster(as.matrix(c(min(vS1),max(vS1))))
								plot(projections_list_1[[g]][[h]][[i]][[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
								if (h == 1) mtext("Historical reconstruction", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
								if (h == 2) mtext("Counterfactual baseline", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
								mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.0, cex=0.50, col="gray30")
								plot(rast, legend.only=T, add=T, col=colourScales[[h]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
									 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
									 at=c(0.2,0.4,0.6,0.8)), alpha=1, side=3)
							}
					}
				for (i in 1:length(pastPeriods))
					{
						difference_obsclim_counterclim = projections_list_1[[g]][[1]][[i]][[j]]
						difference_obsclim_counterclim[] = difference_obsclim_counterclim[]-projections_list_1[[g]][[2]][[i]][[j]][]
						index1 = (((min(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
						index2 = (((max(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
						plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[3]][index1:index2]; rast = raster(as.matrix(c(minVS2,maxVS2)))
						plot(difference_obsclim_counterclim, col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
						mtext("Historical - counterfactual", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
						mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.5, cex=0.50, col="gray30")
						plot(rast, legend.only=T, add=T, col=colourScales[[3]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
							 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
							 at=c(-0.6,-0.3,0,0.3,0.6)), alpha=1, side=3)
					}
				dev.off()
			}
	}

# 7. Computation and mapping of the ESI and SRI metrics for the different time periods

SIppcs = read.csv("All_SIppc_values.csv", head=T)
ESI_list_1 = list(); SRI_list_1 = list()
for (g in 1:length(models_isimip3a))
	{
		ESI_list_2 = list(); SRI_list_2 = list()
		for (h in 1:length(scenarios))
			{
				ESI_list_3 = list(); SRI_list_3 = list()
				for (i in 1:length(pastPeriods))
					{
						counterRasters = list(); bufferRasters = list()
						for (j in 1:dim(species)[1])
							{
								column = paste0("OptTres_",gsub("-","\\.",models_isimip3a_names[g]))
								cutOff = as.numeric(unlist(strsplit(SIppcs[which(SIppcs[,1]==species[j,1]),column]," \\["))[1])
								c = projections_list_1[[g]][[h]][[i]][[j]]
								c[c[]<cutOff] = 0; c[c[]>cutOff] = 1; counterRasters[[j]] = c
								bufferRasters[[j]] = projections_list_1[[g]][[h]][[i]][[j]]
							}
						ESI_list_3[[i]] = mean(stack(bufferRasters)); SRI_list_3[[i]] = sum(stack(counterRasters))
					}
				ESI_list_2[[h]] = ESI_list_3; SRI_list_2[[h]] = SRI_list_3
			}
		ESI_list_1[[g]] = ESI_list_2; SRI_list_1[[g]] = SRI_list_2
	}
r = envVariables[[1]]; r[!is.na(r[])] = 1
contour = rasterToPolygons(r>0, dissolve=T); colourScales = list()
colourScales[[1]] = rev(hcl.colors(100,palette="terrain2")[1:100])
colourScales[[2]] = rev(hcl.colors(100,palette="terrain2")[1:100])
colourScales[[3]] = colorRampPalette(brewer.pal(11,"RdBu"))(100)
for (g in 1:length(models_isimip3a))
	{
		pdf(paste0("All_the_figures_&_SI/ESI_evolutions_m",g,".pdf"), width=8, height=5.8)
		par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); vS1 = c(); vS2 = c()
		for (h in 1:length(scenarios))
			{
				for (i in 1:length(pastPeriods)) vS1 = c(vS1, ESI_list_1[[g]][[h]][[i]][])
			}
		for (i in 1:length(pastPeriods))
			{
				vS2 = c(vS2, ESI_list_1[[g]][[1]][[i]][]-ESI_list_1[[g]][[2]][[i]][])
			}
		vS1 = vS1[which(!is.na(vS1))]; vS2 = vS2[which(!is.na(vS2))]
		minVS2 = min(vS2); maxVS2 = max(vS2)
		if (abs(minVS2) < abs(maxVS2)) minVS2 = -maxVS2
		if (abs(maxVS2) < abs(minVS2)) maxVS2 = -minVS2
		if (maxVS2 < 0.1) { minVS2 = -0.1; maxVS2 = 0.1 }
		for (h in 1:length(scenarios))
			{
				for (i in 1:length(pastPeriods))
					{				
						index1 = (((min(ESI_list_1[[g]][[h]][[i]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
						index2 = (((max(ESI_list_1[[g]][[h]][[i]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
						par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
						plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[h]][index1:index2]; rast = raster(as.matrix(c(min(vS1),max(vS1))))
						plot(projections_list_1[[g]][[h]][[i]][[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
						if (h == 1) mtext("Historical reconstruction", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
						if (h == 2) mtext("Counterfactual baseline", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
						mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.9, at=5.0, cex=0.50, col="gray30")
						par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
						plot(rast, legend.only=T, add=T, col=colourScales[[h]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
							 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
							 at=c(0.2,0.4,0.6,0.8)), alpha=1, side=3)
					}
			}
		for (i in 1:length(pastPeriods))
			{
				difference_obsclim_counterclim = ESI_list_1[[g]][[1]][[i]]
				difference_obsclim_counterclim[] = difference_obsclim_counterclim[]-ESI_list_1[[g]][[2]][[i]][]
				index1 = (((min(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
				index2 = (((max(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
				par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[3]][index1:index2]; rast = raster(as.matrix(c(minVS2,maxVS2)))
				plot(difference_obsclim_counterclim, col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Historical - counterfactual", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
				mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.9, at=5.5, cex=0.50, col="gray30")
				par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
				plot(rast, legend.only=T, add=T, col=colourScales[[3]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
					 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
					 at=c(-0.2,-0.1,0,0.1,0.2)), alpha=1, side=3)
			}
		dev.off()
	}
for (g in 1:length(models_isimip3a))
	{
		pdf(paste0("All_the_figures_&_SI/SRI_evolutions_m",g,".pdf"), width=8, height=5.8); vS1 = c(); vS2 = c()
		par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
		for (h in 1:length(scenarios))
			{
				for (i in 1:length(pastPeriods)) vS1 = c(vS1, SRI_list_1[[g]][[h]][[i]][])
			}
		for (i in 1:length(pastPeriods))
			{
				vS2 = c(vS2, SRI_list_1[[g]][[1]][[i]][]-SRI_list_1[[g]][[2]][[i]][])
			}
		vS1 = vS1[which(!is.na(vS1))]; vS2 = vS2[which(!is.na(vS2))]
		minVS2 = min(vS2); maxVS2 = max(vS2)
		if (abs(minVS2) < abs(maxVS2)) minVS2 = -maxVS2
		if (abs(maxVS2) < abs(minVS2)) maxVS2 = -minVS2
		for (h in 1:length(scenarios))
			{
				for (i in 1:length(pastPeriods))
					{				
						index1 = (((min(SRI_list_1[[g]][[h]][[i]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
						index2 = (((max(SRI_list_1[[g]][[h]][[i]][],na.rm=T)-min(vS1))/(max(vS1)-min(vS1)))*100)+1
						par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
						plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[h]][index1:index2]; rast = raster(as.matrix(c(min(vS1),max(vS1))))
						plot(projections_list_1[[g]][[h]][[i]][[j]], col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
						if (h == 1) mtext("Historical reconstruction", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
						if (h == 2) mtext("Counterfactual baseline", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
						mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.0, cex=0.50, col="gray30")
						par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
						plot(rast, legend.only=T, add=T, col=colourScales[[h]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
							 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
							 at=c(0,10,20,30,40,50)), alpha=1, side=3)
					}
			}
		for (i in 1:length(pastPeriods))
			{
				difference_obsclim_counterclim = SRI_list_1[[g]][[1]][[i]]
				difference_obsclim_counterclim[] = difference_obsclim_counterclim[]-SRI_list_1[[g]][[2]][[i]][]
				index1 = (((min(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
				index2 = (((max(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
				par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
				plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[3]][index1:index2]; rast = raster(as.matrix(c(minVS2,maxVS2)))
				plot(difference_obsclim_counterclim, col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
				mtext("Historical - counterfactual", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
				mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.5, cex=0.50, col="gray30")
				par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
				plot(rast, legend.only=T, add=T, col=colourScales[[3]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
					 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
					 at=c(-20,-10,0,10,20)), alpha=1, side=3)
			}
		dev.off()
	}

for (g in 1:length(models_isimip3a)) print(round(mean(((ESI_list_1[[g]][[1]][[1]][]-ESI_list_1[[g]][[1]][[length(pastPeriods)]][])/ESI_list_1[[g]][[1]][[1]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(mean(((ESI_list_1[[g]][[2]][[1]][]-ESI_list_1[[g]][[2]][[length(pastPeriods)]][])/ESI_list_1[[g]][[2]][[1]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(mean(((ESI_list_1[[g]][[2]][[length(pastPeriods)]][]-ESI_list_1[[g]][[1]][[length(pastPeriods)]][])/ESI_list_1[[g]][[2]][[length(pastPeriods)]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(max(((ESI_list_1[[g]][[1]][[length(pastPeriods)]][]-ESI_list_1[[g]][[2]][[length(pastPeriods)]][])/ESI_list_1[[g]][[1]][[length(pastPeriods)]][])*100,na.rm=T),2))

for (g in 1:length(models_isimip3a)) print(round(mean(((SRI_list_1[[g]][[1]][[1]][]-SRI_list_1[[g]][[1]][[length(pastPeriods)]][])/SRI_list_1[[g]][[1]][[1]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(mean(((SRI_list_1[[g]][[2]][[1]][]-SRI_list_1[[g]][[2]][[length(pastPeriods)]][])/SRI_list_1[[g]][[2]][[1]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(mean(((SRI_list_1[[g]][[2]][[length(pastPeriods)]][]-SRI_list_1[[g]][[1]][[length(pastPeriods)]][])/SRI_list_1[[g]][[2]][[length(pastPeriods)]][])*100,na.rm=T),2))
for (g in 1:length(models_isimip3a)) print(round(max(((SRI_list_1[[g]][[1]][[length(pastPeriods)]][]-SRI_list_1[[g]][[2]][[length(pastPeriods)]][])/SRI_list_1[[g]][[1]][[length(pastPeriods)]][])*100,na.rm=T),2))

g = 2; vS = (SRI_list_1[[g]][[2]][[length(pastPeriods)]][]-SRI_list_1[[g]][[1]][[length(pastPeriods)]][])/SRI_list_1[[g]][[2]][[length(pastPeriods)]][]
vS = vS[which(!is.na(vS))]; vS = vS[which(is.finite(vS))]; print(round(mean(vS*100),2)) # (*)

	# 4.34% (GSWP3-W5E5), 4.00% (20CRv3), 5.11% (20CRv3-ERA5), 5.81% (20CRv3-W5E5) --> an average loss of 4-6% of local ecological suitability between 2000-2019 and 1900-1919 (obsclim)
	# -0.70% (GSWP3-W5E5), 0.22% (20CRv3), 0.11% (20CRv3-ERA5), 0.04% (20CRv3-W5E5) --> an average gain of ~0% of local ecological suitability between 2000-2019 and 1900-1919 (counterclim)
	# 5.17% (GSWP3-W5E5), 4.09% (20CRv3), 5.36% (20CRv3-ERA5), 6.20% (20CRv3-W5E5) --> an average loss of 4-6% of local ecological suitability solely due to climate change
	# 18.79% (GSWP3-W5E5), 14.69% (20CRv3), 16.53% (20CRv3-ERA5), 17.83% (20CRv3-W5E5) --> an average loss up to 15-19% of local ecological suitability solely due to climate change
	# 6.72 (GSWP3-W5E5), 4.34% (20CRv3), 5.76% (20CRv3-ERA5), 6.71% (20CRv3-W5E5) --> an average loss of 4-7% of local species diversity between 2000-2019 and 1900-1919 (obsclim)
	# -1.79% (GSWP3-W5E5), -0.64% (20CRv3), -1.67% (20CRv3-ERA5), -1.70% (20CRv3-W5E5) --> an average gain of 1-2% of local species diversity between 2000-2019 and 1900-1919 (counterclim)
	# 8.16% (GSWP3-W5E5), 4.84%* (20CRv3), 6.94% (20CRv3-ERA5), 8.30% (20CRv3-W5E5) --> an average loss of 5-8% of local species diversity solely due to climate change
	# 60.00% (GSWP3-W5E5), 100.00% (20CRv3), 66.67% (20CRv3-ERA5), 80.00% (20CRv3-W5E5) --> a loss up to 60-100% of local species diversity solely due to climate change

avg_ES_loss = matrix(nrow=dim(species)[1], ncol=length(models_isimip3a)); max_ES_loss = matrix(nrow=dim(species)[1], ncol=length(models_isimip3a)) # for Bastien
row.names(avg_ES_loss) = species[,1]; colnames(avg_ES_loss) = models_isimip3a_names; row.names(max_ES_loss) = species[,1]; colnames(max_ES_loss) = models_isimip3a_names
for (g in 1:length(models_isimip3a))
	{
		for (i in 1:dim(species)[1])
			{
				avg_ES_loss[i,g] = round(mean(((projections_list_1[[g]][[2]][[length(pastPeriods)]][[i]][]-projections_list_1[[g]][[1]][[length(pastPeriods)]][[i]][])
				/projections_list_1[[g]][[2]][[length(pastPeriods)]][[i]][])*100,na.rm=T),2) # ((mean_ES_counterclim_2000_2019 - mean_ES_obsclim_2000_2019)/mean_ES_counterclim_2000_2019)*100
				max_ES_loss[i,g] = round(max(((projections_list_1[[g]][[2]][[length(pastPeriods)]][[i]][]-projections_list_1[[g]][[1]][[length(pastPeriods)]][[i]][])
				/projections_list_1[[g]][[2]][[length(pastPeriods)]][[i]][])*100,na.rm=T),2) # ((max_ES_obsclim_2000_2019 - max_ES_counterclim_2000_2019)/max_ES_obsclim_2000_2019)*100
		    }
	}
write.csv(avg_ES_loss, "Average_ES_loss.csv", quote=F); write.csv(max_ES_loss, "Maximum_ES_loss.csv", quote=F)

target_countries = list(); country_ESIs = TRUE; country_SRIs = FALSE
target_countries[[1]] = shapefile("Countries_shapefiles/Spain_GADM_0.shp")
target_countries[[2]] = shapefile("Countries_shapefiles/France_GADM_0.shp")
target_countries[[3]] = shapefile("Countries_shapefiles/Sweden_GADM_0.shp")
for (i in 1:length(target_countries))
	{
		polIndex1 = NA; polIndex2 = NA; maxArea = -9999
		for (j in 1:length(target_countries[[i]]@polygons))
			{
				for (k in 1:length(target_countries[[i]]@polygons[[j]]@Polygons))
					{
						if (maxArea < target_countries[[i]]@polygons[[j]]@Polygons[[k]]@area)
							{
								polIndex1 = j; polIndex2 = k; maxArea = target_countries[[i]]@polygons[[j]]@Polygons[[k]]@area
							}
					}
			}
		pol = target_countries[[i]]@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = target_countries[[i]]@proj4string
		for (g in 1:length(models_isimip3a))
			{
				if (country_ESIs)
					{
						r1 = mask(crop(ESI_list_1[[g]][[1]][[length(pastPeriods)]], pol), pol)
						r2 = mask(crop(ESI_list_1[[g]][[2]][[length(pastPeriods)]], pol), pol)						
					}	else	{
						r1 = mask(crop(SRI_list_1[[g]][[1]][[length(pastPeriods)]], pol), pol)
						r2 = mask(crop(SRI_list_1[[g]][[2]][[length(pastPeriods)]], pol), pol)		
					}
				if ((i == 1)&(g == 2))
					{
						vS = (r2[]-r1[])/r2[]; vS = vS[which(!is.na(vS))]; vS = vS[which(is.finite(vS))]
						print(round(c(mean(vS),max((r2[]-r1[])/r2[],na.rm=T))*100,2)) # (*)
					}	else	{
						print(round(c(mean((r2[]-r1[])/r2[],na.rm=T),max((r2[]-r1[])/r2[],na.rm=T))*100,2))
					}
			}
	}
		# ESI, Spain: 11.40% (GSWP3-W5E5), 7.06% (20CRv3), 6.60% (20CRv3-ERA5), 6.85% (20CRv3-W5E5) --> an average loss of 7-11% of local species diversity solely due to climate change
		# ESI, France: 18.13% (GSWP3-W5E5), 11.31% (20CRv3), 18.05% (20CRv3-ERA5), 18.74% (20CRv3-W5E5) --> an average loss of 11-19% of local species diversity solely due to climate change
		# ESI, Sweden: -2.36% (GSWP3-W5E5), -0.03% (20CRv3), -0.92% (20CRv3-ERA5), 2.76% (20CRv3-W5E5) --> an average gain of -2 to 3% of local species diversity solely due to climate change
		# SRI, Spain: 22.87% (GSWP3-W5E5), 12.59* (20CRv3), 2.02% (20CRv3-ERA5), 8.95% (20CRv3-W5E5) --> an average loss of 10-20% of local species diversity solely due to climate change
		# SRI, France: 29.86% (GSWP3-W5E5), 15.93% (20CRv3), 30.19% (20CRv3-ERA5), 30.39% (20CRv3-W5E5) --> an average loss of 16-30% of local species diversity solely due to climate change
		# SRI, Sweden: -6.68% (GSWP3-W5E5), -0.33% (20CRv3), -5.32% (20CRv3-ERA5), -0.28% (20CRv3-W5E5) --> an average gain of 0 to 7% of local species diversity solely due to climate change

selected_regions = c("Alpine","Atlantic","Boreal","Continental","Mediterranean","Pannonian")
if (!file.exists("Biogeo_regions_shp/Biogeo_regions_2016.rds"))
	{
		biogeo_regions = shapefile("Biogeo_regions_shp/Biogeo_regions_2016.shp")
		biogeo_regions = subset(biogeo_regions, code%in%selected_regions)
		saveRDS(biogeo_regions, "Biogeo_regions_shp/Biogeo_regions_2016.rds")
	}	else	{
		biogeo_regions = readRDS("Biogeo_regions_shp/Biogeo_regions_2016.rds")
	}
biogeo_region_ESIs = TRUE; biogeo_region_SRIs = FALSE
for (i in 1:length(selected_regions))
	{
		pol = subset(biogeo_regions, code==selected_regions[i])
		pol = spTransform(pol, crs(ESI_list_1[[1]][[1]][[length(pastPeriods)]]))
		for (g in 1:length(models_isimip3a))
			{
				if (biogeo_region_ESIs)
					{
						r1 = mask(crop(ESI_list_1[[g]][[1]][[length(pastPeriods)]], pol), pol)
						r2 = mask(crop(ESI_list_1[[g]][[2]][[length(pastPeriods)]], pol), pol)						
					}	else	{
						r1 = mask(crop(SRI_list_1[[g]][[1]][[length(pastPeriods)]], pol), pol)
						r2 = mask(crop(SRI_list_1[[g]][[2]][[length(pastPeriods)]], pol), pol)		
					}
				print(round(c(mean((r2[]-r1[])/r2[],na.rm=T),max((r2[]-r1[])/r2[],na.rm=T))*100,2))
			}
	}
		# ESI, Alpine: 0.44% (GSWP3-W5E5), 1.21% (20CRv3), -0.73% (20CRv3-ERA5), 0.53% (20CRv3-W5E5)
		# ESI, Atlantic: 6.20% (GSWP3-W5E5), 6.73% (20CRv3), 7.18% (20CRv3-ERA5), 8.30% (20CRv3-W5E5)
		# ESI, Boreal: 0.02% (GSWP3-W5E5), -1.17% (20CRv3), 0.95% (20CRv3-ERA5), 4.40% (20CRv3-W5E5)
		# ESI, Continental: 7.79% (GSWP3-W5E5), 7.13% (20CRv3), 9.99% (20CRv3-ERA5), 9.05% (20CRv3-W5E5)
		# ESI, Mediterranean: 9.50% (GSWP3-W5E5), 5.14% (20CRv3), 4.58% (20CRv3-ERA5), 4.80% (20CRv3-W5E5)
		# ESI, Pannonian: 16.29% (GSWP3-W5E5), 10.71% (20CRv3), 16.90% (20CRv3-ERA5), 15.81% (20CRv3-W5E5)

pdf(paste0("All_the_figures_&_SI/ESI_&_SRI_GSWP3_NEW.pdf"), width=8, height=5.8); vS2 = c()
par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
for (i in 1:length(pastPeriods))
	{
		vS2 = c(vS2, ESI_list_1[[1]][[1]][[i]][]-ESI_list_1[[1]][[2]][[i]][])
	}
vS2 = vS2[which(!is.na(vS2))]
minVS2 = min(vS2); maxVS2 = max(vS2)
if (abs(minVS2) < abs(maxVS2)) minVS2 = -maxVS2
if (abs(maxVS2) < abs(minVS2)) maxVS2 = -minVS2
for (i in 1:length(pastPeriods))
	{
		difference_obsclim_counterclim = ESI_list_1[[1]][[1]][[i]]
		difference_obsclim_counterclim[] = difference_obsclim_counterclim[]-ESI_list_1[[1]][[2]][[i]][]
		index1 = (((min(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
		index2 = (((max(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
		par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
		plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[3]][index1:index2]; rast = raster(as.matrix(c(minVS2,maxVS2)))
		plot(difference_obsclim_counterclim, col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
		mtext("Historical - counterfactual", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
		mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.5, cex=0.50, col="gray30")
		par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
		plot(rast, legend.only=T, add=T, col=colourScales[[3]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
			 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
			 at=c(-0.2,-0.1,0,0.1,0.2)), alpha=1, side=3)
	}
vS2 = c()
for (i in 1:length(pastPeriods))
	{
		vS2 = c(vS2, SRI_list_1[[1]][[1]][[i]][]-SRI_list_1[[1]][[2]][[i]][])
	}
vS2 = vS2[which(!is.na(vS2))]
minVS2 = min(vS2); maxVS2 = max(vS2)
if (abs(minVS2) < abs(maxVS2)) minVS2 = -maxVS2
if (abs(maxVS2) < abs(minVS2)) maxVS2 = -minVS2
for (i in 1:length(pastPeriods))
	{
		difference_obsclim_counterclim = SRI_list_1[[1]][[1]][[i]]
		difference_obsclim_counterclim[] = difference_obsclim_counterclim[]-SRI_list_1[[1]][[2]][[i]][]
		index1 = (((min(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
		index2 = (((max(difference_obsclim_counterclim[],na.rm=T)-minVS2)/(maxVS2-minVS2))*100)+1
		par(lwd=0.2, col="gray30", col.axis="gray30", fg=NA)
		plot(europe3, lwd=0.1, border=NA, col=NA); cols = colourScales[[3]][index1:index2]; rast = raster(as.matrix(c(minVS2,maxVS2)))
		plot(difference_obsclim_counterclim, col=cols, border=NA, lwd=0.1, add=T, legend=F); plot(contour, lwd=0.4, border="gray50", col=NA, add=T)
		mtext("Historical - counterfactual", side=3, line=-1.1, at=3.5, cex=0.45, col="gray30")
		mtext(paste0("(",pastPeriods[i],")"), side=3, line=-1.8, at=5.5, cex=0.50, col="gray30")
		par(lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
		plot(rast, legend.only=T, add=T, col=colourScales[[3]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.085,0.69,0.89), adj=3,
			 axis.args=list(cex.axis=0.60, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.40,0),
			 at=c(-20,-10,0,10,20)), alpha=1, side=3)
	}
dev.off()

# 8. Computation and mapping of the evolution of the Boyce index (BI) through time

if (!file.exists("Boyce_index_lists.rds"))
	{
		BIs_list_1 = list(); g = 1; h = 1; i = 1; j = 1; k = 1; l = 1
		for (g in 1:length(models_isimip3a))
			{
				BIs_list_2 = list()
				for (h in 1:length(scenarios))
					{
						BIs_list_3 = list()
						for (i in 1:length(pastPeriods))
							{
								BIs_list_4 = list(); timePoints = as.numeric(unlist(strsplit(pastPeriods[i],"-"))); errorInRasterName = FALSE
								rasters_stack = stack(envVariables_list_1[[g]][[h]][[i]]); names(rasters_stack) = envVariableNames; print(c(g,h,i))
								first_brt_model = readRDS(paste0("BRT_projection_files/BRT_models/",species[1,1],"_",models_isimip3a_names[g],"_10_SCV2.rds"))[[1]]
								if ("tprecipitation_spring" %in% colnames(first_brt_model$data$x.order)) errorInRasterName = TRUE
								for (j in 1:length(names(rasters_stack)))
									{
										if (names(rasters_stack)[j] == "primary_non.forest_areas") names(rasters_stack)[j] = "primary_non-forest_areas"
										if (names(rasters_stack)[j] == "secondary_non.forest_areas") names(rasters_stack)[j] = "secondary_non-forest_areas"
										if (errorInRasterName)
											{
												if (names(rasters_stack)[j] == "precipitation_spring") names(rasters_stack)[j] = "tprecipitation_spring"
											}
									}
								j = 1
								for (j in 1:dim(species)[1])
									{
										BIs_list_5 = list()
										brt_models = readRDS(paste0("BRT_projection_files/BRT_models/",species[j,1],"_",models_isimip3a_names[g],"_10_SCV2.rds"))
										observations = read.csv(paste0(directory,"/",species[j,1],".csv"), header=T); data_to_discard = c()
										observations = observations[which((observations[,"year"]>=timePoints[1])&(observations[,"year"]<=timePoints[2])),c("longitude","latitude")]
										for (k in 1:length(rasters_stack@layers)) # to discard the observations that do not fall within the rasters
											{
												data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[k]],observations))))
											}
										data_to_discard = unique(data_to_discard); data_to_discard = data_to_discard[order(data_to_discard)]
										if (length(data_to_discard) > 0) observations = observations[which(!c(1:dim(observations)[1])%in%data_to_discard),]
										cellIDs = unique(cellFromXY(nullRaster, observations)); presence_points = xyFromCell(nullRaster, cellIDs)
										if (dim(presence_points)[1] >= 30)
											{
												for (k in 1:length(brt_models))
													{
														df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df)))
														for (l in 1:dim(df)[2])
															{
																if (colnames(df)[l] == "primary_non.forest_areas") colnames(df)[l] = "primary_non-forest_areas"
																if (colnames(df)[l] == "secondary_non.forest_areas") colnames(df)[l] = "secondary_non-forest_areas"
															}
														newdata = df[not_NA,]; n.trees = brt_models[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
														projection = predict.gbm(brt_models[[k]], newdata, n.trees, type, single.tree)
														rast = rasters_stack[[1]]; rast[!is.na(rast[])] = projection; projection = rast
														presence_values = raster::extract(rast, presence_points)
														BIs_list_6 = list(); l = 1 # 10 replicate for the selection of random points selected within the study area
														for (l in 1:10)
															{		
																random_points = spsample(contour, n=dim(presence_points)[1], type="random")@coords
																# plot(rast); points(presence_points); points(random_points, col="red")
																random_values = raster::extract(rast, random_points)
																intervals1 = seq(0,0.9,0.01); intervals2 = intervals1+0.1; Fs1 = rep(NA, length(intervals1))
																for (m in 1:length(intervals1))
																	{
																		P = sum((presence_values>intervals1[m])&(presence_values<=intervals2[m]))
																		E = sum((random_values>intervals1[m])&(random_values<=intervals2[m]))
																		Fs1[m] = P/E
																	}
																Fs2 = Fs1[which(!is.na(Fs1))]; intervals = c(1:length(intervals1))
																intervals = intervals[which(!is.na(Fs1))] # plot(c(1:length(intervals1)), Fs2)
																BIs_list_6[[l]] = cor(Fs2, intervals, method="spearman")
															}
														BIs_list_5[[k]] = BIs_list_6
													}
											}
										BIs_list_4[[j]] = BIs_list_5
									}
								BIs_list_3[[i]] = BIs_list_4
							}
						BIs_list_2[[h]] = BIs_list_3
					}
				BIs_list_1[[g]] = BIs_list_2
			}
		saveRDS(BIs_list_1, "Boyce_index_lists.rds")
	}	else	{
		BIs_list_1 = readRDS("Boyce_index_lists.rds")
	}
differences_list1 = list(); g = 1; h = 1; i = 1; j = 1; k = 1; l = 1
for (g in 1:length(models_isimip3a))
	{
		differences_list2 = list()
		for (i in 1:length(pastPeriods))
			{
				differences = matrix(nrow=dim(species)[1], ncol=10*10); timePoints = as.numeric(unlist(strsplit(pastPeriods[i],"-")))
				for (j in 1:dim(species)[1])
					{
						observations = read.csv(paste0(directory,"/",species[j,1],".csv"), header=T); data_to_discard = c()
						observations = observations[which((observations[,"year"]>=timePoints[1])&(observations[,"year"]<=timePoints[2])),c("longitude","latitude")]
						for (k in 1:length(rasters_stack@layers)) # to discard the observations that do not fall within the rasters
							{
								data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[k]],observations))))
							}
						data_to_discard = unique(data_to_discard); data_to_discard = data_to_discard[order(data_to_discard)]
						if (length(data_to_discard) > 0) observations = observations[which(!c(1:dim(observations)[1])%in%data_to_discard),]
						cellIDs = unique(cellFromXY(nullRaster, observations)); presence_points = xyFromCell(nullRaster, cellIDs)
						brt_models = readRDS(paste0("BRT_projection_files/BRT_models/",species[j,1],"_",models_isimip3a_names[g],"_10_SCV2.rds"))
						if (dim(presence_points)[1] >= 30)
							{
								for (k in 1:length(brt_models))
									{
										for (l in 1:10)
											{
												differences[j,((k-1)*10)+l] = BIs_list_1[[g]][[1]][[i]][[j]][[k]][[l]]-BIs_list_1[[g]][[2]][[i]][[j]][[k]][[l]]
											}	
									}
							}
					}
				differences_list2[[i]] = differences
			}
		differences_list1[[g]] = differences_list2
	}
yMin = 9999; yMax = -9999
for (g in 1:length(models_isimip3a))
	{
		mean_differences = matrix(nrow=dim(species)[1], ncol=length(pastPeriods))
		colnames(mean_differences) = pastPeriods
		for (i in 1:length(pastPeriods))
			{
				for (j in 1:dim(species)[1])
					{
						mean_differences[j,i] = mean(differences_list1[[g]][[i]][j,], na.rm=T)
					}
			}
		if (yMin > min(mean_differences,na.rm=T)) yMin = min(mean_differences,na.rm=T)
		if (yMax < max(mean_differences,na.rm=T)) yMax = max(mean_differences,na.rm=T)
	}
pdf(paste0("All_the_figures_&_SI/Mean_BI_differences_NEW.pdf"), width=8, height=2.5) # dev.new(width=8, height=2.5)
par(mfrow=c(1,4), oma=c(0,0.75,0,0), mar=c(2.5,3,1,1), lwd=0.4, col="gray30"); boxplots = FALSE
for (g in 1:length(models_isimip3a))
	{
		mean_differences = matrix(nrow=dim(species)[1], ncol=length(pastPeriods))
		colnames(mean_differences) = pastPeriods
		for (i in 1:length(pastPeriods))
			{
				for (j in 1:dim(species)[1])
					{
						mean_differences[j,i] = mean(differences_list1[[g]][[i]][j,], na.rm=T)
					}
			}
		if (boxplots) { dev.new(); boxplot(mean_differences) }
		buffer = c(); col1 = rgb(70,118,187,255,maxColorValue=255); col2 = rgb(70,118,187,130,maxColorValue=255) # blue
		# stripchart(data.frame(mean_differences), method="jitter", vertical=T, pch=16, cex=1, col=col2, ann=F, axes=F, ylim=c(yMin,yMax))
		for (i in 1:dim(mean_differences)[2]) buffer = rbind(buffer, cbind(mean_differences[,i], rep(colnames(mean_differences)[i],dim(mean_differences)[1])))
		buffer = data.frame(buffer); colnames(buffer) = c("mean_BI_difference","period")
		buffer$mean_BI_difference = as.numeric(buffer$mean_BI_difference); buffer$period = as.factor(buffer$period)
		beeswarm(formula=as.formula("mean_BI_difference ~ period"), data=buffer, method="compactswarm", pch=16, cex=0.5, spacing=0.95, col=col2, ann=F, axes=F)
		beeswarm(formula=as.formula("mean_BI_difference ~ period"), data=buffer, method="compactswarm", pch=1, lwd=2, cex=0.5, spacing=0.95, col=col1, ann=F, axes=F, add=T)
		abline(h=0, lty=2, lwd=0.5, col="gray50")
		axis(side=1, lwd.tick=0.4, cex.axis=0.9, lwd=0, tck=-0.040, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.40,0), at=c(1:6), label=pastPeriods)
		axis(side=2, lwd.tick=0.4, cex.axis=0.9, lwd=0, tck=-0.035, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.45,0))
		if (g == 1) title(ylab="Mean BI differences", cex.lab=1, mgp=c(1.8,0,0), col.lab="gray30")
		box(lwd=0.4, col="gray30")
	}
dev.off()
pdf(paste0("All_the_figures_&_SI/Mean_BI_detla_1901_NEW.pdf"), width=8, height=2.5) # dev.new(width=8, height=2.5)
par(mfrow=c(1,4), oma=c(0,0.75,0,0), mar=c(2.5,3,1,1), lwd=0.4, col="gray30"); boxplots = FALSE
for (g in 1:length(models_isimip3a))
	{
		mean_differences = matrix(nrow=dim(species)[1], ncol=length(pastPeriods))
		colnames(mean_differences) = pastPeriods
		for (i in 1:length(pastPeriods))
			{
				for (j in 1:dim(species)[1])
					{
						mean_differences[j,i] = mean(differences_list1[[g]][[i]][j,], na.rm=T)
					}
			}
		if (boxplots) { dev.new(); boxplot(mean_differences) }
		col = rgb(70,118,187,255,maxColorValue=255) # blue
		stripchart(data.frame(mean_differences), method="jitter", vertical=T, pch=16, cex=1, col=NA, ann=F, axes=F, ylim=c(yMin,yMax))
		for (i in 1:dim(mean_differences)[1])
			{
				xS = c(1:dim(mean_differences)[2])[which(!is.na(mean_differences[i,]))]
				yS = mean_differences[i,which(!is.na(mean_differences[i,]))]
				yS = yS-0; lines(xS, yS, col=col, lwd=0.2)	
			}
		abline(h=0, lty=2, lwd=0.5, col="gray50")
		axis(side=1, lwd.tick=0.4, cex.axis=0.9, lwd=0, tck=-0.040, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.40,0), at=c(1:6), label=pastPeriods)
		axis(side=2, lwd.tick=0.4, cex.axis=0.9, lwd=0, tck=-0.035, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.45,0))
		if (g == 1) title(ylab="Mean BI differences", cex.lab=1, mgp=c(1.8,0,0), col.lab="gray30")
		box(lwd=0.4, col="gray30")
	}
dev.off()
pValues_negative_BI_differences = matrix(nrow=dim(species)[1], ncol=length(models_isimip3a)*length(pastPeriods))
pValues_positive_BI_differences = matrix(nrow=dim(species)[1], ncol=length(models_isimip3a)*length(pastPeriods))
for (g in 1:length(models_isimip3a))
	{
		for (i in 1:length(pastPeriods))
			{
				if ((g == 1)&(i == 1)) colNames = c()
				colNames = c(colNames, paste0(models_isimip3a_names[g],"-",pastPeriods[i]))
			}
	}
row.names(pValues_negative_BI_differences) = species[,1]; colnames(pValues_negative_BI_differences) = colNames
row.names(pValues_positive_BI_differences) = species[,1]; colnames(pValues_positive_BI_differences) = colNames
for (g in 1:length(models_isimip3a))
	{
		pValues_WMW_negative_differences = matrix(nrow=dim(species)[1], ncol=length(pastPeriods))
		pValues_WMW_positive_differences = matrix(nrow=dim(species)[1], ncol=length(pastPeriods))
		colnames(pValues_WMW_negative_differences) = pastPeriods
		colnames(pValues_WMW_positive_differences) = pastPeriods
		for (i in 1:length(pastPeriods))
			{
				for (j in 1:dim(species)[1])
					{
						if (sum(is.na(differences_list1[[g]][[i]][j,])) != length(differences_list1[[g]][[i]][j,]))
							{
								pValues_WMW_negative_differences[j,i] = wilcox.test(differences_list1[[g]][[i]][j,], mu=0, alternative="less")$p.value
								pValues_WMW_positive_differences[j,i] = wilcox.test(differences_list1[[g]][[i]][j,], mu=0, alternative="greater")$p.value						
							}
					}
			}
		for (i in 1:length(pastPeriods)) # Benjamini-Hochberg correction
			{
				buffer1 = pValues_WMW_negative_differences[,i]; buffer2 = buffer1; 
				buffer2 = buffer2[order(buffer2)]; n = sum(!is.na(buffer2))
				for (j in 1:dim(pValues_WMW_negative_differences)[1])
					{
						if (!is.na(pValues_WMW_negative_differences[j,i]))
							{
								rank = which(buffer2 == pValues_WMW_negative_differences[j,i])[1]
								pValueText = as.character(round(pValues_WMW_negative_differences[j,i],3))
								if (nchar(pValueText) == 4) pValueText = paste0(pValueText,"0")
								if (nchar(pValueText) == 3) pValueText = paste0(pValueText,"00")
								if (pValueText == "1") pValueText = ">0.999"
								if (pValueText == "0") pValueText = "<0.001"
								corrected_cut_off_value = (rank/n)*0.05
								if (pValues_WMW_negative_differences[j,i] < corrected_cut_off_value)
									{
										pValueText = paste0(pValueText,"*")
									}
								pValues_negative_BI_differences[j,(length(pastPeriods)*(g-1))+i] = pValueText
							}
					}
				buffer1 = pValues_WMW_positive_differences[,i]; buffer2 = buffer1; 
				buffer2 = buffer2[order(buffer2)]; n = sum(!is.na(buffer2))
				for (j in 1:dim(pValues_WMW_positive_differences)[1])
					{
						if (!is.na(pValues_WMW_positive_differences[j,i]))
							{
								rank = which(buffer2 == pValues_WMW_positive_differences[j,i])[1]
								pValueText = as.character(round(pValues_WMW_positive_differences[j,i],3))
								if (nchar(pValueText) == 4) pValueText = paste0(pValueText,"0")
								if (nchar(pValueText) == 3) pValueText = paste0(pValueText,"00")
								if (pValueText == "1") pValueText = ">0.999"
								if (pValueText == "0") pValueText = "<0.001"
								corrected_cut_off_value = (rank/n)*0.05
								if (pValues_WMW_positive_differences[j,i] < corrected_cut_off_value)
									{
										pValueText = paste0(pValueText,"*")
									}
								pValues_positive_BI_differences[j,(length(pastPeriods)*(g-1))+i] = pValueText
							}
					}
			}
	}
write.csv(pValues_negative_BI_differences, "Negative_BI_diffs.csv", quote=F)
write.csv(pValues_positive_BI_differences, "Positive_BI_diffs.csv", quote=F)

