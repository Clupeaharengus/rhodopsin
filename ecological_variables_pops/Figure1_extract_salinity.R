library(data.table)
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(viridis)
library(scatterpie)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(ncdf4)

# import data -------------------------------------------------------------

r.her.coords <- read.csv("data/rh1_poolAF_secchi_sal.csv", header = TRUE)

#rasters from:
#http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=GLOBAL_ANALYSIS_FORECAST_PHY_001_024

# plot over geotiff -------------------------------------------------------



#small region:
rel <- raster("data/environmental_variables/salinity_ESA/global-analysis-forecast-phy-001-024_1548881869772.nc", varname = "so")

crs(rel)

sr <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projected_raster <- projectRaster(rel, crs = sr)
rel <- projected_raster
salt.rel <- rel
#crop raster
b <- as(extent(lon.ext,lat.ext), 'SpatialPolygons')
crs(b) <- crs(rel)
rb <- crop(rel, b)
rel <- rb

rel_spdf <- as(rel, "SpatialPixelsDataFrame")
rel <- as.data.frame(rel_spdf)
head(rel)
names(rel) <- c(variable, "x", "y")

#world <- map_data('world')
wrld_map <- readOGR("data/10m-admin-0-countries/10m_admin_0_countries.shp")

df_wrld_map <- fortify(wrld_map)

myPalette <- colorRampPalette(brewer.pal(11, "RdBu"))

d=data.frame(x1=c(15), x2=c(30.2), y1=c(52), y2=c(66))

q <-ggplot() +
  geom_raster(data = rel, aes_string(x = "x", y = "y", fill = variable)) +
  scale_fill_gradientn(colours = myPalette(11)) +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_map(data = df_wrld_map, map = df_wrld_map, aes(map_id=piece), fill="grey90", color="black") +
  xlim(lon.ext) +
  ylim(lat.ext) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(text = element_text(size=40), 
        axis.text=element_text(colour="black"),
        axis.ticks.length=unit(.65, "cm"),
        legend.position = c(0.95, 0.95)) #coord_map("bonne", lat0 = 3, lat1 = 30.5)
q

pdf("output/salinity_plotting/salinity_base.pdf", width=19.69333333333333, height=24.746666666666666, useDingbats = FALSE)
q
dev.off()

pdf("output/salinity_plotting/salinity_scale.pdf", width=5, height=7, useDingbats = FALSE)
q
dev.off()

# extract values ----------------------------------------------------------

#set extant
lon.ext <- c(3, 30.5)
lat.ext <- c(52, 68)
variable <- "Salinity"

pac <- raster("data/environmental_variables/salinity_ESA/global-analysis-forecast-phy-001-024_1548931193550.nc", varname = "so")
atl <- raster("data/environmental_variables/salinity_ESA/global-analysis-forecast-phy-001-024_1548930859285.nc", varname = "so")

i.af.her.coords <- r.her.coords %>% dplyr::select(-variable)

sp.her.coords <- i.af.her.coords
coordinates(sp.her.coords) <- ~ lon + lat
rasValue.atl <- extract(atl, sp.her.coords, buffer = 200000, fun = mean) #go in radius of 50km and average all cells found 
rasValue.pac <- extract(pac, sp.her.coords, buffer = 200000, fun = mean) #go in radius of 50km and average all cells found 

ras.her.coords <- cbind(i.af.her.coords,rasValue.atl)
ras.her.coords2 <- cbind(ras.her.coords,rasValue.pac)

o.ras.her.coords <- ras.her.coords2 %>% mutate (Salinity = ifelse(is.na(rasValue.pac), rasValue.atl, rasValue.pac)) %>% dplyr::select(-rasValue.pac, -rasValue.atl)


write.csv(o.ras.her.coords, "output/rh1_poolAF_secchi_sal_20190327.csv", row.names = FALSE)




