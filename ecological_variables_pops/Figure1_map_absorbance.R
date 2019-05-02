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
library(ggsn)


# Data input description --------------------------------------------------

#Generally:
#https://oceancolor.gsfc.nasa.gov/l3/
#Go here, select:
#Product status = Provisional
#Sensor = MODIS-Aqua
#Product = Total absorbtion 412nm GIOP
#Period = Entire mission composite (2002-2019)
#Resolution = 4km

#For extracting values, I had to get average values within 30km for points (~20%?) that were NA

# import data -------------------------------------------------------------

r.her.coords <- read.csv("output/rh1_poolAF_secchi_sal_20190327.csv", header = TRUE)
r.her.coords <- r.her.coords %>% rename(AF = F261Y)
pop.order <- fread("data/pop_order_df.txt", header = TRUE) %>% rename(num_pop = V1)

r.her.coords <- left_join(r.her.coords, pop.order, by=c("pop" = "pop_order")) %>% dplyr::select(-X)

# plot over geotiff -------------------------------------------------------

#get outline
wrld_map <- readOGR("data/10m-admin-0-countries/10m_admin_0_countries.shp")
df_wrld_map <- fortify(wrld_map)

#set extant
lon.ext <- c(3, 30.5)
lat.ext <- c(52, 68)
variable <- "absorbance"


#### VERY CONFUSING GEOTIFF STUFF 
# when I load the geotiff exported from QGIS (from .nc input from NASA website), the values range from 0,255 (should be 0,1)
# Literally the only useful resource I found on this: http://www.timassal.com/?p=859
# The source of the small function below
  
rel<- raster("data/environmental_variables/NASA_color/4Feb19_full_412_from_MODIS.tif")

rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

#run the fxn
rescaled_rel<-rasterRescale(rel)
values(rescaled_rel)[values(rescaled_rel) == 0] = NA

no_rescale_rel <- rel
rel <- rescaled_rel
crs(rel)
e <- extent(lon.ext,lat.ext)
rel <- crop(rel,e)
rel_spdf <- as(rel, "SpatialPixelsDataFrame")
rel <- as.data.frame(rel_spdf)

head(rel)
names(rel) <- c(variable, "x", "y")
rel <- rel %>% mutate(absorbance = ifelse(absorbance == 0, NA, absorbance))

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), bias=.5)

q <-ggplot() +
  geom_raster(data = rel, aes_string(x = "x", y = "y", fill = "absorbance")) +
  scale_fill_gradientn(colours = myPalette(11), na.value="white") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_map(data = df_wrld_map, map = df_wrld_map, aes(map_id=piece), fill="grey90", color="black") +
  scalebar(df_wrld_map, dist = 5, dd2km = TRUE, model = 'WGS84') +
  xlim(lon.ext) +
  ylim(lat.ext) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(text = element_text(size=40), 
        axis.text=element_text(colour="black"),
        axis.ticks.length=unit(.65, "cm"),
        legend.position = c(0.95, 0.95))

ggplot() + xlim(lon.ext) + ylim(lat.ext) + 
  ggsn:::scalebar(df_wrld_map, dist = 500,  model = 'WGS84', location="bottomright") +
  geom_map(data = df_wrld_map, map = df_wrld_map, aes(map_id=piece), fill="grey90", color="black") 

pdf("output/absorbance_plotting/abso_base_test.pdf", width=19.69333333333333, height=24.746666666666666, useDingbats = FALSE)
q
dev.off()

pdf("output/absorbance_plotting/abso_scale.pdf", width=5, height=7, useDingbats = FALSE)
q
dev.off()

#the remainder of figure generation was done in Adobe illustrator by combining this with the output of Figure1_map_pi.R

# extract values ----------------------------------------------------------
# with full earth
rel2 <- raster("data/environmental_variables/NASA_color/4Feb19_full_412_from_MODIS.tif")
rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

#run the fxn
rescaled_rel2<-rasterRescale(rel2)
values(rescaled_rel2)[values(rescaled_rel2) == 0] = NA

i.af.her.coords <- r.her.coords #%>% select(-variable)
sp.her.coords <- i.af.her.coords
coordinates(sp.her.coords) <- ~ lon + lat
rasValue <- extract(rescaled_rel2, sp.her.coords, buffer = 100000, fun = mean) #go in radius of 30km and average all cells found 
ras.her.coords <- cbind(i.af.her.coords,rasValue)

o.ras.her.coords <- ras.her.coords %>% rename(Absorbance_412nm = rasValue)
hist(o.ras.her.coords$Absorbance_412nm)
plot(o.ras.her.coords$Absorbance_412nm, o.ras.her.coords$AF)

#precise location for pacific sample is unknown
o.ras.her.coords <- o.ras.her.coords %>% mutate(Salinity = ifelse(name == "Pacific", NA, Salinity)) %>%
  mutate(secchi = ifelse(name == "Pacific", NA, secchi)) %>%
  rename(F261Y = AF)

write.csv(o.ras.her.coords, "output/rh1_poolAF_secchi_sal_abs_20190327_100km.csv", row.names = FALSE)
