
# Dependencies
library(maps)
library(raster)
library(maptools)
library(dplyr)
library(RColorBrewer)
library(sf) # Simple feature geometries
library(rnaturalearth) # National boundaries from https://www.naturalearthdata.com/
library(rnaturalearthdata) 
library(cowplot) # Added simplicity for ggplot themes and export
library(ggrepel) # Avoid overlap when placing ggplot labels


# organize Site coordinates
subsites<-read.csv("data/LOOKUPS/subsites.csv") 
sites <- filter(subsites, site_name=="ALEXFIORD"| site_name=="ADVENTDALEN"| site_name=="ATQASUK"|
site_name=="BARROW"|  site_name=="DARING"|  site_name=="ENDALEN" |site_name=="FINSE"| site_name=="FAROE"|
site_name=="GAVIAPASS"| site_name=="HEALY"| site_name=="IMNAVAIT"|site_name=="JAKOBSHORN"| site_name=="KANGER"|
site_name=="LATNJA"| site_name=="NIWOT"| site_name=="TOOLIK"| site_name=="VALBERCLA"| site_name=="WHITEMTNS")
sites<-dplyr::select(sites, site_name, lat,long, subsite)%>%group_by(site_name)%>%
  mutate(subsite_no=n_distinct(subsite))%>%mutate(lat=mean(lat, na.rm=T), long=mean(long, na.rm=T))%>%
  dplyr::select(-subsite)%>%distinct(.)%>%arrange(site_name, long, lat, subsite_no)

sites$site_name_long <- sort(c(
  "Alexandra Fiord", 
  "Endalen", 
  "Adventdalen",   
  "Barrow-Utqiagvik",
  "Atqasuk", 
  "Toolik Lake", 
  "Imnavait Creek", 
  "Latnjajaure", 
  "Kangerlussuaq", 
  "Daring Lake", 
  "Healy", 
  "Faroe Islands", 
  "Finse", 
  "Jakobshorn", 
  "Val Bercla", 
  "Gavia Pass", 
  "Niwot Ridge", 
  "White Mountains"))

sites$site_name_long[sites$site_name == "BARROW"] <- "Utqiagvik"

# Specify a lambert conformal conical projection, which works nicely
# for the area covered by the sites.
# I use a WKT specification here as this is most easy to understand
# You can change the "Central_Meridian" parameter to shift the central
# meridian of the map: -96 for North America, -40 For Greenland, 
# 0 for Greenwich

lambert_conformic <- 'PROJCS["North_America_Lambert_Conformal_Conic",
                            GEOGCS["GCS_North_American_1983",
                                   DATUM["North_American_Datum_1983",
                                         SPHEROID["GRS_1980",6378137,298.257222101]],
                                   PRIMEM["Greenwich",0],
                                   UNIT["Degree",0.017453292519943295]],
                            PROJECTION["Lambert_Conformal_Conic_2SP"],
                            PARAMETER["False_Easting",0],
                            PARAMETER["False_Northing",0],
                            PARAMETER["Central_Meridian",-40], 
                            PARAMETER["Standard_Parallel_1",20],
                            PARAMETER["Standard_Parallel_2",60],
                            PARAMETER["Latitude_Of_Origin",40],
                            UNIT["Meter",1],
                            AUTHORITY["EPSG","102009"]]'

# Generate political boundaries as sf and transform from lat long to lambert
world <- ne_countries(scale = "medium", returnclass = "sf")
world_lambert <- st_transform(world, crs = st_crs(lambert_conformic))

# Generate sf object for site point sources and transform
sites_sf <- st_as_sf(sites, coords = c("long", "lat"), crs = 4326)
sites_sf <- st_transform(sites_sf , crs = st_crs(lambert_conformic))

# Next we need to define the area we want to plot
# For this we get the minimum and maximum coordinates and 
# add a good buffer, I started the buffer with a guess and then
# modified till it looked good
map_extent <- sites_sf %>% st_coordinates() %>%
  as.data.frame() %>%
  summarise(xmin = min(X) - 1000000,
            xmax = max(X) + 4000000,
            ymin = min(Y) - 2000000,
            ymax = max(Y) + 1000000)

# Plot data using ggplot
otc_site_map <- ggplot() + 
  # geom_sf to plot the country boundaries
  geom_sf(data = world_lambert,
          fill = "#d3d3d3FF", 
          colour = "white",
          size = 0.25) + 
  # geom_sf to plot the sites (it automatically knows they're points)
  geom_sf(data = sites_sf, colour = "#498CB8FF",
          fill = "#ffffff00", shape = 21, size = 1.5,
          stroke = 1.5) +  
  # geom_label_repel with st_coordinates as stats to plot labels 
  # that are not overlapping
  geom_label_repel(data = sites_sf, aes(label = site_name_long,
                                        geometry = geometry),
                   stat = "sf_coordinates",
                   size = 3,
                   min.segment.length = 0,
                   segment.colour = "black", # Just in case you would like to modify these
                   colour = "black") + # Just in case you would like to modify these
  labs(x = "", 
       y = "") + 
  # set extent to be plotted using previously determined coordinates
  coord_sf(xlim = c(map_extent$xmin,
                    map_extent$xmax),
           ylim = c(map_extent$ymin,
                    map_extent$ymax), 
           expand = F) +
  # Map relevant modifiers for the theme
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        axis.title = element_blank(),
        panel.grid.minor = element_line(colour = "white"), # Just in case you would like to modify these
        panel.grid.major = element_line(colour = "white"), # Just in case you would like to modify these
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = element_rect(fill = "white"), # Just in case you would like to modify these
        panel.border = element_rect(colour = "lightgrey", fill=NA, size=0.5)
  ) 

# A quick plot on the native graphics device
otc_site_map


