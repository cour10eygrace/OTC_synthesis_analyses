
# Dependencies
library(maps)
library(raster)
library(maptools)
library(dplyr)
library(RColorBrewer)

# organize Site coordinates
subsites<-read.csv("data/LOOKUPS/subsites.csv") 
sites <- filter(subsites, site_name=="ALEXFIORD"| site_name=="ADVENTDALEN"| site_name=="ATQASUK"|
site_name=="BARROW"|  site_name=="DARING"|  site_name=="ENDALEN" |site_name=="FINSE"| site_name=="FAROE"|
site_name=="GAVIAPASS"| site_name=="HEALY"| site_name=="IMNAVAIT"|site_name=="JAKOBSHORN"| site_name=="KANGER"|
site_name=="LATNJA"| site_name=="NIWOT"| site_name=="TOOLIK"| site_name=="VALBERCLA"| site_name=="WHITEMTNS")
sites<-dplyr::select(sites, site_name, lat,long, subsite)%>%group_by(site_name)%>%
  mutate(subsite_no=n_distinct(subsite))%>%mutate(lat=mean(lat, na.rm=T), long=mean(long, na.rm=T))%>%
  dplyr::select(-subsite)%>%distinct(.)%>%arrange(site_name, long, lat, subsite_no)

sites_map <- SpatialPointsDataFrame(coords = sites[,c(3,2)], data=sites[,c(1,4)], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Create world map using maps
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(35, 80), xlim=c(-180, 30))
# Transform map to SpatialLines object
wrld_sp <- map2SpatialLines(wrld)
# Set CRS to lat lon (conversion does not carry it across)
proj4string(wrld_sp) <- CRS("+proj=longlat")

# Project map and sites into Canadian Centred polar projection.
#laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=aea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
#laea_sites <- spTransform(sites_map,CRS("+proj=laea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") )

# Plot map
plot(wrld_sp, col = 'black', lwd = 0.3, cex= 1)
# Plot Sites onto map
plot(sites_map, cex = 1, pch=21,
     col = "black",  
     #red and blue pallettes for siteT-not in order 
     #bg=c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026", 
     #     "#F7FBFF" ,"#DEEBF7" ,"#C6DBEF" ,"#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C" ,"#08306B"),
     bg=c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D" ,"#666666","#4575B4"),
          #"#324D5CFF", 
          # "#46B29DFF", 
          # "#C2A33EFF",
          # "#E37B40FF"),
     add = T, color=sites_map$siteT)
# Add labels
text(sites_map, c("ADVENTDALEN", "ALEXFIORD",
                   "ATQASUK", "UTQUIAVIK",
                   "DARING", "ENDALEN",
                   "FAROE", "FINSE" ,"GAVIAPASS", 
        "HEALY", "IMNAVAIT",  "JAKOBSHORN", "KANGER",
        "LATNJA", "NIWOT", "TOOLIK" , "VALBERCLA", "WHITEMTNS"),
                    font = 2, cex = 0.75,
     col = c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D" ,"#666666","#4575B4"),
      #c("#324D5CFF", 
            # "#46B29DFF", 
            # "#C2A33EFF",
            # "#E37B40FF",),
     pos = c(1,4,4,2,1,2,1,4,4,4,4,3,4,4,1,2,2,1), offset = 0.5, halo = F)


## From Jakob: Suggestions for maps added 15 September 2020 ----
# Here a simple example using simple feautres (sf) and boundaries from
# natural earth. Simple features easily integrate into ggplot.
# I have left some comments, but it's hopefully pretty self explanatory.

# Additional dependencies
library(sf) # Simple feature geometries

library(rnaturalearth) # National boundaries from https://www.naturalearthdata.com/
library(rnaturalearthdata) 

library(cowplot) # Added simplicity for ggplot themes and export
library(ggrepel) # Avoid overlap when placing ggplot labels

# Specify a lamberd conformal conic projection, which works nicely
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

# Generate politcial boundaries as sf and transform from lat long to lambert
world <- ne_countries(scale = "medium", returnclass = "sf")
world_lambert <- st_transform(world, crs = st_crs(lambert_conformic))

# Generate sf object for site point sources and transform
sites_sf <- st_as_sf(sites, coords = c("long", "lat"), crs = 4326)
sites_sf <- st_transform(sites_sf , crs = st_crs(lambert_conformic))
#fix names
sites_sf<-mutate(sites_sf, site_name2=case_when(site_name=="ADVENTDALEN"~"ADVENTDALEN", 
                                               site_name=="ALEXFIORD"~"ALEXANDRA FIORD", 
                                               site_name=="ATQASUK"~"ATQASUK",
                                               site_name=="BARROW"~"UTQIAGVIK",
                                               site_name=="DARING"~"DARING LAKE",
                                               site_name=="ENDALEN"~"ENDALEN",
                                               site_name=="FAROE"~"FAROE ISLANDS",
                                               site_name=="FINSE"~"FINSE",
                                               site_name=="GAVIAPASS"~"GAVIA PASS",
                                               site_name=="HEALY"~"HEALY",
                                               site_name=="TOOLIK"~"TOOLIK LAKE",
                                               site_name=="IMNAVAIT"~"IMNAVAIT CREEK",
                                               site_name=="JAKOBSHORN"~"JAKOBSHORN",
                                               site_name=="KANGER"~"KANGERLUSSUAQ",
                                               site_name=="KANGER"~"KANGERLUSSUAQ",
                                               site_name=="LATNJA"~"LATNJAJAUARE",
                                               site_name=="NIWOT"~"NIWOT RIDGE",
                                               site_name=="VALBERCLA"~"VAL BERCLA",
                                               site_name=="WHITEMTNS"~"WHITE MTNS"))
                                               
                  

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
          fill = "#d3d3d366", 
          size = 0.25) + 
  # geom_sf to plot the sites (it automatically knows they're points)
  geom_sf(data = sites_sf, colour = "red",
          fill = "#ffffff00", shape = 21, size = 2,
          stroke = 2) +  
  # geom_label_repel with st_coordinates as stats to plot labels 
  # that are not overlapping
  geom_label_repel(data = sites_sf, aes(label = site_name2,
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
        panel.grid.minor = element_line(colour = "white"), # Just in case you would like to modify these
        panel.grid.major = element_line(colour = "white"), # Just in case you would like to modify these
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = element_rect(fill = "white") # Just in case you would like to modify these
        ) 

# A quick plot on the native graphics device
otc_site_map

# Export as file using save_plot and the aspect ratio deterimined from the
# min and max coordinates.
save_plot(otc_site_map, filename = "plots/maps/OTC_map_ver2.png",
          base_aspect_ratio = 
            (map_extent$xmax - map_extent$xmin) /
            (map_extent$ymax - map_extent$ymin),
          #base_height = 4 # Just in case you want to modify the base height
          ) 

