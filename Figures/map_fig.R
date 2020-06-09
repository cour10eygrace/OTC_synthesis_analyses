
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
                    font = 2, cex = 0.5,
     col = c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D" ,"#666666","#4575B4"),
      #c("#324D5CFF", 
            # "#46B29DFF", 
            # "#C2A33EFF",
            # "#E37B40FF",),
     pos = c(1,4,4,2,1,2,1,4,4,4,4,3,4,4,1,2,2,1), offset = 0.5, halo = F)

