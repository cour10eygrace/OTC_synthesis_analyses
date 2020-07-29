
# Dependencies
library(maps)
library(raster)
library(maptools)
library(dplyr)
library(RColorBrewer)

# organize Site coordinates
subsites<-read.csv("Data/LOOKUPS/subsites.csv") 
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
     bg=c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D" ,"#666666","#4575B4"),
          add = T, color=sites_map$siteT)
# Add labels
text(sites_map, c("ADVENTDALEN", "ALEXANDRA FIORD",
                   "ATQASUK", "UTQIAGVIK",
                   "DARING LAKE", "ENDALEN",
                   "FAROE ISLANDS", "FINSE" ,"GAVIA PASS", 
        "HEALY", "IMNAVAIT CREEK",  "JAKOBSHORN", "KANGERLUSSUAQ",
        "LATNJAJAURE", "NIWOT RIDGE", "TOOLIK LAKE" , "VALBERCLA", "WHITE MOUNTAINS"),
                    font = 2, cex = 0.75,
     col = c("#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D" ,"#666666","#4575B4"),
      pos = c(1,4,4,2,1,2,1,4,4,4,4,3,4,4,1,2,2,1), offset = 0.5, halo = F)

