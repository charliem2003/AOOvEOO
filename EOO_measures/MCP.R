####################################################################################################
##########                          Minimum convex polygon (MCP)                          ##########
####################################################################################################
### Script for calculating Extent of Occurrence (EOO)) using the minimum convex polygons (MCP)
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charliem2003@gmail.com)
####################################################################################################

### libraries required
library(shapefiles)
library(tripack)
library(rgeos)

### read in point data (shapefile) for all species
fernsShp <- read.shapefile("Data/Ferns/Ferns")
coords <- fernsShp$dbf$dbf

### remove coords outside of americas
coords <- coords[coords$LONGITUDE > -126, ]
coords <- coords[coords$LONGITUDE <  -34, ]
coords <- coords[coords$LATITUDE  >  -45, ]
coords <- coords[coords$LATITUDE  <   52, ]

### species list
species <- unique(coords$BINOMIAL)

### data frame for storing results using different alpha values
mcp <- data.frame(species = species,
                  mcp = NA)
mcp$species <- gsub(" ", "_", mcp$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  
  if(nrow(coordsSp) > 2) {
    ### generate mcp
    hull <- chull(coordsSp[, 1], coordsSp[, 2])
    coordsHull <- SpatialPoints(coordsSp[hull, ],
                                proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsHull <- spTransform(coordsHull, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    
    ### convert to SpatialPolygons object and calculate area
    hull <- SpatialPolygons(list(Polygons(list(Polygon(coordsHull)), ID=1)))
    alpha[sp, "mcp"] <- sum(area(hull))
  }
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(mcp, "All_mcp.csv", quote = F, row.names = F)
