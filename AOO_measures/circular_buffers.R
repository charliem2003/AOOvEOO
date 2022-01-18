####################################################################################################
##########                                Circular buffers                                ##########
####################################################################################################
### Script for calculating Area of Occupancy (AOO) using circular buffers
### Breiner & Bergamini (2018) Biod & Cons, 27:2443-2448
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charliem2003@gmail.com)
####################################################################################################

### libraries required
library(shapefiles)
library(raster)
library(tripack)
library(rgeos)
library(fossil)
library(ape)

### read in point data (shapefile) for all species
fernsShp <- read.shapefile("Data/Ferns/Ferns")
coords <- fernsShp$dbf$dbf

### species list
species <- unique(coords$BINOMIAL)

### data frame for storing results
buffs <- data.frame(species = species, buffer = NA)
buffs$species <- gsub(" ", "_", buffs$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  
  if(nrow(coordsSp) > 2) {
    coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    
    ### draw 2km buffers around
    buffers <- gBuffer(coordsSp, width = sqrt(4 / pi))
    buffs[sp, "buffer"] <- area(buffers)
  }
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(buffs, "Results/circular_buffers.csv", quote = FALSE, row.names = FALSE)
