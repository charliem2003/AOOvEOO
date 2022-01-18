####################################################################################################
##########                                Variable grain                                  ##########
####################################################################################################
### Script for calculating Area of Occupancy (AOO) using the variable grain method of
### Willis et al. 2003. Biod. and Cons., 12:1537–1552.
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charlie.marsh@mailbox.org)
####################################################################################################

### libraries required
library(shapefiles)
library(raster)

### read in point data (shapefile) for all species
fernsShp <- read.shapefile("Data/Ferns/Ferns")
coords <- fernsShp$dbf$dbf

### species list
species <- unique(coords$BINOMIAL)

### data frame for storing results at different cell widths
aooVar <- data.frame(species = species,
                     aoo.var = NA,
                     varWidth = NA)
aooVar$species <- gsub(" ", "_", aooVar$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  
  if(nrow(coordsSp) > 2) {
    coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    
    ### if only a single point then AOO is 2 km2
    if(length(coordsSp) == 1) {
      aooVar$aoo.var[sp] <- 2
    }
    
    if(length(coordsSp) > 1) {
      ### calculate 1/10th largest pairwise interpoint distances
      dists <- pointDistance(coordsSp, lonlat = FALSE)
      cellWidth <- max(dists) / 10
      aooVar$varWidth[sp] <- cellWidth
      
      ### rasterize points at this cell width
      r <- raster(xmn = floor(extent(coordsSp)[1] / cellWidth) * cellWidth,
                  xmx = ceiling(extent(coordsSp)[2] / cellWidth) * cellWidth,
                  ymn = floor(extent(coordsSp)[3] / cellWidth) * cellWidth,
                  ymx = ceiling(extent(coordsSp)[4] / cellWidth) * cellWidth,
                  resolution = cellWidth)
      aooRaster <- rasterize(coordsSp, r)
      vals <- values(aooRaster)
      vals[vals > 1] <- 1
      vals[is.na(vals)] <- 0
      aooVar$aoo.var[sp] <- sum(vals) * (cellWidth ^ 2)
    }
  }
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(aooVar, "Results/variable_grain.csv", quote = FALSE, row.names = FALSE)
