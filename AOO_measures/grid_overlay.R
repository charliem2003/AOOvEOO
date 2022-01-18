####################################################################################################
##########                                 Grid overlay                                   ##########
####################################################################################################
### Script for calculating Area of Occupancy (AOO) using Cartographic method by conglomerates (CMC)
### Hernández & Navarro (2007) Biod & Cons, 16:2457-2470
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
aoo <- data.frame(species = species,
                  aoo1 = NA,
                  aoo2 = NA,
                  aoo4 = NA,
                  aoo8 = NA,
                  aoo16 = NA,
                  aoo32 = NA,
                  aoo64 = NA,
                  aoo128 = NA)
aoo$species <- gsub(" ", "_", aoo$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  
  if(nrow(coordsSp) > 2) {
    coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    
    ### Loop through generating grids at the different cell widths (km)
    cellWidths <- c(1, 2, 4, 8, 16, 32, 64, 128)
    
    for(cellWidth in cellWidths) {
      ### if only a single point then AOO is area of single cell
      if(length(coordsSp) == 1) {
        aoo[sp, paste0("aoo", cellWidth)] <- cellWidth ^ 2
      }
      
      ### if more than one point then rasterize and sum area of occupied cells
      if(length(coordsSp)) {
        r <- raster(xmn =   floor(extent(coordsSp)[1] / cellWidth) * cellWidth,
                    xmx = ceiling(extent(coordsSp)[2] / cellWidth) * cellWidth,
                    ymn =   floor(extent(coordsSp)[3] / cellWidth) * cellWidth,
                    ymx = ceiling(extent(coordsSp)[4] / cellWidth) * cellWidth,
                    resolution = cellWidth)
        aooRaster <- rasterize(coordsSp, r)
        vals <- values(aooRaster)
        vals[vals > 1] <- 1
        vals[is.na(vals)] <- 0
        aoo[sp, paste0("aoo", cellWidth)] <- sum(vals) * (cellWidth ^ 2)
      }
    }
  }
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(aoo, "Results/grid_overlay.csv", quote = FALSE, row.names = FALSE)
