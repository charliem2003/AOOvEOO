####################################################################################################
##########                           LoCoH (a, r and k methods)                           ##########
####################################################################################################
### Script for calculating Area of Occupancy (AOO) using local convex hulls (LoCoH).
### Getz et al. (2007) PLoS ONE 2:e207; Getz & Wilmers (2004) Ecography, 27:489-505
###
### Three methods for constructing the LoCoHs are investigated:
### r-LoCoH - uses a fixed radius around each reference point
### k-LoCoH - kernels constructed from k-1 nearest points
### a-LoCoH - uses all points within radius a where sum of distances sum to <= a
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charliem2003@gmail.com)
####################################################################################################

### libraries required
library(adehabitatHR) # contains functions for LoCoH methods
library(shapefiles)
library(raster)
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

### data frame for storing results
locoh <- data.frame(species = species,
                    k.locoh = NA, k.val = NA,
                    r.locoh = NA, r.val = NA,
                    a.locoh = NA, a.val = NA)
locoh$species <- gsub(" ", "_", locoh$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  # print(paste(sp, nrow(coordsSp), sep = " - "))
  
  if(nrow(coordsSp) > 2) {
    coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    
    ### k-LoCoH - k = sqrt(n) rounded up to nearest integer
    areas <- NULL
    k.val <- ceiling(sqrt(length(coordsSp)))
    locoh[sp, "k.val"] <- k.val
    try(areas <- LoCoH.k(coordsSp, k = k.val, unin = "km", unout = "km2"),
        silent = TRUE)
    if(!is.null(areas)) {
      k.val <- floor(sqrt(length(coordsSp)))
      locoh[sp, "k.val"] <- k.val
      try(areas <- LoCoH.k(coordsSp, k = k.val, unin = "km", unout = "km2"),
          silent = TRUE)
    }
    if(!is.null(areas)) {
      locoh[sp, "k.locoh"] <- areas@data[nrow(areas@data), "area"]
    }
    rm(areas)
    
    ### r-LoCoH - r = half max nearest neighbour distance
    areas <- NULL
    distances <- as.matrix(dist(coordsSp@coords))
    distances[distances == 0] <- NA
    nnb.dist <- apply(distances, 2, min, na.rm = T)
    r.val <- max(nnb.dist) / 2
    locoh[sp, "r.val"] <- r.val
    try(areas <- LoCoH.r(coordsSp, r = r.val, unin = "km", unout = "km2"),
        silent = TRUE)
    if(!is.null(areas)) {
      locoh[sp, "r.locoh"] <- areas@data[nrow(areas@data), "area"]
    }
    rm(areas)
    
    ### a-LoCoH - a = largest interpoint distance
    areas <- NULL
    distances <- dist(coordsSp@coords)
    a.val <- max(distances)
    locoh[sp, "a.val"] <- a.val
    try(areas <- LoCoH.a(coordsSp, a =  a.val, unin = "km", unout = "km2"),
        silent = TRUE)
    if(!is.null(areas)) {
      locoh[sp, "a.locoh"] <- areas@data[nrow(areas@data), "area"]
    }
    rm(areas)
  }
  removeTmpFiles(h = 0.25)
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(locoh, "Results/All_locoh.csv", quote = FALSE, row.names = FALSE)
