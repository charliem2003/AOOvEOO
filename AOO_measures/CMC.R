  ####################################################################################################
  ##########                                      CMC                                       ##########
  ####################################################################################################
  ### Script for calculating Area of Occupancy (AOO) using Cartographic method by conglomerates (CMC)
  ### Hernández & Navarro (2007) Biod & Cons, 16:2457-2470
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
  cmc <- data.frame(species = species, cmc = NA)
  cmc$species <- gsub(" ", "_", cmc$species)
  
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
      
      ### generate minimum spanning tree
      distances <- dist(coordsSp@coords)
      tree <- mst(distances)
  
      ### calculate mean distances of mst
      distances <- pointDistance(coordsSp)
      distances[tree == 0] <- NA
      mean.mst <- mean(distances, na.rm = TRUE)
      
      ### buffer of mean mst around each point
      buffers <- gBuffer(coordsSp, width = mean.mst)
      buffers <- disaggregate(buffers)
      
      ### if you want to plot the MST
      # plot(coordsSp, axes = T, main = species[sp])
      # mstlines(tree, coordsSp@coords)
      # plot(buffers, add = TRUE)
      
      ### separate conglomerates and determine which ponits are in them
      conglomerates <- lapply(1:length(buffers), function(i) buffers[i, ])
      conglomerates <- lapply(conglomerates, function(x) over(coordsSp, x))
      
      aooCons <- vector()
      for(con in 1:length(conglomerates)) {
        ### extract just the coordinates for conglomerate
        conglomerate <- conglomerates[[con]]
        coordsCon <- coordsSp[!is.na(conglomerate)]
        
        ### if there is only a single unique point
        if(length(coordsCon) == 1) {
          aooCons[con] <- 2
        }
        
        if(length(coordsCon) > 1) {
          ### cell width = largest distance between points / 10
          dists <- pointDistance(coordsCon, lonlat = FALSE)
          cellWidth <- max(dists) / 10
          
          ### create species grid
          r <- raster(xmn =   floor(extent(coordsCon)[1] / cellWidth) * cellWidth,
                      xmx = ceiling(extent(coordsCon)[2] / cellWidth) * cellWidth,
                      ymn =   floor(extent(coordsCon)[3] / cellWidth) * cellWidth,
                      ymx = ceiling(extent(coordsCon)[4] / cellWidth) * cellWidth,
                      resolution = cellWidth)
          aooRaster <- rasterize(coordsCon, r)
          vals <- values(aooRaster)
          vals[vals > 1] <- 1
          vals[is.na(vals)] <- 0
          aooCons[con]<- sum(vals) * (cellWidth ^ 2)
        }
      }
      cmc[sp, "cmc"] <- sum(aooCons)
    }
    removeTmpFiles(h = 0.25)
    setTxtProgressBar(pb, value = sp)
  }
  
  ### if you want to save the results
  write.csv(cmc, "Results/cmc.csv", quote = FALSE, row.names = FALSE)
