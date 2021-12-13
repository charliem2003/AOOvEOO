####################################################################################################
##########                             Occupancy downscaling                              ##########
####################################################################################################
### Script for calculating Area of Occupancy (AOO) using occupancy downscaling
### Kunin (1998) Science, 281:1513-1515; Marsh et al. (2019) Div. & Distr. 25:1832-1845
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charliem2003@gmail.com)
####################################################################################################

### libraries required
library(raster)
library(sp)
library(spatstat)
library(rgdal)
library(maptools)
library(shapefiles)
library(downscale)

### read in point data (shapefile) for all species
fernsShp <- read.shapefile("Data/Ferns/Ferns")
coords <- fernsShp$dbf$dbf

### map of coastlines for masking out ocean areas
coastline <- readOGR(dsn = "Data/coastlines", layer = "continent")
coastline <- spTransform(coastline, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))

### remove coords outside of americas
coords <- coords[coords$LONGITUDE > -126, ]
coords <- coords[coords$LONGITUDE <  -34, ]
coords <- coords[coords$LATITUDE  >  -45, ]
coords <- coords[coords$LATITUDE  <   52, ]

### species list
species <- unique(coords$BINOMIAL)

### cell width (km) for creating base maps
cellWidth <- 2

### data frame for storing results
downscale <- data.frame(species = species,
                        no_points = NA,
                        unique_points = NA,
                        ens.down = NA,
                        pl.down = NA)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
  coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
  
  ### check total number of records
  downscale$no_points[sp] <- length(coordsSp)
  
  ### and number of unique points
  if(length(coordsSp) > 1) {
    coordsSp@coords <- as.matrix(data.frame(LONGITUDE = coordsSp@coords[!duplicated(coordsSp@coords), "LONGITUDE"],
                                            LATITUDE  = coordsSp@coords[!duplicated(coordsSp@coords), "LATITUDE"]))
  }
  downscale$unique_points[sp] <- length(coordsSp)
  
  if(length(coordsSp) > 2) {
    
    ###########################################################
    ####### Step 1. Generate atlases at multiple scales #######
    
    ### generate atlas map at 2 x 2 km cell width
    r2km <- raster(xmn =   floor(extent(coordsSp)[1] / cellWidth) * cellWidth,
                   xmx = ceiling(extent(coordsSp)[2] / cellWidth) * cellWidth,
                   ymn =   floor(extent(coordsSp)[3] / cellWidth) * cellWidth,
                   ymx = ceiling(extent(coordsSp)[4] / cellWidth) * cellWidth,
                   resolution = cellWidth,
                   crs = CRS("+proj=cea +datum=WGS84 +units=km"))
    r2km <- rasterize(coordsSp, r2km, background = 0)
    r2km[r2km@data@values > 1] <- 1
    
    ### set ocean areas as NA
    r2km <- mask(r2km, coastline)
    
    all <- data.frame(extent = NA, cellWidth = NA,
                      area1 = NA, area2 = NA, area3 = NA, # store cell areas
                      occ1  = NA, occ2  = NA, occ3  = NA) # store occupancies
    
    ### loop through potential atlas grain sizes
    aggs <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
    for(i in 1:length(aggs)) {
      up <- NULL
      rAgg <- NULL
      
      ### aggregate 2km raster
      try(rAgg <- aggregate(r2km, aggs[i], fun = max),
          silent = TRUE)
      
      ### upgrain 2 more grain sizes using downscale::upgrain
      if(!is.null(rAgg)) {
        try(up <- upgrain(rAgg,
                          scales = 2,
                          method = "All_Sampled",
                          plot = FALSE),
            silent = TRUE)
        
        ### store results
        if(!is.null(up)) { all[i, ] <- c(up$extent.stand,
                                         sqrt(up$occupancy.stand[1, 1]),
                                         up$occupancy.stand$Cell.area,
                                         up$occupancy.stand$Occupancy)
        }
      }
      rm(rAgg)
      rm(up)
    }
    
    ###########################################################
    ############# Step 2. Select best atlas scale #############
    
    ### remove scales at the lower 2 grains if we reached scale of saturation (occ = 1)
    all$occ1[all$occ1 == 1] <- NA
    all$occ2[all$occ2 == 1] <- NA
    
    ### largest atlas grain size that retains 3 grains
    all <- all[max(which(complete.cases(all))), ]
    
    ### save upgrain info
    occ <- data.frame(Cell.area = as.numeric(all[, c("area1", "area2", "area3")]),
                      Occupancy = as.numeric(all[, c("occ1", "occ2", "occ3")]))
    
    ###########################################################
    ################### Step 2. Downscaling ###################
    if(!is.na(sum(occ))) {
      mod <- ensemble.downscale(occupancies = occ,
                                new.areas = 4,    # ie 2 x 2 km
                                extent = all$extent,
                                cell.width = all$cellWidth,
                                models = c("Nachman", "PL", "Logis", "Poisson", "NB"),
                                plot = FALSE, verbose = FALSE)
      
      ### the ensemble estimates
      downscale$ens.down[sp] <- mod$AOO$Means
      
      ### the power law prediction for comparison (the IUCN-suggested model)
      downscale$pl.down[sp] <- mod$AOO$PL
    }
  }
  setTxtProgressBar(pb, value = sp)
}

### if you want to save the results
write.csv(downscale, "All_downscale.csv", quote = F, col.names = T, row.names = F)



