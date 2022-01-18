####################################################################################################
##########                                   Alpha hull                                   ##########
####################################################################################################
### Script for calculating Extent of Occurrence (EOO)) using alpha-hulls
### Burgman & Fox (2003) Animal Cons., 6:19–28
###
### Part of the methods for the manuscript: CJ Marsh, MM Syfert, E Aletrari, Y Gavish, WE Kunin & 
### N Brummitt. 'The effect of sampling effort and methodology on range size estimates of
### poorly-recorded species for IUCN Red List assessments'. In prep.
###
### Charlie Marsh (charlie.marsh@mailbox.org)
####################################################################################################

### libraries required
library(shapefiles)
library(tripack)
library(rgeos)
library(raster)

### read in point data (shapefile) for all species
fernsShp <- read.shapefile("Data/Ferns/Ferns")
coords <- fernsShp$dbf$dbf

### species list
species <- unique(coords$BINOMIAL)

### data frame for storing results using different alpha values
alpha <- data.frame(species = species,
                    alpha_1 = NA, alpha_2 = NA, alpha_3 = NA,
                    alpha_4 = NA, alpha_5 = NA, alpha_6 = NA)
alpha$species <- gsub(" ", "_", alpha$species)

### loop through species
pb <- txtProgressBar(max = length(species), style = 3)
for(sp in 1:length(species)) {
  ### subset points for selected species and remove duplicates
  coordsSp <- coords[coords$BINOMIAL == as.character(species[sp]), c("LONGITUDE", "LATITUDE")]
  coordsSp <- coordsSp[!duplicated(coordsSp), ]
  
  if(nrow(coordsSp) > 2) {
    coordsSp <- SpatialPoints(coordsSp, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    coordsSp <- spTransform(coordsSp, CRS = CRS("+proj=cea +datum=WGS84 +units=km"))
    coordsSp <- coordsSp@coords
    
    ### generate delaunay triangle
    x <- tri.mesh(coordsSp[, "LONGITUDE"], coordsSp[, "LATITUDE"])
    
    ### extract connections (code extractd from plot.tri function)
    tnabor <- integer(x$tlnew)
    nnabs  <- integer(x$n)
    nptr   <- integer(x$n)
    nptr1  <- integer(x$n)
    nbnos  <- integer(x$n)
    ans    <- .Fortran("troutq", as.integer(x$nc), as.integer(x$lc),
                       as.integer(x$n), as.double(x$x), as.double(x$y), as.integer(x$tlist),
                       as.integer(x$tlptr), as.integer(x$tlend), as.integer(6),
                       nnabs = as.integer(nnabs), nptr = as.integer(nptr), nptr1 = as.integer(nptr1),
                       tnabor = as.integer(tnabor), nbnos = as.integer(nbnos),
                       na = as.integer(0), nb = as.integer(0), nt = as.integer(0),
                       PACKAGE = "tripack")
    
    ### get connection lengths
    connections <- data.frame(coord1_id = NA, coord1_x = NA, coord1_y = NA,
                              coord2_id = NA, coord2_x = NA, coord2_y = NA)
    for (i in 1:x$n) {
      inb <- ans$tnabor[ans$nptr[i]:ans$nptr1[i]]
      coords.i <- data.frame(coord1_id = NA, coord1_x = NA, coord1_y = NA,
                             coord2_id = NA, coord2_x = NA, coord2_y = NA)
      for (j in 1:length(inb)) {
        coords.i[j, ] <- c(i, x$x[i], x$y[i], inb[j], x$x[inb[j]], x$y[inb[j]])
      }
      connections <- rbind(connections, coords.i)
    }
    
    ### find distances between connections and delete duplicates
    connections$dist <- pointDistance(connections[, 2:3], connections[, 5:6], lonlat = FALSE)
    connections <- connections[!duplicated(rowSums(connections)), ]
    connections <- connections[complete.cases(connections), ]
    
    ### remove connections within alpha value
    mean_dist <- mean(connections$dist, na.rm = TRUE)
    connections_1 <- connections[connections$dist < mean_dist, ]
    connections_2 <- connections[connections$dist < (mean_dist * 2), ]
    connections_3 <- connections[connections$dist < (mean_dist * 3), ]
    connections_4 <- connections[connections$dist < (mean_dist * 4), ]
    connections_5 <- connections[connections$dist < (mean_dist * 5), ]
    connections_6 <- connections[connections$dist < (mean_dist * 6), ]
    
    ### if you want t plot the connection
    # plot(x$x, x$y)
    # for (i in 1:x$n) {
    #   inb <- ans$tnabor[ans$nptr[i]:ans$nptr1[i]]
    #   for (j in inb) lines(c(x$x[i], x$x[j]), c(x$y[i], x$y[j]))
    # }
    # apply(connections_3, 1,
    #       function(x) {lines(x[c(2, 5)], x[c(3, 6)], col = "blue")})
    # apply(connections_2, 1,
    #       function(x) {lines(x[c(2, 5)], x[c(3, 6)], col = "orange")})
    # apply(connections_1, 1,
    #       function(x) {lines(x[c(2, 5)], x[c(3, 6)], col = "red")})
    
    if(!is.null(connections_1)) {
      connections_1 <- SpatialLines(list(
        Lines(apply(connections_1, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_1 <- gPolygonize(connections_1)
      if(!is.null(connections_1)) {
        projection(connections_1) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_1"] <- sum(area(connections_1))
      }
    }
    
    if(!is.null(connections_2)) {
      connections_2 <- SpatialLines(list(
        Lines(apply(connections_2, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_2 <- gPolygonize(connections_2)
      if(!is.null(connections_2)) {
        projection(connections_2) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_2"] <- sum(area(connections_2))
      }
    }
    
    if(!is.null(connections_3)) {
      connections_3 <- SpatialLines(list(
        Lines(apply(connections_3, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_3 <- gPolygonize(connections_3)
      if(!is.null(connections_3)) {
        projection(connections_3) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_3"] <- sum(area(connections_3))
      }
    }
    
    if(!is.null(connections_4)) {
      connections_4 <- SpatialLines(list(
        Lines(apply(connections_4, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_4 <- gPolygonize(connections_4)
      if(!is.null(connections_4)) {
        projection(connections_4) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_4"] <- sum(area(connections_4))
      }
    }
    
    if(!is.null(connections_5)) {
      connections_5 <- SpatialLines(list(
        Lines(apply(connections_5, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_5 <- gPolygonize(connections_5)
      if(!is.null(connections_5)) {
        projection(connections_5) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_5"] <- sum(area(connections_5))
      }
    }
    
    if(!is.null(connections_6)) {
      connections_6 <- SpatialLines(list(
        Lines(apply(connections_6, 1,
                    function(x) Line(cbind(c(x[2], x[5]), c(x[3], x[6])))),
              ID = "all")))
      connections_6 <- gPolygonize(connections_6)
      if(!is.null(connections_6)) {
        projection(connections_6) <- crs("+proj=cea +datum=WGS84 +units=km")
        alpha[sp, "alpha_6"] <- sum(area(connections_6))
      }
    }
    setTxtProgressBar(pb, value = sp)
  }
}

### if you want to save the results
write.csv(alpha, "Results/alpha_hull.csv", quote = FALSE, row.names = FALSE)
