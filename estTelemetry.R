#' @export
#' 
#' @name estTelemetry
#' @title Estimates target location from ground-based telemetry.
#' @description 
#'
#' @param fLoc Location of file which the data are to be read from (character string).
#'
#' @param fName Name of either .csv or .txt file containing the data (character string).
#'
#' @param colTrackID Name of column containing subject ID (character string).
#'
#' @param colDate Name of column containing date (character string).
#'
#' @param colAzimth Name of column containing compass azimuth/bearing (character string).
#'
#' @param colX Name of column containing x location either northing or latitude (character string).
#' 
#' @param colY Name of column containing y location either easting or longitude (character string).
#'
#' @param confLevel Confidence region either 0.90 or 0.95 (numeric; default = 0.9).
#'
#' @param method Estimation method: 'MLE' for the Von Misis standard (non-robust) MLE, 
#' 'Andrews' for the Andrews robust procedure, or 'Huber' for the Huber robust procedure. 
#' (character; default = 'Huber').
#' 
#' @param k (optional) Bearing standard error estimated from the data expressed as the von 
#' Mises Consentration parameter (integer; default = NA).
#'
#' @param tc (optional) Tuning Constant. interger; default = 1).
#'
#' @details
#'
#' @note See Lenth, R.V. 1981. On finding the source of a signal.
# Technometrics 23:149-154.
#'
#' @return list of data frames and shapefiles
#'
#' @author Lacey Jeroue
#'
#' @examples
#' \dontrun{
#' ex = ReadData(projectName = ListProjects()$ProjectName[1],
#'               surveyType = NULL,
#'               path = '//chy-file-srv/GeneralFiles/All_Wind_Projects/RCodeForSQL/Testing/QAQC/')
#' }


estTelemetry <- structure( 
  function 
  (fLoc,
   fName,
   colTrackID,
   colDate,
   colAzimth,
   colZone,
   colX,
   colY,
   confLevel = .90,
   method = "Huber",
   k = NA,
   tc = 1
  ) {

    
    errorDuplicate <- NULL
    error1Station <- NULL
    
    # Check for and install packages if needed, and add to library
    for(packages in c("sp","geosphere","plyr","reshape2","rgdal")){
      if(!packages %in% rownames(installed.packages())){
        install.packages(packages)
      }
      require(packages)
    } 
  

    # developer tmp***
    confLevel <- .90
    method <-  "MLE"
    k = NA
    tc = 1
    options(stringsAsFactors = F)
    # fName <- "rawBat.csv"
    # fLoc <- "C:/LACEY/CODE/R/Triangulate/"
    fLoc <- "C:/LACEY/CODE/R/Triangulate/sigloc/"
    fName <- "bear.txt"
    colTrackID = "GID"
    colDate = "Date"
    colBearing ="Azimuth"
    colZone = "Zone"
    colX ="Easting"
    colY = "Northing"
    
    
    # Read in data
    if(grepl(".txt", fName)){
      raw <- read.table(paste0(fLoc,fName), header=TRUE)
      raw$Zone <- 15  # temp
    }
    if(grepl(".csv", fName)){
      raw <- read.csv(paste0(fLoc,fName), header=TRUE)
    }
    
    
    # ERROR HANDELING
    
    if(missing("colTrackID")) stop('argument "colTrackID" missing with no default')
    if(missing("colDate")) stop('argument "colDate" missing with no default')
    if(missing("colBearing")) stop('argument "colBearing" missing with no default') 
    if(missing("colZone")) stop('argument "colZone" missing with no default') 
    if(missing("colX")) stop('argument "colX" missing with no default') 
    if(missing("colY")) stop('argument "colY" missing with no default') 
    
    
    columnNames <- c(colTrackID,
                     colDate,
                     colBearing,
                     colZone,
                     colX,
                     colY)
    
    
    if (!all(columnNames %in% names(raw))) {
      missingCol <- columnNames[! columnNames %in% names(raw)]
      stop(paste0("The column names specified do not exist in your dataset and need to be added: ", 
                  missingCol))
    }
    
    stringArgs <- c(fLoc, fName, columnNames)
    if(!all(is.character(stringArgs))){
      stop(paste0("These arguments must be characters: ", fLoc, fName, columnNames))
    }
    
    
    if ( length(raw[duplicated(raw[,columnNames]), ][,1]) > 0){
      warning("The record(s) above are duplicated and will be removed from analysis. Further qaqc is likely needed. See error report.")
      print(raw[duplicated(raw[,columnNames]), ])
      errorDuplicate <- raw[duplicated(raw[,columnNames]), ]
      errorDuplicate$ErrorMessage <- "Duplicated record"
      raw <- raw[-which(duplicated(raw[,columnNames])), ]
    }
    
    if (confLevel != .9 && confLevel != .95) {
      stop("Confidence level must be either .90 or .95")
    }
    
    method = tolower(method)
    if(!method %in% c("mle", "andrews", "huber")){
      stop("Input for method must be either \'MLE\', \'Andrews\', or \'Huber\'")
    }
    
    
    # ADD PARAMETERS
    
    # Add Estimate ID
    uniqueTrack <- aggregate(formula(paste0(colZone, " ~ ", colTrackID, "+", colDate)),
                             data=raw[,c(colTrackID, colDate, colZone)], length)
    names(uniqueTrack)[3] <- "Signals"
    uniqueTrack$EstimateID <- seq(1,nrow(uniqueTrack))
    raw <- merge(raw, uniqueTrack[,c(colTrackID,colDate,"EstimateID")],
                 by.x = c(colTrackID,colDate), by.y = c(colTrackID,colDate))
    

    # Add Station ID
    raw <- raw[order(raw$EstimateID), ]
    StationID <- c()
    for(i in 1:nrow(uniqueTrack)){
      StationIDx <- seq(1, uniqueTrack$Signals[i])
      StationID <- c(StationID, StationIDx)
    }
    raw$StationID <- StationID
    raw$Station <- paste0("F",raw$EstimateID, "S", raw$StationID) # TODO this is likely not necessary; only used in bearing-intersect calculations
    
    
    # CREATE SHAPEFILE FOR TRANSMISSION LOCATIONS
    
    # Designate a coordinate system
    zone <- as.numeric(gsub("[[:alpha:]]","", raw[1,colZone]))
    crsProj <- paste("+proj=utm +zone=", zone, " +ellps=GRS80 +datum=NAD83 +units=m +no_defs", sep="")
    utm <- sp::CRS(crsProj)

      
    # Create shapefile for Transmission stations
    coords <- cbind(raw[,colX], raw[,colY])
    raw.df <- sp::SpatialPointsDataFrame(coords, raw, proj4string = utm, match.ID = T)
    suppressWarnings(rgdal::writeOGR(raw.df, dsn= paste0(fLoc,"/Shapefiles") ,layer="TransmitterLocations",
             driver="ESRI Shapefile"))
   

    # Remove any targets which have only one azimuth/bearing
    # (handeling here so that this station is included in station shapefile created above)
    if (any(uniqueTrack$Signals == 1)){
      warning("The record(s) above have only 1 bearing records and will be removed from analysis. Further qaqc is likely needed. See error report.")
      print(uniqueTrack[uniqueTrack$Signals == 1, ])
      
      # Save targest with only one azimuth for error report
      error1Station <- raw[raw$EstimateID %in% uniqueTrack$EstimateID[uniqueTrack$Signals == 1], ]
      error1Station$ErrorMessage <- "Only one bearing"
      
      raw <- raw[!raw$EstimateID %in% uniqueTrack$EstimateID[uniqueTrack$Signals == 1], ]
    }
    
    
    # Write all Error report
    errors <- rbind(errorDuplicate[, c(columnNames,"ErrorMessage")], 
                    error1Station[, c(columnNames, "ErrorMessage")])
    write.csv(errors, paste0(fLoc,"/ErrorReport.csv"), row.names = FALSE)

    
    
    #   print(sprintf("BatID: %s on %s - %d records removed from analysis due to duplicate locations", 
    #                 dupLocs$BatID[dupLocs$x > 1], dupLocs$Date[dupLocs$x > 1], dupLocs$x[dupLocs$x > 1]))
    # }
    
    
    ## ESTIMATE TARGET LOCATION ####
    targetLocation <- estLocation(raw, method=method, kappa=k, tuningConstant = tc)
    
    
    estimatedPoints <- data.frame(EstimateID = do.call("rbind", lapply(targetLocation, "[[", 1)),
                             do.call("rbind", lapply(targetLocation, "[[", 2)))
    estimatedPoints <- merge(estimatedPoints, unique(raw[,c(colTrackID, colDate, "EstimateID")]))
    
    targetLocation
      
      
    # Remove any attempts that did not result in an estimate
    # TODO throw out the one x,y,azimuth that is permitting convergence so that we can still get
    # an estimate; single out that record for a qaqc report.
    estimatedPoints <- estimatedPoints[!is.na(estimatedPoints$xHat), ]
    
    
    # Create Shapefile of the estimated triangulated Points
    coords <- cbind(estimatedPoints$xHat, estimatedPoints$yHat)
    spdf <- sp::SpatialPointsDataFrame(coords, estimatedPoints, proj4string = utm, match.ID = T)
    rgdal::writeOGR(spdf, dsn= paste0(fLoc,"/Shapefiles") ,layer="EstimatedLocations",
                    driver="ESRI Shapefile")
    
    
    ## CONSTRUCT CONFIDENCE REGION AND CREATE ELLIPSE SHAPEFILE
    confRegion <- confEllipse(estimatedLocations = targetLocation, 
                              confLevel=confLevel)
    confRegion$area
    confRegion$ellipsePoints
    
    
    # Create Ellipse Shapefile
    SLines <- sp::SpatialLines(confRegion$ellipsePoints[!sapply(confRegion$ellipsePoints, is.null)], 
                           proj4string = utm)
    list_of_Lines <- slot(SLines, "lines")
    elip <- sp::SpatialPolygons(lapply(list_of_Lines, function(x) {
      sp::Polygons(list(sp::Polygon(slot(slot(x, "Lines")[[1]],
                                 "coords"))), ID = slot(x, "ID"))
    }), proj4string = utm)
    
    spdf <- sp::SpatialPolygonsDataFrame(elip, confRegion$area, match.ID = F)
    suppressWarnings(rgdal::writeOGR(spdf, dsn= paste0(fLoc, "/Shapefiles"),
                              layer="Confidence_Region",
                              driver="ESRI Shapefile"))
    
    # return(list(DF, hold))
  # }
  
    
    ## PREPARE DATA AND CREATE BEARING SHAPEFILES
    
    # Get all possible station combinations for bearing intersections
    tmp3 <- aggregate(raw[,2], by=list(EstimateID=raw[,c("EstimateID")]), length)
    combos <- NULL
    for(i in 1:nrow(tmp3)){
      combo1 <- utils::combn(raw$Station[raw$EstimateID==tmp3[i,1]], 2)
      meltdf <- reshape2::melt(combo1, id.vars=c())
      meltdf$StationID <- rep(c("Station_1", "Station_2"), length(meltdf[,1])/2)
      combo2 <- reshape2::dcast(meltdf, Var2 ~ StationID, value.var="value")
      combo3 <- combo2[,-c(1)]
      combos <- rbind(combos, combo2)
    }
    
    # Create Data Frame for holding Intersections and Bearings
    st1 <- raw[,c("EstimateID","Station", colBearing,colX,colY)]
    st2 <- raw[,c("EstimateID","Station", colBearing,colX,colY)]
    names(st1) <- paste0(names(st1),"_1")
    names(st2) <- paste0(names(st2),"_2")
    int <- plyr::join(plyr::join(st1, combos[,c("Station_1", "Station_2")]), st2)
    int <- int[!is.na(int$EstimateID_2), ]
    int$EstimateID_2 <- NULL
    names(int) <- gsub("EstimateID_1", "EstimateID", names(int))
    
    
    int4 <- int
    
    
    
    # Solve for Intersections with Linear Algebra and Add to df
    int4$slope_1 <- (tan((int4$Azimuth_1*pi)/180))^(-1)
    int4$b_1 <- int4$Northing_1 - (int4$Easting_1*int4$slope_1)
    
    int4$slope_2 <- (tan((int4$Azimuth_2*pi)/180))^(-1)
    int4$b_2 <- int4$Northing_2 - (int4$Easting_2*int4$slope_2)
    
    int5 <- int4
    
    Easting_int <- NULL
    Northing_int <- NULL
    data.int5 <- NULL

    for(i in 1:length(int5[,1])){
      tmp6 <- cbind(int5[i,10], -1, -(int5[i,11]), int5[i,12], -1, -(int5[i,13]))
      intersect <- NULL
      DF2 <- NULL
      if( !is.null(tryCatch({solve(rbind(tmp6[,1:2],tmp6[,4:5])) %*% rbind(tmp6[,3], tmp6[,6])
      }, error=function(e){})) ){
        DF2 <- int5[i,1:9]
        intersect <- tryCatch({
          solve(rbind(tmp6[,1:2],tmp6[,4:5])) %*% rbind(tmp6[,3], tmp6[,6])
        }, error=function(e){})
      }
      data.int5 <- rbind(data.int5, DF2)
      Easting_int <- rbind(Easting_int, intersect[1])
      Northing_int <- rbind(Northing_int, intersect[2])
    }
    
    int5.1 <- cbind(data.int5, Easting_int, Northing_int)
    int5 <- int5.1
    # head(int5)
    
    # Reshape dataframe to fields of interest for each POLYLINE
    tmp7 <- cbind(int5[,1:6], int5[,10:11])
    tmp8 <- cbind(int5[,1],int5[,6:9], int5[,2], int5[,10:11])
    names(tmp7) <- c("EstimateID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2", "Easting2", "Northing2")
    names(tmp8) <- c("EstimateID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2", "Easting2", "Northing2")
    int6 <- arrange(rbind(tmp7,tmp8), EstimateID)
    head(int6)
    
    # Reshape dataframe to fields of interest for each POLYLINE
    tmp4.1 <- cbind(int4[,1:6])
    tmp4.2 <- cbind(int4[,1],int4[,6:9], int4[,2])
    names(tmp4.1) <- c("EstimateID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2")
    names(tmp4.2) <- c("EstimateID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2")
    int4.3 <- arrange(rbind(tmp4.1,tmp4.2), EstimateID)
    int4.4 <- int4.3[,c(1:5)]
    int4.5 <- int4.4[!duplicated(int4.4),]
    # head(int4.5)
    
    # Determine which intersections are the wrong direction
    int6$reject <- ifelse(int6$Azimuth <= 90 & 
                            int6$Easting1-int6$Easting2 < 0 & 
                            int6$Northing1 - int6$Northing2 < 0, "accept",
                          ifelse(int6$Azimuth > 90 & 
                                   int6$Azimuth <= 180 & 
                                   int6$Easting1-int6$Easting2 < 0 & 
                                   int6$Northing1 - int6$Northing2 > 0, "accept",
                                 ifelse(int6$Azimuth > 180 & 
                                          int6$Azimuth <= 270 & 
                                          int6$Easting1-int6$Easting2 > 0 & 
                                          int6$Northing1 - int6$Northing2 > 0, "accept",
                                        ifelse(int6$Azimuth > 270 & 
                                                 int6$Easting1-int6$Easting2 > 0 & 
                                                 int6$Northing1 - int6$Northing2 < 0, "accept",
                                               "reject"))))
    
    ##  Create Bearings Shapefile
    int7 <- int6[int6$reject=="accept",-c(9)]
    int7 <- int7[!duplicated(int7),]
    int7.1 <- int7
    names(int7) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2", "Easting1", "Northing1")
    
    # Use dest point for those that do not intersect"
    # head(int4)
    # head(int6)
    bad <- int4.5[!int4.5$StID1 %in% int7$StID1,]
    
    if(length(bad[,1]) > 0){
      bad1 <- bad[!duplicated(bad$StID1),]
      bad1$StID2 <- NA
      bad2 <- bad1
      longlat = CRS("+proj=longlat + ellps=WGS84")
      coordinates(bad2) <- c("Easting1", "Northing1")
      proj4string(bad2) <- utm
      longlat_1 <- spTransform(bad2, longlat)
      ll_1 <- data.frame(longlat_1@coords)
      names(ll_1) <- c("Long1", "Lat1")
      bad3 <- cbind(bad1[c("FixID", "StID1", "Azimuth", "Easting1", "Northing1","StID2")], ll_1)
      badcoords <- cbind(bad3[,7:8])
      coordinates(bad3) <- destPoint(badcoords, bad3$Azimuth, 600)
      proj4string(bad3) <- longlat
      utm_bad <- spTransform(bad3, utm)
      badUTMends <- data.frame(utm_bad@coords)
      names(badUTMends) <- c("Easting1", "Northing1")
      bad5 <- cbind(bad1[c("FixID", "StID1", "Azimuth", "Easting1", "Northing1", "StID2")], badUTMends)
      
      allBear <- arrange(rbind(int7, bad5), FixID) # combine all bearings
      allBear2 <- allBear
      names(allBear) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1", "StID2", 
                          "Easting2", "Northing2")
      names(allBear2) <- c("FixID", "StID1", "Azimuth", "Easting", "Northing", "StID2", 
                           "Easting", "Northing")
      
      hold <- list()
      for(i in 1:length(allBear2[,1])){
        hold[[i]] <- Lines(list(Line(data.frame(rbind(allBear2[i,4:5], allBear2[i,7:8])))), ID=i)
      }
      SL1 <- SpatialLines(hold, proj4string = utm)
      SLDF1 <- SpatialLinesDataFrame(SL1, allBear, match.ID = T)
      writeOGR(SLDF1, dsn= paste0(fLoc, "Shapefiles"), layer="Bearings1",driver="ESRI Shapefile")
      
    } else {
      
      hold <- list()
      for(i in 1:length(int7[,1])){
        hold[[i]] <- Lines(list(Line(data.frame(rbind(int7[i,4:5], int7[i,7:8])))), ID=i)
      }
      SL1 <- SpatialLines(hold, proj4string = utm)
      SLDF1 <- SpatialLinesDataFrame(SL1, int7.1, match.ID = F)
      writeOGR(SLDF1, dsn= paste0(fLoc, "Shapefiles") ,layer="Bearings2",driver="ESRI Shapefile")
    }

})
