Triangulate <- structure( 
  function # Estimates target location and creates triangulation shapefiles
  ###############################################################################
  # File: Trangulate.R
  # Authors: Lacey Jeroue Contact: \email{<ljeroue@west-inc.com}
  # Chris Nations: \email{<cnations@west-inc.com}
  # Created: July 12, 2015
  # Last Edited: 08/12/2015 by Lacey
  ##description<< Estimates target location (potentially with error) using known
  # station locations and bearings and with user defined triangulation methods 
  # based on methods in: Lenth, R.V. 1981.  On finding the source of a signal.
  # Technometrics 23:149-154. Produces associated ArcGIS Shapefiles.
  # Returns:
  # A) A folder in user defined file location including shapefiles: 
  # 1) Triangulation_Stations,
  # 2) Triangulated_pts, 
  # 3) Bearings, and 
  # 4) Confidence_Region.
  # B) User can call a list composed of three elements: 
  # 1) A list contaning the von Mises estimate of location (vm), covariance 
  # matrix (qhat), weights (w), k, n, tuning constant, bearings, outcome, 
  # and method, 
  # 2) a data frame including location estimates and 
  # 3) a data frame including area of confidence region.
  ###############################################################################
  
  (fLoc,
   ### The location of csv file which the data are to be read from (character string).
   fName,
   ### The name of csv file which the data are to be read from (character string).
   # fDest = paste0(fLoc,"Shapefiles"),
   ### (optional) Choose names of the folder to hold shapefiles (character string).
   colTrackID,
   colDate,
   colBearing,
   colZone,
   colEasting,
   colNorthing,
   conf=.90,
   ### (optional) Confidence level. Default is 0.90 but 0.95 is an option (integer).
   k = NA,
   ### (optional) Bearing standard error estimated from the data expressed as the von 
   ### Mises Consentration parameter (integer).
   method = "h",
   ### (optional) Estimation method: 'm' for the Von Misis standard (non-robust) MLE, 
   ### 'a' for the Andrews robust procedure, or 'h' for the Huber robust procedure. 
   ### Default is 'h' (character).
   tc = 1
   ### (optional) Tuning Constant. Default is 1.
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
  

    # development
    options(stringsAsFactors = F)
    fName <- "rawBat.csv"
    fLoc <- "C:/LACEY/CODE/R/Triangulate/"
    raw <- read.csv(paste0(fLoc,fName), header=TRUE)
    colTrackID = "Bat"
    colDate = "Date"
    colBearing ="Azimuth"
    colZone = "Zone"
    colEasting ="Easting"
    colNorthing = "Northing"
    

    columnNames <- c(colTrackID,
                     colDate,
                     colBearing,
                     colZone,
                     colEasting,
                     colNorthing)
    
    
    # Error handeling
    stringArgs <- c(fLoc, fName, colTrackID, colDate, colBearing,
                    colZone, colEasting, colNorthing)
    if(!all(is.character(stringArgs))){
      stop("The arguments fLoc, fName, colTrackID, colDate, colBearing, colZone, colEasting, colNorthing should be character stings and surrounded by quotations")
    }
    
    if (!all(columnNames %in% names(raw))) {
      missingCol <- columnNames[! columnNames %in% names(raw)]
      stop(print(paste0("The column names specified do not exist in your dataset: ", 
                        missingCol)))
    }
    
    if ( length(raw[duplicated(raw[,columnNames]), ][,1]) > 0){
      warning("The record(s) above are duplicated and will be removed from analysis. Further qaqc is likely needed. See error report.")
      print(raw[duplicated(raw[,columnNames]), ])
      errorDuplicate <- raw[duplicated(raw[,columnNames]), ]
      errorDuplicate$ErrorMessage <- "Duplicated record"
      raw <- raw[-which(duplicated(raw[,columnNames])), ]
    }
    
    if (conf != .9 && conf != .95) {
      stop("Confidence level must be either .90 or .95")
    }
    
    method = casefold(method)
    if(!method %in% c("m", "a", "h")){
      stop("Input for method must be either \'m\', \'a\', or \'h\'")
    }
    
    # if( length(k) > 1 ){
    #   if( any(is.na(k)) ){   
    #     stop("Missing standard deviations not allowed")
    #   }
    # }
    # if( length(k) > 1 ){
    #   if( n != length(k) ){
    #     stop("Length of k must equal length of angles")
    #   }
    # }
    # 
    # if( n != nrow(station.locs) ){
    #   stop("Number of rows in station.locs much equal length of angles")
    # }
    
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
    
    
    # Designate a coordinate system
    zone <- as.numeric(gsub("[[:alpha:]]","", raw[1,colZone]))
    crsProj <- paste("+proj=utm +zone=", zone, " +ellps=GRS80 +datum=NAD83 +units=m +no_defs", sep="")
    utm <- sp::CRS(crsProj)
    
    
    # Create shapefile for Triangulation_Stations
    coords <- cbind(raw[,colEasting], raw[,colNorthing])
    raw.df <- sp::SpatialPointsDataFrame(coords, raw, proj4string = utm, match.ID = T)
    suppressWarnings(rgdal::writeOGR(raw.df, dsn= paste0(fLoc,"/Shapefiles") ,layer="Triangulation_Stations",
             driver="ESRI Shapefile"))
   

    # Error handeling
    # handeling here so that this station is included in station shapefile created above
    if (any(uniqueTrack$Signals == 1)){
      warning("The record(s) above have only 1 bearing records and will be removed from analysis. Further qaqc is likely needed. See error report.")
      print(uniqueTrack[uniqueTrack$Signals == 1, ])
      error1Station <- raw[raw$EstimateID %in% uniqueTrack$EstimateID[uniqueTrack$Signals == 1], ]
      error1Station$ErrorMessage <- "Only one bearing"
      raw <- raw[!raw$EstimateID %in% uniqueTrack$EstimateID[uniqueTrack$Signals == 1], ]
    }
    
    
    # Write all Error report
    errors <- rbind(errorDuplicate[, c(columnNames,"ErrorMessage")], 
                    error1Station[, c(columnNames, "ErrorMessage")])
    write.csv(errors,
              paste0(fLoc,"/ErrorReport.csv"), row.names = FALSE)

    
    
    #   print(sprintf("BatID: %s on %s - %d records removed from analysis due to duplicate locations", 
    #                 dupLocs$BatID[dupLocs$x > 1], dupLocs$Date[dupLocs$x > 1], dupLocs$x[dupLocs$x > 1]))
    # }
    
    
    # Helper functions
    # ---- Function to compute condition number for 2x2 matrix.
    #      COND2(X) returns the 2-norm condition number (the ratio of the
    #      largest singular value of X to the smallest).  Large condition
    #      numbers indicate a nearly singular matrix.
    #
    #       Modified to handle only the 2-norm, eliminate some error checking,
    #       and other condition testing.
    #       Chris Nations  8/20/02
    #       further modified by Trent 30sep09, to add eps as agrument
    #       further modified by Lacey Jeroue 6/2015 for inclusion in this function
    
    cond2 <- function(x, eps=.Machine$double.eps){
      sv <- svd(x)$d
      if( any( abs(sv) <= eps )){   # not sure need abs
        ans <- Inf
      } else {
        ans <- max(sv) / min(sv)
      }
      ans
    }
    
    # ---- Function to compute angle of 2D system
    cart2pol <- function(x,y){
      atan2(y,x)
    }
    
    
    ## ESTIMATE TARGET LOCATION ####

    getTriangulatedPoint <- function(data = raw, kappa=k, tuningConstant = tc){
      
      ansList <- list()
      triangulatedPoints <- NULL
      
      for(i in unique(raw$EstimateID)){  # i = 1
        station.locs <- cbind(raw[raw$EstimateID == i, colEasting],
                              raw[raw$EstimateID == i, colNorthing])
        angles <- raw[raw$EstimateID == i, colBearing]
        
        # ---- Main code
        cond.num = 1e15        # For test of singularity.
        dist1 = 1             
        x = station.locs[,1] 
        y = station.locs[,2] 
          
        theta = (90-angles)*pi/180  # Convert bearing to standard radian measure
        theta = as.numeric(theta < -pi) *2*pi + theta 
        
        s = sin(theta)  
        c = cos(theta) 
        z = s*x - c*y 
        sstar = s 
        cstar = c 
        n = length(angles)
        w = rep(1, n) 
        M1 = rbind(sstar, -cstar) %*% cbind(s, -c) 
        converge = FALSE 
        trace=FALSE
        iter = 0
        maxiter = 50
        robust = method %in% c("a", "h")
        
        if( all(is.na(k)) ){
          kest = TRUE
        } else {
          kest = FALSE
        }
        
        failed = c("less than 2 bearings",
                   "negative variance estimates",
                   "solution behind stations",
                   "failed to converge", 
                   "parallel bearings")
        
        if( cond2(M1)<cond.num ){
          M2 = rbind(sstar, -cstar) %*% z 
          xyhat = solve(M1,M2) 
          
          # --- Main loop
          
          while (!converge & (iter<maxiter)){
            iter = iter+1 
            xyold = xyhat 
            
            d = sqrt( (xyhat[1]-x)^2 + (xyhat[2]-y)^2 )
            
            if( trace ){
              cat(paste("--- Iteration:", iter, "Current solution: "))
              cat(xyhat)
              cat("\n")
            }
            
            if( (n>2) & robust){ 
              # Need 3 or more bearings to calculate weights for robust methods
              dxy = c(xyhat)-rbind(x, y)       
              muhat = cart2pol(dxy[1,],dxy[2,]) 
              Cd = cos(theta-muhat) 
              if( kest ){
                Cbar = abs(sum(w*Cd)/sum(w))^(n/(n-2))   # Abs is ad hoc but may avoid temporary numeric problems in iteration
                k = (2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar)
                k = ifelse( abs(k) < .Machine$double.eps, Inf, 1/k)  # k is Inf if bearings all cross at single point. Happens with n=2 bearings.
              }
              tt = sqrt(2*k*(1-Cd)) 
              if( method == "a" ){
                phit = tc*sin(tt/tc)*(abs(tt)<(tc*pi)) 
              } else {
                phit = sign(tt)*apply( cbind(abs(tt), tc), 1, min) 
              }
              same = (Cd==1)      # Take care when estimated & observed bearings
              tt[same] = 1         #   are identical  avoid division by zero.
              phit[same] = 1      #   Usually occurs with just 2 intersecting bearings.
              w = phit/tt 
            }    #  if robust>1
            
            sstar = w*(xyhat[2]-y)/(d^3) 
            cstar = w*(xyhat[1]-x)/(d^3) 
            
            M1 = rbind(sstar, -cstar) %*% cbind(s, -c)
            
            if( ((n-sum(!as.logical(w)))>1) & (cond2(M1)<cond.num)){
              M2 = rbind(sstar,  -cstar) %*% z 
              xyhat = solve(M1, M2) 
              converge = sum(abs(xyhat-xyold)<dist1)==2 
            } else {
              break    # If either condition above occurs, convergence will very
            }          #    likely fail, so break out of while loop.
          }  # while
          
          
          if( converge ){
            dxy = c(xyhat) - rbind(x, y) 
            muhat = cart2pol(dxy[1,],dxy[2,]) 
            Cd = cos(theta-muhat) 
            if( kest & (n>2) ){
              Cbar = (sum(w*Cd)/sum(w))^(n/(n-2))    # Exponent is small sample size correction
              k = (2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar) 
              k = ifelse( abs(k) < .Machine$double.eps, Inf, 1/k )   # k = Inf only if perfect solution w/ n>2.  Only in pathological case
            }
            
            if( sum(Cd*w) > 0){     # Weighted average of cosine differences (check 
              VM = c(xyhat)           #   on bad solutions behind DASARs)
              names(VM) = c("x","y")
              
              if( kest & (n==2) ){     # Cannot estimate Qhat with only 2 bearings
                Qhat = matrix( NA, 2, 2) 
                outcome = "successful  2 bearings  no Kappa" 
              } else if( kest & all(is.infinite(k)) ){   #   Perfect solution -> Zero variance
                Qhat = matrix( 0, 2, 2) 
                outcome = "successful  perfect solution  Zero Kappa estimate" 
              } else {
                cv = -(sum(k*sstar*c) + sum(k*cstar*s))/2 
                M3 = cbind(c(sum(k*sstar*s), cv), c(cv, sum(k*cstar*c))) 
                Qhat = solve(M3) 
                
                if( all(diag(Qhat)>0) ){
                  outcome = "successful"        # Successful solution
                } else {
                  outcome = failed[2]         # Implausible variance estimate(s)
                }
              }
              
            } else {
              outcome = failed[3]           # Bad solution behind DASARs
            } # if all(Cd>0)
            
          } else {
            outcome = failed[4]             # No convergence
          }   # if converge
          
        } else {
          outcome = failed[5]
        } # if cond2
        
        
      # }     # if n<=1
      
      if(outcome %in% failed){
        VM = c(NA,NA) 
        Qhat = matrix( NA, 2,2)
        w = rep(0,n) 
      }
      
      dimnames(Qhat) = list(names(VM), names(VM))
      names(w) = dimnames(station.locs)[[1]]
      
      ansList[[i]] <- list( EstimateID=i, loc.est = VM, qhat=Qhat, angle = angles,  
                            outcome=outcome, w=w, k=k, tc=tc, method=method )
      
      t <- data.frame(cbind(data.frame(VM)[1,], data.frame(VM)[2,]))
      names(t) <- c("estEasting", "estNorthing")
      t$EstimateID <- i
      t$TrackID <- unique(raw[i, colTrackID])
      triangulatedPoints <- data.frame(rbind(triangulatedPoints, t))
    } # End Main loop
    } # End function
    

      # do.call("rbind", lapply(ansList, "[[", ))
    
    # Create Shapefile of the estimated triangulated Points
    triangulatedPoints <- merge(triangulatedPoints, unique(raw[,c("EstimateID",colDate)]))
    triangulatedPoints <- triangulatedPoints[!is.na(triangulatedPoints$estEasting),]
    coords <-cbind(triangulatedPoints$estEasting, triangulatedPoints$estNorthing)
    t.spdf <- sp::SpatialPointsDataFrame(coords, triangulatedPoints, proj4string = utm, match.ID = T)
    rgdal::writeOGR(t.spdf, dsn= paste0(fLoc,"/Shapefiles") ,layer="EstimatedLocations",
                              driver="ESRI Shapefile")
    
    
    ## CONSTRUCT CONFIDENCE REGION AND CREATE ELLIPSE SHAPEFILE
    
    hold <- list()
    DF <- NULL
    
    for(i in 1:length(ansList))
    {
      if(is.na(ansList[[i]]$qhat[1,1]))
      {
        next
      } else {
        n <- ansList[[i]]$n
        id <- ansList[[i]]$FixID[1]
        xc <- ansList[[i]]$loc.est[1] # x center pt
        yc <- ansList[[i]]$loc.est[2] # y center pt
        s <- ansList[[i]]$qhat        # covariance matrix
        eig <- eigen(s)
        area_m2 <- (pi*sqrt(sqrt(det(s)^2))*qchisq(conf, 2)) # area of ellipse in m
        d <- sort(2*sqrt(qchisq(conf, 2)*sqrt(eig$val^2))) # length of both axis
        a <- d[2]       # major axis length
        b <- d[1]       # minor axis length
        v <- eig$vect[,1] # dominant eigenvector
        ang <- atan(v[2]/v[1]) # Absolute orientation of ellipse major axis
        if(ang < 0){    # Returns an angle between 0 and pi following 
          ang <- ang+pi # math convention (0 radians = 90 degrees, East;
        }               # pi radians = 270 degrees, West).
        cth <- cos(ang)
        sth <- sin(ang)
        c <- seq(0,200,1)*2*pi/200
        acc <- a*cos(c)
        bsc <- b*sin(c)
        xp <- acc*cth - bsc*sth + xc # x points on ellipse
        yp <- acc*sth + bsc*cth + yc # y points on ellipse
        m <- cbind(xp[-c(201)], yp[-c(201)], xp[c(-1)], yp[c(-1)])
        names(m) <- c("x", "y", "x", "y")
        df <- data.frame(subset(raw, raw$FixID==id))
        df <- df[1,c(1,6:9)]
        df$A.m2 <- round(area_m2, 1)
        df$A.km2 <- round(df$A.m2 * 0.000001, 4)
        hold[[i]] <- Lines(list(Line(data.frame(cbind(xp,yp)))), ID=id)
        DF <- rbind(DF, df)
      } # End 'if/else'
    } # End 'for-loop'
    
    # Create Ellipse Shapefile
    SLines <- SpatialLines(hold[!sapply(hold, is.null)], proj4string = utm)
    list_of_Lines <- slot(SLines, "lines")
    elip <- SpatialPolygons(lapply(list_of_Lines, function(x) {
      Polygons(list(Polygon(slot(slot(x, "Lines")[[1]],
                                 "coords"))), ID = slot(x, "ID"))
    }), proj4string = utm)
    
    SPDF <- SpatialPolygonsDataFrame(elip, DF, match.ID = F)
    suppressWarnings(writeOGR(SPDF, dsn= fDest ,layer="Confidence_Region",
                              driver="ESRI Shapefile"))
    
    
    ## PREPARE DATA AND CREATE BEARING SHAPEFILES
    
    # Get all possible station combinations for bearing intersections
    tmp3 <- aggregate(rawClean[,2], by=list(FixID=rawClean$FixID), length)
    combos <- NULL
    for(i in 1:length(tmp3[,1])){
      combo1 <- combn(rawClean$StationID[rawClean$FixID==tmp3[i,1]], 2)
      melt <- melt(combo1, id.vars=c())
      melt$StationID <- rep(c("Station_1", "Station_2"), length(melt[,1])/2)
      combo2 <- dcast(melt, Var2 ~ StationID, value.var="value")
      combo3 <- combo2[,-c(1)]
      combos <- rbind(combos, combo2)
    }
    
    # Create Data Frame for holding Intersections and Bearings
    tmp4 <- data.frame("StationID" = combos[,2])
    st1 <- join(tmp4, rawClean, "StationID")
    tmp5 <- data.frame("StationID" = combos[,3])
    st2 <- join(tmp5, rawClean, "StationID")
    int <- cbind(st1[,2], st1[,1], st1[,3], st1[,5:6], st2[,1], st2[,3], st2[,5:6])
    names(int) <- c("FixID", "StID_1", "Azimuth_1",  "Easting_1", "Northing_1",
                    "StID_2", "Azimuth_2",  "Easting_2", "Northing_2")
    
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
    # tmp8 <- cbind(int5[,1:2],int5[,7:9], int5[,6], int5[,14:15])
    tmp8 <- cbind(int5[,1],int5[,6:9], int5[,2], int5[,10:11])
    names(tmp7) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2", "Easting2", "Northing2")
    names(tmp8) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2", "Easting2", "Northing2")
    int6 <- arrange(rbind(tmp7,tmp8), FixID)
    head(int6)
    
    # Reshape dataframe to fields of interest for each POLYLINE
    tmp4.1 <- cbind(int4[,1:6])
    tmp4.2 <- cbind(int4[,1],int4[,6:9], int4[,2])
    names(tmp4.1) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2")
    names(tmp4.2) <- c("FixID", "StID1", "Azimuth", "Easting1", "Northing1",
                     "StID2")
    int4.3 <- arrange(rbind(tmp4.1,tmp4.2), FixID)
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
      writeOGR(SLDF1, dsn= fDest ,layer="Bearings1",driver="ESRI Shapefile")
      
    } else {
      
      hold <- list()
      for(i in 1:length(int7[,1])){
        hold[[i]] <- Lines(list(Line(data.frame(rbind(int7[i,4:5], int7[i,7:8])))), ID=i)
      }
      SL1 <- SpatialLines(hold, proj4string = utm)
      SLDF1 <- SpatialLinesDataFrame(SL1, int7.1, match.ID = F)
      writeOGR(SLDF1, dsn= fDest ,layer="Bearings2",driver="ESRI Shapefile")
    }
    return(list( ansList, t.pts2, DF ))
})
