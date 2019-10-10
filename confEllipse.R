## CONSTRUCT CONFIDENCE REGION AND CREATE ELLIPSE SHAPEFILE

confEllipse <- function(estimatedLocations, confLevel){
  ellipsePoints <- list()
  area <- NULL
  
  for(i in 1:length(estimatedLocations))
  {
    if(is.na(estimatedLocations[[i]]$qhat[1,1]))
    {
      next
    } else {
      n <- length(estimatedLocations[[i]]$angle)
      id <- estimatedLocations[[i]]$EstimateID[1]
      xc <- estimatedLocations[[i]]$loc.est[1] # x center pt
      yc <- estimatedLocations[[i]]$loc.est[2] # y center pt
      s <- estimatedLocations[[i]]$qhat        # covariance matrix
      eig <- eigen(s)
      area_m2 <- (pi*sqrt(sqrt(det(s)^2)) * qchisq(confLevel, 2)) # area of ellipse in m
      d <- sort(2*sqrt(qchisq(confLevel, 2) * sqrt(eig$val^2))) # length of both axis
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
      df <- data.frame(subset(raw, raw$EstimateID==id))
      df <- df[1,c(1,6:9)]
      df$A.m2 <- round(area_m2, 1)
      df$A.km2 <- round(df$A.m2 * 0.000001, 4)
      ellipsePoints[[i]] <- sp::Lines(list(sp::Line(data.frame(cbind(xp,yp)))), ID=id)
      area <- rbind(area, df)
    } # End 'if/else'
  } # End 'for-loop'
  
  return(list(ellipsePoints=ellipsePoints, area=area))
}
