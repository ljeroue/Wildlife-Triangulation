# Helper functions
# ---- Function to compute condition number for 2x2 matrix.
#      COND2(X) returns the 2-norm condition number (the ratio of the
#      largest singular value of X to the smallest).  Large condition
#      numbers indicate a nearly singular matrix.
#
#       Modified to handle only the 2-norm, eliminate some error checking,
#       and other condition testing.


cond2 <- function(x, eps=.Machine$double.eps){
  sv <- svd(x)$d  # singular values of x
  if( any( abs(sv) <= eps )){   # not sure need abs
    ans <- Inf
  } else {
    ans <- max(sv) / min(sv)
  }
  ans
}


## ESTIMATE TARGET LOCATION ####
estLocation <- function(raw, method=method, kappa=k, tuningConstant = tc){
  
  ansList <- list()
  pointEstimate <- NULL
  
  for(i in unique(raw$EstimateID)){  # i = 1
    station.locs <- cbind(raw[raw$EstimateID == i, colX],
                          raw[raw$EstimateID == i, colY])
    angles <- raw[raw$EstimateID == i, colBearing]
    
    # ---- Main code
    cond.num = 1e15        # For test of singularity.
    dist1 = 1             
    x = station.locs[,1] 
    y = station.locs[,2] 
    
    # Convert bearing to standard radian measure
    theta = (90-angles)*pi/180 
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
    trace=T
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
          # muhat = cart2pol(dxy[1,],dxy[2,]) 
          muhat = atan2(dxy[2,],dxy[1,]) # compute angle of 2D system
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
        }    #  if robust
        
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
        muhat = atan2(dxy[2,],dxy[1,]) # compute angle of 2D system
        Cd = cos(theta-muhat) 
        if( kest & (n>2) ){
          Cbar = (sum(w*Cd)/sum(w))^(n/(n-2))    # Exponent is small sample size correction
          k = (2*(1-Cbar)+(1-Cbar)^2*(0.48794-0.82905*Cbar-1.3915*Cbar^2)/Cbar) 
          k = ifelse( abs(k) < .Machine$double.eps, Inf, 1/k )   # k = Inf only if perfect solution w/ n>2.  Only in pathological case
        }
        
        if( sum(Cd*w) > 0){     # Weighted average of cosine differences (check 
          VM = c(xyhat)           #   on bad solutions behind DASARs)
          names(VM) = c("xHat","yHat")
          
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
    
    # t <- data.frame(cbind(data.frame(VM)[1,], data.frame(VM)[2,]))
    # names(t) <- c("estEasting", "estNorthing")
    # t$EstimateID <- i
    # t$TrackID <- unique(raw[i, colTrackID])
    # triangulatedPoints <- data.frame(rbind(triangulatedPoints, t))
  } # End Main loop
  
  # } # End function
  
  
   return(ansList)

} #end