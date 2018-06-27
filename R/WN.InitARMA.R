##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
#
# Returns initial values of ARMAff(), except noiseSD (noiseVar).
# -> c(mean, theta1, theta2, ..., theta_p, phi1, ..., phi_q) 
# If whiteN = TRUE, it also returns the estimated white noise (required
# by ARMAff() and MAff() ).  Modified on 2015/11/04


WN.InitARMA <- function(tsData = NULL, 
                        order  = c(1, 0, 1),
                        whiteN = FALSE, 
                        moreOrder = 0,
                        updateWN  = FALSE) {
  
  if (length(tsData) && !is.data.frame(tsData))
    stop("Data must be a 'data.frame' object.")
  
  if (dim(tsData)[2] > 1)
    stop("'tsData' must be a one-column data frame.")
  
  if (length (order) != 3 || !is.Numeric(order, integer.valued = TRUE)) 
    stop("Invalid input for argument 'order'.") 
  
  if (( moreOrder < 0 ) || !is.Numeric(moreOrder, integer.valued = TRUE))
    stop("Wrong input for argument 'moreOrder'.")
  
  if (!is.logical(whiteN))
    stop("Bad input for argument 'whiteN'.")
  
  if (!is.logical(updateWN))
    stop("Wrong input for argument 'updateWN'.")
  
  arOrd <- order[1]; difOrd <- order[2]; maOrd <- order[3]
  #if (difOrd != 0)
  #  stop("Integrated time series not implemented yet.")
  
  initials  <- vector("list", 2)
  names(initials) <- c("Coeffs" , "WhiteNoise")
  initCoeff <- numeric(0)
  order2    <- arOrd + maOrd
  
  myY <- as.matrix(tsData); nn <- nrow(myY) ; names(tsData) <- c("y")
  newOrd <- if (!maOrd) arOrd else min(4 + moreOrder, 4)
  inifit <- vglm(y ~ 1,
                 ARXff(order = newOrd,
                       type.EIM  = "exact",
                       var.arg   = TRUE, 
                       lARcoeff   = "identitylink",
                       noChecks  = TRUE), 
                 crit = "coe", 
                 trace = FALSE, data = tsData,
                 maxit = 20, epsilon = 1e-5,
                 smart = FALSE, noWarning  = TRUE)

  thePhiEst <- Coef(inifit)[3:(newOrd + 2)]
  
  # mean(myY) * (1 - sum(thePhiEst)) 2016/06/10
  initDri   <- if (!arOrd) mean(myY) else Coef(inifit)[1]
  initCoeff <- c(initCoeff, initDri)
  
  if ( maOrd ) { # First stage: Estimating white noise from the AR model #
    
    whiNo <- matrix(0.0, nrow = nn, ncol = 1)
    m.Sum <- matrix(0.0, nrow = nn, ncol = 1)
    for (ii in 1:newOrd) {
      m.pars <- matrix(0.0, nrow = nn, ncol = 1)
      laG.Y  <- matrix(0.0, nrow = nn, ncol = 1)
      m.pars[-(1:ii), 1] <- matrix(thePhiEst[ii], nrow = nn - ii, ncol = 1)
      laG.Y[-(1:ii), 1]  <- myY[1:(nn -ii), 1]
      m.Sum[, 1] <- m.Sum[, 1] + m.pars * laG.Y
    }
    m.Sum <- m.Sum + Coef(inifit)[1]
    whiNo[1:newOrd, 1] <- myY[1:newOrd, 1] - m.Sum[1:newOrd, 1]
    
    for (ii in 1:newOrd)
      whiNo[(newOrd + 1):nn, ] <- whiNo[(newOrd + 1):nn, ] +
                     thePhiEst[ii] * myY[(newOrd - ii + 1):(nn - ii), ]
    
    whiNo[-(1:newOrd), 1] <- myY[-(1:newOrd), 1] - 
                                 (initDri  + whiNo[-(1:newOrd), 1]) 
    
    if (!maOrd) {
      myYshift <- myY - Coef(inifit)[1]/( 1 - sum( thePhiEst ) )  
    } else {
      myYshift <- myY - Coef(inifit)[1]
    }
    
  } else {
    initials[[2]] <- NA
  }
  
  if (arOrd && maOrd) {
    
    # Regress yt onto yt-1,..., yt-p, Wt-1, ..., Wt-q #
    MaxOrd <- max(arOrd, maOrd)
    ysOnly <- matrix(0.0, nrow = nn - MaxOrd, ncol = arOrd )
    WNOnly <- matrix(0.0, nrow = nn - MaxOrd, ncol = maOrd )
    
    for (kk in 1:arOrd) 
      if (MaxOrd - kk + 1 > 0) {
        ysOnly[, kk] <- myYshift[(MaxOrd - kk + 1):(nn - kk), ]
      } else {
        stop("Insufficient number of observations.")
      } 
    
    for (ji in 1:maOrd)
      if (MaxOrd - ji + 1 > 0) {
        WNOnly[, ji] <- whiNo[(MaxOrd - ji + 1):(nn - ji), ] 
      } else {
        stop("Insufficient number of observations.")
      }
    
    ysAndWN <- cbind(ysOnly, WNOnly)
    initCoeff <-  c(initCoeff, lsfit(ysAndWN, myYshift[-(1:MaxOrd), ], 
                                             intercept = FALSE)$coef)
    whiNo.fin <- whiNo
    
    if (updateWN) { # only initials.
      
      for (ii in 1:4) {
        # Regress again. 
        whiNo2 <- matrix( 0.0 , nrow = nn, ncol = 1 )
        lag.y  <- matrix( 0.0 , nrow = nn, ncol = arOrd )
        w.noi  <- matrix( 0.0 , nrow = nn, ncol = maOrd )
        whiNoo <- matrix(initDri, nrow = nn, ncol = 1)
        
        for (kk in 1:arOrd) {
          lag.y[-(1:kk), kk] <- myY[1:(nn - kk), ]
          whiNoo <- whiNoo + initCoeff[kk + 1] * lag.y[, kk]
        }
        
        for (kk in 1:maOrd) {
          w.noi[-(1:kk), kk] <- whiNo[1:(nn - kk), ]
          whiNoo <- whiNoo + initCoeff[kk + 1 + arOrd] * w.noi[, kk]
        }
        
        whiNo2 <- myY - whiNoo
        WNOnly2 <- matrix( NA, nrow = nn - MaxOrd, ncol = maOrd )
        for ( ji in 1:maOrd )
          WNOnly2[, ji] <- whiNo2[(MaxOrd - ji + 1):(nn - ji), ]
        
        # Model 2 #
        ysAndWN <- cbind(ysOnly, WNOnly2)
        initCoeff[-1] <- c(lsfit(ysAndWN, myYshift[-(1:MaxOrd)], 
                                 intercept = FALSE)$coef)
        
        whiNo <- whiNo2
      }
      colnames(whiNo2) <- NULL
      initials[[2]] <- if (whiteN) whiNo.fin else NA
    } else {
      initials[[2]] <- if (whiteN) whiNo else NA
    }
  }
  
  if (!arOrd && maOrd) {
    
    WNOnly <- matrix(NA_real_, nrow = nn - maOrd, ncol = maOrd)
    for (ji in 1:maOrd)
      WNOnly[, ji] <- whiNo[(maOrd - ji + 1):(nn - ji), ]
    
    initCoeff <- c(initCoeff, lsfit(WNOnly, myYshift[-(1:maOrd),], 
                                    intercept = FALSE)$coef)
    whiNo.fin <- whiNo
    
    if (updateWN) { # initials plus WN
      
      store.WN <- matrix(NA, nrow = nn, ncol = 5)
      store.co <- matrix(NA, nrow = 5, ncol = 1 + maOrd)
      store.WN[, 1] <- whiNo
      store.co[1, ] <- initCoeff
      store.ss <- numeric(4)
      
        for (ii in 1:5) {
          w.noi <- matrix( 0.0 , nrow = nn, ncol = maOrd )
          for (kk in 1:maOrd) 
            w.noi[-(1:kk), kk] <- whiNo[1:(nn - kk), ]
          
          ss.1 <- matrix(initCoeff[2:(1 + maOrd)], 
                         nrow = nn, ncol = maOrd, byrow = TRUE)
          ss.1 <- initDri + rowSums(ss.1 * w.noi) + whiNo
          store.ss[ii] <- sqrt( sum((ss.1 - myY)^2) / (nn - maOrd) )
          if (ii == 5)
            break()
          # Regress again
          whiNo2 <- matrix( 0.0 , nrow = nn, ncol = 1 )
          w.noi  <- matrix( 0.0 , nrow = nn, ncol = maOrd )
          whiNoo <- matrix(initDri, nrow = nn, ncol = 1)
          
          for (kk in 1:maOrd) {
            w.noi[-(1:kk), kk] <- whiNo[1:(nn - kk), ]
            whiNoo <- whiNoo + initCoeff[kk + 1 + arOrd] * w.noi[, kk]
          }
          
          whiNo2 <- myY - whiNoo
          WNOnly2 <- matrix( NA, nrow = nn - maOrd, ncol = maOrd )
          for ( ji in 1:maOrd )
            WNOnly2[, ji] <- whiNo2[(maOrd - ji + 1):(nn - ji), ]
          
          # Model 2 #
          ysAndWN <- cbind(WNOnly2)
          initCoeff[-1] <- c(lsfit(ysAndWN, myYshift[-(1:maOrd)], 
                                   intercept = FALSE)$coef)
          store.WN[, ii + 1] <- whiNo2
          store.co[ii + 1, ] <- initCoeff
          
          whiNo <- whiNo2
        }
      k.coe <- which(store.ss == min(store.ss))[1]
      whiNo2 <- store.WN[, k.coe, drop = FALSE]
      initCoeff <- store.co[k.coe, ]
      colnames(whiNo2) <- NULL
      initials[[2]] <- if (whiteN) whiNo2 else NA
    } else {
      initials[[2]] <- if (whiteN) whiNo else NA
    }
  }
  
  if (arOrd && !maOrd) {
    
    if (whiteN && FALSE)
      stop("No white Noise to estimate in a AR process", 
           " Set 'whiteN = FALSE'.")
    if (updateWN && FALSE)
      stop("No white Noise to update in a AR process.", 
           " Set 'updateWN = FALSE'.")
    
    GammaS <- array(NA, dim = arOrd + 1)
    GammaS[1] <-  cov(myY, myY)  
    for (jj in 1:arOrd) 
      GammaS[jj + 1] <- cov(myY[-(1:jj) , ], myY[-((nn - jj + 1):nn),] )

    myToeplitz <- toeplitz(c(GammaS[1:arOrd]))
    myresp     <- c(GammaS[2:( arOrd + 1)])
    initCoeff  <- c(initCoeff,  solve(myToeplitz, myresp))

    # No white noise to estimate. GammaS is returned #
    names(initials) <- c("Coeffs", "Covariances")
    initials[[2]] <- GammaS
  }
  
  arNames <- if (arOrd) paste("arInit", 1:arOrd, sep ="") else NULL
  maNames <- if (maOrd) paste("maInit", 1:maOrd, sep ="") else NULL
  names(initCoeff) <- c("initDri", arNames, maNames)
  initials[[1]] <- initCoeff
  
  initials  
}

