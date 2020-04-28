##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
### In this function vector'y' is lagged --> yy[1:(n - ii) , 1] ###

  cross.gammas <- function(x, y = NULL, lags = 1) {
    xx  <- matrix(x, ncol = 1)
    nx  <- nrow(xx)
    
    if (lags < 0 || !Is.Numeric(lags, isInteger = TRUE))
      stop("'lags' must be a non-negative integer.")
  
    if (length(y)) {
      yy <- matrix(y, ncol = 1)
      ny <- nrow(yy)
      if (nx != ny)
        stop("Number of rows differs.") else
          n <- nx
    } else {
      yy <- xx
      n  <- nrow(xx)
    }
  
    myD <- numeric(lags + 1)
    myD[1] <- if (length(y)) cov(xx, yy) else cov(xx, xx)  # i.e. var(xx)    
    if (lags > 0)
      for (ii in 1:lags) 
        myD[ii + 1]  <- cov(xx[-(1:ii), 1], yy[1:(n - ii) , 1])
  
    myD
  }
  