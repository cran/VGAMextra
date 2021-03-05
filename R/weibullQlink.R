##########################################################################
# These functions are
# Copyright (C) 2014-2021 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


# 01012021... weilbullRPercentilelink
# renamed as weibperclink().  then weibquanlink 200121 then weibullQlink

weibullQlink <-  function(theta, 
                          percentile = stop("Enter percentiles."),
                          shape = NULL, wrt.param = NULL,
                          bvalue = NULL, inverse = FALSE,
                          deriv = 0, short = TRUE, tag = FALSE) {
 
   if (is.character(theta)) {
    string <- if (short)
      paste("weibquanlink(",theta,"; ", percentile, 
             ", shape)", sep = "") else
        paste("log[",
        as.char.expression(theta),"* [-log (1 - p)]^(1/shape)]")
    
    if (tag) 
      string <- paste("2-parameter Weibull Quantile link:", string)
    return(string)
    
  }
  
  if (!is.Numeric(percentile, positive = TRUE) || (any(percentile > 100)))
    stop("Wrong input for 'percentile'.")
  
  #if (length(percentile) > 1)
  #  stop("Only single values (no vectors) handled by 'percentile' atm.")
  
  if (!length(shape))
    stop("Enter a valid value for 'shape'.")
  
  theta <- cbind(theta)
  mydim <-dim(theta)
  shape <- cbind(shape)
  
  ### Note this:
  percentile <- matrix(percentile/1e2, nrow = nrow(theta), 
                       ncol = length(percentile), byrow = TRUE)
  
  
  
  if (ncol(shape) != ncol(theta)) {
    warning("Unequal number of cols for 'theta' and 'shape'.")
    shape <- matrix(shape, nrow = nrow(theta), ncol(theta), byrow = TRUE)
  }
    
  
  if (nrow(shape) != nrow(theta)) {
    vec2com <- c(nrow(shape), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      shape <- matrix(shape, nrow = max(vec2com), ncol = ncol(shape))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(theta))
    }
    warning("Unequal number of rows for arguments 'theta' and 'shape'.")
  }
  
  ## 20210121
  if ( ncol(theta) != ncol(percentile) )
    theta <- matrix(theta, nrow = nrow(theta),
                    ncol = ncol(percentile), byrow = FALSE)
  
  if ( ncol(shape) != ncol(percentile) )
    shape <- matrix(shape, nrow = nrow(shape),
                    ncol = ncol(percentile), byrow = FALSE)

  
  if ( length(wrt.param) && (!(wrt.param %in% 1:2)  ||
      !is.Numeric(wrt.param, length.arg = 1)) )
    stop("Bad input for argument wrt.param.")
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {

    theta.ret <- switch(deriv + 1,
                        # theta
                        exp(theta) * (-log(1 - percentile))^(-1/shape) ,
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(dscale / d eta1 , dshape/ d eta1 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(dscale / d eta2 , dshape / d eta2)g
                          cbind( matrix(0, nrow(theta), ncol(theta)), 
                                 matrix(0, nrow(theta), ncol(theta)) )
                        },
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(d2scale / deta12 , d2shape / deta12 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(d2scale / deta2 , d2shape / deta22)
                          # d2shape / deta2 Computed by the other link
                          help <- matrix(0, nrow(theta), ncol(theta))
                          cbind(help, NULL)
                        })
    
    eta.ret <- theta.ret
  } else {
    pp <- percentile
    ## Here, 1. eta
    ##       2. d eta / d thetaj, j = wrt = 1, 2 (scale or shape)
    ##       3. d2eta / dthetaj2, j = wrt = 1, 2.
    eta.ret <- switch(deriv + 1,
                      log(theta) + log(-log(1 - pp))/shape,
                      #  wrt = 1 -> wrt scale (= theta)
                      #  wrt = 2 -> wrt shape
                      
                      # d eta1 / dscale (d eta2 / dshape = 0 not returned)
                      if (wrt.param == 1) 1 / theta else 
                        -log(-log(1 - pp)) / shape^2,
                      # ditto
                      if (wrt.param == 1) -1 / theta^2 else 
                        2 * log(-log(1 - pp)) / shape^3)
  }
  
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL

  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim

  eta.ret
  
}