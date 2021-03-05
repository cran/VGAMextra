##########################################################################
# These functions are
# Copyright (C) 2014-2021 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


# 01012021... 
# GammaRMeanlink defined as log(shape / rate) Renamed on Feb 2021
# As gammaRMlink().
gammaRMlink <- function(theta, shape = NULL, wrt.param = NULL,
                           bvalue = NULL, inverse = FALSE,
                           deriv = 0, short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("gammaRMlink(", theta, "; shape)", sep = "") else
        paste("log (shape /", as.char.expression(theta), ")", sep = "")
    
    if (tag)
      string <- paste("2-parameter Gamma mean link: ", string)
    return(string)
  }
  
  theta <- cbind(theta)  
  mydim <- dim(theta)
  shape.mat <- if (length(shape)) cbind(shape) else NULL

  
  if (ncol(shape.mat) != ncol(theta)) {
    stop("# columns of 'theta' and 'shape' do not match. ")
  }
  
  if (nrow(shape.mat) != nrow(theta)) {
    vec2com <- c(nrow(shape.mat), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      shape.mat <- matrix(shape.mat, nrow = max(vec2com), 
                          ncol = ncol(theta))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(shape.mat))
    }
    warning("# rows of 'theta' and 'shape' do not match." )
  }
  
  if (length(wrt.param) && (!(wrt.param %in% 1:2) ||
                            !is.Numeric(wrt.param, length.arg = 1)) )
    stop("argument 'wrt.param' should be 1 or 2")
  
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    theta.ret <- switch(deriv + 1,
                    # theta
                    exp(-theta) * shape.mat,
                    
                    #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                    if (wrt.param == 1) {
                          
                      #  c(d rate / d eta1 , d shape / d eta1 )
                      cbind(-theta, 
                            matrix(0, nrow(theta), ncol(theta)))
                        
                    } else {
                          
                      # c(d rate / d eta2 , d shape / d eta2)
                      cbind(matrix(0, nrow(theta), ncol(theta)),
                            matrix(0, nrow(theta), ncol(theta))) 
                    },
                        
                    #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                    if (wrt.param == 1) {
                      
                      #  c(d2rate / deta12 , d2shape / deta12 )
                      cbind(theta, matrix(0, nrow = nrow(theta),
                                            ncol = ncol(theta)))
                          
                    } else {
                          
                    # c(d2rate / deta2 , d2shape / deta22)
                cbind(matrix(0, nrow = nrow(theta), ncol = ncol(theta)), 
                      matrix(0, nrow = nrow(theta), ncol = ncol(theta)))
                })
    eta.ret <- theta.ret
  } else {
    
    eta.ret <- switch(deriv + 1,
                      log( shape.mat / theta),
                      #  wrt = 1 -> wrt rate (=theta = beta)
                      #  wrt = 2 -> wrt shape (= alpha)
                      
                      # d eta1 / dalpha (d eta2 / dalpha = 0 not returned)
                      if (wrt.param == 1) -1 / theta else 1 / shape.mat,
                      # ditto
                      if (wrt.param == 1) 1 / theta^2 else -1/shape.mat^2)
  }
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
}   # END OF R-function.





# 01012021 weibullRMeanlink, renamed as weibullMTTFlink
# renamed as weibmeanlink then weibullMlink
weibullMlink <-  function(theta, shape = NULL, wrt.param = NULL,
                          bvalue = NULL, inverse = FALSE,
                          deriv = 0, short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("weibullMlink(",theta,"; shape)", sep = "") else
        paste("log(",as.char.expression(theta)," * gamma(1 + 1/shape))")
    
    if (tag) 
      string <- 
      paste("2-parameter Weibull Mean link:",  
            string)
    return(string)
    
  }
  
  if (!length(shape))
    stop("Enter a valid value for 'shape'.")
  
  theta <- cbind(theta)
  mydim <-dim(theta)
  shape <- cbind(shape)

  
  if (nrow(shape) != nrow(theta)) {
    vec2com <- c(nrow(shape), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      shape <- matrix(shape, nrow = max(vec2com), ncol = ncol(shape))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(theta))
    }
    warning("Arguments 'theta' and 'shape' differ in dimension.")
  }
  
  if ( length(wrt.param) && (!(wrt.param %in% 1:2)  ||
                             !is.Numeric(wrt.param, length.arg = 1)) )
    stop("Bad input for argument 'wrt.param'.")
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv + 1,
                        # theta
                        exp(theta - lgamma(1 + 1/shape)),
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(dscale / d eta1 , dshape/ d eta1 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(dscale / d eta2 , dshape / d eta2)
                          cbind( matrix(0, nrow(theta), ncol(theta)), 
                                 matrix(0, nrow(theta), ncol(theta)) )
                        },
                        
                        if (wrt.param == 1) {
                          #  c(d2scale / deta12 , d2shape / deta12 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          help <- matrix(0, nrow(theta), ncol(theta))
                          cbind(help, help)
                        })
    
    eta.ret <- theta.ret
  } else {
    
    eta.ret <- switch(deriv + 1,
                      log(theta) + lgamma(1 + 1/shape),
                      #  wrt = 1 -> wrt scale (= theta)
                      #  wrt = 2 -> wrt shape
                      
                      # d eta1 / dscale (d eta2 / dshape = 0 not returned)
                      if (wrt.param == 1) 1 / theta else
                        -digamma(1 + 1/shape) / shape^2,
                      # ditto
                      if (wrt.param == 1) -1 / theta^2 else
                        -trigamma(1 + 1/shape) / shape^2 + 
                        2 * digamma(1 + 1/shape)/shape^3)
  }
  
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
  
}





## 01022021 Weibull Reliable Life
weibullRLifelink <-  function(theta, shape = NULL, 
                              Reliab = stop("Enter a value for 'Reliab'."),
                                wrt.param = NULL,
                                bvalue = NULL, inverse = FALSE,
                                deriv = 0, short = TRUE, tag = FALSE) {
  if (is.character(theta)) {
    string <- if (short)
  paste("weibullRLifelink(",theta,"; ", Reliab, ", shape)", sep = "") else
        paste("log[",
              as.char.expression(theta),"* [-log ('Reliab')]^(1/shape) ]")
    
    if (tag) 
      string <- paste("2-parameter Weibull Reliable Life link:", string)
    return(string)
    
  }
  
  if (!length(Reliab))
    stop("Enter a valid value for 'Reliab', between 0 and 1.")
  
  if (length(Reliab) && (log(Reliab) >= 0))
    stop("Wrong input for 'Reliab. Must be beetween 0 and 1.")
  
  if (!length(shape))
    stop("Enter a valid value for 'shape'.")
  
  theta <- cbind(theta)
  mydim <-dim(theta)
  shape <- cbind(shape)
  
  if (ncol(shape) != ncol(theta))
    stop("Unequal number of cols for 'theta' and 'shape'.")
  
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
  
  if ( length(wrt.param) && (!(wrt.param %in% 1:2)  ||
                             !is.Numeric(wrt.param, length.arg = 1)) )
    stop("Bad input for argument wrt.param.")
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv + 1,
                        # theta
                        exp(theta) * (-log(Reliab))^(-1/shape) ,
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(dscale / d eta1 , dshape/ d eta1 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(dscale / d eta2 , dshape / d eta2)
                          # dshape / d eta2 (not required), computed by the
                          #  other link function.
                          #help <- -theta * shape * digamma(1 + 1/shape)
                          cbind( matrix(0, nrow(theta), ncol(theta)), 
                                 matrix(0, nrow(theta), ncol(theta)) )
                        },
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(d2scale / deta12 , d2shape / deta12 )
                          cbind(theta, matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(d2scale / deta2 , d2shape / deta22)
                          # d2shape / deta2 (not required). Computed by the
                          #  other link function.
                          # PENDING
                          help <- matrix(0, nrow(theta), ncol(theta))
                          #cbind(theta, NULL)
                          cbind(help, NULL)
                        })
    
    eta.ret <- theta.ret
  } else {
    
    ## Here, 1. eta
    ##       2. d eta / d thetaj, j = wrt = 1, 2 (scale or shape)
    ##       3. d2eta / dthetaj2, j = wrt = 1, 2.
    eta.ret <- switch(deriv + 1,
                      log(theta) + log(-log(Reliab))/shape,
                      #  wrt = 1 -> wrt scale (= theta)
                      #  wrt = 2 -> wrt shape
                      
                      # d eta1 / dscale (d eta2 / dshape = 0 not returned)
                      # else d eta1 / dbeta (d eta2/d beta... depends on
                      # the link. Default is  'loglink')
                      if (wrt.param == 1) 1 / theta else 
                        -log(-log(Reliab)) / shape^2,
                      # ditto
                      if (wrt.param == 1) -1 / theta^2 else 
                        2 * log(-log(Reliab)) / shape^3)
  }
  
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
  
}





# 01012021 The Weibull Failure Rate function -
weibullFRFlink <-  function(theta, shape = NULL, wrt.param = NULL,
                            ycons = NULL, bvalue = NULL, inverse = FALSE,
                             deriv = 0, short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("weibullFRFlink(",theta,"; shape)", sep = "") else
        paste("log [shape * ycons^(shape - 1) *",as.char.expression(theta),
              "^(-shape)]")
    
    if (tag) 
      string <-
        paste("2-parameter Weibull Failure Rate function (FRF) link:", 
              string)
    return(string)
    
  }
  
  if (!length(shape))
    stop("Enter a valid value for 'shape'.")
  
  theta <- cbind(theta)
  mydim <-dim(theta)
  shape <- cbind(shape)
  
  
  if (nrow(shape) != nrow(theta)) {
    vec2com <- c(nrow(shape), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      shape <- matrix(shape, nrow = max(vec2com), ncol = ncol(shape))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(theta))
    }
    warning("Arguments 'theta' and 'shape' differ in dimension.")
  }
  
  if ( length(wrt.param) && (!(wrt.param %in% 1:2)  ||
                             !is.Numeric(wrt.param, length.arg = 1)) )
    stop("Bad input for argument 'wrt.param'.")
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv + 1,
                        # theta
                        #exp(theta - lgamma(1 + 1/shape)),
                        ((shape * (ycons)^(shape - 1))^(1/shape)) *
                          exp(-theta/shape),
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          #  c(dscale / d eta1 , dshape/ d eta1 )
                          cbind(-theta/shape, 
                                matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          # c(dscale / d eta2 , dshape / d eta2)
                          cbind( matrix(0, nrow(theta), ncol(theta)), 
                                 matrix(0, nrow(theta), ncol(theta)) )
                        },
                        
                        if (wrt.param == 1) {
                          #  c(d2scale / deta12 , d2shape / deta12 )
                          cbind(theta / shape^2, 
                                matrix(0, nrow(theta), ncol(theta)))
                          
                        } else {
                          help <- matrix(0, nrow(theta), ncol(theta))
                          cbind(help, help)
                        })
    
    eta.ret <- theta.ret
  } else {
 
    eta.ret <- switch(deriv + 1,
                      #log(theta) + lgamma(1 + 1/shape),
                      log(shape) + (shape - 1) * log(ycons) - 
                        shape * log(theta),
                      #  wrt = 1 -> wrt scale (= theta)
                      #  wrt = 2 -> wrt shape
                      
                      # d eta1 / dscale (d eta2 / dshape = 0 not returned)
                      if (wrt.param == 1) -shape / theta else
                        1/ shape + log (ycons/theta),
                      # ditto
                     if (wrt.param == 1) shape / theta^2 else -1 / shape^2)
  }
  
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
  
}






## --  Yet to do (below)




# 20170221... invweibullMeanlink defined as log(scale / gamma(1 - 1/shape))
invweibullMeanlink <- function(theta, shape = NULL, wrt.param = NULL,
                               bvalue = NULL, inverse = FALSE,
                               deriv = 0, short = TRUE, tag = FALSE,
                               extra = list(dkdeta2 = NULL,
                                            d2kdeta22 = NULL)) {
  
  if (is.character(theta)) {
    ret.cha <- if (short) 
      paste("invweibullMeanlink(", as.character(theta),
            "; shape)", sep = "") else 
              paste("log(", as.character(theta),
                    " * Gamma(1 - 1/shape))", sep = "")
    
    if (tag)
      ret.cha <- paste("2-parameter Inverse Weibull mean link: ", ret.cha)
    return(ret.cha)
  }
  
  mydim <- dim(theta)
  shape.mat <- cbind(shape)
  theta <- cbind(theta)
  
  if (ncol(shape.mat) != ncol(theta)) {
    stop("Arguments 'theta' and 'shape' do not have ",
         "an equal number of cols")
  }
  
  if (nrow(shape.mat) != nrow(theta)) {
    vec2com <- c(nrow(shape.mat), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      shape.mat <- matrix(shape.mat, nrow = max(vec2com), ncol = ncol(theta))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(shape.mat))
    }
    warning("Arguments 'theta' and 'shape' do not have ",
            "an equal number of rows")
  }
  
  if (length(wrt.param) && (!(wrt.param %in% 1:2) ||
                            !is.Numeric(wrt.param, length.arg = 1)) )
    stop("argument 'wrt.param' should be 1 or 2")
  
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv + 1,
                        # theta (scale)
                        exp(theta) / gamma(1 - 1/shape),
                        
                        # der = 1
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d scale / d eta1 , d shape / d eta1 )
                          cbind(theta, matrix(0, nrow = nrow(theta),
                                              ncol = ncol(theta)))
                          
                        } else {
                          
                          # c(d scale / d eta2 , d shape / d eta2)
                          # d shape / d eta2 (not required) computed by the
                          #  other link function.
                          extra   <- extra$dkdeta2
                          to.ret  <- -theta * digamma(1 - 1/shape) * 
                            extra/shape^2
                          
                          cbind(to.ret, NULL)
                        },
                        
                        # der = 2
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d2shape / deta12 , d2rate / deta12 )
                          cbind(theta, matrix(0, nrow = nrow(theta),
                                              ncol = ncol(theta)))
                          
                        } else {
                          
                          extra   <- extra$d2kdeta22
                          to.ret  <- -theta * digamma(1 - 1/shape) * 
                            extra/shape^2
                          
                          # c(d2shape / deta2 , d2rate / deta22)
                          # d2rate / deta2 (not required) computed by the
                          #  other link function.
                          cbind(to.ret, NULL)
                        })
    eta.ret <- theta.ret
  } else {
    
    ##-- FROM HERE PENDING>>>>
    
    ## Here, 1. eta
    ##       2. d eta / d thetaj, j = wrt = 1, 2 (scale or shape)
    ##       3. d2eta / dthetaj2, j = wrt = 1, 2.
    eta.ret <- switch(deriv + 1,
                      log(theta * gamma(1 - 1 / shape)),
                      #  wrt = 1 -> wrt scale (= lambda = l)
                      #  wrt = 2 -> wrt shape (= k)
                      
                      # d eta1 / dscale (d eta2 / dscale = 0 not returned)
                      # else d eta1 / dshape (d eta2/d shape use loge or so)
                      if (wrt.param == 1) 1 / theta else 
                        - digamma(1 - 1/shape) / shape^2,
                      # ditto
                      if (wrt.param == 1) -1 / theta^2 else 1 / shape^2)
  }
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
}   # END OF R-function.
