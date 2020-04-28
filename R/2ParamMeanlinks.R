##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


# 20170221... GammaRMeanlink defined as log(shape / rate)
gammaRMeanlink <- function(theta, rate = NULL, wrt.param = NULL,
                           bvalue = NULL, inverse = FALSE,
                           deriv = 0, short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("gammaRMeanlink(", theta, "; rate)", sep = "") else
        paste("log (", as.char.expression(theta), " / rate)", sep = "")
    
    if (tag)
      string <- paste("2-parameter Inverse Gamma mean link: ", string)
    return(string)
  }
  
  mydim <- dim(theta)
  rate.mat <- if (length(rate)) cbind(rate) else NULL
  theta <- cbind(theta)
  
  if (ncol(rate.mat) != ncol(theta)) {
    stop("Arguments 'theta' and 'rate' do not have ",
         "an equal number of cols")
  }
  
  if (nrow(rate.mat) != nrow(theta)) {
    vec2com <- c(nrow(rate.mat), nrow(theta))
    mymax <- which(vec2com == max(vec2com))[1]
    if (mymax == 2) {
      rate.mat <- matrix(rate.mat, nrow = max(vec2com), ncol = ncol(theta))
    } else{
      theta <- matrix(theta, nrow = max(vec2com), ncol = ncol(rate.mat))
    }
    warning("Arguments 'theta' and 'rate' do not have ",
            "an equal number of rows")
  }
  
  if (length(wrt.param) && (!(wrt.param %in% 1:2) ||
                            !is.Numeric(wrt.param, length.arg = 1)) )
    stop("argument 'wrt.param' should be 1 or 2")
  
  
  if (!inverse)
    theta[theta <= 0.0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv + 1,
                        # theta
                        exp(theta) * rate.mat,
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d shape / d eta1 , d rate / d eta1 )
                          cbind(theta, matrix(0, nrow = nrow(theta),
                                              ncol = ncol(theta)))
                          
                        } else {
                          
                          # c(d shape / d eta2 , d rate / d eta2)
                          # d rate / d eta2 (not required) computed by the
                          #  other link function.
                          cbind(theta, NULL)
                        },
                        
                        #  wrt = 1 -> wrt eta1 ELSE wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d2shape / deta12 , d2rate / deta12 )
                          cbind(theta, matrix(0, nrow = nrow(theta),
                                              ncol = ncol(theta)))
                          
                        } else {
                          
                          # c(d2shape / deta2 , d2rate / deta22)
                          # d2rate / deta2 (not required) computed by the
                          #  other link function.
                          cbind(theta, NULL)
                        })
    eta.ret <- theta.ret
  } else {
    
    ## Here, 1. eta
    ##       2. d eta / d thetaj, j = wrt = 1, 2 (shape or rate)
    ##       3. d2eta / dthetaj2, j = wrt = 1, 2.
    eta.ret <- switch(deriv + 1,
                      log(theta / rate.mat),
                      #  wrt = 1 -> wrt shape (= alpha = theta)
                      #  wrt = 2 -> wrt beta (rate)
                      
                      # d eta1 / dalpha (d eta2 / dalpha = 0 not returned)
                      # else d eta1 / dbeta (d eta2/d beta... depends on
                      # the link. Default is  'loglink')
                      if (wrt.param == 1) 1 / theta else -1 / rate.mat,
                      # ditto
                      if (wrt.param == 1) -1 / theta^2 else 1 / rate.mat^2)
  }
  
  if (is.matrix(eta.ret))
    colnames(eta.ret) <- NULL else
      names(eta.ret) <- NULL
  
  if (is.null(mydim) || (length(mydim) > 2))
    dim(eta.ret) <- mydim
  
  eta.ret
}   # END OF R-function.








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
