##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.



gamma1Qlink <- function(theta,
                        p = stop("Argument 'p' must be specified."),
                        bvalue = NULL, inverse = FALSE, 
                        deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Wrong input for argument 'p'. Must lie between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    g.string <- 
      if (short) paste("gamma1Qlink(", theta, "; ", p, ")", sep = "") else
        paste("log qgamma(p = ", p,  ", shape = ", theta, ")", sep = "")
    
    if (tag)
      g.string <- paste("1-parameter Gamma quantile link: ", 
                        g.string, sep = "")
    return(g.string)
  }
  
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }
  
  my.help <- if (!is.null(dim(theta))) dim(theta) else NULL
  etas.tm <- theta
  dim(theta) <- dim(p) <- NULL
  
  if (deriv) {
    
    der1.qg <- (qgamma(p = p, shape = theta + 1e-4) - 
                  qgamma(p = p, shape = theta - 1e-4)) / (2*1e-4)
    
    der2.qg <- (qgamma(p = p, shape = theta + 1e-4) - 
                  2 * qgamma(p = p, shape = theta) + 
                  qgamma(p = p, shape = theta - 1e-4)) / (1e-4^2)
    
    d2Eta.d2t <- if (deriv == 2)  (qgamma(p = p, shape = theta) * 
             der2.qg - der1.qg^2 ) / qgamma(p = p, shape = theta)^2 else 0
  }
  
  if (inverse) {
    
    if (!deriv) {
      
      my.dim  <- if (!is.null(dim(etas.tm))) dim(etas.tm) else NULL
      dim(etas.tm) <- NULL
      etas.tm[etas.tm > 9] <- Inf; gamma.ret <- 0
      to.ret  <- numeric(length(etas.tm))
      
      ets.na  <- which(is.na(etas.tm) & !is.nan(etas.tm))
      ets.nan <- which(is.nan(etas.tm))
      big.vec <- c(ets.na, ets.nan)
      etas.tm <- if (!length(big.vec)) etas.tm else etas.tm[-big.vec] 
      
      if (length(etas.tm)) {
        a <- rep(1e-1, length(etas.tm))  
        b <- rep(1e5,  length(etas.tm))
        gamma.ret <-  newtonRaphson.basic(f = function(x, p, eta) {
          log(qgamma(p = p, shape = x)) - sign(eta) * abs(eta)  },
          fprime = function(x, p, eta) {
            der1.qg <- (qgamma(p = p, shape = x + 1e-4) - 
                          qgamma(p = p, shape = x - 1e-4)) / (2*1e-4)
            der1.qg / qgamma(p = p, shape = x) }, 
          a = a, b = b, n.Seq = 2e2, nmax = 30,
          p = p, eta = etas.tm)
      }
    }
    
    gamma.ret <- switch(deriv + 1,  gamma.ret, 
                   qgamma(p = p, shape = theta) / der1.qg,
                  -(qgamma(p = p, shape = theta) / der1.qg)^3 * d2Eta.d2t)
    
    if (!deriv) {
      to.ret[ets.na ] <- if (length(ets.na ))  NA else NULL
      to.ret[ets.nan] <- if (length(ets.nan)) NaN else NULL
      
      if (length(to.ret[-big.vec])) {
        to.ret[-big.vec] <- gamma.ret 
      } else {
        to.ret[to.ret == 0] <- gamma.ret
      }
      
      dim(to.ret) <- my.dim
      return(to.ret)
      
    } else {
      dim(gamma.ret) <- my.help
      return(gamma.ret)
    }
    
  } else {
    my.ret <- switch(deriv + 1,
                     log(qgamma(p = p, shape = theta)),
                     der1.qg / qgamma(p = p, shape = theta),
                     (qgamma(p = p, shape = theta) * der2.qg - der1.qg^2) / 
                       qgamma(p = p, shape = theta)^2)
    dim(my.ret) <- my.help
    return(my.ret)
  }
}


