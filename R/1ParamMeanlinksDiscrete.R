##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
# 20170510. Mean link functions for DISCRETE distributions.

# 20161026.
borel.tannerMeanlink <- function(theta, Qsize = 1, bvalue = NULL, 
                                 inverse = FALSE, deriv = 0, short = TRUE,
                                 tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("borel.tannerlink(", theta, ")", sep = "") else
        paste("-log(Q^(-1) - ", 
              as.char.expression(theta),
              " * Q^(-1))", sep = "")
    if (tag)
      string <- paste("Borel-Tanner distribution mean link:", string)
    return(string)
  }
  
  if (!is.Numeric(x = Qsize, integer.valued = TRUE, positive = TRUE))
    stop("Wrong input for argument 'Qsize'.")
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 4)
    stop("Argument 'deriv' must be 0, 1, or 2.")
  
  if (length(bvalue)) {
    
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- if (!inverse || (inverse && deriv != 0)) 
      1 - bvalue else theta[theta >= 1.0]
  } else {
    
    theta[theta <= 0.0] <- NaN 
    
    if (!inverse || (inverse && deriv != 0))
      theta[theta >= 1] <- NaN
  }
  
  if (inverse) {
    
    switch(deriv + 1, 
           1 - Qsize * exp(-theta),
           exp(log1p(-theta)), #(theta -1)/theta
           -exp(log1p(-theta)),
           stop("argument 'deriv' unmatched"))
    
  } else {
    
    switch(deriv + 1,
           log(Qsize) - log1p(-theta),
           exp(-log1p(-theta)),
           exp(-2 * log1p(-theta)),
           stop("argument 'deriv' unmatched"))
  }
}






# 20161228
geometricffMeanlink <- function(theta, bvalue = NULL,
                                inverse = FALSE, deriv = 0, 
                                short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("geometricfflink(",  theta, ")", sep = "") else
        paste("-logit(", as.char.expression(theta),")", sep = "")
    if (tag)
      string <- paste("geometric distribution mean link:",
                      string, sep = "")
    return(string)
  }
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 3)
    stop("Invalid value for argument 'deriv'. Must be 0, 1 or 2. ")
  
  if (length(bvalue)) {
    theta[theta <= 0.0] <- bvalue
    theta[theta >= 1.0] <- bvalue 
  } else {
    if (!inverse) {
      theta[theta <= 0.0] <- NaN
      theta[theta >= 1.0] <- NaN 
    }
  }
  
  if (inverse) {
    switch(deriv + 1,
           logitlink(theta = -theta, inverse = inverse),
           -(theta - theta^2),
           -(1 - 2 * theta) * theta * (1 - theta),
           stop("argument 'deriv' unmatched"))
  } else {
    switch(deriv + 1,
           -logitlink(theta = theta),
           -logitlink(theta = theta, deriv = deriv),
           -logitlink(theta = theta, deriv = deriv),
           stop("argument 'deriv' unmatched")) 
  }
}






# 20160909.
logffMeanlink <- function(theta, bvalue = NULL, 
                          alg.roots = c("Newton-Raphson", "bisection")[1],
                          inverse = FALSE, deriv = 0, 
                          short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
    paste("logfflink(", theta, ")", sep = "") else
      paste("logitlink(",  theta, ") - clogloglink(", theta, ")", sep = "")
    if (tag)
      string <- paste("Logarithmic distribution mean link:", string)
    return(string)
  }
  
  alg.roots <- match.arg(alg.roots,
                         c("Newton-Raphson", "bisection"))[1]
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 3)
    stop("Argument 'deriv' must be 0, 1, or 2.")
  
  if (length(bvalue)) {
    
    theta[theta <= 0] <- bvalue
    theta[theta >= 1] <- if (deriv != 0 || (deriv == 0 & !inverse)) 
      1 - bvalue else theta[theta >= 1]
  } else{
    
    theta[theta <= 0] <- NaN
    
    if (deriv != 0 || (deriv == 0 & !inverse) )
      theta[theta >= 1] <- NaN
  }
  
  if (inverse) {
    
    sPrime <- if (deriv != 0) (logitlink(theta, deriv = 1) - 
                            clogloglink(theta, deriv = 1))^(-1) else NULL
    switch(deriv + 1,
           logfflink.inv.deriv0(etas = theta, 
                                alg.roots = alg.roots),
           sPrime,
           (-1) * ( sPrime^3 ) * (logitlink(theta = theta, deriv = 2) - 
                                    clogloglink(theta = theta, deriv = 2)),
           stop("Argument 'deriv' unmatched"))
    
  } else {
    return(logitlink(theta = theta, inverse = FALSE, deriv = deriv) - 
             clogloglink(theta = theta, inverse = FALSE, deriv = deriv))
  }
}





logfflink.inv.deriv0 <- function(etas, eps = 1e-8, 
                                 alg.roots = c("Newton-Raphson", 
                                               "bisection")[1]) {
  
  my.dim    <- if (!is.null(dim(etas))) dim(etas) else NULL
  dim(etas) <- NULL 
  to.ret    <- numeric(length(etas))
  inv.0.ret <- 0  # length(inv.0.r) always TRUE below. Avoid using NULL
  
  alg.roots <- match.arg(alg.roots, 
                         c("Newton-Raphson", "bisection"))[1]
  
  etas[etas > 15] <- exp(log(15))
  ets.na  <- which(is.na(etas) & !is.nan(etas))
  ets.nan <- which(is.nan(etas))
  big.vec <- c(ets.na, ets.nan)
  etas <- if (!length(big.vec)) etas else etas[-big.vec] 
  
  if (length(etas)) 
    if(alg.roots == "Newton-Raphson") {
      
      inv.0.ret <- 
        newtonRaphson.basic(f = function(x, etas) {
          log(x /( (-1) * (1 - x) * log1p(-x) )) - etas
        },
        fprime = function(x, etas) {
          logitlink(theta = x, inverse = FALSE, deriv = 1) - 
            clogloglink(theta = x, inverse = FALSE, deriv = 1)
        },
        a = rep(0 + eps, length(etas)), 
        b = rep(1 - eps, length(etas)), 
        tol = eps, etas = etas)
      
    } else {
      inv.0.ret <- bisection.basic(f = function(x, etas) {
        log(x /( (-1) * (1 - x) * log1p(-x) )) - etas
      }, 
      a = rep(0 + eps, length(etas)), 
      b = rep(1 - eps, length(etas)), 
      etas = etas)
    }
  
  to.ret[ets.nan] <- if (length(ets.nan)) NaN else NULL
  to.ret[ets.na ] <- if (length(ets.na ))  NA else NULL
  
  if (length(inv.0.ret) && length(big.vec)) {
    to.ret[-big.vec] <- inv.0.ret
  } else {
    to.ret <- inv.0.ret
  }
  
  dim(to.ret) <- my.dim
  
  to.ret
}








# 20161026.
posPoiMeanlink <- function(theta, bvalue = NULL,
                           alg.roots = c("Newton-Raphson", "bisection")[1],
                           inverse = FALSE, deriv = 0, 
                           short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("pospoissonlink(", theta, ")", sep = "") else
        paste("-log(", as.char.expression(theta),"^(-1)", 
              " - ", as.char.expression(theta),"^(-1) * exp(-",
              as.char.expression(theta),") )", sep = "")
    if (tag)
      string <- paste("Positive Poisson distribution mean link:", string)
    return(string)
  }
  
  alg.roots <- match.arg(alg.roots, c("Newton-Raphson", "bisection"))[1]
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 4)
    stop("Invalid value for argument 'deriv'. Must be 0, 1, or 2.")
  
  if (length(bvalue)) {
    theta[theta <= 0.0] <- 0.0 + bvalue
  } else {
    theta[theta <= 0.0] <- NaN
  }
  
  if (inverse) {
    
    if (deriv == 2) {
      der.the.et  <- ( exp(-log(theta)) - 1 / expm1(theta) )^(-1)
      der2.et.the <- exp(theta)/expm1(theta)^2 - exp(-2 * log(theta))
      deriv.2invT <- ( (-1) * der.the.et^3 * der2.et.the)
      rm(der.the.et, der2.et.the)
    }
    
    switch(deriv + 1,
           pospoilink.inv.deriv0(etas = theta, alg.roots = alg.roots),
           ( exp(-log(theta)) - 1 / expm1(theta) )^(-1),
           deriv.2invT,
           stop("argument 'deriv' unmatched"))
    
  } else {
    
    switch(deriv + 1,
           log(theta) - log1p(-exp(-theta)),
           exp(-log(theta)) - 1 / expm1(theta),
           exp(theta)/expm1(theta)^2 - exp(-2 * log(theta)),
           stop("argument 'deriv' unmatched"))
  }
}




pospoilink.inv.deriv0 <- function(etas, eps = 1e-8, 
                                  alg.roots = c("Newton-Raphson", 
                                                "bisection")[1]) {
  
  my.dim    <- if (!is.null(dim(etas))) dim(etas) else NULL
  dim(etas) <- inv.0.ret <- NULL 
  to.ret    <- numeric(length(etas))
  inv.ret   <- 
    
    alg.roots <- match.arg(alg.roots, c("Newton-Raphson", "bisection"))[1]
  
  to.igno <- which(is.infinite(exp(etas)))
  ets.na  <- which(is.na(etas) & !is.nan(etas))
  ets.nan <- which(is.nan(etas))
  big.vec <- c(to.igno, ets.na, ets.nan)
  etas <- if (!length(big.vec)) etas else etas[-big.vec] 
  
  if (length(etas)) {
    
    a.eps <- rep(0 + eps, length(etas))
    b.eps <- rep(10 - eps, length(etas))
    
    if(any(etas > 2)) {
      a.eps[etas > 2.0] <- exp(etas[etas > 2.0]) - 2e0
      b.eps[etas > 2.0] <- exp(etas[etas > 2.0]) + 2e0
    }
    
    if (alg.roots == "Newton-Raphson") {
      inv.0.ret <- newtonRaphson.basic(f = function(x, etas = 0) {
        log(x) - log1p(-exp(-x)) - etas },
        fprime = function(x, etas) {
          exp(-log(x)) - 1 / (exp(x) - 1) },
        a = a.eps, b = b.eps, tol = eps,
        n.Seq = 10, etas = etas, nmax = 15) 
    } else {
      inv.0.ret <- bisection.basic(f = function(x, etas =0) {
        log(x) - log1p(-exp(-x)) - etas },
        a = a.eps, b = b.eps, etas = etas)
    }
  } 
  
  to.ret[to.igno] <- if (length(to.igno)) Inf else NULL
  to.ret[ets.nan] <- if (length(ets.nan)) NaN else NULL
  to.ret[ets.na ] <- if (length(ets.na ))  NA else NULL
  
  if (length(big.vec)) {
    to.ret[-big.vec] <- inv.0.ret  
  } else {
    to.ret <- inv.0.ret
  }
  
  dim(to.ret) <- my.dim
  to.ret
}







# 20161002
yulesimonMeanlink <- function(theta, bvalue = NULL,
                              inverse = FALSE, deriv = 0, 
                              short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("yulesimonlink(",  theta, ")", sep = "") else
        paste("-log(1 - ", as.char.expression(theta),"^(-1))", sep = "")
    if (tag)
      string <- paste("Yule-Simon distribution mean link:", string)
    return(string)
  }
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 3)
    stop("Invalid value for argument 'deriv'. Must be 0, 1 or 2. ")
  
  if (length(bvalue)) {
    
    theta[theta <= 0.0] <- bvalue
    
    if (!inverse && (deriv == 0))  # Assure r > 1 when inv = F, d = 0.
      theta[theta >= 0 & theta <= 1.0] <- 1.0 + bvalue
    
  } else {
    
    theta[theta <= 0.0] <- NaN
    
    if (!inverse && (deriv == 0))  # Assure r > 1 when inv = F, d = 0.
      theta[theta >= 0 & theta <= 1.0] <- NaN
  }
  
  if (inverse) {
    my.der <- switch(deriv + 1,
                     exp(theta) / expm1(theta),
                     theta * (1 - theta),
                     theta * (theta - 1) * (2 * theta - 1),
                     stop("argument 'deriv' unmatched"))
    
    return(my.der)
    
  } else {
    
    pre.the <- switch(deriv + 1,
                      -log1p(-theta^(-1)),
                      ( theta * (1 - theta) )^(-1),
                      (2 * theta - 1)/(theta * (1 - theta))^2,
                      stop("argument 'deriv' unmatched")) 
    return(pre.the)
  }
}







# 20160923.
zetaffMeanlink <- function(theta, bvalue = NULL, 
                           alg.roots = c("Newton-Raphson", "bisection")[1],
                           inverse = FALSE, deriv = 0, 
                           short = TRUE, tag = FALSE) {
  
  if (is.character(theta)) {
    string <- if (short)
      paste("zetafflink(", theta, ")", sep = "") else
        paste("log(zeta(", theta,") / zeta(", 
              as.char.expression(theta), " + 1))", sep = "")
    if (tag)
      string <- paste("Zeta distribution mean link:", string)
    return(string)
  }
  
  alg.roots <- match.arg(alg.roots, 
                         c("Newton-Raphson", "bisection"))[1]
  
  if (!is.Numeric(x = deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv >= 3)
    stop("Invalid value for argument 'deriv'. Must be 0, 1, or 2.")
  
  if (length(bvalue)) {
    
    theta[theta <= 0.0] <- bvalue
    
    if (!inverse && (deriv == 0))  # Assure p > 1 when inv = F, d = 0.
      theta[theta >= 0 & theta <= 1.0] <- 1.0 + bvalue
    
  } else {
    
    theta[theta <= 0.0] <- NaN
    
    if (!inverse && (deriv == 0)) # Assure p > 1 when inv = F, d = 0.
      theta[theta >= 0 & theta <= 1.0] <- NaN
  }
  
  # remove NaN, < 0 and NAs before calling 'zeta()', it returns ERROR.
  my.dim  <- if (!is.null(dim(theta))) dim(theta) else NULL
  dim(theta) <- NULL; my.ret  <- 0
  to.ret  <- numeric(length(theta))
  ets.na  <- which(is.na(theta) & !is.nan(theta))
  ets.nan <- which(is.nan(theta))
  big.vec <- c(ets.na, ets.nan)
  theta   <- if (!length(big.vec)) theta else theta[-big.vec] 
  
  if (length(theta)) {
    
    off.theta <- theta + 1 # In concordance with VGAM.
    
    if (deriv != 0) {
      pEta.dthe <- (zeta(x = theta, deriv = 1)/zeta(x = theta) - 
                      zeta(x = off.theta, deriv = 1)/zeta(x = off.theta))
      der1.ratio <- (zeta(x = theta) * zeta(x = theta, deriv = 2)- 
                       zeta(x = theta, deriv = 1)^2) / zeta(x = theta)^2
      
      der2.ratio <- (zeta(x = off.theta) * zeta(x = off.theta, deriv = 2) -
                       zeta(x = off.theta, deriv = 1)^2) / zeta(x = off.theta)^2
    }
    
    if (inverse) {
      if (deriv == 0)
        theta[theta > 17] <- exp(log(17))
      my.ret <- switch(deriv + 1,   
                       zetafflink.inv.deriv0(etas = theta, 
                                             alg.roots = alg.roots), 
                       pEta.dthe^(-1), 
                       (-1) * (pEta.dthe^(-3)) * (der1.ratio - der2.ratio))
      
    } else {
      if (deriv == 0) {
        theta[theta > 1e10] <- exp(log(1e10))
        off.theta[theta > 1e10] <- exp(log(1e10))
      }
      my.ret <- switch (deriv + 1,
                        log( zeta(x = theta) / zeta(x = off.theta) ),
                        pEta.dthe,
                        der1.ratio - der2.ratio)
    }
  }
  
  to.ret[ets.nan] <- if (length(ets.nan)) NaN else NULL
  to.ret[ets.na ] <- if (length(ets.na ))  NA else NULL
  
  if (length(my.ret) && length(big.vec)) {
    to.ret[-big.vec] <- my.ret
  } else {
    to.ret <- my.ret
  }
  dim(to.ret) <- my.dim
  to.ret
}




zetafflink.inv.deriv0 <- function(etas, eps = 1e-8,
                                  alg.roots = c("Newton-Raphson",
                                                "bisection")[1]) {
  
  alg.roots <- match.arg(alg.roots,
                         c("Newton-Raphson", "bisection"))[1]
  
  if (any(etas < 1e-9)) 
    alg.roots <- "bisection"
  
  if (alg.roots == "Newton-Raphson") {
    newtonRaphson.basic(f = function(x, etas){
      log(zeta(x)/zeta(x + 1)) - etas
    }, 
    fprime = function(x, etas){
      zeta(x = x, deriv = 1)/zeta(x = x) - 
        zeta(x = exp(log1p(x)), deriv = 1)/ 
        zeta(x = exp(log1p(x))) 
    },
    a = rep(1 + eps, length(etas)),
    b = rep(60, length(etas)), nmax = 20,
    n.Seq = 59, tol = eps, etas = etas)
  } else {
    bisection.basic(f = function(x, etas){
      log(zeta(x)/zeta(x + 1)) - etas
    }, 
    a = rep(1 + eps, length(etas)), 
    b = rep(60 - eps, length(etas)), 
    etas = etas, nmax = 150)
  }
}
