###########################################################################
# These functions are
# Copyright (c) 2014-2018 V. Miranda and T. W. Yee, University of Auckland.
# All rights reserved.

dgen.betaII <- function(x, 
                        scale    = 1.0, 
                        shape1.a = 1.0, 
                        shape2.p = 1.0, 
                        shape3.q = 1.0, 
                        log      = FALSE)  {
  
  if (!length(shape2.p) || !length(shape3.q))
    stop("Argument 'shape2.p' or 'shape3.q' is missing, with no default.")
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("Bad input for argument 'log'.")
  b  <- scale ; a  <- shape1.a ; p  <- shape2.p ; q  <- shape3.q
  rm(scale, shape1.a, shape2.p, shape3.q, log)

  mydata <- log(a) + (a * p - 1) * 
                 log(ifelse(is.nan(x), NaN, 
                            ifelse(abs(x) == Inf, abs(1/x), 
                                   ifelse(x < 0, 1/Inf, x)))) -
             a * p * log(b) - lbeta(p, q) - 
             (p + q) * log1p(ifelse(is.nan(x), NaN, 
                                    (ifelse(x >= 0, x, 1/Inf)/b)^a))
  
  if (log.arg) mydata else exp(mydata)   
}

pgen.betaII <- function(q, 
                        scale    = 1.0, 
                        shape1.a = 1.0,  
                        shape2.p = 1.0, 
                        shape3.q = 1.0, 
                        lower.tail = TRUE, 
                        log.p    = FALSE) {
  
  if (!length(shape2.p) || !length(shape3.q))
    stop("Argument 'shape2.p' or 'shape3.q' is missing, with no default.")
  if (!is.logical(log.p) || length(log.p) != 1)  
    stop("Bad input for argument 'log.p'.")
  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("Bad input for argument 'lower.tail'.")
  b   <- scale ; a   <- shape1.a ; p   <- shape2.p ; qsh <- shape3.q
  rm(scale, shape1.a, shape2.p, shape3.q)
  i <- which(abs(q) != Inf, arr.ind = TRUE)
  q[i] <- ifelse(q[i] < 0, 0.0, (q[i]^a)/(q[i]^a + b^a))
  mydata <- pbeta(q = q, shape1 = p, shape2 = qsh, ncp = 0, 
                  lower.tail = lower.tail, log.p = log.p)
  
  mydata
}

qgen.betaII <-function(p, 
                       scale    = 1.0, 
                       shape1.a = 1.0,
                       shape2.p = 1.0, 
                       shape3.q = 1.0,
                       lower.tail = TRUE, 
                       log.p    = FALSE) {
  
  if (!length(shape2.p) || !length(shape3.q))
    stop("Argument 'shape2.p' or 'shape3.q' is missing, with no default.")
  if (!is.logical(log.p) || length(log.p) != 1)  
    stop("Bad input for argument 'log.p'.")
  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("Bad input for argument 'lower.tail'.")
  b   <- scale ; a   <- shape1.a ; psh <- shape2.p ; q   <- shape3.q
  rm(scale, shape1.a, shape2.p, shape3.q)
  if (log.p) {
    p <- qbeta(p = ifelse(is.nan(p), NaN,
                     ifelse(p <= 0.0, p, 0.0)), 
          shape1 = psh, shape2 = q, ncp = 0,
          lower.tail = lower.tail, log.p = log.p)
    
    b * (p / (1 - p))^(1/a)
  } else {
    p <- qbeta(p = ifelse(is.nan(p), NaN,
                     ifelse(p <= 0.0, 0.0,
                            ifelse(p >= 1.0, 1.0, p))),
          shape1 = psh, shape2 = q, ncp = 0,
          lower.tail = lower.tail, log.p = log.p)
    
    b * (p / (1 - p))^(1/a)
  }
}

rgen.betaII <- function(n, 
                        scale    = 1.0, 
                        shape1.a = 1.0, 
                        shape2.p = 1.0, 
                        shape3.q = 1.0) {
  
  if (!length(shape2.p) || !length(shape3.q))
    stop("Argument 'shape2.p' or 'shape3.q' is missing, with no default.")
  b   <- scale ; a   <- shape1.a ; p   <- shape2.p ; q   <- shape3.q
  rm(scale, shape1.a, shape2.p, shape3.q)
  mydata <- rbeta(n = n, shape1 = p, shape2 = q, ncp = 0)
  
  b * (mydata / (1 - mydata))^(1/a)
}

