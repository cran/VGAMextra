##########################################################################
# These functions are
# Copyright (C) 2014-2018 V. Miranda a& T. W. Yee, University of Auckland.
# All rights reserved.

dinvgamma <- function (x, 
                       scale = 1/rate, 
                       shape, 
                       rate  = 1, 
                       log   = FALSE) {
  
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("Bad input for argument 'log'") 
  
  rm(log) 
  i <- which(x < .Machine$double.eps)
  if (log.arg) {
    data <- dgamma(x = 1/x, shape = shape, scale = 1/scale, 
                   log = log.arg) - 2*log(ifelse(x > 0, x, Inf))
  } else {
    data <- dgamma(x = 1/x, shape = shape, scale = 1/scale, 
                   log = log.arg)/(x^2)
  }
  if(length(i)) 
    data[i] <- dgamma(x = 1/x, shape = ifelse(shape == 1, 2*shape, shape), 
                      scale = 1/scale, log = log.arg)[i]
  data
}


pinvgamma <- function (q, 
                       scale = 1/rate, 
                       shape, 
                       rate  = 1, 
                       lower.tail = TRUE, 
                       log.p = FALSE) {
  
     i    <- which(q < 0)
  data    <- pgamma(1/q, shape = shape, scale = 1/scale, 
                 lower.tail = !lower.tail, log.p = log.p)
  if(length(i))
  data[i] <- pgamma(1/q, shape = shape, scale = 1/scale,
                    lower.tail = lower.tail, log.p = log.p)[i]
  return(data)
  
}


qinvgamma <- function (p, 
                       scale = 1/rate, 
                       shape, 
                       rate  = 1,
                       lower.tail = TRUE,
                       log.p = FALSE) 
   1/qgamma(p = p, shape = shape, scale = 1/scale, 
            lower.tail = !lower.tail, log.p = log.p)


rinvgamma <- function (n, 
                       scale = 1/rate, 
                       shape, 
                       rate  = 1) 
  1/rgamma(n = n, shape = shape, scale = 1/scale)

