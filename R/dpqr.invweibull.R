##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

dinvweibull <- function (x, 
                         scale = 1, 
                         shape, 
                         log   = FALSE) {
  
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("Bad input for argument 'log'") 
  rm(log) 
  
  i  <- which(x < .Machine$double.eps)
  if (log.arg) {
    data <- dweibull(x = 1/x, shape = shape, scale = 1/scale, 
                     log = log.arg) - 2*log(ifelse(x > 0, x, Inf))
  } else {
    data <- dweibull(x = 1/x, shape = shape, scale = 1/scale, 
                     log = log.arg)/(x^2)
  }
  if (length(i))
    data[i] <- dweibull(x = 1/x, shape = ifelse(shape == 1.0, 2*shape, shape),
                        scale = 1/scale, log = log.arg)[i]
  data
}


pinvweibull <- function (q, 
                         scale  = 1, 
                         shape, 
                         lower.tail = TRUE, 
                         log.p  = FALSE) {
  
  i   <- which(q < 0)
  data    <- pweibull(1/q, shape  = shape, scale = 1/scale, 
                      lower.tail = !lower.tail, log.p = log.p)
  if(length(i))
  data[i] <- pweibull(1/q, shape = shape, scale = 1/scale,
                      lower.tail = lower.tail, log.p = log.p)[i]
  return(data)
  
}


qinvweibull <- function (p, 
                         scale  = 1, 
                         shape, 
                         lower.tail = TRUE, 
                         log.p  = FALSE) 
  1/qweibull(p = p, shape = shape, scale = 1/scale, 
             lower.tail = !lower.tail, log.p = log.p)


rinvweibull <- function (n, 
                         scale = 1, 
                         shape) 
  1/rweibull(n = n, shape = shape, scale = 1/scale)

