##########################################################################
# These functions are
# Copyright (C) 2014-2024 V. Miranda & T. Yee & V Liu
# Auckland University of Technology & University of Auckland
# All rights reserved.

dtrunclnorm <- function(x, meanlog = 0, sdlog = 1,
                          min.support = 0, max.support = Inf,
                          log = FALSE) {

  if (!is.logical(log))
    stop("Bad input for argument log.")
  
  if (log) {
    ans <- dnorm(log(x), meanlog, sdlog, log = TRUE) - log(x) -
             log(pnorm(log(max.support), meanlog, sdlog) - 
                 pnorm(log(min.support), meanlog, sdlog) )
    
  } else {
     ans <- (1 / x) * dnorm(log(x), meanlog, sdlog) / 
                (pnorm(log(max.support), meanlog, sdlog) - 
                         pnorm(log(min.support), meanlog, sdlog))
    
  }
  
  ans[x > max.support] <- 0
  ans[x < min.support] <- 0
  ans[min.support < 0] <- NaN
  ans[min.support >= max.support]  <- NaN
  ans[sdlog <= 0] <- NaN
  
  ans
}




ptrunclnorm <- function(q, meanlog = 0, sdlog = 1,
                           min.support = 0, max.support = Inf) {
  
  ans <- ( pnorm(log(q), meanlog, sdlog) - 
              pnorm(log(min.support), meanlog, sdlog) ) / 
                    ( pnorm(log(max.support), meanlog, sdlog) - 
                          pnorm(log(min.support), meanlog, sdlog) )
  
  ans[q >= max.support] <- 1
  ans[q <= min.support] <- 0
  ans[min.support < 0] <- NaN
  ans[min.support >= max.support]  <- NaN
  ans[sdlog <= 0] <- NaN
  
  ans
  
  
}



qtrunclnorm <- function(p, meanlog = 0, sdlog  = 1,
                        min.support = 0, max.support = Inf,
                        log.p = FALSE) {
  
  
  if (!is.Numeric(p, positive = TRUE)) 
    stop("bad input for argument 'p'")
  
  if (!is.logical(log.p))
    stop("Wrong inout for argument log.p")
  
  if (log.p)
    p <- exp(p)
  
  if (max(p) > 1) 
    stop("Argument 'p' or 'log(p)' must be in (0, 1)")
  
  ans <- ( pnorm(log(max.support), meanlog, sdlog) - 
             pnorm(log(min.support), meanlog, sdlog) ) * p +
    pnorm(log(min.support), meanlog, sdlog)
  
  
  ans[min.support >= max.support]  <- NaN
  ans[sdlog <= 0] <- NaN
  
  exp(qnorm(ans, meanlog, sdlog))
  
  
}



rtrunclnorm <- function(n, meanlog = 0, sdlog  = 1,
                        min.support = 0, max.support = Inf) {
  
  ### Inverse truncation limits - CDF
  nn <- n; rm(n)
  inv_LL <- plnorm(min.support, meanlog, sdlog)
  inv_UL <- plnorm(max.support, meanlog, sdlog)
  
  ## Unif in boundaries + inverse CDF
  uu <- runif(nn, inv_LL, inv_UL)
  
  ## Truncated-lognormal data
  ans <- qlnorm(uu, meanlog, sdlog)
  ans[min.support >= max.support]  <- NaN
  ans[sdlog <= 0] <- NaN
  
  ans
  
  
}


