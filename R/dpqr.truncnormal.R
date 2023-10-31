##########################################################################
# These functions are
# Copyright (C) 2014-2024 V. Miranda & T. Yee & V Liu
# Auckland University of Technology & University of Auckland
# All rights reserved.

dtruncnorm <- function(x, mean = 0, sd  = 1, 
                         min.support = -Inf, max.support = Inf,
                         log = FALSE) {

  if (!is.logical(log))
    stop("Bad input for argument log.")
  
  if (log) {
    
    ans <- dnorm(x, mean, sd, log = TRUE) - 
                 log(pnorm((max.support - mean)/ sd) - 
                          pnorm((min.support - mean) / sd))
    ## Same as
    #  ans <- dnorm( (x - mean) / sd, log = TRUE) - 
    #             log(pnorm(max.support, mean, sd) -
    #                 pnorm(min.support, mean, sd)) - log(sd)
    
  } else {
    ans <- dnorm(x, mean, sd) / (pnorm(max.support, mean, sd) - 
                                        pnorm(min.support, mean, sd))
    
    ## Same as
    # ans <- (1 / sd) * dnorm( (x - mu) / sd ) / 
    #    (pnorm(max.support, mean, sd) - pnorm(min.support, mean, sd))
    
  }
  
  ans[x > max.support] <- 0
  ans[x < min.support] <- 0
  ans[min.support >= max.support]  <- NaN
  ans[sd <= 0] <- NaN
  
  ans
}




ptruncnorm <- function(q, mean = 0, sd  = 1,
                         min.support = -Inf, max.support = Inf) {

  ans <- ( pnorm(q, mean, sd) - pnorm(q = min.support, mean, sd) ) / 
                    ( pnorm(max.support, mean, sd) - 
                             pnorm(min.support, mean, sd) )
  
  ans[q >= max.support] <- 1
  ans[q <= min.support] <- 0
  ans[min.support >= max.support]  <- NaN
  ans[sd <= 0] <- NaN
  
  ans
  
}




qtruncnorm <- function(p, mean = 0, sd  = 1, 
                         min.support = -Inf, max.support = Inf,
                         log.p = FALSE) {
  
  if (!is.Numeric(p, positive = TRUE)) 
    stop("bad input for argument 'p'")
  
  if (!is.logical(log.p))
    stop("Wrong inout for argument log.p")
  
  if (log.p)
    p <- exp(p)
  
  if (max(p) > 1) 
    stop("Argument 'p' or 'log(p)' must be in (0, 1)")
  
  ans <- ( pnorm(max.support, mean, sd) - 
                   pnorm(min.support, mean, sd) ) * p +
                                 pnorm(min.support, mean, sd)
  
  
  ans[min.support >= max.support]  <- NaN
  ans[sd <= 0] <- NaN
  
  qnorm(ans, mean, sd)
}





rtruncnorm <- function(n, mean = 0, sd  = 1,
                         min.support = -Inf, max.support = Inf) {
  
  ### Inverse truncation limits - CDF
  nn <- n; rm(n)
  inv_LL <- pnorm(min.support, mean, sd)
  inv_UL <- pnorm(max.support, mean, sd)
  
  ## Unif in boundaries + inverse CDF
  uu <- runif(nn, inv_LL, inv_UL)
  
  ## Truncated-normal data
  ans <- qnorm(uu, mean, sd)
  
  ans[min.support >= max.support]  <- NaN
  ans[sd <= 0] <- NaN
  
  ans
  
  
}

