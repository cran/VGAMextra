##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.

# 20170102
# Is ncp an argument of the inverse chisq? So far, set to zero. 

dinv.chisq <- function(x, df, log = FALSE) {
  if (log) {
    dchisq(x = 1/x, df = df, ncp = 0, log = TRUE ) - 2 * log(x)
  } else {
    dchisq(x = 1/x, df = df, ncp = 0, log = FALSE ) * (1 / x^2)
  }
}


pinv.chisq <- function(q, df, lower.tail = TRUE, log.p = FALSE)
  pchisq(q = 1/q, ncp = 0, df = df, lower.tail = !lower.tail, log.p = log.p)


qinv.chisq <- function(p, df, lower.tail = TRUE, log.p = FALSE)
     1 / qchisq(p = p, df = df, ncp = 0, 
                   lower.tail = !lower.tail, log.p = log.p)


rinv.chisq <- function(n, df) 
     1 / rchisq(n = n, df = df, ncp = 0)