##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
# 20161212

Q.reg <- function(y, pvector = NULL, length.arg = NULL) {
  
  NOS <- if (!length(dim(y))) exp(0) else dim(y)[2]
  length.arg <- if (length(pvector)) length(pvector) else length.arg
  if (is.null(length.arg))
    stop("Enter a prototype vector of quantiles")

  kronecker(rbind(rep(1, NOS * length.arg)), Y = y)[, 
          interleave.VGAM(NOS * length.arg, M1 = length.arg)]
}