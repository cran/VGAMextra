##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0
# Modified on 26/11/2018, to be included in the paper.

maxwellQlink <- function(theta,
                         p = stop("Argument 'p' must be specified."),
                         bvalue = NULL, inverse = FALSE,
                         deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'. Must lie between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    m.string <- 
      if (short) paste("maxwellQlink(", 
                       theta, "; ", p, ")", sep = "") else
                         paste("(1/2) * log[2 * qgamma(", p, ", 1.5)/ ", 
                               as.char.expression(theta), 
                                "]", sep = "")
    
    if (tag)
      m.string <- paste("Maxwell quantile link: ", m.string, sep = "")
    
    return(m.string)
  }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  #if (length(p) > 1)
  #  p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE)
  
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }
  
  
  if (inverse) {
    
    switch(deriv + 1,
           2 * qgamma(p = p, shape = 1.5, log.p = FALSE) / exp(2 * theta),
           -2 * theta,
           4 * theta )
    
  } else {
    
    switch(deriv + 1, 
           (1 / 2) * log(2 * qgamma(p = p, shape = 1.5) / theta),
           -1 / (2 * theta),
           1 / (2 * theta^2) )
    
  }
}