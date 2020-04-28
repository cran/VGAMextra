##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0



### Link function for any quantile  ###
### 22/Nov/2016 
### Modified on 26/Nov/2018, to be included in the paper.


toppleQlink <- function(theta,
                        p = stop("Argument 'p' must be specified."),
                        bvalue = NULL, inverse = FALSE,
                        deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for 'p'. Must lie between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    t.string <- 
      if (short) paste("toppleQlink(", 
                       as.char.expression(theta), "; ", p,
                       ")", sep = "") else
        paste("logit[1 - sqrt( 1 - ", p, "^(1/", 
              as.char.expression(theta), 
              "))]", sep = "")
    
    if (tag)
      t.string <- paste("Topp-Leone quantile link: ", t.string, sep = "")
    
    return(t.string)
  }
  
  ###  No need to include '!inverse' since both require theta in (0, 1) ###
  ### 27 / Nov / 2018.. Above is no longer true. Taking the logit of the
  ### quantile function implies eta (inverse = TRUE) to be any real number.
  #if (length(bvalue)) {  
  #  theta[theta <= 0] <- bvalue
  #  theta[theta >= 1] <- 1 - bvalue
  #}
  if (!inverse) {
        theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
        theta[theta >= 1] <- if (length(bvalue)) 1 - bvalue else NaN
  }
  
  ## 5/12/18
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }

  hh <- 1e-5
  second.deriv.noInv <-  
    (logitlink(1 - sqrt(1 - p^(1/(theta + hh))), inverse = FALSE) -
    2 *  logitlink(1 - sqrt(1 - p^(1/theta)), inverse = FALSE) +
    logitlink(1 - sqrt(1 - p^(1/(theta - hh))), inverse = FALSE)) / hh^2
  
  second.deriv.Inv <-  
    ( log(p) / log1p(-(1 - logitlink(theta + hh, inverse = TRUE))^2) -
       2 *   log(p) / log1p(-(1 - logitlink(theta, inverse = TRUE))^2) +
      log(p) / log1p(-(1 - logitlink(theta - hh, inverse = TRUE))^2) )/hh^2

  if (inverse) {
    my.k <- 2 * ( theta^2 ) * sqrt(1 - p^(1/theta))
    switch(deriv + 1, 
           log(p) / log1p(-(1 - logitlink(theta, inverse = TRUE))^2),
           -my.k / (logitlink(1 - sqrt(1 - p^(1/theta)),
                          inverse = FALSE, deriv = 1) *
                      (p^(1 / theta)) * log(p)),
           second.deriv.Inv)
    
  } else {
    my.k <- 2 * ( theta^2 ) * sqrt(1 - p^(1/theta))
    switch(deriv + 1,
      logitlink(1 - sqrt(1 - p^(1/theta)), inverse = FALSE),
       -logitlink(1 - sqrt(1 - p^(1/theta)), inverse = FALSE, deriv = 1) *
             (p^(1 / theta)) * log(p) / my.k,
           second.deriv.noInv)
    #-(log(p) * p^(1 / theta)) /(2 * ((1 - p^(1/theta))^(-1/2) - 1) *
    #  theta^2 * (1 - p^(1 / theta))^(1.5))
  }
}



rayleighQlink <- function(theta,
                          p = stop("Argument 'p' must be specified."),
                          bvalue = NULL, inverse = FALSE,
                          deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid value for 'p'. Must be between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    r.string <- 
      if (short) paste("rayleighQlink(", theta, "; ", p, ")", sep = "") else
        paste("log(", as.char.expression(theta), ")",
              " + (1/2) * loglog[(1 - ", p, ")^(-2)]", sep = "")
    
    if (tag)
      r.string <- paste("Rayleigh quantile link: ", r.string, sep = "")
    
    return(r.string)
  }
  
  ## 5/12/18
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }

  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN

  
  if (inverse) {
    
    switch(deriv + 1, 
           exp(theta) / sqrt(-2 * log1p(-p)),
           theta,
           theta )
    
  } else {
    
    switch(deriv + 1,
           log(theta) + (1/2) * logloglink((1 - p)^(-2)),
           1/theta ,
           -1 / theta^2 )
    
  }
}






benini1Qlink <- function(theta, 
                         p = stop("Argument 'p' must be entered."),
                         y0 = stop("Argument 'y0' must be specified."),
                         bvalue = NULL, inverse = FALSE,
                         deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'.")
  
  if (length(y0) & !is.Numeric(y0, positive = TRUE, length.arg = 1) )
    stop("Invalid input for argument 'y0'.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    b.string <- 
      if (short) paste("beniniQlink(", theta, "; ", p, ")", sep = "") else
        paste("log(", y0, ") + sqrt(-log(1 - ", p, ") / ", 
              as.char.expression(theta), ")", sep = "")
    
    if (tag)
      b.string <- paste("Benini quantile link:", b.string)
    
    return(b.string)
    
  }
  
  ## 5/12/18
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  ## 5/12/18
  if (inverse)
    theta[theta <= log(y0)] <- 
      if(!length(bvalue)) log(y0) + 1e-1 else log(y0) + NaN
                                                    
  
  if (inverse) {
    m.min <- log(y0)
    if (!deriv & any(theta <= m.min)) {
      theta[theta <= m.min] <- NaN
      warning("Values of 'eta' out of range. These were replaced by NaN.")
    }
    
    switch(deriv + 1,
           -log1p(-p) / ( log( exp(theta) / y0 ) )^2,
           -2 * theta^(1.5) / sqrt( -log1p(-p) ),
           -6 * theta^2 / log1p(-p))
    #log( y0 * sqrt(-log1p(-p) / theta))
  } else {
    
    switch(deriv + 1,
           log(y0) + sqrt(-log1p(-p) / theta ),
           -sqrt( -log1p(-p) ) / (2 * theta^(1.5)) ,
           3 * sqrt( -log1p(-p) ) / (4 * theta^(2.5)) )
    
  }
}



