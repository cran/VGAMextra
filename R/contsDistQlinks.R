##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.

### Link function for any quantile  ###
### 22/Nov/2016


toppleQlink <- function(theta,
                        p = stop("Argument 'p' must be specified."),
                        bvalue = NULL, inverse = FALSE,
                        deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for 'p'. Must lie between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    m.help <- paste("Max( 1 - sqrt[ 1 - ", p, "^(1/theta) ] ) ", sep = "")
    t.string <- 
      if (short) paste("toppleQlink(", theta, "; ", p, ")" , sep = "") else
      paste("( 1 - sqrt[ 1 - ", p, "^(1/", 
            as.char.expression(theta), 
            ") ] ) / ", m.help, sep = "")
    
    if (tag)
      t.string <- paste("Topp-Leone quantile link: ", t.string, sep = "")
    
    return(t.string)
  }
  
  ###  No need to include '!inverse' since both require theta in (0, 1) ###
  if (length(bvalue)) {  
    theta[theta <= 0] <- bvalue
    theta[theta >= 1] <- 1 - bvalue
  }
  
  m.thets <- ppoints(10^3)
  m.pthe  <- exp(log1p(-p^(1/m.thets))); rm(m.thets) # max(1 - p^(1/theta))
  m.max   <- max( exp(log1p(-sqrt(m.pthe))) )[1]
  
  m.pthe  <- exp(log1p(-p^(1/theta)))              ## 1 - p^(1/theta)
  der.thp <- -( p^(1/theta) ) * log(p) / theta^2   ## dp^(1/theta)/d theta
  d1.dthe <- der.thp / (2 * sqrt(m.pthe))          ## d eta/d theta
  d.prod  <- 2 * theta * sqrt(m.pthe) - theta^2 * (d1.dthe) 
  d2e.dt2 <- -(log(p)/2) * (1/m.max) * ( (sqrt(m.pthe) * der.thp * theta^2- 
                  d.prod * p^(1/theta) ) / ( sqrt(m.pthe) * theta^2)^2 ) 
  
  if (inverse) {
    
    switch(deriv + 1, 
           log(p) / log1p(-(1 - theta * m.max)^2),
           m.max / d1.dthe,
           -(m.max / d1.dthe)^3 * (d2e.dt2))
    
  } else {
    
    switch(deriv + 1,
           exp(log1p(-sqrt(m.pthe))) / m.max,
           d1.dthe / m.max,
           d2e.dt2)
    
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
        paste(as.char.expression(theta),
              " * sqrt( -2 * log(1 - ", p, ") )", sep = "")
    
    if (tag)
      r.string <- paste("Rayleigh quantile link: ", r.string, sep = "")
    
    return(r.string)
  }
  
  if (length(bvalue)) 
    theta[theta <= 0] <- bvalue
  
  if (any(theta <= 0)) 
    theta[theta <= 0] <- NaN

  if (inverse) {
    
    switch(deriv + 1, 
           theta / sqrt(-2 * log1p(-p)),
           1 / sqrt(-2 * log1p(-p)),
           log(1))
    
  } else {
    
    switch(deriv + 1,
           theta * sqrt(-2 * log1p(-p)),
           sqrt(-2 * log1p(-p)),
           log(1))
    
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
  
  if (!inverse & length(bvalue))
    theta[theta <= 0] <- bvalue
  
  if (inverse & length(bvalue))
    theta[theta <= log(y0)] <- log(y0) + bvalue
  
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
    
  } else {
    
    switch(deriv + 1,
           log(y0) + sqrt(-log1p(-p) / theta ),
           -sqrt( -log1p(-p) ) / (2 * theta^(1.5)) ,
           3 * sqrt( -log1p(-p) ) / (4 * theta^(2.5)) )
    
  }
}



