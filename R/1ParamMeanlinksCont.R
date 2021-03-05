##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

# All links renamed in Feb2021 Mean to 'M'
expMlink <- function(theta, location = 0, 
                     bvalue = NULL, inverse = FALSE,
                     deriv = 0, short = TRUE, tag = FALSE) {
  
  if (!is.Numeric(deriv, length.arg = 1,
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  Aloc <- location; rm(location)
  
  if (is.character(theta)) {
    e.ret <- if (short) paste("expMlink(", 
                              Aloc, " + ", as.char.expression(theta), 
                              ")", sep = "" ) else
      paste("log(", Aloc, " + ", as.char.expression(theta), ")", sep = "")
    
    if (tag)
      e.ret <- paste("Exponential mean link: ", e.ret, sep = "")
    return(e.ret)
  }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    switch(deriv + 1,
           (exp(theta) - Aloc)^(-1),
           -(Aloc * theta^2 + theta),
           theta^2 * (Aloc + theta^(-1)) * (2 * Aloc * theta + 1))
    
  } else {
    
    switch(deriv + 1,
           log(Aloc + theta^(-1)),
           -1 /( theta^2 * (Aloc + theta^(-1)) ),
           (2 * Aloc * theta + 1) / (theta^2 * (Aloc + theta^(-1)) )^2 )
    
  }
}






inv.chisqMlink <- function(theta, bvalue = NULL, inverse = FALSE, 
                            deriv = 0, short = TRUE, tag = FALSE) {
  
  
  if (!is.numeric(deriv) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
      inv.ret <- if (short) paste("inv.chisqMlink(",
                                  theta, ")", sep = "") else
      paste("log( 1 / (", as.char.expression(theta), " - 2)", sep = "")
    
    if (tag)
      inv.ret <- paste("Inverse chi-squared mean link: ",
                       inv.ret, sep = "")
    
    return(inv.ret)
  }
  
  if (!inverse) {
    theta[theta <= 2] <- if (length(bvalue)) bvalue else NaN
  }
  
  if (inverse) {
    switch(deriv + 1, 
           exp(-theta) + 2,
           -(theta - 2),
           (theta -2))
  } else {
    switch(deriv + 1, 
           -log(theta - 2), 
           -exp(-log(theta - 2)),
           1/(theta - 2)^2)
    
  }
}






ffff.help <- function(x) 1 - 4^x * gamma(1 + x)^2/gamma(2 + 2* x)


toppleMlink <- function(theta, bvalue = NULL, inverse = FALSE,
                        deriv = 0, short = TRUE, tag = FALSE) {
  
  if (!is.Numeric(deriv, length.arg = 1,
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    if (short) {
      tp.ret <- paste("toppleMlink(", theta, ")", sep = "" ) 
      tp.ret <- if (!tag) tp.ret else
                    paste("Toople mean link: ", tp.ret, sep = "")
      return(tp.ret)
    } else {
      tp.ret <- c("logit(mean.topple / max(mean.topple))") 
      tp.ret <- if(!tag) tp.ret  else
        paste("Toople mean link: ", tp.ret, sep = "")
      return(cat(tp.ret, ".\n Here: ",
                 c("'mean.toople = 1 - 4^(", 
                   as.character(theta),") * gamma(1 +   ",
                   as.character(theta), ")^2 / " , "gamma(2 + 2 * ",
                   as.character(theta), ")'\n")))
    }
  }
  
  if (!inverse) {
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
    theta[theta >= 1] <- if (length(bvalue)) bvalue else NaN
  }
    
  
  t.max  <- max(ffff.help(ppoints(1e4))); hh     <- 1e-3
  der.1  <- if(!deriv) NA else
      (1/t.max) * (ffff.help(theta + hh) - ffff.help(theta - hh)) / (2*hh)
  der.2  <- if (!deriv) NA else
            (1/t.max) * (ffff.help(theta + hh) - 2 * ffff.help(theta) + 
                             ffff.help(theta - hh)) / hh^2
  
  if (inverse) {
    
    if (!deriv) {
      
      m.ret <- NULL
      m.ins <- inspectVGAMextra(x = c(theta), b.valueG = 7.5,
                                inf.s = FALSE)
      the.t <- m.ins[[1]]
      
      if (length(the.t))
        m.ret <- newtonRaphson.basic(f = function(x, h, eta) {
          logitlink(theta = ffff.help(x)/t.max) - sign(eta) * abs( eta ) },
                                     fprime = function(x, h, eta) {
          int.d1 <- (1/t.max)*(ffff.help(x + h) - ffff.help(x - h)) / (2*h)
            logitlink(theta = ffff.help(x)/t.max, deriv = 1) * int.d1 }, 
        a = rep(1e-6, length(the.t)), 
        b = rep(1 - hh, length(the.t)),
        eta = the.t, h = hh)
      
      # Passing down 'theta' rather than c(theta). Dimension required.
      return(inspectVGAMextra(x = theta, the.NR = m.ret,
                              inf.s = exp(0), neginf = 0,
                              inverse = TRUE, b.valueG = exp(0),
                              b.valueL = 0, extra = m.ins[[2]]))
    } else {
      switch(deriv,
      (logitlink(theta = ffff.help(theta)/t.max, deriv = 1) * der.1)^(-1),
       -(logitlink(ffff.help(theta)/t.max, deriv = 1) * der.1)^(-3) * 
          logitlink(theta = ffff.help(theta)/t.max, deriv = 1) * der.2 + 
          (der.1^2) * logitlink(theta = ffff.help(theta)/t.max, deriv = 2))
      
    }
  } else {
    switch(deriv + 1, 
           logitlink(theta = ffff.help(theta)/t.max), 
           logitlink(theta = ffff.help(theta)/t.max, deriv = 1) * der.1,
           logitlink(theta = ffff.help(theta)/t.max, deriv = 1) * der.2 + 
          (der.1^2) * logitlink(theta = ffff.help(theta)/t.max, deriv = 2))
  }
}






rayleighMlink <- function(theta, bvalue = NULL, inverse = FALSE,
                          deriv = 0, short = TRUE, tag = FALSE) {
  
  if (!is.Numeric(deriv, length.arg = 1,
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
   ry.ret <- if (short) paste("rayleighMlink(",theta, ")", sep = "") else
      paste("log(", as.char.expression(theta), 
             " * gamma(0.5) / sqrt(2) )", sep = "")
    
    if (tag)
      ry.ret <- paste("Rayleigh mean link: ", ry.ret, sep = "")
    return(ry.ret)
  }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    switch(deriv + 1,
           exp(theta) * sqrt(2) / gamma(0.5),
           theta, theta)
    
  } else {
    
    switch(deriv + 1,
           log(theta * sqrt(pi / 2)),
           1/theta, -1/theta^2 )
    
  }
}








maxwellMlink <- function(theta, bvalue = NULL, inverse = FALSE,
                         deriv = 0, short = TRUE, tag = FALSE) {
  
  if (!is.Numeric(deriv, length.arg = 1,
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    mx.ret <- if (short) paste("MaxwellMlink(", theta, ")", sep = "") else
      paste("log(", as.char.expression(theta), 
            " * 2 / sqrt(pi/2) )", sep = "")
    
    if (tag)
      mx.ret <- paste("Maxwell mean link: ", mx.ret, sep = "")
    return(mx.ret)
  }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    switch(deriv + 1,
           8 * exp(-2 * theta) / pi ,
           -2 * theta, 4 * theta) 
    
  } else {
    
    switch(deriv + 1,
           log(sqrt(8 / pi) * theta^(-0.5)),
           -1 / (2 * theta), 1/(2 * theta^2))
    
  }
}
