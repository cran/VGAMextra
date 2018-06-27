##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
# 20161212

normal1sdff <- function(zero = NULL, link = "loge", fixed.mean = 0,
                        p.quant = NULL,  var.arg = FALSE){
  
  
  if (length(p.quant)) {
    #p.quant <- rep(p.quant, 10)[1:3]
    if (any(p.quant <= 0) || any(p.quant >= 1)) 
      stop("Invalid quantiles. Must lie between 0 and 1.")
    
    if (identical(link, "normal1sdQlink") ||
                          identical(link, normal1sdQlink)) {
      link <- as.list(substitute(link(p = p.quant, mean = fixed.mean)))
    } else {
      stop("Invalid link for quantile regression with this ",
           "family function.")
    }
  } else {
    if (identical(link, "normal1sdQlink") ||
                          identical(link, normal1sdQlink)) {
      stop("Invalid link to model the standard deviation.")
    } else {
      link <- as.list(substitute(link)) 
    }
    
  }
  
  if (!is.logical(var.arg))
    stop("Wrong input for argument 'var.arg'")
  
  earg <- link2list(link)
  link <- attr(earg, "function.name")
  
  mynames.of <-  if ( var.arg) {
    namesof("var", link, earg = earg, tag = FALSE)
  } else {
    namesof("sd", link, earg = earg, tag = FALSE)
  }
  
  
  new("vglmff",
      blurb = c("1-parameter Normal distribution (sigma^2) \n\n",
                "Link:     ", mynames.of, "\n", 
                "Mean:     FIXED (must be entered)\n",
                "Variance: sigma^2"),
      
  
      constraints = eval(substitute(expression({
        M1 <- 1
        constraints <- cm.zero.VGAM(constraints, x = x, zero = .zero , 
              M = M, predictors.names = mynames1, M1 = M1)
      }), list( .zero = zero ))),
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 1,
             Q1 = 1,
             fixed.mean = .fixed.mean ,
             p.quant    = .p.quant ,
             zero = .zero )
      }, list( .zero = zero , .fixed.mean = fixed.mean , 
               .p.quant = p.quant ))),
      
      
      
      initialize = eval(substitute(expression({
        my.temp <- w.y.check(w = w, y = y,
                             Is.positive.y = FALSE,
                             ncol.w.max = Inf,
                             ncol.y.max = Inf,
                             out.wy = TRUE,
                             colsyperw = 1,
                             maximize = TRUE)
        w <- my.temp$w
        y <- my.temp$y
        n <- nrow(y); NOS <- ncol(y)
        
        if (length( .p.quant ) && (ncol(y) != length( .p.quant ))) 
          warning("\n\nSeems like the quantile vectors entered does not",
               " match the number of responses. It was recycled.\n")
        
      
        M  <- if (is.matrix(y)) ncol(y) else 1
        M1 <- ncol(y)
        
        mynames1 <- param.names(if ( .var.arg ) "var" else "sd", M)
        predictors.names <- namesof(mynames1, .link , 
                                    earg = .earg , short = TRUE)
        #my.start <- matrix( sqrt(sum( (y - .fixed.mean)^2 ) / n), 
        #                    nrow = n, ncol = M1)
        my.start <- matrix(apply(y, 2, if (.var.arg ) var else sd),
                           nrow = n, ncol = M1, byrow = TRUE) 
        if (!length(etastart))
          etastart <- cbind(theta2eta(my.start, .link , earg = .earg ))
        
      }), list( .link = link, .earg = earg , .var.arg = var.arg ,
                .fixed.mean = fixed.mean , .p.quant = p.quant ))), 
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta, .link , earg = .earg )
      }),list( .link = link, .earg = earg , .fixed.mean = fixed.mean )),
      
      
      
      last = eval(substitute(expression({
        n  <- nrow(y)
        misc$link <- rep_len( .link , M)
        names(misc$link) <- mynames1
        
        misc$earg <- vector("list", M)
        names(misc$earg) <- names(misc$link)
        for (ii in 1:M)
          misc$earg[[ii]] <- .earg
        
        misc$expected <- TRUE
        misc$multipleResponses <- TRUE
        misc$M1 <- M1
        misc$var.arg <- .var.arg
        misc$p.quant <- .p.quant
        misc$fixed.mean <- .fixed.mean
      }), list( .link = link, .earg = earg , .fixed.mean = fixed.mean , 
                .p.quant = p.quant , .var.arg = var.arg ))),
      
      
       linkfun = eval(substitute(function(mu, extra = NULL) {
        theta2eta(mu, .link , earg = .earg )
      }), list( .link = link, .earg = earg )),
      
      
      
      loglikelihood = eval(substitute(
                 function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL, summation = TRUE) {
                   
        if ( .var.arg ) {
          my.sd <- sqrt(abs(eta2theta(eta, .link , .earg )))
        } else {
          my.sd <- eta2theta(eta, .link , .earg )
        }
        
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dnorm(x = y, mean = .fixed.mean, 
                                    sd = my.sd , log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
      }}, list( .fixed.mean = fixed.mean , .p.quant = p.quant, 
                .link = link , .earg = earg , .var.arg = var.arg ))),
      
      
      vfamily = c("normal1ffsd"),
      
      simslot = eval(substitute(
        function(object, nsim) {
          pwts <- if (length(pwts <- object@prior.weights) > 0)
           pwts else weights(object, type = "prior")
          if (any(pwts != 1)) 
            warning("ignoring prior weights")
          mu <- fitted(object)

          rnorm(nsim * length(my.sd), mean = .fixed.mean,  sd = my.sd)
        }, list( .link = link, .earg = earg , 
                 .fixed.mean = fixed.mean ))),
      
      
      
      deriv = eval(substitute(expression({
        n <- nrow(y)
         if ( .var.arg ) {
           my.var  <- eta2theta(eta, .link , .earg )
           dl.dvar <- (y - .fixed.mean )^2 / (2 * my.var^2) -
                                                  1 / (2 * my.var)
           dvar.deta <- dtheta.deta(my.var, .link , earg = .earg )
           ans <- c(w) * dl.dvar * dvar.deta
         } else {
           my.sd  <- eta2theta(eta, .link , .earg )
           dl.dsd <- (y - .fixed.mean)^2 / my.sd^3 - 1 / my.sd
           dsd.deta <- dtheta.deta(my.sd, .link , earg = .earg )
           ans <- c(w) * dl.dsd * dsd.deta
         }
        
        ans
        
      }), list( .link = link, .earg = earg ,  .var.arg = var.arg ,
                .fixed.mean = fixed.mean , .p.quant = p.quant ))),
      
      
      weight = eval(substitute(expression({
         if ( .var.arg ) {
           ned2l.dvar <- 1 / (2 * my.var^2)
           wz <- ned2l.dvar * dvar.deta^2
         } else {
           ned2l.dsd <- 2 /  my.sd^2
           wz <- ned2l.dsd * dsd.deta^2
         }
        
        c(w) * wz
      }), list( .var.arg = var.arg ))))
}