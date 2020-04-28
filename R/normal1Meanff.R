##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0
# 2018/11/19

normal1Meanff <- function(zero = NULL, link = "loglink",
                        p.quant = NULL, fixed.sd = 1){
  
  var.arg <-  FALSE
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")
  
  if (FALSE)
  if (length(p.quant)) {
    #p.quant <- rep(p.quant, 10)[1:3]
    if (any(p.quant <= 0) || any(p.quant >= 1)) 
      stop("Invalid quantiles. Must lie between 0 and 1.")
    
    if (identical(link, "normal1MeanQlink") ||
                          identical(link, normal1MeanQlink)) {
      link <- as.list(substitute(link(p = p.quant, sd.cons = fixed.sd)))
    } else {
      stop("Invalid link for quantile regression with this ",
           "family function.")
    }
  } else {
    if (identical(link, "normal1MeanQlink") ||
                          identical(link, normal1MeanQlink)) {
      stop("Invalid link to model the Mean.")
    } else {
      link <- as.list(substitute(link)) 
    }
    
  }
  
  #fixed.mean <- fixed.sd  # Just to reproduce the code from normal1sdff()
  if (!is.logical(var.arg))
    stop("Wrong input for argument 'var.arg'")
  
  mynames.of <-  if ( var.arg) {
    namesof("var", link, earg = earg, tag = FALSE)
  } else {
    namesof("sd", link, earg = earg, tag = FALSE)
  }
  
  
  new("vglmff",
      blurb = c("1-parameter Normal distribution (estimates the Mean) \n\n",
                "Link:     ", mynames.of, "\n", 
                "SD:     FIXED (must be entered)\n",
                "Mean: Mean"),
      
  
      constraints = eval(substitute(expression({
        M1 <- 1
        constraints <- cm.zero.VGAM(constraints, x = x, zero = .zero , 
              M = M, predictors.names = mynames1, M1 = M1)
      }), list( .zero = zero ))),
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 1,
             Q1 = 1,
             fixed.sd = .fixed.sd ,
             p.quant    = .p.quant ,
             zero = .zero )
      }, list( .zero = zero , .fixed.sd = fixed.sd , 
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
        
        mynames1 <- param.names("Mean", M)
        predictors.names <- namesof(mynames1, .link , 
                                    earg = .earg , short = TRUE)
        #my.start <- matrix( sqrt(sum( (y - .fixed.mean)^2 ) / n), 
        #                    nrow = n, ncol = M1)
        my.start <- matrix(apply(y, 2, mean),
                           nrow = n, ncol = M1, byrow = TRUE) 
        if (!length(etastart))
          etastart <- cbind(theta2eta(my.start, .link , earg = .earg ))
        
      }), list( .link = link, .earg = earg , .var.arg = var.arg ,
                .fixed.sd = fixed.sd , .p.quant = p.quant ))), 
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        eta2theta(eta, .link , earg = .earg )
      }),list( .link = link, .earg = earg )),
      
      
      
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
        misc$fixed.sd <- .fixed.sd
      }), list( .link = link, .earg = earg , .fixed.sd = fixed.sd , 
                .p.quant = p.quant , .var.arg = var.arg ))),
      
      
       linkfun = eval(substitute(function(mu, extra = NULL) {
        theta2eta(mu, .link , earg = .earg )
      }), list( .link = link, .earg = earg )),
      
      
      
      loglikelihood = eval(substitute(
                 function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL, summation = TRUE) {
      
          my.mean <- eta2theta(eta, .link , .earg )
        
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dnorm(x = y, mean = my.mean , 
                                    sd = .fixed.sd , log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
      }}, list( .fixed.sd = fixed.sd , .p.quant = p.quant, 
                .link = link , .earg = earg , .var.arg = var.arg ))),
      
      
      vfamily = c("normal1Meanff"),
      
      simslot = eval(substitute(
        function(object, nsim) {
          pwts <- if (length(pwts <- object@prior.weights) > 0)
           pwts else weights(object, type = "prior")
          if (any(pwts != 1)) 
            warning("ignoring prior weights")
          mu <- fitted(object)

          rnorm(nsim * length(my.mean), mean = my.mean,  sd = .fixed.sd )
        }, list( .link = link, .earg = earg , 
                 .fixed.sd = fixed.sd ))),
      
      
      
      deriv = eval(substitute(expression({
        
        sd.fix <- matrix( .fixed.sd , nrow(y), ncol(y), byrow = TRUE)
        n <- nrow(y)
        my.mean  <- eta2theta(eta, .link , .earg )
        dl.dmean <- (y - my.mean) /  sd.fix^2
        dmean.deta <- dtheta.deta(my.mean, .link , earg = .earg )
        ans <- c(w) * dl.dmean * dmean.deta 
        ans
        
      }), list( .link = link, .earg = earg ,
                .fixed.sd = fixed.sd , .p.quant = p.quant ))),
      
      
      weight = eval(substitute(expression({
      
        ned2l.dmean <- 1/sd.fix^2
        wz <- ned2l.dmean * dmean.deta^2
        
        c(w) * wz
      }), list( .fixed.sd = fixed.sd ))))
}