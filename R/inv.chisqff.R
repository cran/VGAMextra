##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
# 20170102

inv.chisq <- function(link = "loge", zero = NULL){
  
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  validpar <- (identical(link, "inv.chisqMeanlink") || 
                               identical(link, inv.chisqMeanlink))
  
  
  new("vglmff",
      blurb = c("Inverse chi-squared distribution\n\n",
                "Link:     ",
                namesof("nu", link, earg = earg, tag = FALSE), "\n", 
                "Mean:     1/(nu - 2), nu > 2 \ n",
                "Variance: 2/(nu - 2)^2 * (nu - 4), nu > 4"),
      
      
      constraints = eval(substitute(expression({
        M1 <- 1
        constraints <- cm.zero.VGAM(constraints, x = x, zero = .zero , 
                                    predictors.names =  predictors.names, 
                                    M1 = M1, M = M)
      }), list( .zero = zero ))),
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 1,
             Q1 = 1,
             zero = .zero )
      }, list( .zero = zero ))),
      
      
      
      initialize = eval(substitute(expression({
        my.temp <- w.y.check(w = w, y = y,
                             Is.positive.y = FALSE,
                             ncol.w.max = Inf,
                             ncol.y.max = Inf,
                             out.wy = TRUE,
                             colsyperw = 1,
                             maximize = TRUE)
        w  <- my.temp$w
        y  <- my.temp$y
        n  <- nrow(y); NOS <- ncol(y)
        M  <- if (is.matrix(y)) ncol(y) else 1
        M1 <- ncol(y)
        
        dofnames <- param.names("dof", M)
        predictors.names <- namesof(dofnames, .link , 
                                    earg = .earg , short = TRUE)
        
        vec.init <- matrix(colMeans(y), nrow = n, 
                              ncol = M1, byrow = TRUE) + 2
        if (!length(etastart))
          etastart <- cbind(theta2eta(vec.init, .link , earg = .earg ))
        
      }), list( .link = link, .earg = earg  ))), 
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        1 / (eta2theta(eta, link = .link , earg = .earg ) - 2)
      }, list( .link = link, .earg = earg ) )),
      
      
      
      #validparams = eval(substitute(function(eta, y, extra = NULL) {
      #  print(.validpar)
      #  dofs2 <- eta2theta(eta, link = .link , earg = .earg)
      #  alright <- if ( .validpar ) all(dofs2 > 2) else TRUE
      #  alright <- all(dofs2 > 2)
      #  print(head(dofs2))
      #  alright
      
      #}, list( .link = link , .earg = earg , .validpar = validpar ) )),
      
      
      
      last = eval(substitute(expression({
        n  <- nrow(y)
        misc$link <- rep_len( .link , M)
        names(misc$link) <- dofnames
        
        misc$earg <- vector("list", M)
        names(misc$earg) <- names(misc$link)
        for (ii in 1:M)
          misc$earg[[ii]] <- .earg
        
        misc$expected <- TRUE
        misc$multipleResponses <- TRUE
        misc$M1 <- M1
      }), list( .link = link, .earg = earg ))),
      
      
      #linkfun = eval(substitute(function(mu, extra = NULL) {
      #  theta2eta(mu, .link , earg = .earg )
      #}), list( .link = link, .earg = earg )),
      
      
      loglikelihood = eval(substitute( function(mu, y, w, 
                                       residuals = FALSE, eta,
                                       extra = NULL, summation = TRUE) {
        
          dofs <- eta2theta(eta, link = .link , earg = .earg )
          
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.chis <- c(w) * dinv.chisq(x = y, df = dofs, log = TRUE)
            if (summation) {
              sum(ll.chis)
            } else {
              ll.chis
            }
          }}, list( .link = link , .earg = earg ))),
      
      
      vfamily = c("InvChisq"),
      
      simslot = eval(substitute(
        function(object, nsim) {
          pwts <- if (length(pwts <- object@prior.weights) > 0)
            pwts else weights(object, type = "prior")
          if (any(pwts != 1)) 
            warning("ignoring prior weights")
          mu <- fitted(object)
          
          1 / rnorm(nsim * length(dofs), dof = dofs )
        }, list( .link = link, .earg = earg ))),
      
      
      
      deriv = eval(substitute(expression({
        nn   <- nrow(y)
        dofs <- eta2theta(eta, link = .link , earg = .earg )
        
        dl.dnu <- -(0.5) * (digamma(dofs/2) + log(y) + log(2))
        dnu.de <- dtheta.deta(dofs, .link , earg = .earg )
        
         c(w) * dl.dnu * dnu.de
        
      }), list( .link = link, .earg = earg  ))),
      
      
      
      weight = expression({
        ned2l.dnu <- trigamma(dofs/2)/4
        wz <- ned2l.dnu * dnu.de^2
        c(w) * wz
      }))
}