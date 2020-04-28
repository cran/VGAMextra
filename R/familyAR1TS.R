##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


### Exact EIms for the AR1 model.
# Here, theta = c(m*, serrors^2, theta)
AR1EIM.G2 <- function(y, drift, sdError, AR1coeff, var.arg = TRUE,
                      order = 1, nodrift = FALSE) {
  
  y <- cbind(y); nn <- nrow(y)
  RMat <- diag(c(sdError^2))
  dRMat <- cbind( if(!nodrift) matrix(0, nrow = nn, ncol = 1) else NULL, 
                  if (var.arg) matrix(1, nrow = nn, ncol = 1) else
                    2 * sdError,  
                  matrix(0, nrow = nn, ncol = 1))
  
  dmt0 <- matrix(1, nrow = nn, ncol = 1)
  dmt  <- cbind(if (!nodrift) dmt0 else NULL,
                rep_len(0, nn), c(1 / (1 - AR1coeff[1])^2, y[-nn]))
  
  M <- 3  # For AR1
  try.comb <- combVGAMextra(1:M, nodrift = nodrift)
  
  finMat <- apply(try.comb, 1, function(x) {
    kl <- x
    invRMat <- solve(RMat)
    dRdthel <- diag(dRMat[, kl[1]])
    dRdthek <- diag(dRMat[, kl[2]])
    term.1 <- (0.5) * invRMat %*% dRdthel %*% invRMat %*% dRdthek
    term.2 <- diag(dmt[, kl[1]]) %*%
      invRMat %*% diag(dmt[, kl[2]])
    diag(term.1 + term.2)
  })
  finMat[ , ncol(finMat)] <- (0.995) * finMat[ , ncol(finMat)]
  finMat
  
}





### Density fo the AR1 model. dARff() may also work.
dAR1extra <- function(x,
                      drift = 0,  # Stationarity is the default
                      var.error = 1, ARcoef1 = 0.0,
                      type.likelihood = c("exact", "conditional"),
                      log = FALSE) {
  
  type.likelihood <- match.arg(type.likelihood,
                               c("exact", "conditional"))[1]
  
  is.vector.x <- is.vector(x)
  
  x <- as.matrix(x)
  drift <- as.matrix(drift)
  var.error <- as.matrix(var.error)
  ARcoef1 <- as.matrix(ARcoef1)
  LLL <- max(nrow(x), nrow(drift), nrow(var.error), nrow(ARcoef1))
  UUU <- max(ncol(x), ncol(drift), ncol(var.error), ncol(ARcoef1))
  x          <- matrix(x,         LLL, UUU)
  drift    <- matrix(drift,   LLL, UUU)
  var.error  <- matrix(var.error, LLL, UUU)
  rho        <- matrix(ARcoef1,   LLL, UUU)
  
  if (any(abs(rho) > 1))
    warning("Values of argument 'ARcoef1' are greater ",
            "than 1 in absolute value")
  
  if (!is.logical(log.arg <- log) || length(log) != 1) 
    stop("Bad input for argument 'log'.")
  rm(log)
  
  ans <- matrix(0.0, LLL, UUU)
  
  var.noise <- var.error / (1 - rho^2)
  
  ans[ 1, ] <- dnorm(x    = x[1, ],
                     mean = drift[ 1, ] / (1 - rho[1, ]), 
                     sd   = sqrt(var.noise[1, ]), log = log.arg)
  ans[-1, ] <- dnorm(x    = x[-1, ],
                     mean = drift[-1, ] + rho[-1, ] * x[-nrow(x), ],
                     sd   = sqrt(var.error[-1, ]), log = log.arg)
  
  if (type.likelihood == "conditional")
    ans[1, ] <- NA
  
  if (is.vector.x) as.vector(ans) else ans
}






if (FALSE)
AR1extra.control <- function(epsilon  = 1e-6,
                        maxit    = 30,
                        stepsize = 1,...){
  list(epsilon  = epsilon,
       maxit    = maxit,
       stepsize = stepsize,
       ...)
}





## Family function AR1extra() to fit order-1 Autoregressive models
AR1extra <-
  function(zero = c(if (var.arg) "var" else "sd", "rho"),
           type.EIM  = c("exact", "approximate")[1],
           var.arg = TRUE,  
           nodrift = FALSE,
           ldrift = "identitylink",
           lsd  = "loglink",
           lvar = "loglink",
           lrho = "rhobitlink",
           idrift  = NULL,
           isd  = NULL,
           ivar = NULL,
           irho = NULL,
           print.EIM = FALSE) {
    
    if (length(isd) && !is.Numeric(isd, positive = TRUE))
      stop("Bad input for argument 'isd'")
    
    if (length(ivar) && !is.Numeric(ivar, positive = TRUE))
      stop("Bad input for argument 'ivar'")
    
    if (length(irho) &&
          (!is.Numeric(irho) || any(abs(irho) > 1.0)))
      stop("Bad input for argument 'irho'")
    
    type.EIM <- match.arg(type.EIM, c("exact", "approximate"))[1]
    poratM   <- (type.EIM == "exact")
    imethod = 1
    if (!is.logical(nodrift) ||
          length(nodrift) != 1)
      stop("Argument 'nodrift' must be a single logical")
    
    if (!is.logical(var.arg) ||
          length(var.arg) != 1)
      stop("Argument 'var.arg' must be a single logical")
    
    if (!is.logical(print.EIM))
      stop("Invalid 'print.EIM'.")
    
    type.likelihood <- "exact"
    
    ismn <- idrift
    lsmn <- as.list(substitute(ldrift))
    esmn <- link2list(lsmn)
    lsmn <- attr(esmn, "function.name")     
    
    lsdv <- as.list(substitute(lsd))
    esdv <- link2list(lsdv)
    lsdv <- attr(esdv, "function.name")
    
    lvar  <- as.list(substitute(lvar))
    evar  <- link2list(lvar)
    lvar  <- attr(evar, "function.name")
    
    lrho <- as.list(substitute(lrho))
    erho <- link2list(lrho)
    lrho <- attr(erho, "function.name")     
    
    n.sc <- if (var.arg) "var" else "sd"
    l.sc <- if (var.arg) lvar else lsdv
    e.sc <- if (var.arg) evar else esdv
    
    new("vglmff", 
        blurb = c(ifelse(nodrift, "Two", "Three"),
                  "-parameter autoregressive process of order-1\n\n",
                  "Links:       ",
                  if (nodrift) "" else
                    paste(namesof("drift", lsmn, earg = esmn), ", ", 
                          sep = ""),
                  namesof(n.sc , l.sc, earg = e.sc), ", ",
                  namesof("rho", lrho, earg = erho), "\n",
                  "Model:       Y_t = drift + rho * Y_{t-1} + error_{t},",
                  "\n",
                  "             where 'error_{2:n}' ~ N(0, sigma^2) ",
                  "independently",
                  if (nodrift) ", and drift = 0" else "",
                  "\n",
                  "Mean:        drift / (1 - rho)", "\n",
                  "Correlation: rho = ARcoef1", "\n",
                  "Variance:    sd^2 / (1 - rho^2)"),
        
        
        
        
     constraints = eval(substitute(expression({
       M1 <- 3 - .nodrift
       dotzero <- .zero
       constraints <- 
         cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                      predictors.names = parameter.names,
                      M1 = M1)
        }), list( .zero = zero, .nodrift = nodrift ))),
     
     
     
     
     infos = eval(substitute(function(...) {
       list(M1 = 3 - .nodrift , 
            Q1 = 1, 
            expected = TRUE, 
            multipleResponse = TRUE,
            type.likelihood = .type.likelihood ,
            ldrift = if ( .nodrift ) NULL else .lsmn ,
            edrift = if ( .nodrift ) NULL else .esmn ,
            lvar = .lvar ,
            lsd  = .lsdv ,
            evar = .evar ,
            esd  = .esdv ,
            lrho = .lrho ,
            erho = .erho ,
            zero = .zero )
     }, list( .lsmn = lsmn, .lvar = lvar, .lsdv = lsdv, .lrho = lrho,
              .esmn = esmn, .evar = evar, .esdv = esdv, .erho = erho,
              .type.likelihood = type.likelihood,
              .nodrift = nodrift, .zero = zero))),
     
     
     
     
     initialize = eval(substitute(expression({
       extra$M1 <- M1 <- 3 - .nodrift
       check <- w.y.check(w = w, y = y,
                          Is.positive.y = FALSE,
                          ncol.w.max = Inf,
                          ncol.y.max = Inf,
                          out.wy = TRUE,
                          colsyperw = 1, 
                          maximize = TRUE)
       w <- check$w
       y <- check$y
       if ( .type.likelihood == "conditional") {
         w[1, ] <- 1.0e-6
       } else {
         w[1, ] <- 1.0e-6
       }
         
       NOS <- ncoly <- ncol(y)
       n <- nrow(y)
       M <- M1 * NOS
       extra$y <- y
       extra$print.EIM <- FALSE
       
       var.names <- param.names("var", NOS)
       sdv.names <- param.names("sd",  NOS)
       smn.names <- if ( .nodrift ) NULL else 
                                param.names("drift",   NOS)
       rho.names <- param.names("rho", NOS)
       
       parameter.names <- c(smn.names, 
                            if ( .var.arg ) var.names else sdv.names,
                            rho.names)
       parameter.names <- 
         parameter.names[interleave.VGAM(M, M1 = M1)]
       
       predictors.names <-
         c(if ( .nodrift ) NULL else
           namesof(smn.names, .lsmn , earg = .esmn , tag = FALSE),
           if ( .var.arg ) 
             namesof(var.names, .lvar , earg = .evar , tag = FALSE) else
               namesof(sdv.names, .lsdv , earg = .esdv , tag = FALSE),
           namesof(rho.names, .lrho , earg = .erho , tag = FALSE))
       
       predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]
       
        if (!length(etastart)) {
         init.smn <- matrix( if (length( .ismn )) .ismn else 0,
                             nrow = n, ncol = NOS, byrow = TRUE)    
         init.rho <- matrix(if (length( .irho )) .irho else 0.1, 
                            n, NOS, byrow = TRUE)           
         init.sdv <- matrix(if (length( .isdv )) .isdv else 1.0,  
                            n, NOS, byrow = TRUE)          
         init.var <- matrix(if (length( .ivar )) .ivar else 1.0,
                            n, NOS, byrow = TRUE)           
         for (jay in 1:NOS) {
           mycor <- cov(y[-1, jay], y[-n, jay]) / 
                       apply(y[, jay, drop = FALSE], 2, var)
           init.smn[ , jay] <- mean(y[, jay]) * (1 - mycor) 
           if (!length( .irho ))
             init.rho[, jay] <-  sign(mycor) * min(0.95, abs(mycor))
           if (!length( .ivar ))
             init.var[, jay] <- var(y[, jay]) * (1 -  mycor^2)
           if (!length( .isdv ))
             init.sdv[, jay] <- sqrt(init.var[, jay])
         }  # for
         
         etastart <-
           cbind(if ( .nodrift ) NULL else
             theta2eta(init.smn, .lsmn , earg = .esmn ),
             if ( .var.arg ) 
               theta2eta(init.var, .lvar , earg = .evar ) else
                 theta2eta(init.sdv, .lsdv , earg = .esdv ),
             theta2eta(init.rho, .lrho , earg = .erho ))

         etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
       }
       
     }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
               .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
               .ismn = ismn, .irho = irho, .isdv = isd , .ivar = ivar,
               .type.likelihood = type.likelihood, 
               .var.arg = var.arg, .nodrift = nodrift ))),
     
     
     
     
     linkinv = eval(substitute(function(eta, extra = NULL) {
       n   <- nrow(eta)
       M1  <- 3 - .nodrift
       NOS <- ncol(eta)/M1
       
       ar.smn <- if ( .nodrift ) 0 else
         eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                   .lsmn , earg = .esmn )
       ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                           .lrho , earg = .erho )
       y.lag <- matrix(0, nrow = n, ncol = NOS)
       y.lag[-1, ] <- extra$y[-n, ]
       ar.smn + ar.rho * y.lag
       
     }, list ( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
               .var.arg = var.arg, .type.likelihood = type.likelihood,
               .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ,
               .nodrift = nodrift ))),
     
     
     
     
     last = eval(substitute(expression({
       if (any(abs(ar.rho) > 1)) 
         warning("Regularity conditions are violated at the final",
                 "IRLS iteration, since 'abs(rho) > 1")
       
       M1 <- extra$M1
       
       temp.names <- parameter.names
       temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]
       
       misc$link <- rep( .lrho , length = M1 * ncoly)
       misc$earg <- vector("list", M1 * ncoly)
       names(misc$link) <- names(misc$earg) <- temp.names
       for (ii in 1:ncoly) {
         if ( !( .nodrift ))
           misc$link[ M1*ii-2 ] <- .lsmn
         
         misc$link[ M1*ii-1 ] <- if ( .var.arg ) .lvar else .lsdv
         misc$link[ M1*ii   ] <- .lrho
         
         if ( !( .nodrift ))
           misc$earg[[M1*ii-2]] <- .esmn
         misc$earg[[M1*ii-1]] <- if ( .var.arg ) .evar else .esdv
         misc$earg[[M1*ii  ]] <- .erho
       }
       
      #if (( .poratM ) && any(flag.1) ) 
      #   warning("\nExact EIM approach currently implemented for",
      #           " AR(1) processes \n with stationary noise. Shifting",
      #           " to the 'approximate' EIM approach.")
       
       misc$M1 <- M1
       misc$var.arg   <- .var.arg
       misc$expected  <- TRUE
       misc$nodrift   <- .nodrift
       misc$poratM    <- .poratM
       misc$print.EIM <- FALSE
       misc$type.likelihood   <- .type.likelihood
       misc$multipleResponses <- TRUE
       
     }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
               .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
               .irho = irho, .isdv = isd , .ivar = ivar,
               .nodrift = nodrift, .poratM = poratM, 
               .var.arg = var.arg, .type.likelihood = type.likelihood ))),
     
     
     
     
     loglikelihood = eval(substitute(
       function(mu, y, w, residuals= FALSE, eta, 
                extra = NULL, summation = TRUE) {
         
         M1  <- 3 - .nodrift
         NOS <- ncol(eta)/M1
         
         if ( .var.arg ) {
           ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                               .lvar , earg = .evar )
           ar.sdv <- sqrt(ar.var)
         } else {
           ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                               .lsdv , earg = .esdv )
           ar.var <- ar.sdv^2
         }  
         ar.smn <- if ( .nodrift ) 0 else
           eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                     .lsmn , earg = .esmn )
         ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                             .lrho , earg = .erho )
         
         if (residuals) {
           stop("Loglikelihood not implemented yet to handle",
                "residuals.")
         } else {
           loglik.terms <- 
             c(w) * dAR1extra(x = y,
                         drift = ar.smn ,
                         var.error = ar.var,
                         type.likelihood = .type.likelihood ,
                         ARcoef1 = ar.rho, log = TRUE)
           loglik.terms <- as.matrix(loglik.terms)
           
           if (summation) {
             sum(if ( .type.likelihood == "exact") loglik.terms else
               loglik.terms[-1, ] )
           } else {
             loglik.terms
           }
         }
      }, list( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
               .var.arg = var.arg, .type.likelihood = type.likelihood,
               .nodrift = nodrift,
               .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),
     
     
     vfamily = c("AR1"),
     
     
     validparams = eval(substitute(function(eta, y, extra = NULL) {
       M1    <- 3 - .nodrift
       n     <- nrow(eta)
       NOS   <- ncol(eta)/M1
       ncoly <- ncol(as.matrix(y))
       
       if ( .var.arg ) {
         ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                             .lvar , earg = .evar )
         ar.sdv <- sqrt(ar.var)
       } else {
         ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                             .lsdv , earg = .esdv )
         ar.var <- ar.sdv^2
       }
       ar.smn <- if ( .nodrift ) matrix(0, n, NOS) else
         eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                   .lsmn , earg = .esmn )
       
       ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                           .lrho , earg = .erho )
       okay1 <- all(is.finite(ar.sdv)) && all(0 < ar.sdv) &&
         all(is.finite(ar.smn)) &&  all(is.finite(ar.rho))
       okay1
     }, list( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
              .var.arg = var.arg, .type.likelihood = type.likelihood,
              .nodrift = nodrift,
              .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),
     
     
      
     
     simslot = eval(substitute(function(object, nsim) {
         
         pwts <- if (length(pwts <- object@prior.weights) > 0)
           pwts else weights(object, type = "prior")
         if (any(pwts != 1))
           warning("ignoring prior weights")
         eta <- predict(object)
         fva <- fitted(object)      
         M1  <- 3 - .nodrift
         NOS <- ncol(eta)/M1
         
         if ( .var.arg ) {
           ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                               .lvar , earg = .evar )
           ar.sdv <- sqrt(ar.var)
         } else {
           ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                               .lsdv , earg = .esdv )
           ar.var <- ar.sdv^2
         }  
         ar.smn <- if ( .nodrift ) matrix(0, n, NOS) else
           eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                     .lsmn , earg = .esmn )
         ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                             .lrho , earg = .erho )
         
         ans <- array(0, c(nrow(eta), NOS, nsim))
         for (jay in 1:NOS) {
           ans[1, jay, ] <- ar.smn[1, jay] +
             rnorm(nsim, sd = sqrt(ar.var[1, jay]))
           for (ii in 2:nrow(eta))
             ans[ii, jay, ] <- ar.smn[ii, jay] +
               ar.rho[ii, jay] * ans[ii-1, jay, ] +
               rnorm(nsim, sd = sqrt(ar.var[ii, jay]))
         }
         ans <- matrix(c(ans), c(nrow(eta) * NOS, nsim))
         ans
         
     }, list( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
              .var.arg = var.arg, .nodrift = nodrift,
              .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),
     
     
     
     
     deriv = eval(substitute(expression({
       M1    <- 3 - .nodrift
       NOS   <- ncol(eta)/M1
       ncoly <- ncol(as.matrix(y))
       
       if ( .var.arg ) {
         ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                             .lvar , earg = .evar )
         ar.sdv <- sqrt(ar.var)
       } else {
         ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                             .lsdv , earg = .esdv )
         ar.var <- ar.sdv^2
       }
       
       ar.smn <- if ( .nodrift ) matrix(0, n, NOS) else
             eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                       .lsmn , earg = .esmn )
       
       ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                           .lrho , earg = .erho )
       
       if (any(abs(ar.rho) < 1e-2))
         warning("Estimated values of 'rho' are too close to zero.")
       y.lags  <- apply(y, 2, function(x) WN.lags(y = cbind(x), lags = 1))
       y.means <- y - (ar.smn + ar.rho * y.lags)
       
       dl.dsmn      <- y.means / ar.var
       if ( .var.arg ) {
         dl.dvarSD <-  y.means^2 / ( 2 * ar.var^2) - 1 / (2 * ar.var)
       } else {
         dl.dvarSD <- y.means^2 / ar.sdv^3 - 1 / ar.sdv
       }
       dl.drho <- y.means * y.lags / ar.var
       
       dsmn.deta <- dtheta.deta(ar.smn, .lsmn , earg = .esmn )
       drho.deta <- dtheta.deta(ar.rho, .lrho , earg = .erho )
       if ( .var.arg ) {
         dvarSD.deta <- dtheta.deta(ar.var, .lvar , earg = .evar )
       } else {
         dvarSD.deta <- dtheta.deta(ar.sdv, .lsdv , earg = .esdv )
       }
       
       myderiv <-
         c(w) * cbind(if ( .nodrift ) NULL else dl.dsmn * dsmn.deta,
                      dl.dvarSD * dvarSD.deta, 
                      dl.drho * drho.deta)
       
       myderiv <- myderiv[, interleave.VGAM(M, M1 = M1)]
       
       myderiv
       
     }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
               .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
               .nodrift = nodrift , .var.arg = var.arg, 
               .type.likelihood = type.likelihood ))),
    
     
     
     
     weight = eval(substitute(expression({
       helpPor <- .poratM
       
       ### The EXACT EIMs
       M.fin <- if ( .nodrift ) M1 + (M1 - 1) else
         M1 + (M1 - 1) + (M1 - 2)
       pre.wz <- array(NA_real_, dim = c(n, M.fin , NOS))
       
       for (jj in 1:NOS) {
         pre.wz[, , jj] <- AR1EIM.G2(y  = y[, jj, drop = FALSE],
                                     drift    = ar.smn[, jj, drop = FALSE],
                                     sdError  = ar.sdv[, jj, drop = FALSE],
                                     AR1coeff = ar.rho[, jj, drop = FALSE],
                                     order = 1, nodrift = .nodrift ,
                                     var.arg = .var.arg )
       }
       
       if ( !(.nodrift) ) {
         dTHE.dETA <- cbind(dsmn.deta^2, dvarSD.deta^2, drho.deta^2,
                            matrix(0, nrow = n, ncol = NOS),
                            matrix(0, nrow = n, ncol = NOS),
                            dsmn.deta * drho.deta)
       } else {
         dTHE.dETA <- cbind(dvarSD.deta^2, drho.deta^2,
                            matrix(0, nrow = n, ncol = NOS))
       }
       
       wzExact <- arwzTS(wz = pre.wz, w = w,
                         M1 = M1, dTHE.dETA = dTHE.dETA)
       
       if ( .print.EIM )
         wzPrint1 <- arwzTS(wz = pre.wz, w = w, M1 = M1,
                            dTHE.dETA = dTHE.dETA, print.EIM = TRUE)
       
       ### Approximate EIMs
       pre.wz <- matrix(0, nrow = n,
                        ncol = ifelse( .nodrift, 2 * M - 1, 3 * M - 3))
       
       gamma0 <- ar.var / (1 - ar.rho^2)
       ned2l.dsmn   <- 1 / ar.var 
       ned2l.dvarSD <- if ( .var.arg ) 1 / (2 * ar.var^2) else 2 / ar.var
       ned2l.drho   <- gamma0 / ar.var
       
       if (!( .nodrift ))
         pre.wz[, M1*(1:NOS) - 2] <- ned2l.dsmn * dsmn.deta^2
       
       pre.wz[, M1*(1:NOS) - 1] <- ned2l.dvarSD * dvarSD.deta^2
       pre.wz[, M1*(1:NOS)    ] <- ned2l.drho * drho.deta^2
       wzApp <- w.wz.merge(w = w, wz = pre.wz, n = n,
                           M = ncol(pre.wz), ndepy = NOS)
       
       wz <- if (helpPor) wzExact else wzApp
       
       if ( .print.EIM ) {
         wzEx1 <- matrix(NA_real_, nrow = n, ncol = NOS)
         wzEx2 <- matrix(NA_real_, nrow = n, ncol = NOS)
         wzApp <- wzApp[, 1:M, drop = FALSE]
         wzApp <- array(wzApp, dim = c(n, M / NOS, NOS))
         for (jj in 1:NOS) {
           wzEx1[, jj] <- if (NOS == 1) rowSums(wzPrint1) else 
                          rowSums(wzPrint1[, , jj, drop = FALSE])
           wzEx2[, jj] <- rowSums(wzApp[, , jj, drop = FALSE])
         }
         print.Mat <- cbind(wzEx1, wzEx2)
         colnames(print.Mat) <- c(paste("Exact", 1:NOS), 
                                  paste("Approximate", 1:NOS))
         rownames(print.Mat) <- paste("", 1:n)
         print(head(print.Mat))
       } else {
         rm(wzExact, wzApp)
       }
       
       wz
       
     }), list( .var.arg = var.arg, .type.likelihood = type.likelihood,
               .nodrift = nodrift, .poratM  = poratM, 
               .print.EIM = print.EIM  ))) )
 }
