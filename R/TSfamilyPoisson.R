##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.




poissonTSff.control <- function(criterion = "loglikelihood",
                                save.weights = TRUE,
                                summary.HDEtest = FALSE,
                                ...) {
  list(criterion = criterion, 
       save.weights = save.weights,
       summary.HDEtest = summary.HDEtest, ...)
}





poissonTSff <- function(Order = c(1, 1),
                        link = "loge",
                        lagged.fixed.obs   = NULL,
                        lagged.fixed.means = NULL,
                        interventions = list(),
                        init.p.ARMA = NULL,
                        f.transform.Y = NULL,
                        transform.lambda = FALSE,
                        imethod = 1,
                        dispersion = 1) {
  
  ### Default from poissonff, April-2017
  earg.link   <- FALSE
  type.fitted <- c("mean", "quantiles")[1]
  percentiles <- c(25, 50, 75)
  imethod   <- imethod + 1
  init.p    <- init.p.ARMA; rm(init.p.ARMA)
  parallel  <- FALSE
  onedpar   <- FALSE
  bred <- FALSE
  imu  <- NULL
  zero <- NULL
  link <- match.arg(link, c("identitylink", "negloge",
                            "reciprocal", "loge"))[1]
  
  if (length(f.transform.Y) && !is.function(f.transform.Y))
    stop("Wrong input for argument 'f.transform.Y'. Must be a function.",
         " Enter NULL for the identity function.")
  
  if (!is.logical(transform.lambda))
    stop("Argument 'transform.lambda' must be logical.")
  
  fy   <- f.transform.Y; rm(f.transform.Y)
  flam <- transform.lambda; rm(transform.lambda)
  
  if (length(init.p) && !Is.Numeric(init.p, isInteger = TRUE, 
                                    length.arg = 1))
    stop("Bad input for argument 'init.p'")
  
  if (length(interventions)) {
    
    if (any(interventions$tau == 0))
      stop("Bad input for 'tau' in argument 'interventions'.")
    
    if (!is.list(interventions) || (length(interventions) != 3))
      stop("Bad input for argument 'interventions'.")
    
    if (any(interventions$delta < 0) || any(interventions$delta > 1))
      stop("Bad input for 'delta' in argument 'interventions'.")
  }
  
  if (!is.logical(bred) || length(bred) > 1)
    stop("Argument 'bred' must be a single logical")
  
  if (!is.Numeric(Order, length.arg = 2))
    stop("Bad input for argument 'Order'.")
  
  estimated.dispersion <- (dispersion == 0)
  interv     <- interventions; rm(interventions)
  fixed.obs  <- unique(lagged.fixed.obs); rm(lagged.fixed.obs)
  fixed.mean <- unique(lagged.fixed.means); rm(lagged.fixed.means)
  
  ord1 <- max(fixed.obs, Order[1])
  ord2 <- max(fixed.mean, Order[2])
  
  if (length(fixed.obs) && !is.vector(fixed.obs))
    stop("Wrong input for argument 'lagged.fixed.obs'.")
  
  if (length(fixed.mean) && !is.vector(fixed.mean))
    stop("Wrong input for argument 'lagged.fixed.means'.")
  
  if (length(fixed.obs) && 
       !Is.Numeric(fixed.obs, Nnegative = TRUE, isInteger = TRUE))
    stop("Wrong input for argument 'lagged.fixed.obs'")
  
  if (length(fixed.mean) &&
      !Is.Numeric(fixed.mean, Nnegative = TRUE, isInteger = TRUE))
    stop("Wrong input for argument 'lagged.fixed.mean'")
  
  if (earg.link) {
    earg <- link
  } else {
    link <- as.list(substitute(link))
    earg <- link2list(link)
  }
  link <- attr(earg, "function.name")
  
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("Argument 'imethod' must be 1 or 2.")
  if (length(imu) &&
      !is.Numeric(imu, positive = TRUE))
    stop("bad input for argument 'imu'")
  
  
  alo <- new("vgltsmff",
      blurb = c("VGLM-INGARCH Poisson TS model of ",
                "order - (", Order[1],",",Order[2],")\n\n",
                "Link:     ",
                paste(namesof("lambda", link, earg = earg),
                      "(t)", sep = ""), "\n",
                "Uncond. Variance: lambda(t)"),
      
      
    first = eval(substitute(expression({
      
      x.matrix <- x
      counts <- ncol(x)
      my.ord <- .Order
      # intercept.only <- (ncol(x) == 1 && colnames(x) == "(Intercept)")
      #if (!intercept.only)
      #  stop("Currently, this family function only handles ",
      #       "intercept-only models.")
      
      if (NCOL(y) > 2)
        stop("Currently, only univariate time series handled.")
      
      iniOrd <- if (length( .init.p )) .init.p else
        max(25, .ord1 + .ord2 + 1)  # Dec 2017, Old is init.p = 10
      nn <- NROW(y)
      
      temp.pd <- data.frame(y = y, WN.lags(y = cbind(y), 
                                           lags = iniOrd ))
      colnames(temp.pd) <- c("y", paste("x", 2:( iniOrd + 1), 
                                        sep = ""))
      
      myform <- character(0)
      for (ii in 2:(iniOrd + 1)) {
        pre <- paste(paste("x", ii, sep = ""),
                     if (ii == iniOrd + 1) "" else " + ", sep ="")
        myform <- paste(myform, pre, sep = "")
      }
      
      myform <- paste("y ~ 1 +", myform, sep = " ")
      vglm.temp <- vglm(as.formula(myform),
                        poissonff(link = .link , imethod = 2),
                        trace = FALSE, data = temp.pd, smart = FALSE)
      a.help  <- fitted.values(vglm.temp) # Lambda hat
      
      if (FALSE) {
        a.help <- lsfit(x = WN.lags(y = cbind(y), lags = iniOrd),
                        y = y, intercept = TRUE)
        a.help <- cbind(rep(1, nn), WN.lags(y = cbind(y),
                                        lags = iniOrd)) %*% coef(a.help)
      }
      
      x1.mat <- x2.mat <- x3.mat <- tau <- NULL
      
      
      if ( .ord1 ) {
        x1.mat <- WN.lags(y = cbind(y), lags = .ord1 ,
                          to.complete = rep( round(0* mean(y)), .ord1 ))
        colnames(x1.mat) <- paste("Ylag", 1:( .ord1 ), sep = "")
        
        if ( length( .fixed.obs ) )
          x1.mat <- x1.mat[, if (!my.ord[1]) .fixed.obs else 
                                unique( c(1:(my.ord[1]) , .fixed.obs )) ,
                           drop = FALSE]
        
        if (length( .fy )) {
          fy <- .fy
          x1.mat <- if (identical(fy, log)) log1p(x1.mat) else fy(x1.mat)
          colnames(x1.mat) <- 
            if (identical(fy, log) || identical(fy, log1p))
               paste("Log(Ylag + 1)", 1:( .ord1 ), sep = "")  else 
                    paste("f(Ylag)", 1:( .ord1 ), sep = "")
        }
        counts <- length(if (!my.ord[1]) .fixed.obs else 
             unique( c(1:(my.ord[1]) , .fixed.obs ))) + counts
      }
      
      
      
      if ( .ord2 ) {
        x2.mat <- WN.lags(y = cbind(a.help), lags = .ord2 ,
                      to.complete =rep(round(0 * mean(a.help)), .ord2 ))
        
        colnames(x2.mat) <- paste("lambLag", 1:( .ord2 ), sep = "")
        
        if ( length( .fixed.mean ) )
          x2.mat <- x2.mat[, if (!my.ord[2])  .fixed.mean else
                               unique(c(1:(my.ord[2]), .fixed.mean )) ,
                           drop = FALSE]
        
        if ( .flam ) {
          x2.mat <- if (identical( .link, "loge")) 
                     theta2eta(x2.mat + 1, .link , .earg) else
                       theta2eta(x2.mat, .link , .earg)
          if (!identical( .link, "identitylink")) 
            colnames(x2.mat) <- paste(attributes( .earg )$function.name,
                        paste("lambLag", 1:NCOL(x2.mat), sep = ""))
        }
        counts <- length(if (!my.ord[2])  .fixed.mean else
            unique(c(1:(my.ord[2]), .fixed.mean ))) + counts
      } 
      
      
      
      if (length( .interv )) {
        tau <- .interv$tau
        del <- .interv$delta
        noI <- .interv$No.Inter
        sum.ad <- 0
        x3.mat <- matrix(0, nn, ncol = length(tau))
        
        for (kk in seq_len(length(tau))) 
          x3.mat[tau[kk]:nn, kk] <- del[kk]^(0:(nn - tau[kk]))
        
        colnames(x3.mat) <- paste("Interv.", 1:length(tau), sep="")
        
        if (!noI) {
          here.first <- which(del > 0 & del < 1)
          if (!length(here.first))
            stop("Only exponential intervention effects handled ",
                 "currently.")
          
          if (length(here.first) == 1) {
            warning("No interactions terms needed. Only a single ",
                    "'delta' value lies in (0, 1).")
          } else {
            here.first <-rbind(combVGAMextra(x =here.first)[
                                            -(1:length(here.first)), ])
            tt.mat <- matrix(NA_real_, nrow = nn, ncol = NROW(here.first))
            tt.nams <- character(0)
            for (ii in 1:NROW(here.first)) {
              tt.mat[, ii] <- x3.mat[, here.first[ii, ][1]] *
                         x3.mat[, here.first[ii, ][2]]
              tt.nams <- c(tt.nams, 
                           paste(paste("I", here.first[ii, ][1], sep = ""),
                           paste("I", here.first[ii, ][2], sep = ""),
                           sep = ""))
            }
            colnames(tt.mat) <- tt.nams
            sum.ad <- length(tt.nams); rm(tt.nams)
            x3.mat <- cbind(x3.mat, tt.mat); rm(tt.mat)
          }
        } 
      } else {
        sum.ad <- 0
      }
      
      x.matrix <- cbind(x.matrix, x1.mat, x2.mat, x3.mat)
      list.names <- vector("list", counts + length(tau) + sum.ad)
      names(list.names) <- colnames(x.matrix)
      
      for (ii in 1:(counts + length(tau) + sum.ad))
        list.names[[ii]] <- ii
      attr(x.matrix, "assign") <- list.names
      x <- x.matrix
      
      
    }), list( .ord1 = ord1 , .ord2 = ord2 , .Order = Order , 
              .fixed.mean = fixed.mean , .fixed.obs = fixed.obs ,
              .link = link , .earg = earg , .interv = interv ,
              .init.p = init.p , .fy = fy , .flam = flam ))),
      
      
      
      
      
   constraints = eval(substitute(expression({
     #constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
     #                      bool = .parallel , 
     #                     constraints = constraints)
     constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                 predictors.names = predictors.names,
                                 M1 = 1)
     
   }), list( .parallel = parallel, .zero = zero ))),
      
      
      
   infos = eval(substitute(function(...) {
     list(M1 = 1,
          Q1 = 1,
          expected = TRUE,
          multipleResponses = TRUE,
          parameters.names = c("lambda"),
          type.fitted = .type.fitted ,
          percentiles = .percentiles ,
          bred = .bred ,
          zero = .zero )
     
   }, list( .zero = zero,
               .type.fitted = type.fitted,
               .percentiles = percentiles,
               .bred = bred ))),
      
      
   #   deviance = eval(substitute(
  #      function(mu, y, w, residuals = FALSE, eta, extra = NULL,
  #               summation = TRUE) {
  #        mupo <- eta2theta(eta, link = .link , earg = .earg )
  #        nz <- (y > 0)
  #        devi <-  -(y - mupo)
  #        devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mupo[nz])
  #        if (residuals) {
  #          sign(y - mupo) * sqrt(2 * abs(devi) * c(w))
  #        } else {
  #          dev.elts <- 2 * c(w) * devi
  #          if (summation) {
  #            sum(dev.elts)
  #          } else {
  #            dev.elts
  #          }
  #        }
  #      }, list( .link = link, .earg = earg ))),
  
  
      
   initialize = eval(substitute(expression({
     
     temp5 <-
       w.y.check(w = w, y = y,
                 Is.nonnegative.y = TRUE,
                 ncol.w.max = Inf,
                 ncol.y.max = Inf,
                 out.wy = TRUE,
                 colsyperw = 1,
                 maximize = TRUE)
     w <- temp5$w
     y <- temp5$y
     
     
     M <- ncoly <- ncol(y)
     
     assign("CQO.FastAlgorithm", ( .link == "loge"), envir = VGAMenv)
     
     old.name <- "mu"
     new.name <- "lambda"
     dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
     dn2 <- if (length(dn2)) {
       paste("E[", dn2, "]", sep = "")
     } else {
       paste(new.name, 1:M, sep = "")
     }
     predictors.names <-
       namesof(if (M > 1) dn2 else new.name, # was "mu" == old.name
               .link ,
               earg = .earg , short = TRUE)
     
     
     if ( .bred ) {
       if ( !control$save.weights ) {
         save.weights <- control$save.weights <- TRUE
       }
     }
     
     extra$type.fitted <- .type.fitted
     extra$colnames.y  <- colnames(y)
     extra$percentiles <- .percentiles
     
     
     if (!length(etastart)) {
       mu.init <- pmax(y, 1/8)
       for (iii in 1:ncol(y)) {
         if ( .imethod == 2) {
           mu.init[, iii] <- weighted.mean(y[, iii], w[, iii]) + 1/8
         } else if ( .imethod == 3) {
           mu.init[, iii] <- median(y[, iii]) + 1/8
         }
       }
       if (length( .imu ))
         mu.init <- matrix( .imu , n, ncoly, byrow = TRUE)
       etastart <- theta2eta(mu.init, link = .link , earg = .earg )
       
     }
     
     
   }), list( .link = link, .estimated.dispersion = estimated.dispersion,
             .type.fitted = type.fitted,
             .percentiles = percentiles,
             .bred = bred,
             .imethod = imethod, .imu = imu, .earg = earg))),
  
  
   linkinv = eval(substitute(function(eta, extra = NULL) {
     mupo <- eta2theta(eta, link = .link , earg = .earg )
     
     type.fitted <-
       if (length(extra$type.fitted)) {
         extra$type.fitted
       } else {
         warning("cannot find 'type.fitted'. Returning 'mean'.")
         "mean"
       }
     
     type.fitted <- match.arg(type.fitted,
                              c("mean", "quantiles"))[1]
     
     if (type.fitted == "mean") {
       return(label.cols.y(mupo, colnames.y = extra$colnames.y,
                           NOS = NOS))
     }
     
     percvec <- extra$percentiles
     lenperc <- length(percvec)
     NOS <- NCOL(eta) / c(M1 = 1)
     jvec <- lenperc * (0:(NOS - 1))
     ans <- matrix(0, NROW(eta), lenperc * NOS)
     for (kay in 1:lenperc)
       ans[, jvec + kay] <-
       qpois(0.01 * percvec[kay], lambda = mupo)
     
     rownames(ans) <- rownames(eta)
     label.cols.y(ans, colnames.y = extra$colnames.y,
                  NOS = NOS, percentiles = percvec,
                  one.on.one = FALSE)
     
   }, list( .link = link, .earg = earg))),
      
  
  
   last = eval(substitute(expression({
     
     if (exists("CQO.FastAlgorithm", envir = VGAMenv))
       rm("CQO.FastAlgorithm", envir = VGAMenv)
     dpar <- .dispersion
     if (!dpar) {
       temp87 <- (y-mu)^2 *
         wz / (dtheta.deta(mu, link = .link , earg = .earg )^2)  # w cancel
       if (M > 1 && ! .onedpar ) {
         dpar <- rep_len(NA_real_, M)
         temp87 <- cbind(temp87)
         nrow.mu <- if (is.matrix(mu)) nrow(mu) else length(mu)
         for (ii in 1:M)
           dpar[ii] <- sum(temp87[, ii]) / (nrow.mu - ncol(x))
         if (is.matrix(y) && length(dimnames(y)[[2]]) == length(dpar))
           names(dpar) <- dimnames(y)[[2]]
       } else {
         dpar <- sum(temp87) / (length(mu) - ncol(x))
       }
     }
     misc$dispersion <- dpar
     misc$default.dispersion <- 1
     misc$estimated.dispersion <- .estimated.dispersion
     
     misc$imethod <- .imethod
     
     
     misc$link <- rep_len( .link , M)
     names(misc$link) <- if (M > 1) dn2 else new.name  # Was old.name=="mu"
     
     misc$earg <- vector("list", M)
     names(misc$earg) <- names(misc$link)
     for (ii in 1:M)
       misc$earg[[ii]] <- .earg
     
   }), list( .dispersion = dispersion, .imethod = imethod,
             .estimated.dispersion = estimated.dispersion,
             .bred = bred,
             .onedpar = onedpar, .link = link, .earg = earg))),
      
    
  
  
  
  linkfun = eval(substitute( function(mu, extra = NULL) {
        theta2eta(mu, link = .link , earg = .earg )
  }, list( .link = link, .earg = earg))),
      
  
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL, summation = TRUE) {
      
      mupo <- eta2theta(eta, link = .link , earg = .earg )
      if (residuals) {
        c(w) * (y / mupo - 1)
      } else {
        ll.elts <- c(w) * dpois(x = y, lambda = mupo, log = TRUE)
        if (summation) {
          sum(ll.elts)
        } else {
          ll.elts
        }
      }
    
  }, list( .link = link, .earg = earg ))),
  
  
   
  vfamily = c("poissonTSff", "vgltsmff", "VGLMINGARCH"),
  
  
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    mupo <- eta2theta(eta, link = .link , earg = .earg )
    okay1 <- all(is.finite(mupo)) && all(0 < mupo)
    okay1
    
  }, list( .link = link, .earg = earg ))),
      
    
      
      
   simslot = function(object, nsim) {
     pwts <- if (length(pwts <- object@prior.weights) > 0)
       pwts else weights(object, type = "prior")
     if (any(pwts != 1))
       warning("ignoring prior weights")
     ftd <- fitted(object)
     rpois(nsim * length(ftd), ftd)
     
    },
      
      
      
   deriv = eval(substitute(expression({
     mupo <- eta2theta(eta, link = .link , earg = .earg )
     yBRED <- if ( .bred ) {
       Hvector <- hatvaluesbasic(X.vlm = X.vlm.save,
                                 diagWm = c(t(c(w) * mupo)))  # Handles M>1
       
       
       varY <- mupo  # Is a matrix if M>1.
       d1.BRED <-   dtheta.deta(mupo, .link , earg = .earg )
       d2.BRED <- d2theta.deta2(mupo, .link , earg = .earg )
       y + matrix(Hvector, n, M, byrow = TRUE) *
         varY * d2.BRED / (2 * d1.BRED^2)
     } else {
       y
     }
     
     
     answer <- if ( .link == "loge" && (any(mupo < .Machine$double.eps))) {
       c(w) * (yBRED - mupo)
     } else {
       lambda <- mupo
       dl.dlambda <- (yBRED - lambda) / lambda
       dlambda.deta <- dtheta.deta(theta = lambda,
                                   link = .link , earg = .earg )
       c(w) * dl.dlambda * dlambda.deta
     }
     
     answer
     
   }), list( .link = link, .earg = earg, .bred = bred))),
      
  
  
  
   weight = eval(substitute(expression({
     if ( .link == "loge" && (any(mupo < .Machine$double.eps))) {
       tmp600 <- mupo
       tmp600[tmp600 < .Machine$double.eps] <- .Machine$double.eps
       c(w) * tmp600
     } else {
       ned2l.dlambda2 <- 1 / lambda
       ned2lambda.deta2 <- d2theta.deta2(theta = lambda,
                                         link = .link , earg = .earg )
       c(w) * dlambda.deta^2 * ned2l.dlambda2
     }
     
   }), list( .link = link, .earg = earg))))
  
  slot(alo, "typeTS") <- "poisson"
  
  alo
}