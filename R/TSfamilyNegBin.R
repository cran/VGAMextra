##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.



NegBinomTSff <- function(Order = c(1, 1),
                         link = "identitylink",
                         lagged.fixed.obs   = NULL,
                         lagged.fixed.means = NULL,
                         interventions = list(),
                         init.p.ARMA = NULL,
                         f.transform.Y = NULL,
                         transform.lambda = FALSE,
                         imethod = 1) {
  
  init.p    <- init.p.ARMA; rm(init.p.ARMA)
  if (length(init.p) && !Is.Numeric(init.p, isInteger = TRUE, 
                                    length.arg = 1))
    stop("Bad input for argument 'init.p'")
  
  if (length(interventions)) {
    if (!is.list(interventions) || (length(interventions) != 3))
      stop("Bad input for argument 'interventions'.")
    
    if (any(interventions$tau == 0))
      stop("Bad input for 'tau' in argument 'interventions'.")
    
    if (any(interventions$delta < 0) || any(interventions$delta > 1))
      stop("Bad input for 'delta' in argument 'interventions'.")
    
    if (length(interventions$No.Inter) &&
         !is.logical(interventions$No.Inter))
      stop("Wrong input for argument 'No.Inter'.")
    
  }
  
  if (length(f.transform.Y) && !is.function(f.transform.Y))
    stop("Wrong input for argument 'f.transform.Y'. Must be a function.",
         " Enter NULL for the identity function.")
  
  if (!is.logical(transform.lambda))
    stop("Argument 'transform.lambda' must be logical.")
  
  interv   <- interventions; rm(interventions)
  fixed.obs  <- unique(lagged.fixed.obs); rm(lagged.fixed.obs)
  fixed.mean <- unique(lagged.fixed.means); rm(lagged.fixed.means)
  fy   <- f.transform.Y; rm(f.transform.Y)
  flam <- transform.lambda; rm(transform.lambda)

  
  ord1 <- max(fixed.obs, Order[1])
  ord2 <- max(fixed.mean, Order[2])
  
  if (length(fixed.mean) && !is.vector(fixed.mean))
    stop("Wrong input for argument 'lagged.fixed.means'.")
  
  # Defaults, 17-07
  zero <- "size"
  parallel <- FALSE
  deviance.arg <- FALSE
  type.fitted <- c("mean", "quantiles")[1]
  percentiles <- c(25, 50, 75)
  mds.min <- 1e-3
  nsimEIM <- 500
  cutoff.prob <- 0.999  # Maxiter = 5000,
  eps.trig <- 1e-7
  max.support <- 4000
  max.chunk.MB <- 30  # max.memory = Inf is allowed
  imu <- NULL
  iprobs.y <- NULL  # 0.35,
  gprobs.y <- ppoints(6)
  isize <- NULL
  gsize.mux <- exp(c(-30, -20, -15, -10, -6:3))
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  
  
  if (!is.logical( deviance.arg ) || length( deviance.arg ) != 1)
    stop("argument 'deviance.arg' must be TRUE or FALSE")
  
  type.fitted <- match.arg(type.fitted,
                           c("mean", "quantiles"))[1]
  
  lmu   <- link
  lsize <- "identitylink"

  lmunb <- as.list(substitute(lmu))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")
  
  imunb <- imu
  
  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")
  
  
  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 1e-5)
    stop("argument 'eps.trig' must be positive and smaller in value")
  
  if (length(imunb) && !is.Numeric(imunb, positive = TRUE))
    stop("bad input for argument 'imu'")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("bad input for argument 'isize'")
  
  if (!is.Numeric(cutoff.prob, length.arg = 1) ||
      cutoff.prob < 0.95 || cutoff.prob >= 1)
    stop("range error in the argument 'cutoff.prob'; ",
         "a value in [0.95, 1) is needed")
  
  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
    stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 10)
    warning("argument 'nsimEIM' should be an integer ",
            "greater than 10, say")
  
  
  if (is.logical( parallel ) && parallel  && length(zero))
    stop("need to set 'zero = NULL' when parallel = TRUE")
  
  
  
  ans <-
    new("vgltsmff",
        blurb = c("VGLM-INGARCH Negative binomial TS model of ",
                  "order - (", Order[1],",",Order[2],")\n\n",
                  "Links:    ",
                  namesof("lambda",   lmunb, earg = emunb), ", ",
                  namesof("size", lsize, earg = esize), "\n\n",
                  "Unc. Mean:     lambda",  
                  ", where lambda = lambda(t)", "\n",
                  "Unc.  Var: lambda * (1 + lambda / size) for NB-2"),
        
        
        
        
    first = eval(substitute(expression({
      
      x.matrix <- x
      counts <- ncol(x)
      my.ord <- .Order
      
      if (NCOL(y) > 2)
        stop("Currently, only univariate time series handled.")
      
      iniOrd <- if (length( .init.p )) .init.p else
        max(8, .ord1 + .ord2 + 1)  # Dec 2017, Old is init.p = 10
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
                        negbinomial(lmu = .lmu ),
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
                     to.complete = rep(round(0 * mean(a.help)), .ord2 ))
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
      }  else {
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
              .link = lmu , .interv = interv , .init.p = init.p ,
              .earg = emunb , .fy = fy , .flam = flam , .lmu = lmu ))),
        
        
        
    constraints = eval(substitute(expression({
      #constraints <- cm.VGAM(matrix(1, M, 1), x = x,
      #                       bool = .parallel ,
      #                       constraints = constraints)
      
      constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                  predictors.names = predictors.names,
                                  M1 = 2)
      
    }), list( .parallel = parallel, .zero = zero ))),
        
        
        
     infos = eval(substitute(function(...) {
       list(M1    = 2,
            Q1    = 1,
            expected = TRUE,
            imethod = .imethod ,
            mds.min = .mds.min ,
            multipleResponses = TRUE,
            parameters.names = c("mu", "size"),
            type.fitted = .type.fitted ,
            percentiles = .percentiles ,
            lmu   = .lmunb ,
            lsize = .lsize ,
            nsimEIM = .nsimEIM ,
            eps.trig = .eps.trig ,
            zero  = .zero ,
            max.chunk.MB = .max.chunk.MB ,
            cutoff.prob = .cutoff.prob)
       
     }, list(.zero = zero, .lsize = lsize, .lmunb = lmunb,
             .type.fitted = type.fitted,
             .percentiles = percentiles ,
             .eps.trig = eps.trig,
             .imethod = imethod,
             .mds.min = mds.min,
             .cutoff.prob = cutoff.prob,
             .max.chunk.MB = max.chunk.MB,
             .nsimEIM = nsimEIM ))),
      
    
    
    initialize = eval(substitute(expression({
      M1 <- 2
      
      temp12 <-
        w.y.check(w = w, y = y,
                  Is.nonnegative.y = TRUE,
                  Is.integer.y = TRUE,
                  ncol.w.max = Inf,
                  ncol.y.max = Inf,
                  out.wy = TRUE,
                  colsyperw = 1, maximize = TRUE)
      w <- temp12$w
      y <- temp12$y
      
      
      assign("CQO.FastAlgorithm",
             ( .lmunb == "loge") && ( .lsize == "loge"),
             envir = VGAMenv)
      
      if (any(function.name == c("cqo", "cao")) &&
          ((is.Numeric( .zero , length.arg = 1) && .zero != -2) ||
           (is.character( .zero ) && .zero != "size")))
        stop("argument zero = 'size' or zero = -2 is required")
      
      
      extra$type.fitted <- .type.fitted
      extra$percentiles <- .percentiles
      extra$colnames.y  <- colnames(y)
      NOS <- ncoly <- ncol(y)  # Number of species
      M <- M1 * NOS
      predictors.names <-
        c(namesof(param.names("mu",   NOS),
                  .lmunb , earg = .emunb , tag = FALSE),
          namesof(param.names("size", NOS),
                  .lsize , earg = .esize , tag = FALSE))
      predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]
      
      gprobs.y <- .gprobs.y
      imunb <- .imunb  # Default is NULL
      if (length(imunb))
        imunb <- matrix(imunb, n, NOS, byrow = TRUE)
      
      if (!length(etastart)) {
        munb.init <-
          size.init <- matrix(NA_real_, n, NOS)
        if (length( .iprobs.y ))
          gprobs.y <-  .iprobs.y
        gsize.mux <- .gsize.mux  # gsize.mux is on a relative scale
        
        for (jay in 1:NOS) {  # For each response 'y_jay'... do:
          wm.yj <- weighted.mean(y[, jay], w = w[, jay])
          munb.init.jay <- if ( .imethod == 1 ) {
            negbinomial.initialize.yj(y[, jay], w[, jay],
                                      gprobs.y = gprobs.y,
                                      wm.yj = wm.yj)
          } else {
            wm.yj
          }
          if (length(imunb))
            munb.init.jay <- imunb[, jay]
          
          
          gsize <- gsize.mux * 0.5 * (mean(munb.init.jay) +
                                        wm.yj)
          if (length( .isize ))
            gsize <- .isize  # isize is on an absolute scale
          
          
          try.this <-
            grid.search2(munb.init.jay, gsize,
                         objfun = NBD.Loglikfun2,
                         y = y[, jay], w = w[, jay],
                         ret.objfun = TRUE)  # Last value is the loglik
          
          munb.init[, jay] <- try.this["Value1"]
          size.init[, jay] <- try.this["Value2"]
        }  # for (jay ...)
        
        
        
        newemu <- .emunb
        if ( .lmunb == "nbcanlink") {
          newemu$size <- size.init
        }
        
        etastart <-
          cbind(theta2eta(munb.init, link = .lmunb , earg = newemu ),
                theta2eta(size.init, link = .lsize , earg = .esize ))
        
        
        if ( .lmunb == "nbcanlink") {
          if (any(cond1 <- is.na(etastart[, c(TRUE, FALSE)])) ||
              any(cond2 <- etastart[, c(TRUE, FALSE)] >= 0))
            etastart[c(cond1) || c(cond2), c(TRUE, FALSE)] <- -0.1
        }
        
        
        if (M > M1)
          etastart <-
          etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
      }
      
    }), list(.lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize,
             .imunb = imunb,
             .gprobs.y = gprobs.y, .gsize.mux = gsize.mux,
             .deviance.arg = deviance.arg,
             .isize = isize, .iprobs.y = iprobs.y,
             .nsimEIM = nsimEIM,
             .zero = zero, .imethod = imethod,
             .type.fitted = type.fitted,
             .percentiles = percentiles ))),
    
    
    linkinv = eval(substitute(function(eta, extra = NULL) {
      NOS <- ncol(eta) / c(M1 = 2)
      kmat <- NULL
      
      munb <- if ( .lmunb == "nbcanlink") {
        
        newemu <- .emunb
        kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                          .lsize , earg = .esize )
        
        
        munb <- kmat / expm1(-eta[, c(TRUE, FALSE), drop = FALSE])
        munb
      } else {
        eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                  .lmunb , earg = .emunb )
      }
      
      type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
        warning("cannot find 'type.fitted'. ",
                "Returning the 'mean'.")
        "mean"
      }
      type.fitted <- match.arg(type.fitted,
                               c("mean", "quantiles"))[1]
      if (type.fitted == "mean") {
        return(label.cols.y(munb, colnames.y = extra$colnames.y,
                            NOS = NOS))
      }
      
      
      if (is.null(kmat))
        kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                          .lsize , earg = .esize )
      percvec <- extra$percentiles
      lenperc <- length(percvec)
      jvec <- lenperc * (0:(NOS - 1))
      ans <- matrix(0, nrow(eta), lenperc * NOS)
      for (kay in 1:lenperc)
        ans[, jvec + kay] <-
        qnbinom(0.01 * percvec[kay], mu = munb, size = kmat)
      
      rownames(ans) <- rownames(eta)
      
      
      label.cols.y(ans, colnames.y = extra$colnames.y,
                   NOS = NOS, percentiles = percvec,
                   one.on.one = FALSE)
      
    }, list( .lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize))),
    
    
    
    last = eval(substitute(expression({
    
      if (exists("CQO.FastAlgorithm", envir = VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAMenv)
      
      
      if (function.name == "cao")
        ind2 <- FALSE
      
      
      save.weights <- control$save.weights <- !all(ind2)
      
      
      temp0303 <- c(rep_len( .lmunb , NOS),
                    rep_len( .lsize , NOS))
      names(temp0303) <- c(param.names("mu",   NOS),
                           param.names("size", NOS))
      misc$link <- temp0303[interleave.VGAM(M, M1 = M1)] # Already named
      
      misc$earg <- vector("list", M)
      names(misc$earg) <- names(misc$link)
      for (ii in 1:NOS) {
        misc$earg[[M1*ii-1]] <- newemu
        misc$earg[[M1*ii  ]] <- .esize
      }
      
      misc$additional.coe <- kmat # fit$coefficients[2]
      
    }), list( .lmunb = lmunb, .lsize = lsize,
              .emunb = emunb, .esize = esize ))),
    
    
        
    linkfun = eval(substitute(function(mu, extra = NULL) {
      M1 <- 2
      
      newemu <- .emunb
      
      
      if ( .lmunb == "nbcanlink") {
        newemu$size <- eta2theta(eta.size, .lsize , earg = .esize )
      }
      
      
      
      eta.munb <- theta2eta(mu, .lmunb , earg = newemu)
      eta.size <- theta2eta(if (is.numeric( .isize )) .isize else 1.0,
                            .lsize , earg = .esize )
      eta.size <- 0 * eta.munb + eta.size  # Right dimension now.
      
      
      eta.temp <- cbind(eta.munb, eta.size)
      eta.temp[, interleave.VGAM(ncol(eta.temp), M1 = M1), drop = FALSE]
      
    }, list( .lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize,
             .isize = isize ))),
    
    
    
    loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, 
               eta, extra = NULL, summation = TRUE) {
        
    vecTF <-  c(TRUE, FALSE)
    kmat <- eta2theta(eta[, !vecTF, drop=FALSE], .lsize , earg = .esize )
         
        
    munb <- if ( .lmunb == "nbcanlink") {
          
          munb <- kmat / expm1(-eta[, vecTF, drop = FALSE])
          
          if (min(munb) <= 0) {
            munb[munb <= 0] <- median(munb[munb > 0])  # 0.1
            warning("'munb' has some negative values. ",
                    "Using a temporary fix.")
          }
          
          munb
        } else {
          eta2theta(eta[, vecTF, drop = FALSE], .lmunb , earg = .emunb )
        }
        
        
        if (residuals) {
          stop("loglikelihood residuals not implemented yet")
        } else {
          ll.elts <- c(w) * dnbinom(x = y, mu = munb, size = kmat,
                                    log = TRUE)
          if (summation) {
            sum(ll.elts)
          } else {
            ll.elts
          }
        }
        
      
    }, list( .lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize))),
        
        
        
    vfamily = c("NegBinomTSff", "vgltsmff", "VGLMINGARCH"),
        
        
        
    simslot = eval(substitute(function(object, nsim) {
      
      pwts <- if (length(pwts <- object@prior.weights) > 0)
        pwts else weights(object, type = "prior")
      if (any(pwts != 1))
        warning("ignoring prior weights")
      eta <- predict(object)
      vecTF <- c(TRUE, FALSE)
      munb <- cbind(eta2theta(eta[,  vecTF], .lmunb , earg = .emunb ))
      size <- cbind(eta2theta(eta[, !vecTF], .lsize , earg = .esize ))
      rnbinom(nsim * length(munb), mu = munb, size = size)
      
    }, list( .lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize ))),
        
        
    
    
    
    validparams = eval(substitute(function(eta, y, extra = NULL) {
      size <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                        .lsize , earg = .esize )
      
      munb <- if ( .lmunb == "nbcanlink") {
        munb <- size / expm1(-eta[, c(TRUE, FALSE), drop = FALSE])
        munb
      } else {
        eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                  .lmunb , earg = .emunb )
      }
      
      smallval <- .mds.min  # .munb.div.size
      okay1 <- all(is.finite(munb)) && all(0 < munb) &&
        all(is.finite(size)) && all(0 < size)
      
      
      okay0 <- if ( .lmunb == "nbcanlink")
        all(eta[, c(TRUE, FALSE)] < 0) else TRUE
      
      
      overdispersion <- if (okay1) all(smallval < munb / size) else FALSE
      if (!overdispersion)
        warning("parameter 'size' has very large values relative ",
                "to 'mu'; ",
                "try fitting a quasi-Poisson ",
                "model instead.")
      okay1 && overdispersion && okay0
    }, list( .lmunb = lmunb, .emunb = emunb,
             .lsize = lsize, .esize = esize,
             .mds.min = mds.min))),
        
        
        
    deriv = eval(substitute(expression({
      
      if (iter == 1 && .deviance.arg ) {
        if (control$criterion != "coefficients" &&
            control$half.step)
       warning("Argument 'criterion' should be 'coefficients' or ",
              "'half.step' should be 'FALSE' when 'deviance.arg = TRUE'")
        
        
        low.index <- ifelse(names(constraints)[1] == "(Intercept)", 2, 1)
        if (low.index <= length(constraints))
          for (iii in low.index:length(constraints)) {
            conmat <- constraints[[iii]]
            if (any(conmat[c(FALSE, TRUE), ] != 0))
              stop("argument 'deviance.arg' should only be TRUE for NB-2 ",
                   "models; ",
                   "non-zero elements detected for the 'size' parameter." )
          }  # for iii
      }  # (iter == 1 && .deviance.arg )

      vecTF <- c(TRUE, FALSE)
      M1 <- 2
      NOS <- ncol(eta) / M1
      kmat <- eta2theta(eta[, !vecTF, drop = FALSE],
                        .lsize , earg = .esize )

      munb <- if ( .lmunb == "nbcanlink") {
        
        munb <- kmat / expm1(-eta[, vecTF, drop = FALSE])
        
        if (iter <= 2 && min(munb) <= 0) {
          munb[munb <= 0] <- median(munb[munb > 0])
          warning("'munb' has some negative values. ",
                  "Using a temporary fix.")
        }
        
        munb
      } else {
        eta2theta(eta[, vecTF, drop = FALSE], .lmunb , earg = .emunb )
      }
      smallval <- .mds.min  # Something like this is needed
      if (any(big.size <- (munb / kmat < smallval))) {
        kmat[big.size] <- munb[big.size] / smallval
      }
      
      
      dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat)
      dl.dsize <- digamma(y + kmat) - digamma(kmat) +
        log1p(-munb / (kmat + munb)) - (y - munb) / (munb + kmat)
      if (any(big.size)) {
        dl.dsize[big.size] <- 1e-8  # A small number
      }
      
      
      dsize.deta2 <- dtheta.deta(kmat, .lsize , earg = .esize )
      
      
      dmunb.deta1 <- if ( .lmunb == "nbcanlink") {
        
        dl.dsize <- digamma(y + kmat) - digamma(kmat) +
          log1p(-munb / (kmat + munb))
        
        dmunb.deta1 <- nbcanlink(munb, size = kmat, wrt.param = 1,
                                 inverse = TRUE, deriv = 1)
        dmunb.deta1
      } else {
        dtheta.deta(munb, .lmunb , earg = .emunb )
      }
      
      myderiv <- c(w) * cbind(dl.dmunb * dmunb.deta1,
                              dl.dsize * dsize.deta2)
      if (M > M1)
        myderiv <- myderiv[, interleave.VGAM(M, M1 = M1)]
      myderiv
      
    }), list( .lmunb = lmunb, .lsize = lsize,
              .emunb = emunb, .esize = esize,
              .deviance.arg = deviance.arg, .mds.min = mds.min ))),
     
    
      
     
   weight = eval(substitute(expression({
     wz <- matrix(NA_real_, n, M)
     
     max.support <- .max.support
     max.chunk.MB <- .max.chunk.MB
     
     
     ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
     for (jay in 1:NOS) {
       eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
       Q.mins <- 0
       Q.maxs <- round(qnbinom(p = eff.p[2],
                               mu = munb[, jay],
                               size = kmat[, jay]) * 1.1) + 30
       
       eps.trig <- .eps.trig
       Q.MAXS <- if ( .lsize == "loge")
         pmax(10, ceiling(kmat[, jay] / sqrt(eps.trig))) else Inf
       Q.maxs <- pmin(Q.maxs, Q.MAXS)
       
       ind1 <- if (max.chunk.MB > 0)
         (Q.maxs - Q.mins < max.support) else FALSE
       if ((NN <- sum(ind1)) > 0) {
         Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
         n.chunks <- if (intercept.only) 1 else
           max(1, ceiling( Object.Size / max.chunk.MB))
         chunk.rows <- ceiling(NN / n.chunks)
         ind2[, jay] <- ind1  # Save this
         wind2 <- which(ind1)
         
         upr.ptr <- 0
         lwr.ptr <- upr.ptr + 1
         while (lwr.ptr <= NN) {
           upr.ptr <- min(upr.ptr + chunk.rows, NN)
           sind2 <- wind2[lwr.ptr:upr.ptr]
           wz[sind2, M1*jay] <-
             EIM.NB.specialp(mu          = munb[sind2, jay],
                             size        = kmat[sind2, jay],
                             y.max = max(Q.maxs[sind2]),
                             cutoff.prob = .cutoff.prob ,
                             intercept.only = intercept.only)
           
           
           if (any(eim.kk.TF <- wz[sind2, M1*jay] <= 0)) {
             ind2[sind2[eim.kk.TF], jay] <- FALSE
           }
           
           lwr.ptr <- upr.ptr + 1
         }  # while
       }  # if
     }  # end of for (jay in 1:NOS)
     
     for (jay in 1:NOS) {
       run.varcov <- 0
       ii.TF <- !ind2[, jay]  # Not assigned above
       if (any(ii.TF)) {
         kkvec <- kmat[ii.TF, jay]
         muvec <-   munb[ii.TF, jay]
         for (ii in 1:( .nsimEIM )) {
           ysim <- rnbinom(sum(ii.TF), mu = muvec, size = kkvec)
           dl.dsize <- digamma(ysim + kkvec) - digamma(kkvec) -
             (ysim - muvec) / (muvec + kkvec) +
             log1p( -muvec / (kkvec + muvec))
           run.varcov <- run.varcov + dl.dsize^2
         }  # for ii
         
         run.varcov <- c(run.varcov / .nsimEIM )
         ned2l.dsize2 <- if (intercept.only)
           mean(run.varcov) else run.varcov
         
         wz[ii.TF, M1*jay] <- ned2l.dsize2
       }  # any ii.TF
     }  # for jay
     
     save.weights <- !all(ind2)
     
     
     ned2l.dmunb2 <- 1 / munb - 1 / (munb + kmat)
     wz[, c(TRUE, FALSE)] <- ned2l.dmunb2 * dmunb.deta1^2
     
     
     if ( .lmunb == "nbcanlink") {
       wz[, !vecTF] <- wz[, !vecTF] + 1 / kmat - 1 / (kmat + munb)
     }
     wz[, !vecTF] <- wz[, !vecTF] * dsize.deta2^2
     
     if ( .lmunb == "nbcanlink") {
       ned2l.dmunb.dsize <- 1 / (munb + kmat)
       wzoffdiag <- ned2l.dmunb.dsize * dmunb.deta1 * dsize.deta2
       wz <- if (M > M1) {
         wzoffdiag <- kronecker(wzoffdiag, cbind(1, 0))
         cbind(wz, wzoffdiag[, -ncol(wzoffdiag)])
       } else {
         cbind(wz, wzoffdiag)
       }
     }
     
     w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
     
   }), list( .cutoff.prob = cutoff.prob,
             .max.support = max.support,
             .max.chunk.MB = max.chunk.MB,
             .lmunb = lmunb, .lsize = lsize,
             .eps.trig = eps.trig,
             .nsimEIM = nsimEIM ))))
  
  
  
  
  
  if (deviance.arg) {
    ans@deviance <- eval(substitute(
      function(mu, y, w, residuals = FALSE, eta, extra = NULL,
               summation = TRUE) {
        
        
        
        
        
        
        size <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                          .lsize , earg = .esize )
        
        if (residuals) {
          stop("this part of the function has not been written yet.")
        } else {
          dev.elts <- 2 * c(w) *
            (y * log(pmax(1, y) / mu) -
               (y + size) * log((y + size) / (mu + size)))
          if (summation) {
            sum(dev.elts)
          } else {
            dev.elts
          }
        }
      }, list( .lsize = lsize, .esize = esize,
               .lmunb = lmunb, .emunb = emunb )))

  }
  
  slot(ans, "typeTS") <- "negbinomial"
  ans
}


NegBinomTSff.control <- function(save.weights = TRUE, 
                                 summary.HDEtest = FALSE,...) {
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}