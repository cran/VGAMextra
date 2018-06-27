##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.


yulesimonTSff <- function(Order = c(1, 1),
                          link = "loge",
                          lagged.fixed.obs   = NULL,
                          lagged.fixed.means = NULL,
                          interventions = list(),
                          init.p.ARMA = NULL,
                          f.transform.Y = NULL,
                          transform.lambda = FALSE,
                          nsimEIM = 200) {
  
  init.p <- init.p.ARMA; rm(init.p.ARMA)
  
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
  
  if (!is.Numeric(Order, length.arg = 2))
    stop("Bad input for argument 'Order'.")
  
  if (length(f.transform.Y) && !is.function(f.transform.Y))
    stop("Wrong input for argument 'f.transform.Y'. Must be a function.",
         " Enter NULL for the identity function.")
  
  if (!is.logical(transform.lambda))
    stop("Argument 'transform.lambda' must be logical.")
  
  interv     <- interventions; rm(interventions)
  fixed.obs  <- unique(lagged.fixed.obs); rm(lagged.fixed.obs)
  fixed.mean <- unique(lagged.fixed.means); rm(lagged.fixed.means)
  
  ord1 <- max(fixed.obs, Order[1])
  ord2 <- max(fixed.mean, Order[2])
  fy   <- f.transform.Y; rm(f.transform.Y)
  flam <- transform.lambda; rm(transform.lambda)
  
  
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

  zero   <- NULL
  ishape <- NULL
  lshape <- link; rm(link)
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  
  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")
  
  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")
  
  
  
  alo <-
    new("vgltsmff",
        blurb = c("VGLM-INGARCH Yule-Simon TS model of ",
                "order - (", Order[1],",",Order[2],")\n\n",
                "shape > 0, y = 1, 2,..\n\n",
                "Link:    ",
                 paste(namesof("lambda", lshape, earg = eshape),
                       sep = ""), ", with lambda = lambda(t).", "\n\n",
                "Uncond. Mean: lambda / (lambda - 1), provided lambda > 1\n",
                "Uncond.  Var: lambda^2 / ((lambda - 1)^2 * (lambda - 2)), ",
                "provided lambda > 2"),
      
      first = eval(substitute(expression({
        
        x.matrix <- x
        counts <- ncol(x)
        my.ord <- .Order
        # intercept.only <- (ncol(x) == 1 && colnames(x) == "(Intercept)")
        #if (!intercept.only)
        #  stop("Currently, this family function only handles ",
        #       "intercept-only models.")
        
        if (any(y < 1) || any(!is.Numeric(y, integer.valued = TRUE)))
          stop("The response 'y' has values out of range, i.e., < 1.")
        
        if (NCOL(y) > 1)
          stop("Currently, only univariate time series handled.")
        
        # Differs from other VGLTSM for count data.
        iniOrd <- if (length( .init.p )) .init.p else
              min(5, .Order[1] + .Order[2] + 1) # Dec 2017. Didn't change
        
        nn <- NROW(y)
        
        temp.pd <- data.frame(y = y, WN.lags(y = cbind(y),
                              to.complete = rep(mean(y), iniOrd),
                                             lags = iniOrd))
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
                          yulesimon(lshape = .link , zero = NULL),
                          trace = FALSE, data = temp.pd, smart = FALSE)
        
        a.help  <- fitted.values(vglm.temp) # Lambda hat
        
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
              tt.mat <- matrix(NA_real_, nrow = nn,
                               ncol = NROW(here.first))
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
                .link = lshape , .interv = interv , .init.p = init.p ,
                .earg = eshape , .fy = fy , .flam = flam ))),
      
      
      constraints = eval(substitute(expression({
        constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                    predictors.names = predictors.names,
                                    M1 = 1)
      }), list( .zero = zero ))),
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 1,
             Q1 = 1,
             expected = TRUE,
             multipleResponses = TRUE,
             nsimEIM = .nsimEIM ,
             parameters.names = c("lambda"),
             zero = .zero )
      }, list( .zero = zero,
               .nsimEIM = nsimEIM ))),
      
      
      initialize = eval(substitute(expression({
        
        
        temp5 <-
          w.y.check(w = w, y = y,
                    Is.positive.y = TRUE,
                    ncol.w.max = Inf,
                    ncol.y.max = Inf,
                    Is.integer.y = TRUE,
                    out.wy = TRUE,
                    colsyperw = 1,
                    maximize = TRUE)
        w <- temp5$w
        y <- temp5$y
        
        
        
        ncoly <- ncol(y)
        
        M1 <- 1
        extra$ncoly <- ncoly
        extra$M1 <- M1
        M <- M1 * ncoly
        
        
        mynames1  <- param.names("lambda", ncoly)
        predictors.names <-
          namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)
        
        if (!length(etastart)) {
          wmeany <- colSums(y * w) / colSums(w) + 1/8
          
          shape.init <- wmeany / (wmeany - 1)
          shape.init <- matrix(if (length( .ishape )) .ishape else
            shape.init, n, M, byrow = TRUE)
          etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
        }
      }), list( .lshape = lshape, .eshape = eshape, .ishape = ishape ))),
      linkinv = eval(substitute(function(eta, extra = NULL) {
        ans <- shape <- eta2theta(eta, .lshape , earg = .eshape )
        
        ans[shape >  1] <- shape / (shape - 1)
        ans[shape <= 1] <- NA
        ans
      }, list( .lshape = lshape, .eshape = eshape ))),
      
      
      last = eval(substitute(expression({
        M1 <- extra$M1
        misc$link <- c(rep_len( .lshape , ncoly))
        names(misc$link) <- mynames1
        
        misc$earg <- vector("list", M)
        names(misc$earg) <- mynames1
        for (ii in 1:ncoly) {
          misc$earg[[ii]] <- .eshape
        }
        
        misc$M1 <- M1
        misc$ishape <- .ishape
        misc$nsimEIM <- .nsimEIM
        misc$predictors <-  eta2theta(eta, .lshape , earg = .eshape )
      }), list( .lshape = lshape, .eshape = eshape, .nsimEIM = nsimEIM,
                .ishape = ishape ))),
      loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL,
                 summation = TRUE) {
          shape <- eta2theta(eta, .lshape , earg = .eshape )
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dyules(x = y, shape = shape, log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list( .lshape = lshape, .eshape = eshape ))),
      
      
      vfamily = c("yulesimon","vgltsmff", "VGLMINGARCH"),
      
      
      validparams = eval(substitute(function(eta, y, extra = NULL) {
        shape <- eta2theta(eta, .lshape , earg = .eshape )
        okay1 <- all(is.finite(shape)) && all(0 < shape)
        okay1
      }, list( .lshape = lshape, .eshape = eshape ))),
      
      
      
      
      
      simslot = eval(substitute(
        function(object, nsim) {
          
          pwts <- if (length(pwts <- object@prior.weights) > 0)
            pwts else weights(object, type = "prior")
          if (any(pwts != 1))
            warning("ignoring prior weights")
          eta <- predict(object)
          shape <- eta2theta(eta, .lshape , earg = .eshape )
          ryules(nsim * length(shape), shape = shape)
        }, list( .lshape = lshape, .eshape = eshape ))),
      
      
      
      
      
      
      
      deriv = eval(substitute(expression({
        M1 <- 1
        shape <- eta2theta(eta, .lshape , earg = .eshape )
        dl.dshape <- 1/shape + digamma(1+shape) - digamma(1+shape+y)
        dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
        c(w) * dl.dshape * dshape.deta
      }), list( .lshape = lshape, .eshape = eshape ))),
      weight = eval(substitute(expression({
        
        run.var <- 0
        for (ii in 1:( .nsimEIM )) {
          ysim <- ryules(n, shape <- shape)
          dl.dshape <- 1/shape + digamma(1+shape) - digamma(1+shape+ysim)
          rm(ysim)
          temp3 <- dl.dshape
          run.var <- ((ii-1) * run.var + temp3^2) / ii
        }
        wz <- if (intercept.only)
          matrix(colMeans(cbind(run.var)),
                 n, M, byrow = TRUE) else cbind(run.var)
        
        wz <- wz * dshape.deta^2
        
        
        c(w) * wz
      }), list( .nsimEIM = nsimEIM ))))
  
  slot(alo, "typeTS") <- "yulesimon"
  
  alo
}



yulesimonTSff.control <- function(save.weights = TRUE, 
                                  summary.HDEtest = FALSE,
                                  ...) {
  
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}