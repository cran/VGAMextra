##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0


ARIMAX.errors.ff <- function(order = c(1, 1, 1), 
                             zero  = "var",         # optionally, "mean".
                             order.trend = 0,
                             include.int = TRUE,
                             diffCovs  = TRUE, 
                             xLag = 0,
                             include.currentX = TRUE,
                             lvar = "loglink",
                             lmean = "identitylink") {
  
  
  nodrift <-  FALSE
  imethod <- 1
  apply.parint <- FALSE
  var.arg <- TRUE
  lsd <- "loglink"
  isd <- NULL
  #zero <- "var",  # optionally, "mean".
  parallel <- FALSE
  smallno <- 1.0e-5
  toc <- !(imethod) # Set toc = TRUE for method two. 
  include.current <- include.currentX
  rm(include.currentX)
  
  if (!is.logical(include.current))
    stop("'include.current' must be logical.")
  
  if (!is.logical(diffCovs))
    stop("'diffCovs' must be logical.")
  
  if (!Is.Numeric(order.trend, isInteger = TRUE, Nnegative = TRUE))
    stop("Bad input for argument 'order.trend'.")
  
  if (!is.logical(include.int))
    stop("'include.int' must be logical.")
  
  if (any(xLag < 0))
    stop("Wrong input for argument 'xLag'")
  
  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")
  
  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")
  
  lvare <- as.list(substitute(lvar))
  evare <- link2list(lvare)
  lvare <- attr(evare, "function.name")

  ord <- order; rm(order)
  
  if (!is.Numeric(smallno, length.arg = 1,
                  positive = TRUE))
    stop("argument 'smallno' must be positive and close to 0")
  if (smallno > 0.1) {
    warning("replacing argument 'smallno' with 0.1")
    smallno <- 0.1
  }
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1 or 2 or 3 or 4")
  
  if (!is.logical(var.arg) ||
      length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")
  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")
  
  
  if (is.logical(parallel) && parallel && length(zero))
    stop("set 'zero = NULL' if 'parallel = TRUE'")
  
  
  new("vglmff",
      blurb = c("Generalized dynamic regression with non-seasonal ARIMA errors \n\n",
                "Links:    ",
                namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
                if (var.arg)
                  namesof("var",  lvare, earg = evare, tag = TRUE) else
                    namesof("sd" ,  lsdev, earg = esdev, tag = TRUE),
                "\n",
                if (var.arg) "Variance: var" else "Variance: sd^2"),
      
      
      
      constraints = eval(substitute(expression({
        
        constraints <- 
                cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                             predictors.names = predictors.names, M1 = 2)
        
        ord.int <- .ord
        i.int <- .include.int
        p.int <- ord.int[1]
        d.int <- ord.int[2]
        q.int <- ord.int[3]
        int.zero <- .zero
        col.x <- c(grep("coeff", colnames(x)))
        constra.temp <- colnames(x)[c(col.x)]
        
        
        if (length(int.zero)) {
          
          if ((length(int.zero) == 1) && (int.zero == "mean") && 
                                                      length(col.x)) {
            for (ii in constra.temp)
              constraints[[ii]] <- cbind(c(1, 0))
          }
          
          if ((int.zero == "mean") &&  length(col.x) && d.int) 
            warning("Differenced covariates have been included",
                    " in the variance model.")
          
        } else {
          
          constra.temp <- colnames(x)[col.x]
          
          for (ii in constra.temp)  
            constraints[[ii]] <- cbind(c(1, 0))
        }
        
        if (!i.int) 
          constraints[["(Intercept)"]] <- cbind(c(0, 1))
        
      }), list( .zero = zero, .ord = ord , 
                .include.int = include.int ,
                .parallel = parallel, .apply.parint = apply.parint ))),
      
      
      
      
      first = eval(substitute(expression({
        
        ord.int <- .ord
        p.int <- ord.int[1]
        d.int <- ord.int[2]
        q.int <- ord.int[3]
        extra$y.r <- cbind(y)
        nn  <- NROW(y)
        NOS <- NCOL(y)
        x.errU <- z.Resid <- NULL
        pq.max <- max(c(p.int, q.int))
        
        if (NCOL(y) > 1)
          stop("Only univariate responses handled at the moment.")
        
        if ((NCOL(x) == 1) && (colnames(x) == "(Intercept)"))
          stop("No covariates entered. Refer to ARIMAXff() to fit", 
               " ARIMA models with no covariates.")
        
        if (!d.int) {
          
          if ( !( .include.current ) && ( any( .xLag == 0 ) ) )
            stop("Contradictory entries 'include.current' and 'xLag' ")
          
          x.int <- x[, -1, drop = FALSE]
          int.lags <- rep( .xLag , NCOL(x.int))[1:NCOL(x.int)]
          
          temp4 <- NULL
          temp5 <- which(int.lags > 0)
          
          if (length(temp5)) {
            sub.x <- x.int[, temp5, drop = FALSE]
            sub.lag <- int.lags[temp5]
            
            for (ii in 1:NCOL(sub.x)) {
              save.this <- WN.lags(cbind(sub.x[, ii]), 
                                   lags = sub.lag[ii])
              colnames(save.this) <- 
                paste(colnames(x[, ii +1, drop = FALSE]), "Lag",
                      1:sub.lag[ii], sep ="")
              
              temp4 <- cbind(temp4, save.this)
            }
          }
          
          if ( .include.current ) {
            x.int <- cbind(x, temp4)  
          } else {
            x.int <- cbind(x[, 1], temp4)
            colnames(x.int) <- c("(Intercept)", colnames(x.int)[-1])
          }
          
          
          i.int <- ( .include.int ) 
          # 'x.int' already has the intercept column
          u.res <- 
            if (i.int) residuals(lm(y ~ x.int[, -1, drop = FALSE])) else
                          residuals(lm(y ~ x.int[, -1, drop = FALSE] - 1)) 
          
          extra$my.coefs <- coef(lm(y ~ x.int[, -1, drop = FALSE] ))
          
          if (p.int) {
            x.errU <- WN.lags(cbind(u.res), p.int)[-c(1:pq.max), ,
                                                   drop = FALSE]
            colnames(x.errU) <- paste("ARcoeff", 1:p.int, sep = "")
          }
          
          if (q.int) {
            LagsX <- max(12, sum(ord.int))  # May be an argument later 20180316
           # u.res <- scale(u.res, scale = FALSE)
            if ( !( .toc ) ) {
              u.resL   <- WN.lags(cbind(u.res), lags = LagsX)
              z.Resid  <- residuals(lm(u.res ~  -1 + u.resL))
              
              z.Resid  <- WN.lags(cbind(z.Resid), 
                              lags = q.int)[-c(1:pq.max), , drop = FALSE]
              rm(u.resL)
              colnames(z.Resid) <- paste("MAcoeff", 1:q.int, sep = "")
            }
          }
          
          x.matrix <- x.test <- cbind(x.int[-c(1:pq.max),  , drop = FALSE],
                                      x.errU, z.Resid)
          i.trend <- .order.trend
          
          if ( i.trend ) {
            rem.this <- NULL
            for (ii in 1:i.trend) 
              rem.this <- cbind(rem.this, 
                                c((1 + c(1:pq.max)):nn)^(ii))
            x.matrix <- cbind(x.matrix[, 1],
                              rem.this, 
                              x.matrix[, -1, drop = FALSE])
            rm(rem.this)
            colnames(x.matrix) <- c(colnames(x.test)[1], 
                                    paste("trend", 1:i.trend, sep = ""),
                                    colnames(x.test)[-1])
          } else {
            rm(x.test)
          }
          
          list.names <- vector("list", NCOL(x.matrix) )
          names(list.names) <- colnames(x.matrix)
          for (ii in 1:NCOL(x.matrix)) 
            list.names[[ii]] <- ii
          attr(x.matrix, "assign") <- list.names
          
          x <- x.matrix
          y <- extra$y <- y[-c(1:pq.max)] #- mean(y[-c(1:pq.max)])
          w <- w[-c(1:pq.max)]
          
          
      } else {
          
          y <- cbind(y)
          
          mat.save <- matrix(0, nn, d.int)
          
          for (ii in 1:NOS) {
            yaux1 <- y[, ii , drop = FALSE]
            for (jj in 1:( d.int )) {
              mat.save[c(1:(nn - jj)), jj] <-  diff(yaux1)
              yaux1 <- mat.save[c(1:(nn - jj)), jj, drop = FALSE] 
            }
          }
          extra$mat.save <- mat.save
          
          if ( !( .include.current ) && ( any( .xLag == 0 ) ) )
            stop("Contradictory entries 'include.current' and 'xLag' ")
          
          y.diff <-  apply(cbind(y), 2, function(x) {
            diff(x, differences = d.int)
          })
          
          x.int <- x[, -1, drop = FALSE]
          difxs <- .diffCovs
          
          if (difxs) {
            x.diff <-  apply(cbind(x.int), 2, function(x) {
              diff(x, differences = d.int)
            })
            
            colnames(x.diff) <- paste("Diff", colnames(x.int), sep = "")
          } else{
            x.diff <- x.int[ c(1:NROW(y.diff)), , drop = FALSE]
          }
          
          u.res <- residuals(lm(y.diff ~ -1 + x.diff))
          extra$u.res <- u.res
          if (p.int) {
            x.errU <- WN.lags(cbind(u.res), p.int)
            colnames(x.errU) <- paste("ARcoeff", 1:p.int, sep = "")
          }
          
          int.lags <- rep( .xLag , NCOL(x.int))[1:NCOL(x.int)]
          temp4 <- NULL
          temp5 <- which(int.lags > 0)
          
          if (length(temp5)) {
            sub.x <- x.diff[, temp5, drop = FALSE]
            sub.lag <- int.lags[temp5]
            
            for (ii in 1:NCOL(sub.x)) {
              save.this <- WN.lags(cbind(sub.x[, ii]), 
                                   lags = sub.lag[ii])
              colnames(save.this) <- 
                paste(if (difxs) "Diff" else "",
                      colnames(x[, ii +1, drop = FALSE]), "Lag",
                      1:sub.lag[ii], sep ="")
              temp4 <- cbind(temp4, save.this)
            }
          }
          
          x.int <- cbind(1, if ( .include.current ) x.diff else NULL, temp4)
          save.this <- x.int # cbind(1, x.diff)
          colnames(save.this) <- c("(Intercept)", 
                                  if ( .include.current ) colnames(x.diff) 
                                         else NULL,   colnames(temp4))
          
          if (q.int) {
            
            x.diff.temp <- x.diff
            for (ii in 1:NCOL(x.diff))
              x.diff.temp <- cbind(x.diff.temp, WN.lags(cbind(x.diff[, ii]),
                                                        lags = q.int))
            LagsX <- p.int + q.int # May be an argument later 20180316
            u.resL   <- WN.lags(cbind(u.res), lags = LagsX )
            
            z.Resid  <- residuals(lm(u.res ~  -1 + x.diff +  #x.int +
                              WN.lags(cbind(y[-c(1:d.int)]), lags = 9) +
                         x.diff.temp + u.resL))
              
            z.Resid  <- WN.lags(cbind(z.Resid), lags = q.int)
            rm(u.resL, x.diff.temp)
            
            colnames(z.Resid) <- paste("MAcoeff", 1:q.int, sep = "")
          }
          
          
          x.matrix <- cbind(save.this, x.errU, z.Resid)
          
          list.names <- vector("list", NCOL(x.matrix))
          names(list.names) <- colnames(x.matrix)
          for (ii in 1:NCOL(x.matrix)) 
            list.names[[ii]] <- ii
          attr(x.matrix, "assign") <- list.names
          x <- x.matrix
          
          y <- extra$y <- y.diff 
          w <- w[1:length(y.diff)]
          
        }
        
        
      }), list( .ord = ord , .diffCovs = diffCovs,
                .order.trend = order.trend ,
                .nodrift = nodrift , .xLag = xLag , .toc = toc ,
                .include.current = include.current,
                .include.int = include.int ))),
      
      
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 2,
             Q1 = 1,
             expected = TRUE,
             multipleResponses = TRUE,
             parameters.names = c("mean", if ( .var.arg ) "var" else "sd"),
             var.arg = .var.arg ,
             parallel = .parallel ,
             zero = .zero )
      }, list( .zero = zero ,
               .parallel = parallel ,
               .var.arg = var.arg ))),
      
      
      
      initialize = eval(substitute(expression({
        orig.y <- y
        
        if (length(attr(orig.y, "Prior.Weights"))) {
          if (any(c(w) != 1))
            warning("replacing the 'weights' argument by the 'Prior.Weights'",
                    "attribute of the response (probably due to Qvar()")
          
          
          w <- attr(orig.y, "Prior.Weights")
          
          
          extra$attributes.y <- attributes(orig.y)
          
        } else {
        }
        
        
        temp5 <-
          w.y.check(w = w, y = y,
                    ncol.w.max = Inf,
                    ncol.y.max = Inf,
                    out.wy = TRUE,
                    colsyperw = 1,
                    maximize = TRUE)
        w <- temp5$w
        y <- temp5$y
        
        ncoly <- ncol(y)
        M1 <- 2
        extra$ncoly <- ncoly
        extra$M1 <- M1
        M <- M1 * ncoly
        
        
        
        mynames1 <- param.names("mean", ncoly)
        mynames2 <- param.names(if ( .var.arg ) "var" else "sd", ncoly)
        predictors.names <-
          c(namesof(mynames1, .lmean , earg = .emean , tag = FALSE),
            if ( .var.arg )
              namesof(mynames2, .lvare , earg = .evare , tag = FALSE) else
                namesof(mynames2, .lsdev , earg = .esdev , tag = FALSE))
        predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]
        extra$predictors.names <- predictors.names
        
        
        if (!length(etastart)) {
          sdev.init <- mean.init <- matrix(0, n, ncoly)
          for (jay in 1:ncoly) {
            jfit <- lm.wfit(x = x,  y = y[, jay], w = w[, jay])
            mean.init[, jay] <- if ( .lmean == "loglink")
              pmax(1/1024, y[, jay]) else
          if ( .imethod == 1) median(y[, jay]) else
            if ( .imethod == 2) weighted.mean(y[, jay], w = w[, jay]) else
               if ( .imethod == 3) weighted.mean(y[, jay], w = w[, jay]) *
              0.5 + y[, jay] * 0.5 else
                mean(jfit$fitted)
            
            sdev.init[, jay] <-
              if ( .imethod == 1) {
                sqrt( sum(w[, jay] *
                        (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
              } else if ( .imethod == 2) {
                if (jfit$df.resid > 0)
                  sqrt( sum(w[, jay] * jfit$resid^2) / jfit$df.resid ) else
                    sqrt( sum(w[, jay] * jfit$resid^2) / sum(w[, jay]) )
              } else if ( .imethod == 3) {
                sqrt( sum(w[, jay] *
                      (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
              } else {
                sqrt( sum(w[, jay] * abs(y[, jay] -
                                     mean.init[, jay])) / sum(w[, jay]) )
              }
            
            if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
              sdev.init[, jay] <- 1.01
            
          }
          
          
          if (length( .isdev )) {
            sdev.init <- matrix( .isdev , n, ncoly, byrow = TRUE)
          }
          
          
          etastart <-
            cbind(theta2eta(mean.init,   .lmean , earg = .emean ),
                  if ( .var.arg )
                    theta2eta(sdev.init^2, .lvare , earg = .evare ) else
                      theta2eta(sdev.init  , .lsdev , earg = .esdev ))
          etastart <-
            etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
          
          colnames(etastart) <- predictors.names
        }
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .isdev = isd,
                .var.arg = var.arg, .imethod = imethod ))),
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        
        M1 <- extra$M1
        ncoly <- NOS <- extra$ncoly
        pOrd <- .ord[1]
        dOrd <- .ord[2]
        qOrd <- .ord[3]
        
        if ( .lmean == "explink") {
          if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
            warning("turning some columns of 'eta' positive in @linkinv")
            for (ii in 1:ncoly)
              eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
          }
        }
        
        y.est <- eta2theta(eta[, M1*(1:ncoly) - 1],
                                .lmean , earg = .emean ) #+ extra$u.res
        
        if (dOrd) {
          mat.save <- extra$mat.save
          nn <- nrow(mat.save)
          sum.mat <- cbind(extra$y.r,
                           extra$mat.save[, -c(dOrd), drop = FALSE],
                           c(y.est, rep(0, dOrd )))
          
          t.rev <- rev(1:NCOL(sum.mat))
          
          for (ii in 1:dOrd) {
            yaux1 <- sum.mat[c(1:(nn - dOrd + ii - 1)), t.rev[ii], drop = FALSE]
            yaux2 <- sum.mat[c(1:(nn - dOrd + ii - 1)),
                             t.rev[ii] - 1, drop = FALSE]
            yaux3 <- yaux1 + yaux2
            yaux3 <- rbind(sum.mat[1, t.rev[ii] - 1, drop = FALSE], yaux3)
            sum.mat[ , t.rev[ii] - 1] <- c(yaux3, rep(0, dOrd - ii))
          }
          y.est <- sum.mat[, 1] 
          
        } else {
          
          nn.r <- nrow(extra$y.r); nnOrd <- exp(0)
          my.mx <- max(c(pOrd, qOrd))
          y.est.ret <- matrix(NA_real_, nn.r, 1)
          y.est.ret[c(1:my.mx)] <- (nnOrd + 0.10)*extra$y.r[c(1:my.mx)]
          y.est.ret[-c(1:my.mx)] <- y.est
          y.est <- c(y.est.ret)
        }
        
        y.est
        
      }, list( .lmean = lmean, .ord = ord ,
               .emean = emean, .esdev = esdev , .evare = evare,
               .smallno = smallno ))),
      
      last = eval(substitute(expression({
        M1 <- extra$M1
        nn <- nrow(eta)
        temp.names <- c(mynames1, mynames2)
        temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]
        misc$link <- rep_len( .lmean , M1 * ncoly)
        misc$earg <- vector("list", M1 * ncoly)
        names(misc$link) <- names(misc$earg) <- temp.names
        for (ii in 1:ncoly) {
          misc$link[ M1*ii-1 ] <- .lmean
          misc$link[ M1*ii   ] <- if ( .var.arg ) .lvare else .lsdev
          misc$earg[[M1*ii-1]] <- .emean
          misc$earg[[M1*ii  ]] <- if ( .var.arg ) .evare else .esdev
        }
        
        misc$var.arg <- .var.arg
        misc$M1 <- M1
        misc$expected <- TRUE
        misc$imethod <- .imethod
        misc$multipleResponses <- TRUE
        misc$parallel <- .parallel
        misc$apply.parint <- .apply.parint
        misc$smallno <- .smallno
        y <- extra$y.r
        
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .parallel = parallel, .apply.parint = apply.parint,
                .smallno = smallno, .ord = ord ,
                .var.arg = var.arg, .imethod = imethod ))),
      
      loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL,
                 summation = TRUE) {
          ncoly <- extra$ncoly
          M1 <- extra$M1
          mu <-  eta2theta(eta[, M1*(1:ncoly) - 1],
                           .lmean , earg = .emean ) #+ extra$u.res 
          if ( .lmean == "explink") {
            if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
              warning("turning some columns of 'eta' positive in @loglikelihood")
              for (ii in 1:ncoly)
                eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
            }
          }
          
          if ( .var.arg ) {
            Varm <- eta2theta(eta[, M1*(1:ncoly)], .lvare , earg = .evare )
            sdev <- sqrt(Varm)
          } else {
            sdev <- eta2theta(eta[, M1*(1:ncoly)], .lsdev , earg = .esdev )
          }
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dnorm(y, m = mu, sd = sdev, log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list( .lsdev = lsdev, .lvare = lvare,
                 .esdev = esdev, .evare = evare,
                 .lmean = lmean, .emean = emean ,
                 .smallno = smallno,
                 .var.arg = var.arg ))),
      vfamily = c("ARIMAX.errors.ff", "vgltsmff"),
      validparams = eval(substitute(function(eta, y, extra = NULL) {
        M1 <- 2
        ncoly <- NCOL(y)
        
        mymu <- eta2theta(  eta[, M1*(1:ncoly) - 1], .lmean , earg = .emean )
        if ( .var.arg ) {
          Varm <- eta2theta(eta[, M1*(1:ncoly)    ], .lvare , earg = .evare )
          sdev <- 111
        } else {
          sdev <- eta2theta(eta[, M1*(1:ncoly)    ], .lsdev , earg = .esdev )
          Varm <- 111
        }
        okay1 <- all(is.finite(mymu)) &&
          all(is.finite(sdev)) && all(0 < sdev) &&
          all(is.finite(Varm)) && all(0 < Varm)
        okay2 <- TRUE
        if ( .lmean == "explink") {
          okay2 <- all(0 < eta[, M1*(1:ncoly) - 1])
        }
        okay1 && okay2
      }, list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
               .emean = emean, .esdev = esdev, .evare = evare,
               .smallno = smallno,
               .var.arg = var.arg ))),
      
      
      
      
      simslot = eval(substitute(
        function(object, nsim) {
          
          pwts <- if (length(pwts <- object@prior.weights) > 0)
            pwts else weights(object, type = "prior")
          if (any(pwts != 1))
            warning("ignoring prior weights")
          mymu <- fitted(object)
          eta <- predict(object)
          if ( .var.arg ) {
            Varm <- eta2theta(eta[, c(FALSE, TRUE)], .lvare , earg = .evare )
            sdev <- sqrt(Varm)
          } else {
            sdev <- eta2theta(eta[, c(FALSE, TRUE)], .lsdev , earg = .esdev )
          }
          rnorm(nsim * length(mymu), mean = mymu, sd = sdev)
        }, list( .lsdev = lsdev, .lvare = lvare,
                 .esdev = esdev, .evare = evare,
                 .lmean = lmean,
                 .smallno = smallno,
                 .var.arg = var.arg ))),
      
      
      
      
      deriv = eval(substitute(expression({
        ncoly <- extra$ncoly
        M1 <- extra$M1
        y <- extra$y
        
        if ( .lmean == "explink") {
          if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
            warning("turning some columns of 'eta' positive in @deriv")
            for (ii in 1:ncoly)
              eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
          }
        }
        
        
        mymu <- eta2theta(  eta[, M1*(1:ncoly) - 1], .lmean , earg = .emean )
        if ( .var.arg ) {
          Varm <- eta2theta(eta[, M1*(1:ncoly)    ], .lvare , earg = .evare )
          sdev <- sqrt(Varm)
        } else {
          sdev <- eta2theta(eta[, M1*(1:ncoly)    ], .lsdev , earg = .esdev )
        }
        
        dl.dmu <- (y - mymu) / sdev^2
        if ( .var.arg ) {
          dl.dva <- -0.5 / Varm + 0.5 * (y - mymu)^2 / sdev^4
        } else {
          dl.dsd <- -1.0 / sdev +       (y - mymu)^2 / sdev^3
        }
        
        dmu.deta <- dtheta.deta(mymu,   .lmean , earg = .emean )
        if ( .var.arg ) {
          dva.deta <- dtheta.deta(Varm, .lvare , earg = .evare )
        } else {
          dsd.deta <- dtheta.deta(sdev, .lsdev , earg = .esdev )
        }
        
        ans <- c(w) *
          cbind(dl.dmu * dmu.deta,
                if ( .var.arg ) dl.dva * dva.deta else
                  dl.dsd * dsd.deta)
        ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
        
        
        
        
        
        
        ans
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .smallno = smallno,
                .var.arg = var.arg ))),
      weight = eval(substitute(expression({
        wz <- matrix(NA_real_, n, M)  # Diagonal matrix
        
        
        ned2l.dmu2 <- 1 / sdev^2
        
        if ( .var.arg ) {
          ned2l.dva2 <- 0.5 / Varm^2
        } else {
          ned2l.dsd2 <- 2 / sdev^2
        }
        
        wz[, M1*(1:ncoly) - 1] <- ned2l.dmu2 * dmu.deta^2
        wz[, M1*(1:ncoly)    ] <- if ( .var.arg ) {
          ned2l.dva2 * dva.deta^2
        } else {
          ned2l.dsd2 * dsd.deta^2
        }
        
        w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
      }), list( .var.arg = var.arg ))))
}



ARIMAX.errors.ff.control <- function(save.weights = TRUE, 
                                     summary.HDEtest = FALSE,...) {
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}