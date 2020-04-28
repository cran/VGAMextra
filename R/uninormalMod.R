#############################################################################
# These functions are Coyright (C) 2014 - 2019
# V. Miranda-Soberanis, Auckland University of Technology
# T. Yee, University of Auckland


uninormalMod <- function(zero = "sd", tau.arg = NULL,
                         lmean = "identitylink", lsd = "loglink",
                         parallel.qr = TRUE,
                         imethod = 1, isd = NULL, parallel = FALSE,
                         smallno = 1.0e-5) {
  
  apply.parint <- FALSE
  var.arg <- FALSE
  
  p.vector <- c(tau.arg); rm(tau.arg)
  if (length(p.vector) & (!is.Numeric(p.vector, positive = TRUE) ||
                          any(p.vector >= 1) ))
    stop("Invalid input for argument 'p.vector'.")
  
  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")
  
  myflag <- (lmean == "uninormalQlink")
  
  
  if (myflag  && !length(p.vector))
    stop("For vector quantile regression 'p.vector' must be entered.")
  
  lvar <- lsd <- "loglink"
  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")
  
  lvare <- as.list(substitute(lvar))
  evare <- link2list(lvare)
  lvare <- attr(evare, "function.name")
  
  if (!is.Numeric(smallno, length.arg = 1, positive = TRUE))
    stop("argument 'smallno' must be positive and close to 0")
  
  if (smallno > 0.1) {
    warning("replacing argument 'smallno' with 0.1")
    smallno <- 0.1
  }
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) || imethod > 4)
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
      blurb = c("Univariate normal distribution (uninormalQlink)\n\n",
                "Links:    ",
                namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
                if (var.arg)
                  namesof("var",  lvare, earg = evare, tag = TRUE) else
                    namesof("sd" ,  lsdev, earg = esdev, tag = TRUE),
                "\n",
                if (var.arg) "Variance: var" else "Variance: sd^2"),
      
      
  constraints = eval(substitute(expression({
    
    NOS <- M/M1
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel ,
                          constraints = constraints)

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                               M1 = 2)
    
    }), list( .zero = zero, .parallel.qr = parallel.qr ,
              .parallel = parallel, .apply.parint = apply.parint ))),
  
  
      
      first = eval(substitute(expression({
        if (FALSE) {
          # Later, consider turning Q.reg() as deprecated. - Discussion
          NOS <- NCOL(y)
          fooNOS <- Q.reg(y, length(p.vector))
          
        }
        
      }),list( .zero = zero, .parallel = parallel,
               .apply.parint = apply.parint ))),
  
  
  
      infos = eval(substitute(function(...) {
        
        list(M1 = 2,
             Q1 = 1,
             charfun = TRUE,
             expected = TRUE,
             hadof = FALSE,
             multipleResponses = TRUE,
             #parameters.names = c("mean", if ( .var.arg ) "var" else "sd"),
             var.arg = .var.arg ,
             parallel = .parallel ,
             zero = .zero )
      }, list( .zero = zero , .parallel = parallel , .var.arg = var.arg , 
               .lmean = lmean , .emean = emean ))),
  
      
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
        
        
        
        mynames1 <- param.names("mean", ncoly, skip1 = TRUE)
        if ( .myflag ) {
          mynames1 <- rep("mu", ncoly)
        }
        mynames2 <- param.names(if ( .var.arg ) "var" else "sd",
                                ncoly, skip1 = TRUE)
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
          
          neweQlink <- .emean
          if ( .lmean == "uninormalQlink") {
            neweQlink$sd <- sdev.init
          }
          
          extra$mean.init <- mean.init
          extra$sdev.init <- sdev.init
          
          etastart <-
            cbind(theta2eta(mean.init , .lmean , earg = neweQlink),
                  if ( .var.arg )
                    theta2eta(sdev.init^2, .lvare , earg = .evare ) else
                      theta2eta(sdev.init  , .lsdev , earg = .esdev ))
          etastart <-
            etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
          
          colnames(etastart) <- predictors.names
        }
        
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .isdev = isd, .myflag = myflag ,
                .var.arg = var.arg, .imethod = imethod ))),
  
      
      linkinv = eval(substitute(function(eta, extra = NULL) { 
        M1 <- extra$M1
        ncoly <- extra$ncoly
        
        neweQlink <- .emean
        if ( .myflag ) {
          
          int.sd <- eta2theta(eta[, M1 *(1:ncoly)], .lsdev , earg = .esdev )
          neweQlink$sd <- int.sd
          
        } else {
          
          if ( .lmean == "explink") {
            if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
              warning("turning some columns of 'eta' positive in @linkinv")
              for (ii in 1:ncoly)
                eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
            }
          }
       }
        
        eta2theta(eta[, M1*(1:ncoly) - 1, drop = FALSE], .lmean ,
                  earg = neweQlink )
        
      }, list( .lmean = lmean, .lsdev = lsdev ,
               .emean = emean, .esdev = esdev ,
               .evare = evare, .esdev = esdev,
               .smallno = smallno , .myflag = myflag ))),
  
  
      
      last = eval(substitute(expression({
        
        M1 <- extra$M1
        int.sd <- eta2theta(eta[, M1 *(1:ncoly)], .lsdev , earg = .esdev )
        
        neweQlink <- .emean
        if ( .myflag ) {
          # Check this. zzz
          #int.sd <- coef(fit)[interleave.VGAM(M, M1 = 2, inverse = TRUE)]
          #int.sd <- int.sd[-(1:NCOL(y))]
          #matrix(int.sd, NROW(y), NCOL(y), byrow = TRUE)
          neweQlink$sd <- int.sd
        }
        ###  10/08/2019 - Let Thomas know...
        #if ( .myflag ) {
        #  fit$fitted.values <-  uninormalQlink(theta = mymu,
        #                         sd = sdev, p = .p.vector , inverse = FALSE) 
        #}
        
        mymu <- eta2theta(eta[, M1*(1:ncoly) - 1, drop = FALSE], .lmean ,
                          earg = neweQlink )
        
        temp.names <- c(mynames1, mynames2)
        temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]
        misc$link <- rep_len( .lmean , M1 * ncoly)
        misc$earg <- vector("list", M1 * ncoly)
        names(misc$link) <- names(misc$earg) <- temp.names
        for (ii in 1:ncoly) {
          misc$link[ M1*ii-1 ] <- .lmean
          misc$link[ M1*ii   ] <- if ( .var.arg ) .lvare else .lsdev
          misc$earg[[M1*ii-1]] <- neweQlink #.emean
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
        
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .parallel = parallel, .apply.parint = apply.parint,
                .smallno = smallno, .myflag = myflag ,
                .p.vector = p.vector ,
                .var.arg = var.arg, .imethod = imethod ))),
  
      
      loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL,
                 summation = TRUE) {
          
          ncoly <- extra$ncoly
          M1 <- extra$M1
          
          if ( .var.arg ) {
            Varm <- eta2theta(eta[, M1*(1:ncoly)], .lvare , earg = .evare )
            sdev <- sqrt(Varm)
          } else {
            sdev <- eta2theta(eta[, M1*(1:ncoly)], .lsdev , earg = .esdev )
          }
          
          neweQlink <- .emean
          if ( .myflag ) {
            
            neweQlink$sd <- sdev
            
          } else {
            
            if ( .lmean == "explink") {
              if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
                warning("turning some columns of 'eta' positive in @loglikelihood")
                for (ii in 1:ncoly)
                  eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
              }
            }
            
          }
          
          mymean <- eta2theta(eta[, M1 * (1:ncoly) - 1, drop = FALSE], 
                              .lmean , earg = neweQlink )
          
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dnorm(y, mean = mymean, sd = sdev, log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list( .lsdev = lsdev, .lvare = lvare, .lmean = lmean ,
                 .esdev = esdev, .evare = evare, .emean = emean ,
                 .myflag = myflag ,
                 .smallno = smallno, .var.arg = var.arg ))),
  
      vfamily = c("uninormalMod"),
      
################################################################
################################################################

      deriv = eval(substitute(expression({
        NOS <- ncoly <- extra$ncoly
        M1 <- extra$M1
        
        if ( .var.arg ) {
          Varm <- eta2theta(eta[, M1*(1:ncoly) ], .lvare , earg = .evare )
          sdev <- sqrt(Varm)
        } else {
          sdev <- eta2theta(eta[, M1*(1:ncoly) ], .lsdev , earg = .esdev )
        }
        
        neweQlink <- .emean
        if ( .myflag ) {
          
          neweQlink$sd  <- sdev 
          mymu <- eta2theta(eta[, M1*(1:ncoly) - 1, drop = FALSE],
                            .lmean , earg = neweQlink )
          neweQlink$inverse <- TRUE
          neweQlink$deriv   <- 1
        
          neweQlink$wrt.param <- 1
          # Returns a two-column matrix -> dmu/deta1 , dsd/deta1
          dmusd.deta1 <-  dtheta.deta(mymu, .lmean , earg = neweQlink)
          
          dmu.deta1 <- dmusd.deta1[, 1:NOS, drop = FALSE]
          dsd.deta1 <- dmusd.deta1[, -(1:NOS), drop = FALSE]
          
          neweQlink$wrt.param <- 2
          dmu.deta2 <-  dtheta.deta(mymu, .lmean ,  
                                 earg = neweQlink)[, 1:NOS, drop = FALSE]
          dsd.deta2 <-  dtheta.deta(sdev, .lsdev , earg = .esdev )
          
          ## dl.dmu and dl.dsd... 
          myv <- matrix( .p.vector , NROW(eta), NCOL(y), byrow = TRUE)
          dmu.dsd <- (-1) * sqrt(2) * erf(2 * myv - 1, inverse = TRUE) 
          
          dl.dmu <- (y - mymu) / sdev^2
          dl.dsd <- -1/sdev + dmu.dsd * (y - mymu)/sdev^2  + 
                                              (y - mymu)^2 / sdev^3
          # dl.deta
          dl.deta1 <- dl.dmu * dmu.deta1 + dl.dsd * dsd.deta1
          dl.deta2 <- dl.dmu * dmu.deta2 + dl.dsd * dsd.deta2
          
          ans <- c(w) * cbind(dl.deta1, dl.deta2)
          ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
          
        } else {
          
          if ( .lmean == "explink") {
            if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
              warning("turning some columns of 'eta' positive in @deriv")
              for (ii in 1:ncoly)
                eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
            }
          }
          
          mymu <- eta2theta(eta[, M1*(1:ncoly) - 1, drop = FALSE],
                            .lmean , earg = neweQlink )
          
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
          
          ans <- c(w) * cbind(dl.dmu * dmu.deta,
                  if ( .var.arg ) dl.dva * dva.deta else
                                    dl.dsd * dsd.deta)
          ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
          
        }
        ans
      }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
                .emean = emean, .esdev = esdev, .evare = evare,
                .smallno = smallno, .myflag = myflag ,
                .p.vector = p.vector, .var.arg = var.arg ))),
  
  
      weight = eval(substitute(expression({
        
        wz <- matrix(NA_real_, n, M)  # Diagonal matrix
        
        if ( .myflag ) {
          wz <- matrix(0, n, M + (M - 1))  # NOT diagonal matrix
          
          ned2l.dmu2 <- 1 / sdev^2
          ned2l.dsd2 <- 2 / sdev^2 + (1 / sdev^2) * ( dmu.dsd^2)
          ned2l.dmudsd2 <- dmu.dsd * (1 / sdev^2)
          
          
          wz[, M1*(1:NOS) - 1] <- ned2l.dmu2 * dmu.deta1^2 + 
            2  *  ned2l.dmudsd2 * dmu.deta1 * dsd.deta1 +
            ned2l.dsd2 * dsd.deta1^2
          
          wz[, M1*(1:NOS) ] <- ned2l.dmu2 * dmu.deta2^2 + 
            2  *  ned2l.dmudsd2 * dmu.deta2 * dsd.deta2 +
            ned2l.dsd2 * dsd.deta2^2
          
  
          wz[, M1*(1:NOS) + (M - 1)] <-
                     ned2l.dmu2 * dmu.deta1 * dmu.deta2  +
                       ned2l.dmudsd2 *(dmu.deta1 * dsd.deta2 +
                                         dmu.deta2 * dsd.deta1) +
                       ned2l.dsd2 * dsd.deta1 * dsd.deta2

        } else {
          
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
          
       
        }
        w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
   
      }), list( .var.arg = var.arg  , .myflag = myflag ))))
}
