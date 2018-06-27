# These functions are
# Copyright (C) 2014-2018 V. Miranda and T. W. Yee, University of Auckland
# 20170221

gammaRMean <- function(zero = "rate", lmu = "gammaRMeanlink",
                       lrate = "loge", lshape = NULL,
                       irate = NULL,   ishape = NULL,
                       lss = TRUE) {
  
  
  expected <- TRUE  # FALSE does not work well
  
  iratee <- irate
  
  if (length(lmu)) {
    lshape <- lmu
  }
  
  if (length(lmu) && !identical(lmu, "gammaRMeanlink"))
    stop("Wrong input for 'lmu'.")
  
  lratee <- as.list(substitute(lrate))
  eratee <- link2list(lratee)
  lratee <- attr(eratee, "function.name")
  
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  
  
  if (length( iratee) && !is.Numeric(iratee, positive = TRUE))
    stop("bad input for argument 'irate'")
  if (length( ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  
  
  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")
  
  
  ratee.TF <- if (lss) c(TRUE, FALSE) else c(FALSE, TRUE)
  scale.12 <- if (lss) 1:2 else 2:1
  blurb.vec <- c(namesof("rate",  lratee, earg = eratee),
                 namesof("shape", lshape, earg = eshape))
  blurb.vec <- blurb.vec[scale.12]
  
  
  
  new("vglmff",
      blurb = c("2-parameter Gamma distribution\n",
                "Links:    ",
                blurb.vec[1], ", ",
                blurb.vec[2], "\n",
                "Mean:     mu = shape/rate\n",
                "Variance: (mu^2)/shape = shape/rate^2"),
      
      
      
      
      constraints = eval(substitute(expression({
        
        # ADVISE THOMAS about the next change.
        # Otherwise, grep() will detect 'rate' twice therefore no cols
        # predictors.names = parameters.names instead of parameters.names
        constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                    predictors.names = parameters.names,
                                    M1 = 2)
      }), list( .zero = zero ))),
      
      
      
      
      infos = eval(substitute(function(...) {
        list(M1 = 2,
             Q1 = 1,
             expected = .expected ,
             multipleResponses = TRUE,
             zero = .zero )
      }, list( .zero = zero, .scale.12 = scale.12, .ratee.TF = ratee.TF,
               .expected = expected
      ))),
      
      
      
      
      initialize = eval(substitute(expression({
        
        temp5 <-
          w.y.check(w = w, y = y,
                    Is.positive.y = TRUE,
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
        
        
        if ( .lss ) {
          mynames1 <- param.names("rate",  ncoly)
          mynames2 <- param.names("shape", ncoly)
          predictors.names <-
            c(namesof(mynames1, .lratee , earg = .eratee , tag = FALSE),
              namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))
          
        } else {
          mynames1 <- param.names("shape", ncoly)
          mynames2 <- param.names("rate",  ncoly)
          predictors.names <-
            c(namesof(mynames1, .lshape , earg = .eshape , tag = FALSE),
              namesof(mynames2, .lratee , earg = .eratee , tag = FALSE))
        }
        
        parameters.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
        predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]
        
        
        Ratee.init <- matrix(if (length( .iratee )) .iratee else 0 + NA,
                             n, ncoly, byrow = TRUE)
        Shape.init <- matrix(if (length( .ishape )) .iscale else 0 + NA,
                             n, ncoly, byrow = TRUE)
        
        if (!length(etastart)) {
          mymu <- y + 0.167 * (y == 0)
          
          
          for (ilocal in 1:ncoly) {
            junk <- lsfit(x, y[, ilocal], wt = w[, ilocal],
                          intercept = FALSE)
            var.y.est <- sum(c(w[, ilocal]) * junk$resid^2) / (nrow(x) -
                                                        length(junk$coef))
            
            if (!is.Numeric(Shape.init[, ilocal]))
              Shape.init[, ilocal] <- (mymu[, ilocal])^2 / var.y.est
            
            if (!is.Numeric(Ratee.init[, ilocal]))
              Ratee.init[, ilocal] <- Shape.init[, ilocal] / mymu[, ilocal]
          }
          
          if ( .lshape == "loglog")
            Shape.init[Shape.init <= 1] <- 3.1  # Hopefully value's big enough
          
          neweshape <- .eshape
          if ( .lshape == "gammaRMeanlink") {
            neweshape$rate <- Ratee.init
          }
          
          # Allow proper initial values
          mycheck <- which((Shape.init / Ratee.init) < 1)
          if (length(mycheck))
            Shape.init[mycheck] <- Ratee.init[mycheck] + 1e-1
          
          etastart <- if ( .lss )
            cbind(theta2eta(Ratee.init, .lratee , earg = .eratee ),
                  theta2eta(Shape.init, .lshape , earg = neweshape))[,
                        interleave.VGAM(M, M1 = M1)] else
              cbind(theta2eta(Shape.init, .lshape , earg = neweshape),
                    theta2eta(Ratee.init, .lratee , earg = .eratee ))[,
                        interleave.VGAM(M, M1 = M1)]
        }
        
      }), list( .lratee = lratee, .lshape = lshape,
                .iratee = iratee, .ishape = ishape,
                .eratee = eratee, .eshape = eshape,
                .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        
       Ratee <- eta2theta(eta[,    .ratee.TF  ], .lratee , earg = .eratee )
        
        if ( .lshape == "gammaRMeanlink") {
          neweshape <- .eshape
          neweshape$rate <- Ratee
          Shape <- eta2theta(eta[, !( .ratee.TF )], .lshape ,
                             earg = neweshape)
        } else {
          Shape <- eta2theta(eta[, !( .ratee.TF )], .lshape ,
                             earg = .eshape )
        }
       
       Shape / Ratee
        
      }, list( .lratee = lratee, .lshape = lshape,
               .eratee = eratee, .eshape = eshape,
               .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      
      last = eval(substitute(expression({
        
        misc$multipleResponses <- TRUE
        
        M1 <- extra$M1
        avector <- if ( .lss ) c(rep_len( .lratee , ncoly),
                                 rep_len( .lshape , ncoly)) else
                                   c(rep_len( .lshape , ncoly),
                                     rep_len( .lratee , ncoly))
        
        misc$link <- avector[interleave.VGAM(M, M1 = M1)]
        temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
        names(misc$link) <- temp.names
        
        misc$earg <- vector("list", M)
        names(misc$earg) <- temp.names
        
        neweshape <- .eshape
        if ( .lshape == "gammaRMeanlink") {
          neweshape$rate <- if ( .lss )
             eta2theta(coef(fit)[1], .lratee , earg = .eratee) else 
               eta2theta(coef(fit)[2], .lratee , earg = .eratee)
        }
        
        for (ii in 1:ncoly) {
          misc$earg[[M1*ii-1]] <- if ( .lss ) .eratee else neweshape
          misc$earg[[M1*ii  ]] <- if ( .lss ) neweshape else .eratee
        }
        
        misc$M1 <- M1
        
        # Advise Thomas.
        if ( .lshape == "gammaRMeanlink") {
          fit$fitted.values <- 
            gammaRMeanlink(theta = Shape, rate = Ratee, inverse = FALSE) 
        }
        
      }), list( .lratee = lratee, .lshape = lshape,
                .eratee = eratee, .eshape = eshape,
                .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      
      loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL, summation = TRUE) {
         
          Ratee <- eta2theta(eta[,    .ratee.TF  ], .lratee ,
                             earg = .eratee )
          
          if ( .lshape == "gammaRMeanlink") {
            neweshape <- .eshape
            neweshape$rate <- Ratee
            Shape <- eta2theta(eta[, !( .ratee.TF )], .lshape ,
                               earg = neweshape)
          } else {
            Shape <- eta2theta(eta[, !( .ratee.TF )], .lshape ,
                               earg = .eshape)
          }
          
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * dgamma(x=y, shape = Shape,
                                     rate = Ratee, log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list( .lratee = lratee, .lshape = lshape,
                 .eratee = eratee, .eshape = eshape,
                 .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      vfamily = c("gammaR"),
      
      
      
      simslot = eval(substitute(
        function(object, nsim) {
          
          pwts <- if (length(pwts <- object@prior.weights) > 0)
            pwts else weights(object, type = "prior")
          if (any(pwts != 1))
            warning("ignoring prior weights")
          eta <- predict(object)
          
          # All works with drop = TRUE as well.
          Ratee <- eta2theta(eta[,    .ratee.TF , drop = FALSE ],
                             .lratee , earg = .eratee )
          
          if ( .lshape == "gammaRMeanlink") {
            neweshape <- .eshape
            neweshape$rate <- Ratee
            Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE],
                               .lshape , earg = neweshape)
          } else {
            Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE],
                               .lshape , earg = .eshape )
          }
          
          rgamma(nsim * length(Shape), shape = Shape, rate = Ratee)
          
        }, list( .lratee = lratee, .lshape = lshape,
                 .eratee = eratee, .eshape = eshape,
                 .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      
      validparams = eval(substitute(function(eta, y, extra = NULL) {
        
        Ratee <- eta2theta(eta[,    .ratee.TF , drop = FALSE ],
                           .lratee , earg = .eratee )
        
        if ( .lshape == "gammaRMeanlink") {
          neweshape <- .eshape
          neweshape$rate <- Ratee
          Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE], .lshape ,
                             earg = neweshape)
          #ratio.SR <- Shape / Ratee 
        } else {
          Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE],
                             .lshape , earg = .eshape )
        }
        
        okay1 <- all(is.finite(Shape)) && all(0 < Shape) &&
          all(is.finite(Ratee)) && all(0 < Ratee) 
        
        okay1
        
      }, list( .lshape = lshape, .lratee = lratee,
               .eshape = eshape, .eratee = eratee,
               .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss))),
      
      
      
      
      deriv = eval(substitute(expression({
        
        M1 <- 2
        NOS <- extra$ncoly
        Ratee <- eta2theta(eta[,    .ratee.TF , drop = FALSE ],
                           .lratee , earg = .eratee )
        
        if ( .lshape == "gammaRMeanlink") {
          neweshape <- .eshape
          neweshape$rate <- Ratee
          Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE], .lshape ,
                                 earg = neweshape)
        } else {
          Shape <- eta2theta(eta[, !( .ratee.TF ), drop = FALSE],
                             .lshape , earg = .eshape )
        }
        
        dl.dratee <- mu - y
        dl.dshape <- log(y * Ratee) - digamma(Shape)
        
        
        if ( .lshape == "gammaRMeanlink") {
          neweshape$wrt.param <- 1
          neweshape$inverse   <- TRUE
          neweshape$deriv <- 1
          
          dshape.deta1 <- dtheta.deta(Shape, .lshape ,
                                   earg = neweshape)[, 1:NOS, drop = FALSE]
          dratee.deta1 <- dtheta.deta(Shape, .lshape ,
                                earg = neweshape)[, -(1:NOS), drop = FALSE]

          neweshape$wrt.param <- 2
          dshape.deta2 <-  dtheta.deta(Shape, .lshape ,
                                  earg = neweshape)[, 1:NOS, drop = FALSE]
          dratee.deta2 <- dtheta.deta(Ratee, .lratee , earg = .eratee )
          
          # Remove later.
          #dratee.deta <- dtheta.deta(Ratee, .lratee , earg = .eratee )
          #dshape.deta <- dtheta.deta(Shape, .lshape , earg = neweshape )
          
          myderiv <- if ( .lss )
            c(w) * cbind(dl.dshape * dshape.deta2 +
                                 dl.dratee * dratee.deta2,
                         dl.dshape * dshape.deta1 +
                                    dl.dratee * dratee.deta1) else
              c(w) * cbind(dl.dshape * dshape.deta1 +
                                 dl.dratee * dratee.deta1,
                           dl.dshape * dshape.deta2 +
                                 dl.dratee * dratee.deta2)
        } else {
          dratee.deta <- dtheta.deta(Ratee, .lratee , earg = .eratee )
          dshape.deta <- dtheta.deta(Shape, .lshape , earg = .eshape )
          
          myderiv <- if ( .lss )
            c(w) * cbind(dl.dratee * dratee.deta,
                         dl.dshape * dshape.deta) else
                           c(w) * cbind(dl.dshape * dshape.deta,
                                        dl.dratee * dratee.deta)
        }
        
        myderiv[, interleave.VGAM(M, M1 = M1)]
      }), list( .lratee = lratee, .lshape = lshape,
                .eratee = eratee, .eshape = eshape,
                .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss ))),
      
      
      
      
      weight = eval(substitute(expression({
        
        ned2l.dratee2 <- Shape / (Ratee^2)
        ned2l.drateeshape <- -1/Ratee
        ned2l.dshape2 <- trigamma(Shape)
        
        if ( .expected) {
          ratee.adjustment <-  0
          shape.adjustment <-  0
        } else {
          d2ratee.deta2 <- d2theta.deta2(Ratee, .lratee , earg = .eratee )
          d2shape.deta2 <- d2theta.deta2(Shape, .lshape , earg = .eshape )
          ratee.adjustment <- dl.dratee * d2ratee.deta2
          shape.adjustment <- dl.dshape * d2shape.deta2
        }
        
        
        if ( .lshape == "gammaRMeanlink") {
          
          term1 <- ned2l.dshape2 * dshape.deta2^2 + 
                     2 * ned2l.drateeshape * dshape.deta2 * dratee.deta2 +
                        ned2l.dratee2 * dratee.deta2^2
          
          term2 <- ned2l.dshape2 * dshape.deta1^2 + 
                      2 * ned2l.drateeshape * dshape.deta1 * dratee.deta1 +
                         ned2l.dratee2 * dratee.deta1^2
          
          term3 <- ned2l.dshape2 * dshape.deta1 * dshape.deta2 +
                      ned2l.drateeshape * ( dshape.deta1 * dratee.deta2 +
                                  dshape.deta2 * dratee.deta1) +
                            ned2l.dratee2 * dratee.deta1 * dratee.deta2
          
          if ( .lss ) {
            wz <- array(c(c(w) * term1, c(w) * term2, c(w) * term3),
                        dim = c(n, M / M1, 3))
          } else {
            wz <- array(c(c(w) * term2, c(w) * term1, c(w) * term3),
                        dim = c(n, M / M1, 3))
          }
          
        } else {
          wz <- if ( .lss )
      array(c(c(w) * (ned2l.dratee2 * dratee.deta^2 - ratee.adjustment),
            c(w) * (ned2l.dshape2 * dshape.deta^2 - shape.adjustment),
            c(w) * (ned2l.drateeshape * dratee.deta * dshape.deta)),
            dim = c(n, M / M1, 3)) else
         array(c(c(w) * (ned2l.dshape2 * dshape.deta^2 - shape.adjustment),
               c(w) * (ned2l.dratee2 * dratee.deta^2 - ratee.adjustment),
               c(w) * (ned2l.drateeshape * dratee.deta * dshape.deta)),
               dim = c(n, M / M1, 3))
        }
        
        wz <- arwz2wz(wz, M = M, M1 = M1)
        wz
        
      }), list( .lratee = lratee, .lshape = lshape,
                .eratee = eratee, .eshape = eshape, .expected = expected,
                .scale.12 = scale.12, .ratee.TF = ratee.TF, .lss = lss  ))))
}