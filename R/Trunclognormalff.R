####################################################################
# These functions are
# Copyright (C) 1998-2024
# T.W. Yee, University of Auckland.
#  V. Miranda, Siqi Liu, Auckland University of Technology
# All rights reserved.


trunclognormal <-
  function(lmeanlog = "identitylink", lsdlog = "loglink",
           min.support = 1e-6, max.support = Inf, 
           zero = "sdlog") {
    
    
    lmulog <- as.list(substitute(lmeanlog))
    emulog <- link2list(lmulog)
    lmulog <- attr(emulog, "function.name")
    
    lsdlog <- as.list(substitute(lsdlog))
    esdlog <- link2list(lsdlog)
    lsdlog <- attr(esdlog, "function.name")
    
    # Both must be positive.
    if (min.support < 1e-6)
      stop("Bad input for 'min.support'.")
    if (max.support < 1e-6)
      stop("Bad input for 'max.support'.")
    
    
    new("vglmff",
        blurb = c("Doubly Truncated lognormal ",
                  "distribution\n\n",
                  "Links:    ",
                  namesof("meanlog", lmulog, earg = emulog, tag = TRUE),
                  ", ",
                  namesof("sdlog",   lsdlog, earg = esdlog, tag = TRUE),
                  "\n",
                  paste("Truncation point(s): ", 
                        min.support, ", ", max.support, ".")),
        constraints = eval(substitute(expression({
          constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                      predictors.names = predictors.names,
                                      M1 = 2)
        }), list( .zero = zero ))),
        
        
        infos = eval(substitute(function(...) {
          list(M1 = 2,
               Q1 = 1,
               dpqrfun = "lnorm",
               lmeanlog = .lmeanlog ,
               lsdlog   = .lsdlog ,
               expected = TRUE,
               min.support = .min.support,
               max.support = .max.support,
               multipleResponses = FALSE,
               parameters.names = c("meanlog", "sdlog"),
               zero = .zero )
        }, list( .zero = zero,
                 .lmeanlog = lmeanlog,
                 .lsdlog   = lsdlog,
                 .min.support = min.support,
                 .max.support = max.support
        ))),
        
        
        first = eval(substitute(expression({
          
          
          int.min.support <- .min.support
          int.max.support <- .max.support 
          
          extra$truncation <- FALSE # Ordinary lognormal
          if ( (int.min.support > 1e-6) || !is.infinite(int.max.support))
            #if (!is.infinite(int.max.support) || ()) 
            extra$truncation <- TRUE
          
          
        }), list( .min.support = min.support ,
                  .max.support = max.support ))),
        
        
        
        initialize = eval(substitute(expression({
          
          w.y.check(w = w, y = y,
                    Is.positive.y = TRUE)
          
          
          ncoly <- 1 # ncol(yy) so far, one response.
          M1 <- 2
          
          extra$min.support <- matrix( .min.support , n, ncoly, 
                                       byrow =  TRUE)
          if (any(y < extra$min.support))
            stop("\n Some response values smaller ",
                 "than argument 'min.support'.")
          

          extra$max.support <- matrix( .max.support , n, ncoly, 
                                       byrow =  TRUE)
          if (any(y > extra$max.support))
            stop("\n Some response values greater ",
                 "than argument 'max.support'.")
          
          predictors.names <-
            c(namesof("meanlog", .lmulog , .emulog , tag = FALSE),
              namesof("sdlog",   .lsdlog , .esdlog , tag = FALSE))
          
          if (!length(etastart)) {
            mylm <- lm.wfit(x = x, y = c(log(y)), w = c(w))
            sdlog.y.est <- sqrt( sum(c(w) * mylm$resid^2)
                                 / mylm$df.residual )

            etastart <- cbind(
              meanlog = rep_len(theta2eta(fitted(mylm) , .lmulog ,
                                          earg = .emulog ), n),
              sdlog   = rep_len(theta2eta(sdlog.y.est , .lsdlog ,
                                          earg = .esdlog ), n))
          }

          
        }), list( .lmulog = lmulog, .lsdlog = lsdlog,
                  .emulog = emulog, .esdlog = esdlog,
                  .min.support = min.support,
                  .max.support = max.support ))),
        
        
        
        
        linkinv = eval(substitute(function(eta, extra = NULL) {
          
          mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
          sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
          
          qdiffmax <- (log(extra$max.support) - mulog) / sdlog - sdlog
          qdiffmin <- (log(extra$min.support) - mulog) / sdlog - sdlog
          
          
          numdiff <- pnorm(q = qdiffmax) - pnorm(q = qdiffmin)
          dendiff <- pnorm(q = qdiffmax + sdlog) -
            pnorm(q = qdiffmin + sdlog)
          
          ## Mean of a truncated lognormal in [min.support, max.support]
          exp(mulog + 0.5 * sdlog^2) * ( numdiff / dendiff )
          
        }, list( .lmulog = lmulog, .lsdlog = lsdlog,
                 .emulog = emulog, .esdlog = esdlog ))),
        
        
        
        
        last = eval(substitute(expression({
          
          
          misc$link <-    c("meanlog" = .lmulog , "sdlog" = .lsdlog )
          
          misc$earg <- list("meanlog" = .emulog , "sdlog" = .esdlog )
          
          misc$expected <- TRUE
          
          if (!extra$truncation)
            cat("\n",
                "No truncation limits entered - Fitting an ordinary",
                "lognormal distribution.")
          
        }), list( .lmulog = lmulog, .lsdlog = lsdlog,
                  .emulog = emulog, .esdlog = esdlog,
                  .min.support = min.support,
                  .max.support = max.support ))),
        
        
        
        
        
        loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta,
                   extra = NULL,
                   summation = TRUE) {
            
            mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
            sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
            
            
            if (is.infinite( .max.support )) {
              Bstar <- (log( .Machine$double.xmax ) - mulog) / sdlog
            } else {
              Bstar <- (log(extra$max.support) - mulog) / sdlog
            }
            
            Astar <- (log(extra$min.support) - mulog) / sdlog
            
            if (residuals) {
              stop("loglikelihood residuals not implemented yet")
            } else {
        #      ll.elts <- c(w) * ( dlnorm(y, mulog, sdlog = sdlog, 
        #                                 log = TRUE) -
        #                    log( pnorm(q = Bstar, lower.tail = TRUE) -
        #                         pnorm(q = Astar , lower.tail = TRUE) ))  
              
            ll.elts <- c(w) * ( dtrunclnorm(x = y, 
                          meanlog = mulog, sdlog = sdlog,
                          min.support = extra$min.support,
                          max.support = extra$max.support, 
                          log = TRUE))
              
              
              if (summation) {
                sum(ll.elts)
              } else {
                ll.elts
              }
            }
          }, list( .lmulog = lmulog, .lsdlog = lsdlog,
                   .emulog = emulog, .esdlog = esdlog,
                   .min.support = min.support,
                   .max.support = max.support ))),
        
        
        vfamily = c("trunclognormal"),
        
        
        
        validparams = eval(substitute(function(eta, y, extra = NULL) {
          
          mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
          sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
          okay1 <- all(is.finite(mulog)) &&
            all(is.finite(sdlog)) && all(0 < sdlog)
          okay1
        }, list( .lmulog = lmulog, .lsdlog = lsdlog,
                 .emulog = emulog, .esdlog = esdlog,
                 .min.support = min.support,
                 .max.support = max.support ))),
        
        
        
        
        #  simslot = eval(substitute(
        #  function(object, nsim) {#
        
        #pwts <- if (length(pwts <- object@prior.weights) > 0)
        #          pwts else weights(object, type = "prior")
        #if (any(pwts != 1))
        #  warning("ignoring prior weights")
        #eta <- predict(object)
        #mulog <- eta2theta(eta[, c(TRUE, FALSE)], .lmulog , .emulog )
        #sdlog <- eta2theta(eta[, c(FALSE, TRUE)], .lsdlog , .esdlog )
        #rlnorm(nsim * length(mulog),
        #       meanlog = mulog, sdlog = sdlog)
        #}, list( .lmulog = lmulog, .lsdlog = lsdlog,
        #         .emulog = emulog, .esdlog = esdlog ))),
        
        
        
        deriv = eval(substitute(expression({
          
          mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
          sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
          
          dmulog.deta <- dtheta.deta(mulog, .lmulog , earg = .emulog )
          dsdlog.deta <- dtheta.deta(sdlog, .lsdlog ,   earg = .esdlog )
          
          ##################################################################
          
          if (extra$truncation) {
            
            if (is.infinite( .max.support )) {
              Bstar <- (log( .Machine$double.xmax ) - mulog) / sdlog
            } else {
              Bstar <- (log(extra$max.support) - mulog) / sdlog
            }
            
            Astar <- (log(extra$min.support) - mulog) / sdlog
            
            
            dnum <- dnorm(x = Bstar) - dnorm(x = Astar)
            dden <- pnorm(q = Bstar) - pnorm(q = Astar)
            
            dl.dmulog <- ((log(y) - mulog) / sdlog^2) + 
                                 (1 / sdlog)* (dnum / dden)
            dl.dsdlog <- -1 / sdlog + ((log(y) - mulog)^2) / sdlog^3 + 
              (1 / sdlog) * ( (Bstar * dnorm(x = Bstar) - 
                                 Astar * dnorm(x = Astar) )/ dden )
            
            
          } else {
            # Lognormal untruncated
            dl.dmulog <- (log(y) - mulog)/sdlog^2
            dl.dsdlog <- -1 / sdlog + (log(y) - mulog)^2 / sdlog^3
            
          }
          
          
          ##################################################################
          c(w) * cbind(dl.dmulog * dmulog.deta,
                       dl.dsdlog * dsdlog.deta)
          
        }), list( .lmulog = lmulog, .lsdlog = lsdlog,
                  .emulog = emulog, .esdlog = esdlog,
                  .min.support = min.support,
                  .max.support = max.support ))),
        
        
        
        
        weight = expression({
          
          if (extra$truncation) {
            
            wz <- matrix(0, n, 3)  # Non - Diagonal!
        #comp2 <- ( Bstar * dnorm(x = Bstar) - Astar * dnorm(x = Astar) ) /
              #( pnorm(q = Bstar ) - pnorm(q = Astar) )
            
            # dnum2 <-  dnorm(x = Bstar) - dnorm(x = Astar)
            #  dden2 <- pnorm(q = Bstar) - pnorm(q = Astar)
            comp3 <- ( dnum / dden )^2
            
            comp11 <- ( Bstar * dnorm(x = Bstar) - 
                              Astar * dnorm(x = Astar)) / 
                           (pnorm(q = Bstar) - pnorm(q = Astar))
                    
            comp12 <- ((Astar * dnorm(x = Astar) - 
                                Bstar * dnorm(x = Bstar)) /
                         (pnorm(q = Bstar) - pnorm(q = Astar)))
            
            comp21 <- (dnorm(x = Bstar) * Bstar * (1 - Bstar^2) - 
                          dnorm(x = Astar) * Astar * (1 - Astar^2)) /
              ( pnorm(q = Bstar) - pnorm(q = Astar) ) 
            
            comp22 <- (dnorm(x = Astar) - dnorm(x = Bstar)) / 
              (pnorm(q = Bstar) - pnorm(q = Astar))
            
            comp23 <- (dnorm(x = Bstar) * (1 - Bstar^2) - 
                         dnorm(x = Astar) * (1 - Astar^2)) /
                              (pnorm(q = Bstar) - pnorm(q = Astar))
            
            comp24 <- ( (dnorm(x = Bstar) - dnorm(x = Astar)) * 
                (dnorm(x = Bstar) * Bstar - dnorm(x = Astar) * Astar)) /
              (pnorm(q = Bstar) - pnorm(q = Astar))^2
            
            ned2l.dmulog2 <- 1 / sdlog^2 - (1 / sdlog^2) * (comp11 + comp3)
            ned2l.dsdlog2 <- 2 / sdlog^2 + (3 / sdlog^2) * comp12 + 
              (1 / sdlog^2) * (comp21 + comp11^2) + (1/ sdlog^2) * comp11
            
            ned2l.dsdlogdmulog <- (2 / sdlog^2 ) * comp22 + 
                                      (1 / sdlog^2) * comp23 - 
                                         (1 / sdlog^2) * comp24
            
            wz[, iam(1, 1, M)] <- ned2l.dmulog2 * dmulog.deta^2
            wz[, iam(2, 2, M)] <- ned2l.dsdlog2 * dsdlog.deta^2
            wz[, iam(2, 1, M)] <- ned2l.dsdlogdmulog *
                                    dmulog.deta * dsdlog.deta
            
            
            wz <- c(w) * wz
          
            
          } else {
            
            # Lognormal untruncated
            wz <- matrix(NA_real_, n, 2)  #  Diagonal!
            ned2l.dmulog2 <- 1/sdlog^2
            ned2l.dsdlog2 <- 2 * ned2l.dmulog2
            wz[, iam(1, 1, M)] <- ned2l.dmulog2 * dmulog.deta^2
            wz[, iam(2, 2, M)] <- ned2l.dsdlog2 * dsdlog.deta^2
            
            wz <- c(w) * wz
            
          }
          
          wz
          
        }))
  }  # trunclognormal


## ARMAXff control function.
trunclognormal.control <- function(save.weights = TRUE,
                                   summary.HDEtest = FALSE,
                                   epsilon = 1e-7,...) { 
  #criterion <- "loglikelihood"
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest,
       epsilon = 1e-7,...)
}

