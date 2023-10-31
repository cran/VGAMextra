####################################################################
# These functions are
# Copyright (C) 1998-2023
# T.W. Yee, University of Auckland.
# V. Miranda, Siqi Liu, Auckland University of Technology
# All rights reserved.


truncnormal <-
  function(lmean = "identitylink", lsd = "loglink",
           min.support = -Inf, max.support = Inf, 
           zero = "sd") {
    
    
    lmu <- as.list(substitute(lmean))
    emu <- link2list(lmu)
    lmu <- attr(emu, "function.name")
    
    lsd <- as.list(substitute(lsd))
    esd <- link2list(lsd)
    lsd <- attr(esd, "function.name")

    
    if (max.support < min.support)
      stop("Bad input for 'max.support'")
    
    flag.max.support <- FALSE
    
    if (is.infinite(max.support)) {
      flag.max.support <- TRUE
    }
    
    flag.min.support <- FALSE
    
    if (is.infinite(min.support)) {
      flag.min.support <- TRUE
    }

  
    new("vglmff",
        blurb = c("Doubled truncated normal ",
                  "distribution\n\n",
                  "Links:    ",
                  namesof("mean", lmu, earg = emu, tag = TRUE),
                  ", ",
                  namesof("sd",   lsd, earg = esd, tag = TRUE),
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
               lmean = .lmean ,
               lsd   = .lsd ,
               expected = TRUE,
               min.support = .min.support,
               max.support = .max.support,
               multipleResponses = FALSE,
               parameters.names = c("mean", "sd"),
               zero = .zero )
        }, list( .zero = zero,
                 .lmean = lmean,
                 .lsd   = lsd,
                 .min.support = min.support,
                 .max.support = max.support,
                 .flag.max.support = flag.max.support,
                 .flag.min.support = flag.min.support
        ))),
        
        
        first = eval(substitute(expression({
          
          
          int.min.support <- .min.support
          int.max.support <- .max.support 
          
          extra$truncation <- FALSE # Ordinary normal
          if ( (int.min.support > -Inf) || 
                          !is.infinite( .max.support ))
            extra$truncation <- TRUE
          
          extra$colmeans <- numeric()
          
          
        }), list( .min.support = min.support ,
                  .max.support = max.support ,
                  .flag.max.support = flag.max.support,
                  .flag.min.support = flag.min.support))),
        
        
        
        initialize = eval(substitute(expression({
          
          w.y.check(w = w, y = y,
                    Is.positive.y = FALSE)
          
          
          ncoly <- 1 # ncol(yy) so far, one response.
          M1 <- 2
          
          extra$min.support <- matrix( .min.support , n, ncoly, 
                                       byrow =  TRUE)
          
          extra$max.support <- matrix( .max.support , n, ncoly, 
                                       byrow =  TRUE)
          
         
          if (any(y < extra$min.support))
            stop("Some response values less than 'min.support'.")
          
          if (any(y > extra$max.support ))
            stop("Some response values greater than 'max.support'.")
          
          
          my.max.support <- .max.support
          my.flag <- .flag.max.support
          
          if (my.flag) {
            my.max.support <-   .Machine$integer.max 
            extra$max.support <- matrix(my.max.support , n, ncoly, 
                                         byrow =  TRUE)
          }
          
          
          
          my.min.support <- .min.support
          flag.min <- .flag.min.support
          
          if (flag.min) {
            my.min.support <- (-1) *.Machine$integer.max #2.225074e-100
            extra$min.support <- matrix(my.min.support, n, ncoly, 
                                         byrow =  TRUE)
          }
            


          predictors.names <-
            c(namesof("mean", .lmu , .emu , tag = FALSE),
              namesof("sd",   .lsd , .esd , tag = FALSE))
          
          if (!length(etastart)) {
            mylm <- lm.wfit(x = x, y = c(y), w = c(w))
            sd.y.est <- sqrt( sum(c(w) * mylm$resid^2)
                                 / mylm$df.residual )
            
            #print(head(mean(y)))
            #loglike <- function(x, theta) {
            #  sum(-dtruncnormal(x = x, A = 7.4, B = 12.8, 
            #                    mean = 10.09, sd = theta, 
            #                    log = TRUE))
            #}
  
            #tt <- optimize(f = loglike, 
            #         interval = c(1, 4), x = y)$minimum
            
            etastart <- cbind(
              mean = rep_len(theta2eta(fitted(mylm) , .lmu ,
                                          earg = .emu ), n),
              sd   = rep_len(theta2eta(sd.y.est, .lsd ,
                                          earg = .esd ), n))
          }
          
          
        }), list( .lmu = lmu, .lsd = lsd,
                  .emu = emu, .esd = esd,
                  .min.support = min.support,
                  .max.support = max.support ,
                  .flag.max.support = flag.max.support,
                  .flag.min.support = flag.min.support))),
        
        
        linkinv = eval(substitute(function(eta, extra = NULL) {
          
          mu <- eta2theta(eta[, 1], .lmu , earg = .emu )
          sd <- eta2theta(eta[, 2], .lsd , earg = .esd )
          
          
          qdiffmax <- (extra$max.support - mu) / sd - sd
          qdiffmin <- (extra$min.support - mu) / sd - sd
          

          numdiff <- pnorm(q = qdiffmax) - pnorm(q = qdiffmin)
          dendiff <- pnorm(q = qdiffmax + sd) -
            pnorm(q = qdiffmin + sd)
          
          ## Mean of a truncated normal in [min.support, max.support]
          mu + sd * ( numdiff / dendiff )
          
        }, list( .lmu = lmu, .lsd = lsd,
                 .emu = emu, .esd = esd ))),
        
        
        
        
        last = eval(substitute(expression({
          
          
          misc$link <-    c("mean" = .lmu , "sd" = .lsd )
          
          misc$earg <- list("mean" = .emu , "sd" = .esd )
          
          misc$expected <- TRUE
          
          if (!extra$truncation)
            cat("\n",
                "No truncation limits entered - Fitting an ordinary",
                "normal distribution.")
          
        }), list( .lmu = lmu, .lsd = lsd,
                  .emu = emu, .esd = esd,
                  .min.support = min.support,
                  .max.support = max.support ))),
        
   
        
        
        loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta,
                   extra = NULL,
                   summation = TRUE) {
            
            mu <- eta2theta(eta[, 1], .lmu , earg = .emu )
            sd <- eta2theta(eta[, 2], .lsd , earg = .esd )
            
            
            my.max.support <- .max.support
            my.flag <- .flag.max.support
            if (my.flag) {
              Bstar <- extra$max.support
            } else {
              Bstar <- cbind((my.max.support - mu) / sd)
            }
            
            my.min.support <- .min.support
            flag.min <- .flag.min.support
            if(flag.min) {
              Astar <- extra$min.support
            } else {
              Astar <- cbind((my.min.support - mu) /sd)
            }
            
            
            if (residuals) {
              stop("loglikelihood residuals not implemented yet")
            }
            else {
              ll.elts <- c(w) * ( dtruncnorm(x = y,
                          mean = mu, sd = sd,
                          min.support = my.min.support,
                          max.support = my.max.support, log = TRUE))
              
            # Removed on August 2023.
            #ll.elts <- c(w) * ( dnorm(y, mu, sd = sd, log = TRUE) -
            #          log( pnorm(q = Bstar, lower.tail = TRUE) -
            #             pnorm(q = Astar, lower.tail = TRUE) ))
              
              if (summation) {
                sum(ll.elts)
              } else {
                ll.elts
              }
            }
          }, list( .lmu = lmu, .lsd = lsd,
                   .emu = emu, .esd = esd,
                   .min.support = min.support,
                   .max.support = max.support ,
                   .flag.max.support = flag.max.support,
                   .flag.min.support = flag.min.support))),
        
     
  
        vfamily = c("truncnormal"),
        
        
        
        validparams = eval(substitute(function(eta, y, extra = NULL) {
          
          mu <- eta2theta(eta[, 1], .lmu , earg = .emu )
          sd <- eta2theta(eta[, 2], .lsd , earg = .esd )
          okay1 <- all(is.finite(mu)) &&
            all(is.finite(sd)) && all(0 < sd)
          okay1
        }, list( .lmu = lmu, .lsd = lsd,
                 .emu = emu, .esd = esd,
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
          
          mu <- eta2theta(eta[, 1], .lmu , earg = .emu )
          sd <- eta2theta(eta[, 2], .lsd , earg = .esd )
          
          dmu.deta <- dtheta.deta(mu, .lmu , earg = .emu )
          dsd.deta <- dtheta.deta(sd, .lsd ,   earg = .esd )
  
          
          
          ##################################################################
          
          my.max.support <- .max.support
          my.flag <- .flag.max.support
          if (my.flag) {
            Bstar <- extra$max.support
          } else {
            Bstar <- cbind((my.max.support - mu) / sd)
          }
          
          my.min.support <- .min.support
          flag.min <- .flag.min.support
          if(flag.min) {
            Astar <- extra$min.support
          } else {
            Astar <- cbind((my.min.support - mu) /sd)
          }
          
          
          
          if (extra$truncation) {
            
            
            dnum <- dnorm(x = Bstar) - dnorm(x = Astar)
            dden <- pnorm(q = Bstar) - pnorm(q = Astar)
            
            
            dl.dmu <- (y - mu) / sd^2 + 
                                 (1 / sd)* (dnum / dden)
            dl.dsd <- -1 / sd + (y - mu)^2 / sd^3 + 
              (1 / sd) * ( (Bstar * dnorm(x = Bstar) - 
                                 Astar * dnorm(x = Astar) )/ dden )
            
            
          } else {
            
            # normal untruncated
            dl.dmu <- (y - mu) / sd^2
            dl.dsd <- -1 / sd + (y - mu)^2 / sd^3
            
          }
          
          
          ##################################################################
          testtt <- c(w) * cbind(dl.dmu * dmu.deta,
                        dl.dsd * dsd.deta)
          #print(colSums(testtt))
          
          extra$colmeans <- c(extra$colmeans, colMeans(testtt))
          
          c(w) * cbind(dl.dmu * dmu.deta,
                       dl.dsd * dsd.deta)
        }), list( .lmu = lmu, .lsd = lsd,
                  .emu = emu, .esd = esd,
                  .min.support = min.support,
                  .max.support = max.support ,
                  .flag.max.support = flag.max.support,
                  .flag.min.support = flag.min.support))),
        
        
        
        
        weight = expression({

          if (extra$truncation) {
            
            wz <- matrix(0, n, 3)  # Non - Diagonal!
            
            elogyminmu <- sd * (dnorm(x = Astar) - dnorm(x = Bstar)) / (
              pnorm(q = Bstar) - pnorm(q = Astar))
            
            elogyminmu2 <- sd^2 * (1 - (dnorm(x = Bstar) * Bstar - 
                                             dnorm(x = Astar) * Astar) / (
                                    pnorm(q = Bstar) - pnorm(q = Astar)))
            
            
            ned2l.dmu2 <- 1 / sd^2 - 
              (1 / sd^2) * B11(Astar, Bstar, sd)
            
            ned2l.dsd2 <- (-1 / sd^2) + (3 / sd^4) * 
              elogyminmu2 + B22(Astar, Bstar, sd)
            
            
            ned2l.dsddmu <- (2 / sd^3) * elogyminmu +
              B21(Astar, Bstar, sd)
            
            
            wz[, iam(1, 1, M)] <- ned2l.dmu2 * dmu.deta^2
            wz[, iam(2, 2, M)] <- ned2l.dsd2 * dsd.deta^2
            wz[, iam(2, 1, M)] <- ned2l.dsddmu * dmu.deta * dsd.deta
            
            
            wz <- c(w) * wz
        
            
          } else {

            # normal untruncated
            wz <- matrix(NA_real_ , n, 2)  #  Diagonal!
            ned2l.dmu2 <- 1/sd^2
            ned2l.dsd2 <- 2 * ned2l.dmu2
            wz[, iam(1, 1, M)] <- ned2l.dmu2 * dmu.deta^2
            wz[, iam(2, 2, M)] <- ned2l.dsd2 * dsd.deta^2
            
            wz <- c(w) * wz
            
          }
          
          
          wz
          
        }))
  }  # truncnormal


B11 <- function(Astar, Bstar, sd) {
  
  comp11 <- ( Bstar * dnorm(x = Bstar) - Astar * dnorm(x = Astar) ) /(
               pnorm(q = Bstar) - pnorm(q = Astar) )
  
  comp22 <- ( ( dnorm(x = Bstar) - dnorm(x = Astar) ) / (
                pnorm(q = Bstar) - pnorm(q = Astar) ) )^2
  
  comp11 + comp22

}


B21 <- function(Astar, Bstar, sd) {
  
  
  
  comp23 <- ( dnorm(x = Bstar) * (1 - Bstar^2) - 
               dnorm(x = Astar) * (1 - Astar^2) ) / (
                 pnorm(q = Bstar) - pnorm(q = Astar) )
  
  comp24 <- ( ( dnorm(x = Bstar) - dnorm(x = Astar) ) * 
              ( dnorm(x = Bstar) * Bstar - dnorm(x = Astar) * Astar ) ) / (
                  pnorm(q = Bstar) - pnorm(q = Astar) )^2
  
  (1 / sd^2) * comp23 - (1 / sd^2) * comp24
  
}


B22 <- function(Astar, Bstar, sd) {
  
  

 
  comp21 <- ( dnorm(x = Bstar) * Bstar * (1 - Bstar^2) - 
               dnorm(x = Astar) * Astar * (1 - Astar^2) ) /
                  ( pnorm(q = Bstar) - pnorm(q = Astar) ) 
  
  comp11 <- ( Bstar * dnorm(x = Bstar) - 
                Astar * dnorm(x = Astar) ) / 
                  ( pnorm(q = Bstar) - pnorm(q = Astar) )

  
  (1 / sd^2) * (comp21 + comp11^2) + (1/ sd^2) * comp11

}

