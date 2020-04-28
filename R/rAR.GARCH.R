##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

rAR.GARCH <- function(n, nwarmup = 1500, extended = FALSE,
                      ARmodel = list(drift = NULL, ARcoef = NULL),
                      GARCHmodel = list(omega = NULL, ARCHcoef = NULL,
                                       GARCHcoef = NULL, leverage = NULL),
                      type.TS = c("ARCH", "GARCH",
                                  "Taylor-Schwert", "A-GARCH",
                                  "Log-GARCH", "M-GARCH")[1]) {
  
  tyTS <- match.arg(type.TS, c("ARCH", "GARCH", "IGARCH",
                               "Taylor-Schwert", "A-GARCH",
                               "Log-GARCH", "M-GARCH"))[1]
  rm(type.TS)
  newN <- n + nwarmup
  ar.drift  <- ARmodel$drift
  ar.coef   <- ARmodel$ARcoef
  arch.om   <- GARCHmodel$omega
  arch.coe  <- GARCHmodel$ARCHcoef
  garch.coe <- GARCHmodel$GARCHcoef
  leve.coe  <- GARCHmodel$leverage
  epst <- sigmat <- yt <- numeric(0)
  
  if (sum(c(arch.coe, garch.coe, leve.coe)) > 1 - 1e-5)
    stop("The variance model entered is non-stationary.")
    
  if ((tyTS == "ARCH") || (tyTS == "GARCH")) {
    
    m.length  <- max(c(length(arch.coe), length(garch.coe)))[1]
    sigmat[1:m.length] <- sqrt(arch.om)
    epst[1:m.length]   <- rnorm(n = m.length) * sigmat[1:m.length]
    yt[1:m.length]     <- ar.drift + epst[1:m.length]
    
    ar.coef <- if (!length(ARmodel$ARcoef)) rep(0, m.length) else
                                                           ARmodel$ARcoef
    if (length(arch.coe) < m.length) 
      arch.coe <- c(arch.coe, rep(0, length = m.length - length(arch.coe)))
    
    if (length(garch.coe) < m.length)
      garch.coe <- c(garch.coe,
                     rep(0, length = m.length -length(garch.coe)))
    
    if (length(ar.coef) < m.length)
      ar.coef <- c(ar.coef, rep(0, length = m.length - length(ar.coef)))
    
    for (ii in (m.length + 1):newN) {
      sigmat[ii] <- sqrt(sum(arch.coe *
                                (epst[(ii - 1):(ii - m.length)]^2)) +
          sum(garch.coe * (sigmat[(ii - 1):(ii - m.length)]^2)) + arch.om)
    
      epst[ii] <- sigmat[ii] * rnorm(1)
      yt[ii]   <- sum(ar.coef * yt[(ii - 1):(ii - m.length)]) +
                                               epst[ii] + ar.drift
    }
    
  }
    
    
    
   if (tyTS == "Taylor-Schwert") {
      
     m.length  <- max(c(length(arch.coe), length(garch.coe)))[1]
     sigmat[1:m.length] <- arch.om
     epst[1:m.length]   <- rnorm(n = m.length) * sigmat[1:m.length]
     yt[1:m.length]     <- ar.drift + epst[1:m.length]
      
     ar.coef <- if (!length(ARmodel$ARcoef)) rep(0, m.length) else
                                                        ARmodel$ARcoef
      
     if (length(arch.coe) < m.length)
      arch.coe <- c(arch.coe, rep(0, length = m.length - length(arch.coe)))
      
     if (length(garch.coe) < length(arch.coe))
      garch.coe <- c(garch.coe,
                     rep(0, length = m.length -length(garch.coe)))
      
     if (length(ar.coef) < m.length)
       ar.coef <- c(ar.coef, rep(0, length = m.length - length(ar.coef)))
      
     for (ii in (m.length + 1):newN) {
       sigmat[ii] <- sum(arch.coe * abs(epst[(ii - 1):(ii - m.length)])) +
         sum(garch.coe * (sigmat[(ii - 1):(ii - m.length)])) + arch.om
       
       epst[ii] <- sigmat[ii] * rnorm(1)
       yt[ii]   <- sum(ar.coef * yt[(ii - 1):(ii - m.length)]) +
                                                epst[ii] + ar.drift
     }
   }
      
      
      
   if (tyTS == "A-GARCH") {
     
     m.length  <- max(c(length(arch.coe), length(garch.coe),
                         length(leve.coe)))[1]
     sigmat[1:m.length] <- sqrt(arch.om)
     epst[1:m.length]   <- rnorm(n = m.length) * sigmat[1:m.length]
     yt[1:m.length]     <- ar.drift + epst[1:m.length]
    
     ar.coef <- if (!length(ARmodel$ARcoef)) rep(0, m.length) else
                                                        ARmodel$ARcoef
     if (length(arch.coe) < m.length)
       arch.coe <- c(arch.coe, rep(0, length = m.length - length(arch.coe)))
      
     if (length(garch.coe) < m.length )
       garch.coe <- c(garch.coe,
                      rep(0, length = m.length - length(garch.coe)))
      
     if (length(leve.coe) < m.length )
       leve.coe <- c(leve.coe,
                      rep(0, length = m.length - length(leve.coe)))
     
     if (length(ar.coef) < m.length)
       ar.coef <- c(ar.coef, rep(0, length = m.length - length(ar.coef)))
      
      for (ii in (m.length + 1):newN) {
        sigmat[ii] <- sqrt(sum(
                       arch.coe * (epst[(ii - 1):(ii - m.length)]^2) +
                          leve.coe * epst[(ii - 1):(ii - m.length)]) +
        sum(garch.coe * (sigmat[(ii - 1):(ii - m.length)]^2)) + arch.om)
        
        epst[ii] <- sigmat[ii] * rnorm(1)
        yt[ii]   <- sum(ar.coef * yt[(ii - 1):(ii - m.length)]) +
          epst[ii] + ar.drift
      }
   }
  
   if (tyTS == "Log-GARCH") {
     m.length  <- max(c(length(arch.coe), length(garch.coe)))[1]
     sigmat[1:m.length] <- exp(arch.om)
     epst[1:m.length]   <- rnorm(n = m.length) * sigmat[1:m.length]
     yt[1:m.length]     <- ar.drift + epst[1:m.length]
     
     ar.coef <- if (!length(ARmodel$ARcoef)) rep(0, m.length) else
                                                          ARmodel$ARcoef
     
     if (length(arch.coe) < m.length)
      arch.coe <- c(arch.coe, rep(0, length = m.length - length(arch.coe)))
     
     if (length(garch.coe) < length(arch.coe))
       garch.coe <- c(garch.coe,
                    rep(0, length = m.length -length(garch.coe)))
     
     if (length(ar.coef) < m.length)
       ar.coef <- c(ar.coef, rep(0, length = m.length - length(ar.coef)))
     
     for (ii in (m.length + 1):newN) {
      sigmat[ii] <-
        exp(sum(arch.coe * abs(epst[(ii - 1):(ii - m.length)])) +
            sum(garch.coe * (log(sigmat[(ii - 1):(ii - m.length)]))) +
              arch.om)
         
       epst[ii] <- sigmat[ii] * rnorm(1)
       yt[ii]   <- sum(ar.coef * yt[(ii - 1):(ii - m.length)]) +
                                               epst[ii] + ar.drift
     }
   }
  
  
  
  if (tyTS == "M-GARCH") {
    m.length  <- max(c(length(arch.coe), length(garch.coe)))[1]
    sigmat[1:m.length] <- sqrt(exp(arch.om))
    epst[1:m.length]   <- rnorm(n = m.length) * sigmat[1:m.length]
    yt[1:m.length]     <- ar.drift + epst[1:m.length]
    
    ar.coef <- if (!length(ARmodel$ARcoef)) rep(0, m.length) else
      ARmodel$ARcoef
    
    if (length(arch.coe) < m.length)
      arch.coe <- c(arch.coe, rep(0, length = m.length - length(arch.coe)))
    
    if (length(garch.coe) < length(arch.coe))
      garch.coe <- c(garch.coe,
                     rep(0, length = m.length -length(garch.coe)))
    
    if (length(ar.coef) < m.length)
      ar.coef <- c(ar.coef, rep(0, length = m.length - length(ar.coef)))
    
    for (ii in (m.length + 1):newN) {
      sigmat[ii] <- exp(sqrt(
              sum(arch.coe * log(epst[(ii - 1):(ii - m.length)]^2)) +
              sum(garch.coe * (log(sigmat[(ii - 1):(ii - m.length)]^2))) +
              arch.om))
      
      epst[ii] <- sigmat[ii] * rnorm(1)
      yt[ii]   <- sum(ar.coef * yt[(ii - 1):(ii - m.length)]) +
        epst[ii] + ar.drift
    }
  }
  
  ans <- cbind(yt, sigmat, epst)
  colnames(ans) <- c(tyTS, "sigma", "eps")
  ans <- if (!extended) ans[-(1:nwarmup), tyTS, drop = FALSE] else
                              ans[-(1:nwarmup), , drop = FALSE]
  attr(ans, "list") <- list(ARmodel, GARCHmodel)
  ans
}


if (FALSE) {
  test <- rAR.GARCH(n = 400, extended = FALSE,
            ARmodel = list(drift = 1, ARcoef = c(0.1)),
            GARCHmodel = list(omega = 0.5, ARCHcoef = c(0.5, 0.1), GARCHcoef = c(0) ) )
  head(test)
  
  garch.data <- data.frame(y = test)
  colnames(garch.data) <- "y"
  head(garch.data)
  dim(garch.data)
  fit1 <- vglm(y ~ 1, AR.GARCHff(ARorder = 1, GARCHorder = c(2, 0),
                                 type.TS = "ARCH"),
               trace = TRUE, data = garch.data)
  ( mycoef <- coef(fit1, matrix = TRUE))
  c(mycoef[1, 1], mycoef[,2])   # True values are c(0, 0.15, 0.14, 0.31)
  summary(fit1)
  
  
  
  
  
  
  
  
  set.seed(20170623)
  nn <- 550
  ar <- 0    # No autoregressive component.
  omega <- 0.23  # Intercept
  alpha <- c(0.14, 0.31)
  beta  <- 0.19
  
  ## Generate the data using garchSim()  ##
  spec <- garchSpec(model = list(omega = omega, 
                                 alpha = alpha,  beta = beta))
  
  garch.data <- data.frame(y = garchSim(spec, n = nn))
  rownames(garch.data) <- NULL
  colnames(garch.data) <- "y"
  head(garch.data)
  
  
  test <- rAR.GARCH(n = nn, extended = FALSE,
                    ARmodel = list(drift = 1, ARcoef = c(0)),
                    GARCHmodel = list(omega = 0.23,
                                      ARCHcoef = alpha, GARCHcoef = beta ),
                    type.TS = "GARCH")
  head(test)
  
  garch.data <- data.frame(y = test)
  colnames(garch.data) <- "y"
  head(garch.data)
  dim(garch.data)
  
  # Fit GARCH(2, 1) with vglm()
  fit1 <- vglm(y ~ 1, AR.GARCHff(ARorder = 0, GARCHorder = c(2, 1),
                                 type.TS = "GARCH"),
               trace = TRUE, data = garch.data)
  (mycoef <- coef(fit1, matrix = TRUE))
  
  #summary(fit1)
  #constraints(fit1, mat = T)
  
  ## Fit the model with garchFit()
  garch.Fit <- garchFit(formula = ~ arma(0, 0) + garch(2, 1),
                        data = garch.data[, c("y"), drop = FALSE],
                        trace = FALSE, include.mean = TRUE)
  garch.Fit@fit$coef
  c(mycoef[1, 1], mycoef[,2])   # True values are c(0, 0.23, 0.14, 0.31, 0.19)
  
  
}