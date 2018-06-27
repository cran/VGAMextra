##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
# Supports the Wald, score, and lrt tests (20180209)


VGLM.INGARCHff <-  function(Order = c(1, 1),
                            dist.type = c("poisson", "negbinomial",
                                        "logarithmic", "yulesimon")[1],
                          link = c("loge", "identitylink", "negloge",
                                   "reciprocal", "logit")[1],
                          interventions = list(),
                          lagged.fixed.means = NULL,
                          lagged.fixed.obs   = NULL,
                          f.transform.Y = NULL,
                          transform.lambda = FALSE,
                          init.p.ARMA = NULL, ...) {
  
  init.p    <- init.p.ARMA; rm(init.p.ARMA)
  fixed.mean <- lagged.fixed.means; rm(lagged.fixed.means)
  fixed.obs  <- lagged.fixed.obs; rm(lagged.fixed.obs)
  
  #link <- match.arg(link, c("loge", "identitylink", "negloge",
  #                          "reciprocal", "logit"))[1]
  
  dist.type <- match.arg(dist.type, c("poisson", "negbinomial",
                                      "logarithmic", "yulesimon"))[1]
  
  if (length(f.transform.Y) && !is.function(f.transform.Y))
    stop("Wrong input for argument 'f.transform.Y'. Must be a function.",
         " Enter NULL for the identity function.")
  
  if (!is.logical(transform.lambda))
    stop("Argument 'transform.lambda' must be logical.")
  
  if (FALSE) {
    if ( (dist.type == "negbinomial") && !(link == "loge") )
      warning("Special care needed with the 'Negative Binomial' ",
              "distribution using links \n  other than 'loge'.",
              " Try another suitable link if numerical issues arise. ")
    
    if ( (dist.type == "logarithmic") && !(link == "logit") )
      warning("Special care needed with the 'logarithmic' distribution",
              "  using links \n  other than 'logit'.",
              " Try another suitable link if numerical issues arise. ")
    
    if ( (dist.type == "yulesimon") && !(link == "loge") )
      warning("Special care needed with the 'yulesimon' distribution",
              "  using links \n  other than 'loge'.",
              " Try another suitable link if numerical issues arise. ")
  }
  
  dist.type <- if (dist.type == "negbinomial" ) "NegBinom" else dist.type
  
  if (length(interventions)) {
    
    if (any(interventions$tau == 0))
      stop("Bad input for 'tau' in argument 'interventions'.")
    
    if (!length(interventions$No.Inter)) {
      interventions$No.Inter <- TRUE
      #warning("No interactions considered. Change this through the ",
      #        "argument 'No.Inter'.")
    } else {
      if (!is.logical(interventions$No.Inter))
        stop("Wrong input for argument 'No.Inter'.")
    }
        
    if (!is.list(interventions) || (length(interventions) != 3))
      stop("Bad input for argument 'interventions'.")
    
    if (any(interventions$delta < 0) || any(interventions$delta > 1))
      stop("Bad input for 'delta' in argument 'interventions'.")
  }
  
  
  if (length(fixed.mean) && !is.vector(fixed.mean))
    stop("Wrong input for argument 'lagged.fixed.means'.")
  
  if (length(fixed.obs) && !is.vector(fixed.obs))
    stop("Wrong input for argument 'lagged.fixed.obs'.")
  
  if (length(init.p) && !Is.Numeric(init.p, isInteger = TRUE, 
                                    length.arg = 1))
    stop("Bad input for argument 'init.p.ARMA'")
  
  if (!is.character(dist.type))
    stop("Wrong input for argument 'data.type'")
  
  TSfflist <- c("poissonTSff", "NegBinomTSff",
                "logarithmicTSff", "yulesimonTSff")
  
  ts.family <- grep(dist.type, TSfflist)
  ts.family <- get(TSfflist[ts.family], envir = .GlobalEnv )
  
  answer <- ts.family(Order = Order,
                      link = link ,
                      f.transform.Y = f.transform.Y,
                      transform.lambda = transform.lambda,
                      init.p.ARMA = init.p,
                      lagged.fixed.means = fixed.mean,
                      lagged.fixed.obs   = fixed.obs,
                      interventions = interventions, ...)
  
  if (!length(slot(answer, "typeTS"))) 
    slot(answer, "typeTS") <- dist.type
  
  if (!length(slot(answer, "vfamily")))
    slot(answer, "vfamily") <- "VGLM.INGARCHff"
  

  answer
  
}


rINGARCH <- function(n, type.INGARCH = c("INGARCH", "NegLog-INGARCH",
                                         "Recip-INGARCH")[1],
                     Order = c(0, 0), intercept = NULL,
                     alpha.coe = NULL, gamma.coe = NULL,
                     beta.coe = NULL, covs = NULL,
                     burnin = 400) {
  
  nn   <- n + burnin ; rm(n)
  ord1 <- Order[1]
  ord2 <- Order[2]
  maxord <- max(Order)[1]; rm(Order)
  alpha   <- alpha.coe; rm(alpha.coe)
  beta    <- gamma.coe; rm(gamma.coe)
  coeCovs <- beta.coe;  rm(beta.coe)
  beta0   <- intercept; rm(intercept)
  
  type.INGARCH <- match.arg(type.INGARCH, c("INGARCH", "NegLog-INGARCH",
                                            "Recip-INGARCH"))
  
  if (sum(c(alpha, beta, coeCovs)) > 1 - 1e-5)
    stop("The variance model entered is non-stationary.")
  
  if (length(covs) && !is.matrix(covs))
    stop("Argument 'covs' must be a matrix.")
  
  if (length(covs)) {
    covs.2 <- matrix(0, nn, NCOL(covs))
    covs.2[-(1:burnin), ] <- covs; rm(covs)
  } else {
    covs.2 <- matrix(0, nn, 1)
    coeCovs <- exp(0)
  }
  
  if (maxord == 0)
    stop("Only INARCH and INGARCH models handled.",
         "For Poisson regression refer to other VGLM family functions.")
  
  if (!length(alpha) && ord1)
    stop("Enter 'alpha' coefficients ")
  
  if (!length(beta) && ord2)
    stop("Enter 'beta' coefficients ")
  
  if (ord1 == 0) {
    alpha <- rep(0, maxord)
    ord1  <- maxord
  } 
  
  if (ord2 == 0) {
    beta <- rep(0, maxord)
    ord2  <- maxord
  } 
  
  p.data   <- numeric(0)
  lambda.t <- numeric(0)
  
  if (type.INGARCH == "INGARCH") {
    lambda.t[1:(maxord + 1)] <- beta0
    p.data[1:(maxord + 1)]   <- rpois(1, lambda.t[1:(maxord + 1)])
    
    for (ii in (maxord + 2):nn) {
      lambda.t[ii] <- beta0 + sum(alpha * p.data[(ii - 1):(ii - ord1)]) +
        sum(beta  * lambda.t[(ii - 1):(ii - ord2)]) +
        coeCovs * covs.2[ii, ]
      p.data[ii]   <- rpois(1, lambda.t[ii])
    }
  }
  
  if (type.INGARCH == "Recip-INGARCH") {
    lambda.t[1:(maxord + 1)] <- exp(beta0)
    p.data[1:(maxord + 1)]   <- rpois(1, lambda.t[1:(maxord + 1)])
    
    for (ii in (maxord + 2):nn) {
      lambda.t[ii] <- (beta0 + sum(alpha * p.data[(ii - 1):(ii - ord1)]) +
                         sum(beta  * lambda.t[(ii - 1):(ii - ord2)]) +
                         coeCovs * covs.2[ii, ])^(-1)
      p.data[ii]   <- rpois(1, lambda.t[ii])
    }
  }
  
  if (type.INGARCH == "NegLog-INGARCH") {
    lambda.t[1:(maxord + 1)] <- exp(-beta0)
    p.data[1:(maxord + 1)]   <- rpois(1, lambda.t[1:(maxord + 1)])
    
    for (ii in (maxord + 2):nn) {
      lambda.t[ii] <- 
        exp(-( beta0 + sum(alpha * p.data[(ii - 1):(ii - ord1)]) +
                 sum(beta  * lambda.t[(ii - 1):(ii - ord2)]) ) +
              coeCovs * covs.2[ii, ])
      p.data[ii]   <- rpois(1, lambda.t[ii])
    }
  }
  
  # Remove burn-in data
  p.data <- cbind(p.data[-(1:burnin)])
  #my.vch <- timeSequence(timeDate(format(Sys.time(),
  #                                       format = "%Y-%m-%d")) -
  #                         NROW(p.data) * 3600 * 24,
  #                       length.out = NROW(p.data))
  #p.data <- timeSeries(p.data, charvec = my.vch)
  colnames(p.data) <- c("PoissonTS")
  attr(p.data, "lambda") <- lambda.t[-(1:burnin)]
  
  p.data
  
}



VGLM.INGARCHff.control <- function(save.weights = TRUE, 
                                   summary.HDEtest = FALSE,
                                   ...) {
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}
