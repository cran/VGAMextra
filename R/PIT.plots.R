##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

######## Non-randomized PIT histogram.
setGeneric("typeTS" , function(object)
                standardGeneric("typeTS"))
setMethod("typeTS", signature(object = "vglm"),
          function(object) object@family@typeTS)


responseVGAMex <- function(object,...) UseMethod("responseVGAMex", object)
PIT <- function(object,...) UseMethod("PIT", object)


# Extractor function. # vgltsm not a class yet.
responseVGAMex.vglm <- function(object, ...) {
  if (!inherits(object@family, "vgltsmff"))
    stop("Currently, only 'vglm' objects handled.")
  object@y
}


PIT.vglm <- function(object, bins = 10, ...){
  
  distr.type <- typeTS(object)
  if ((distr.type == "poisson") || (distr.type == "negbinomial")) {
    addit.coefs <- object@misc$additional.coe
    pred.shapes <- NULL
  } else {
    addit.coefs <- NULL
    pred.shapes <- object@misc$predictors
  }
    
  
  execute.PIT(y = c(responseVGAMex(object)),
              pred.values = c(fitted.values(object)),
              pred.shapes = pred.shapes,
              distr.type  = typeTS(object), #typeTS(object).. was poisson,
              #object@misc$additional.coe, was NULL
              addit.coefs = addit.coefs,
              bins = bins, ...)
  
}

execute.PIT <- function(y, pred.values, pred.shapes,
                        distr.type = c("poisson", "negbinomial",
                                       "logarithmic", "yulesimon")[1],
                        addit.coefs = NULL,
                        bins = 10, ...){
  
  distr.type <- match.arg(distr.type, c("poisson", "negbinomial",
                                        "logarithmic", "yulesimon"))[1]
  
  nn  <- length(pred.values)
  uu  <- seq(0, 1, length = bins + 1)
  pit <- numeric(length(uu))
  pit.y <- y; rm(y)
  
  ## Check zero values in the response.
  pit.y.min1 <- pit.y - 1
  check.zes <- which(pit.y <= 1e-5)
  
  if (length(check.zes)) 
    pit.y.min1[check.zes] <- 0
  
  qq.t <- switch(distr.type,
                 "poisson"  = ppois(q = pit.y,
                                    lambda = pred.values),
                 "negbinomial" = pnbinom(q = pit.y,
                                         size = addit.coefs,
                                         mu = pred.values),
                 "logarithmic" = plog(q = pit.y,
                                      shape = pred.shapes),
                 "yulesimon" = pyules(q = pit.y,
                                      shape = pred.shapes))
  
  qq.tmin1 <- switch(distr.type,
                     "poisson"  = ppois(q = pit.y.min1,
                                        lambda = pred.values),
                     "negbinomial" = pnbinom(q = pit.y.min1,
                                             size = addit.coefs,
                                             mu = pred.values),
                     "logarithmic" = plog(q = pit.y.min1,
                                          shape = pred.shapes),
                     "yulesimon" = pyules(q = pit.y.min1,
                                          shape = pred.shapes))
  
    
  my.pit <- t(cbind(qq.tmin1, qq.t))
  my.pit <- apply(my.pit, 2, function(x) punif(uu, x[1], x[2])/nn)
  
  my.pit <- rowSums(my.pit)
  histo.out <- list(breaks  = uu,
                    counts  = diff(my.pit) * nn,
                    density = diff(my.pit) * bins,
                    mids    = (uu[-(bins + 1)] + uu[-1]) / 2,
                    xname   = NULL,
                    equidist = TRUE)
  
  class(histo.out) <- "histogram"
  plot_args <- modifyList(list(main = "Non-randomized PIT Histogram",
                               xlab = "Probability Integral Transform",
                               ylab = "Density",
                               freq = FALSE,
                               ylim = range(0, histo.out$density)),
                          list(...))
  
  do.call("plot", args = c(list(x = histo.out), plot_args ))
  abline(h = 1, lty = "dashed", col = "red")
}

