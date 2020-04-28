##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

KPSS.test <- function(x, type.H0 = c("level", "trend")[1],
                      trunc.l = c("short", "large")[1],
                      show.output = TRUE) {
  
  my.class <- match.arg(attributes(x)$class, "ts")
  trunc.l <- match.arg(trunc.l, c("short", "large"))[1]
  trunc.l <- (trunc.l == "short") 
  t.st <- match.arg(type.H0, c("level", "trend"))[1]
  t.st <- which(t.st == c("level", "trend"))[1]
  x.fin <- cbind(x); rm(x)
  
  if (!is.logical(show.output))
    stop("Wrong input for argument 'show.output'")
  
  if (ncol(x.fin) != 1)
    stop("At present, only univariate TSs handled.")
  
  nn <- nrow(x.fin)
  x.fin <- data.frame(t = 1:nn, x = x.fin)
  my.res <- switch(t.st,
                   with(x.fin, residuals(lm(x ~ 1))),
                   with(x.fin, residuals(lm(x ~ t))))
  
  #var.rs <- sum(my.res^2) / nn
  s2l    <- sum(my.res^2) / nn + 
             (2 / nn) * con.est.s2l(resids = my.res,
                        l.short = trunc.l)
  eta.fin <- (sum(cumsum(my.res)^2) / nn^2) / s2l

  crit.val <- if (type.H0 == "level")
    matrix(c(0.347, 0.463, 0.573, 0.739), 1, 4, byrow = TRUE) else
      matrix(c(0.119, 0.146, 0.176, 0.216), 1, 4, byrow = TRUE) 
  colnames(crit.val) <- c("10%", "5%", "2.5%", "1%")
  rownames(crit.val) <- c("Critical value ")
  
  if (show.output) {
    cat("\n",
        "H0:", type.H0 ,"stationary vs. H1: Unit root. 
      Test statistic:", eta.fin, "\n\n")
    
    cat("p-value:", pVal.KPSS.test(x =c(crit.val),
                                   y = c(0.10, 0.05, 0.025, 0.01),
                                   xout = eta.fin))
    cat("\nUpper tail percentiles:\n")
    print(crit.val)
  }
  
  invisible(structure(list(crit.value = crit.val, 
            pvalue =  pVal.KPSS.test(x =c(crit.val),
                             y = c(0.10, 0.05, 0.025, 0.01),
                             xout = eta.fin),
            resids = my.res)))
  
}

# Consistent estimator of s^2(l). Only second term.
# See equation (10) in Kwiatkowski, Philips, Shin (1992)
con.est.s2l <- function(resids, l.short = NULL) {
  
  if (length(l.short) && !is.logical(l.short))
    stop("Bad input for argument 'l.short'.")
  
  nn.i <- length(c(resids))
  bigL <- if (l.short)  ceiling(3 * sqrt(nn.i) /11) else
    ceiling(9 * sqrt(nn.i) /11)
  
  # Term 1 not needed. Given by kpss.test
  # Term 2
  lag.rb <- WN.lags(y = cbind(resids), lags = bigL)
  res.p <- apply(lag.rb, 2, function(x) {
    sum(x * resids)
  })
  
  sum(res.p * (1 - (1:bigL)/(bigL + 1)))

}


## To compute the p--values
pVal.KPSS.test <- function(x, y, xout) {
  first <- spline(x = x, y = y, xout = xout)$y
  vec.neg <- which(first < 0)
  
  if (length(vec.neg))
    first[vec.neg] <- 1e-16
  
  first
}
