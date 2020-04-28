##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

extract.Residuals <- function(object, TSprocess,...) {
  if (!is.character(TSprocess))
    stop("Invalid TS process entry")
  TSprocess <- match.arg(TSprocess, c("ARIMA", "ARMA", "MA", "AR"))[1]
  
  if (!inherits(object, "vglm"))
    stop("Currently, only 'vglm' objects handled.")
  
  object@extra$ResArma
  
}







fittedVGAMextra <- function(object, ...) {
  
  if (!inherits(object, "vglm"))
    stop("Currently, only 'vglm' objects handled.")
  
  object@fitted.values
  
}






weightsVGAMextra <- function(object, type.w = "prior", ...) {
  
  type.w <- match.arg(type.w, c("prior"))[1]
  if (!inherits(object, "vglm"))
    stop("Currently, only 'vglm' objects handled.")
  myw <- if (type.w == "prior") object@prior.weights else object@weights
  if(!length(myw))
    myw <- matrix(weights(object), object@misc$n, object@misc$NOS)
  
  myw
  
}