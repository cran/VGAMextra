##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.


is.FormulaAR <- function(Model = ~ 1, 
                         Resp  = 1) {
  if ( length(Resp) && 
         ( if (length(Model)) (class(Model) == 'formula') else TRUE) )
    TRUE else FALSE
}