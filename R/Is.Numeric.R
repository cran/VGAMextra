##########################################################################
# These functions are
# Copyrigth (C) 2014-2018 V. Miranda and T. W. Yee. University of Auckland
# function Is.Numeric


Is.Numeric <- function(x, isInteger  = FALSE, 
                       length.arg = NULL, 
                       Nnegative  = NULL)  {
  
  if ( all( is.numeric(x) )  &&
       ( if ( isInteger ) all( x%%1 == 0 )  else TRUE )   &&
       ( if ( length( length.arg ) ) 
         (length(x) == length.arg )  else TRUE )  &&
       ( if ( length(Nnegative) && is.logical(Nnegative) ) 
         if (Nnegative) all(x >= 0)  else all(x < 0) else TRUE)) 
    TRUE else FALSE
}
    