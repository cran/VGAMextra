\name{vgltsmff}
\docType{class}
\alias{vgltsff-class}
\alias{vgtsff-class}
\alias{vgltsmff-class}
\alias{ARMAvgltsmff}


\title{Class of Vector Generalized Linear Time Series Models}
\description{ Time series family functions for the \pkg{VGAMextra} package}

\section{Objects from the Class}{
  Objects can be created by calling 
  \code{new("vgltsmff"...)}

}

\section{slots}{

  
Implementation of vector generalized linear time series (TS) family 
functions (\emph{vgltsff}) at \pkg{VGAMextra} is entirely based on 
the structure of family functions of the class 
\code{\link[VGAM:vglmff-class]{vglmff-class}}.

Hence, refer to \code{\link[VGAM:vglmff-class]{vglmff-class}} for a 
thourugh description of slots and features involved when objects 
of class \code{"vgtsff"} are being created.

}

\section{Methods}{

 Thus far, the following methods for objects of class
 \code{"vgltsff-class"} are implemented:
 
 \describe{
 
  \item{\code{summary}}{
  
  Additional information to that displayed by the \code{summary}
  methods from \pkg{VGAM}. That is: 
  
  a) Standard errors  based on the MLEs asymptotic distributions, and
  
  b) Checks on stationarity and/or invertibility via the polynomial roots.
     
  
  
  Currently, summary methods at \pkg{VGAMextra} have been 
  implemented for:
  
      \describe{
          \item{\code{signature(VGAMff = "ARff")}:}{
           For ARX--types family functions.}
          \item{\code{signature(VGAMff = "MAff")}:}{  
          For MAX--types  family functions.}
          \item{\code{signature(VGAMff = "ARMAff")}:}{  
          For ARMAX--like family functions.}
      }
  }
  See 
  \code{\link[VGAMextra:summaryS4VGAMextra]{summaryS4VGAMextra}}
  for further details.
 }
}
\note{


Programmers to write VGAM/VGLM time series family functions are also
allowed to write methods functions either for specific purposes, or to 
extend those current methods to print some extra output required.

In such cases, notice that the class \code{vgltsff-class} is 
\emph{labeled} by an object of class \code{"character"} 
(a character vector) specified at the slot
\code{@vfamily} within the family function. This is, in fact, 
one of the required slots by the class 
\code{\link[VGAM:vglmff-class]{vglmff-class}}.

Additionally, practitioners are encouraged to mantain all previous 
conventions for naming the arguments in Ts family functions
as specified at
\code{\link[VGAM:vglmff-class]{vglmff-class}}, e.g., 
\code{link} is the argument for parameter link functions, etc.

}

\author{
V. Miranda and T.W. Yee.

}


