\name{HKdata}
\alias{HKdata}
\docType{data}
\title{Air pollution and hospital admissions due to respiratory and 
cardiovascular causes, Hong Kong.}
\description{
 Daily air pollution levels and
 hospital admissions due to respiratory and 
cardiovascular causes,
between 1 January 1994 and 31 December 1997,
Hong Kong. 
 
}
\usage{data(HKdata)}
\format{
This is a subset of the data
analyzed in Xia and Tong (2006) which 
stores the following time series:
  \describe{
  \item{no}{
  Time vector, 
  
  }
  \item{cardio, resp}{
  Integer.
  Daily hospital admissions due to respiratory and 
cardiovascular causes, 1 January 1994 and 31 December 1997. 
  }
  \item{no2, so2, rsp, o3}{
  Numeric.
  Daily mean average of \eqn{\textrm{NO}_{2}}{NO2}, 
  \eqn{\textrm{SO}_{2}}{SO2}, respirable
  suspended particles (rsp; that is
  \eqn{\textrm{PM}_{10}}{PM10}), and
  \eqn{\textrm{O}_{3}}{O3}, in parts per billion (ppb). 
  }
  \item{temp, hum}{
  Numeric.
  Daily mean average of temperature (Celsius deg.) and humidity (\%)
  }
  \item{mon, tue, wed, thur, fri, sat}{
  Factors with two levels. Weekdays/weekends indicator variables.
  
  }
    
  }
}
%\details{
%}
\source{
Data set retrieved from 
%\url{https://www.stat.nus.edu.sg/~staxyc/}
\url{https://blog.nus.edu.sg/homepage/research/}

}
\references{

Xia, Y. and Tong, H. (2006)
Cumulative effects of air pollution on public health.
\emph{Statistics in Medicine}.
\bold{25(29)}, 3548-3559.

}
\examples{
data(HKdata)
summary(HKdata[, -1])

}
\keyword{datasets}
