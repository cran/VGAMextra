\name{ap.mx}
\alias{ap.mx}
\docType{data}
\title{Air pollution Data, Mexico City.}
\description{
  Daily air pollution levels in Mexico City, January 2004 -- June 2005.

}
\usage{data(ap.mx)}
\format{
This data frame stores time series vectors with the
following information:
  \describe{
  \item{time}{
  Time vector.
  }
  \item{PM10}{
  24--hr average concentration of
  \eqn{PM_{10}}{PM10}, in micrograms per milliliter.
  }
  \item{O3}{
  Daily maximum 8--hour moving average
  of ozone, in micrograms per milliliter.
  }
  \item{temp}{
  Daily mean average of temperature, in celsius degrees.
  }
  \item{HR}{
  Daily mean average (\%) of relative humidity.
  }
    
  }
}
\details{
These are readings of \eqn{PM_{10}}{PM10},
\eqn{O_{3}}{O3}, temperature and humidity between
1 January 2004 and 30 June 2005 in 
Mexico City Metropolitan Area.
Each observation is the 24--hr mean average (between
00:00 and 23:59 hrs), except for
ozone, where the
maximum over all the sliding 8--hour--windows,
between 00:00 and 23:59 hrs is reported, viz.
the daily maximum 8--hour moving average.



}
\source{
National Institute of Ecology.
Gathers and disseminates the data generated 
by the central air quality monitoring network
in Mexico City. Website:
  https://www.gob.mx/inecc/

}
%\references{
%}
\examples{
data(ap.mx)
summary(ap.mx[, -1])
class(ap.mx[, "PM10"])
\donttest{
layout(matrix(c(1, 1, 2,3), 2, 2, byrow = TRUE))
plot.ts(ts(ap.mx$PM10), ylab = expression(PM[10]~"Series"), 
        col = "brown", xaxt = "n", las = 1)
xtick <- c(1, 92, 183, 275, 367, 457, 518)
xtext <- c("Jan/04", "April/04", "July/04", "Oct/04", "Jan/05",
           "April/05", "June/05")
axis(side = 1, at = xtick, labels = FALSE)
text(x = xtick, par("usr")[3], labels = xtext,
     pos = 1, xpd = TRUE, col = "black")
pacf(ap.mx$PM10, main = "", ylim= c(-0.5, 1), lag.max = 60, las = 1)
acf(ap.mx$PM10, main = "", ylim= c(-0.5, 1), lag.max = 60, las = 1)
}
}
\keyword{datasets}
