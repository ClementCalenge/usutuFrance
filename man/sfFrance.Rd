\name{sfFrance}
\alias{sfFrance}
\docType{data}
\title{
  Map of France
}
\description{
  This object is an object of class \code{sfc} containing the
  coordinates of the contour of France (Lambert II projection).
}
\usage{data("sfFrance")}
\format{
  An object of class \code{sfc} (package sf)
}
\source{
  Institut Geographique National, France.
}
\examples{
data(sfFrance)

library(sf)
plot(sfFrance)
}
\keyword{datasets}
