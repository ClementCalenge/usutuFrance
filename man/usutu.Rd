\name{usutu}
\alias{usutu}
\docType{data}
\title{
Birds Tested for Usutu from the SAGIR Network in France
}
\description{
  This data frame contains the result of the Usutu tests carried out on
  dead birds collected by the SAGIR network in July-August 2018 in France.
}
\usage{data("usutu")}
\format{
  A data frame with 60 observations on the following 5 variables.
  \describe{
    \item{\code{species}}{a character string containing the latin name
      of the bird species}
    \item{\code{date}}{the Date of collection of the bird}
    \item{\code{x}}{the x coordinate of the dead bird at the time of
      collection (Lambert II projection)}
    \item{\code{y}}{the y coordinate of the dead bird at the time of
      collection (Lambert II projection)}
    \item{\code{result}}{the result of the test, a character vector
      containing either \code{"Positive"} or \code{"Negative"}.}
  }
}
\source{
  SAGIR network, Office francais de la biodiversite, Orleans, France.
}
\references{
  Bouchez-Zacria, M., Calenge, C., Villers, A., Lecollinet, S.,
  Gonzalez, G., Quintard, B., Leclerc, A., Baurier, F., Paty, M.-C.,
  Eraud, C. and Decors, A. 2024. Analysis of the Usutu episode of summer
  2018 in birds in France. In prep.
}
\examples{
data(usutu)

head(usutu)
}
\keyword{datasets}
