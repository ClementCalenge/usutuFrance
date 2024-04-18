\name{merula}
\alias{merula}
\docType{data}
\title{
  Results of the ACT Survey Carried Out in France From 1996 to 2022 for the
  Monitoring of the Common Blackbird 
}
\description{
  This dataset contains the results of the count of singing male
  blackbird on 1044 sampling routes distributed over continental France,
  every year from 1996 to 2022.
}
\usage{data("merula")}
\format{
  A data frame with 49266 observations on the following 12 variables.
  \describe{
    \item{\code{tot_MN}}{a numeric vector containing the total number of
      singing male blackbirds counted on the road during the day of
      counting.}
    \item{\code{x}}{the Lambert II X coordinate of the centroid of the sampling
      route}
    \item{\code{y}}{the Lambert II Y coordinate of the centroid of the sampling
      route}
    \item{\code{dateTime}}{a POSIXct containing the date and time of the
    count}
    \item{\code{year}}{a numeric vector containing the year of the count}
    \item{\code{TsSR}}{a numeric vector containing the number of minutes
      since sunrise}
    \item{\code{cd_routeF}}{a factor with 1044 levels containing the
      route ID}
    \item{\code{yearF}}{a factor with 27 levels corresponding to the 27
      years of monitoring}
    \item{\code{jd}}{a numeric vector containing the Julian date of the
      count (number of days elapsed since January 1st. }
    \item{\code{passageN}}{a factor with levels \code{First}
      \code{Second} containing the passage number on the route.}
  }
}
\source{
  Reseau Oiseaux de Passage, Office Francais de la Biodiversite/
  Federation Nationale des Chasseurs / Federation Departementale des
  Chasseurs.
}
\references{
  Boutin, J.-M., Barbier, L. and Roux, D. (2001) Suivi des effectifs
  nicheurs d'alaudidees, colombides et turdidees en France: le
  programme ACT. Alauda, 69, 53 -- 61.
  
  Boutin, J.-M., Roux, D. and Eraud, C. (2003) Breeding bird
  monitoring in France: the ACT Survey. Ornis Hungarica,
  12-13, 1-2.
}
\examples{
data(merula)
head(merula)
}
\keyword{datasets}
