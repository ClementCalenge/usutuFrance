\name{resop}
\alias{resop}
\docType{data}
\title{
  Spatial Distribution of the Sampling Units for the ACT Survey
}
\description{
  This dataset contains the spatial distribution of the sampling routes
  and points used to count singing male blackbirds (see
  \code{data(merula)}) in continental France between 1996 and 2022. 
}
\usage{data("resop")}
\format{
  A data frame with 5217 observations on the following 8 variables.
  \describe{
    \item{\code{id_point}}{a numeric vector containing the ID of the
      sampling point}
    \item{\code{cd_route}}{a character vector containing the ID of the
      corresponding sampling route}
    \item{\code{no_point}}{a numeric vector containing the number of the
      point on the route (there are 5 points per route.}
    \item{\code{cd_route_point}}{a character vector concatenating the ID
      of the route and of the point}
    \item{\code{insee_dep}}{a character vector containing the INSEE code
      of the French department where the route is located.}
    \item{\code{x}}{numeric vectors containing the X Lambert II
      coordinates of the point}
    \item{\code{y}}{numeric vectors containing the Y Lambert II
      coordinates of the point}
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
  12 -- 13, 1 -- 2.
}
\examples{
data(resop)
plot(resop[,c("x","y")], asp=1)
}
\keyword{datasets}
