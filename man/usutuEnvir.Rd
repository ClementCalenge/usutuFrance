\name{usutuEnvir}
\alias{usutuEnvir}
\docType{data}
\title{
  Map of Potential Wetlands and of Log-Human Population Density in France
}
\description{
  This \code{SpatialGridDataFrame} contains the raster maps of two
  environmental variables: the map of potential wetlands and the
  log-human population density (the coordinates are presented in Lambert
  II projection).
}
\usage{data("usutuEnvir")}
\format{
  An object of class \code{SpatialGridDataFrame} (package "sp")
}
\details{

  The map of potential wetlands was derived by the INRA (see source): it
  estimates the probability that each 1km x 1km pixel of the map is a
  wetland and then discretize this probability in three classes (0 = no
  wetlands, 1 = rather strong probability, 2 = strong probability, or 3
  = very strong probability). This map also identifies lakes and
  foreshores in France. We downloaded this map, and redefined pixels
  with lakes and foreshores as "certain wetlands" (new class 4). Then,
  we smoothed this map by calculating for each pixel the sum of this
  ordered factor (0 to 4) for the focus pixel as well as the pixels
  immediately on the left, the right, the top and the bottom of the
  focus pixel.

  The map of human density was calculated by finding the commune in
  which the center of each 1 x 1 km pixel was located and attributing
  the 2005 human density to it (number of inhabitants divided by commune
  area). We then log-transformed this map (transformation: x ->
  log(x+1)).
  
}
\source{
  For the human density, Institut Geographique National (France). The
  density of Humid areas corresponds to the map of potential wetlands
  derived by INRA (Orleans) at the following URL
  \url{http://geowww.agrocampus-ouest.fr/web/?p=1538}. 
}
\examples{
data(usutuEnvir)

## Wetlands:
image(usutuEnvir[1])

## Log-Human population density
image(usutuEnvir[2])


}
\keyword{datasets}
