\name{simKcsr}
\alias{simKcsr}
\docType{data}
\title{
  Simulated Envelope Limits of the K Function for the Usutu Dataset under the Complete Spatial Randomness.
}
\description{
  This dataset contains the K function estimated for 999 simulations of
  complete spatial randomness, to test whether Usutu location data is
  randomly distributed in space (see examples
  for the code used to obtain this object).
}
\usage{data("simKcsr")}
\format{
  A matrix containing the value of the K function for each one of the
  100 values of distance and each one of the 999 simulations of the
  complete spatial randomness (see also the vignette for more context).
}
\references{
  Bouchez-Zacria, M., Calenge, C., Villers, A., Lecollinet, S.,
  Gonzalez, G., Quintard, B., Leclerc, A., Baurier, F., Paty, M.-C.,
  Eraud, C. and Decors, A. 2024. Analysis of the Usutu episode of summer
  2018 in birds in France. In prep.
}
\examples{

## Data frame containing the location of the birds.
data(usutu)

## Only the positive birds
usutup <- usutu[usutu$result=="Positive",]

## get the coordinates in a matrix
xy <- as.matrix(usutup[,c("x","y")])

## contour of the map of France
data(sfFrance)
oo <- sfFrance[[1]][[1]]

## The set of distances for which the function K is calculated
## (between 1 to 250 km)
vs <- seq(1000, 250000, length=100)

## Calculation of khat (division by 1000, since
## it is easier to interpret the function K on a
## scale expressed in kilometers instead of meters)
kh <- splancs::khat(xy/1000, oo/1000, vs/1000)

## Function L
kh <- sqrt(kh/pi)-vs/1000

## Simulate the complete spatial randomness 999 times, and recalculate
## the K function for each simulation.  WARNING: Very slow! (see below
## for results)
\dontrun{
    simKcsr <- splancs::Kenv.csr(nrow(xy),oo/1000,999,vs/1000,quiet=TRUE)
}

## You can avoid this long calculation by loading the resulting dataset
data(simKcsr)

## transformation envelope of K into envelope of L
simKcsr$lower <- sqrt(simKcsr$lower/pi)-vs/1000
simKcsr$upper <- sqrt(simKcsr$upper/pi)-vs/1000

## Display the function L with envelope expected under CSR
plot(vs/1000, kh, ty="l",
     ylim=range(c(kh,unlist(simKcsr))), xlab="distance t", ylab="L(t)")
polygon(cbind(c(vs/1000, rev(vs/1000)),
              c(simKcsr$lower,rev(simKcsr$upper))), col="grey", border = NA)
abline(h=0)


}
\keyword{datasets}
