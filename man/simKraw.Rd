\name{simKraw}
\alias{simKraw}
\docType{data}
\title{
  Simulated Envelope Limits of the K function of the Thomas Process
  fitted on the Raw Usutu Data.
}
\description{
  This dataset contains the K function estimated for 999 simulations of
  the Thomas process fitted to the raw Usutu location data (see examples
  for the code used to obtain this object).
}
\usage{data("simKraw")}
\format{
  A matrix containing the value of the K function for each one of the
  100 values of distance and each one of the 999 simulations of the
  Thomas process (see also the vignette for more context).
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


## Fit the Thomas process (similarly, divide by 1000
## to work in kilometers instead of meters)
## Define a point pattern for the package spatstat
ppo <- spatstat.geom::ppp(xy[,1]/1000, xy[,2]/1000,
                          window =
                          spatstat.geom::owin(poly=oo[nrow(oo):1,]/1000)) 

## Fit the Thomas process
ppro <- spatstat.model::kppm(ppo, ~1, "Thomas")

## Simulate the Thomas process 999 times,
## and recalculate the K function for each simulation.
## WARNING: Very slow!
\dontrun{
simKraw <- sapply(1:999, function(r) {
    cat(r,"\r")
    s <- simulate(ppro,1)[[1]]
    khs <- splancs::khat(cbind(s$x,s$y), oo/1000, vs/1000)
    khs <- sqrt(khs/pi)-vs/1000
    return(khs)
})
}

## You can get the result with
data(simKraw)

## And show the fit with:
rk <- t(apply(simKraw,1,function(x) range(x)))
vsb <- vs/1000

## prepare the plot
plot(vsb,vsb, ty="n", ylim=range(c(c(rk),kh)),
      ylab="Function L", xlab="Distance")

## show the envelope
polygon(c(vsb, rev(vsb)), c(rk[,1], rev(rk[,2])),
        col="grey", border=NA)
## and the function L
lines(vsb, kh)


}
\keyword{datasets}
