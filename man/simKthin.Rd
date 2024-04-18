\name{simKthin}
\alias{simKthin}
\docType{data}
\title{
  Simulated Envelope Limits of the K function of the Thomas Process
  fitted on the Raw Usutu Data after Thinning.
}
\description{
  This dataset contains the K function estimated for 999 simulations of
  the Thomas process fitted to the Usutu location data after removal of
  all dead birds located at less than 5 km from another bird (see
  examples for the code used to obtain this object).

}
\usage{data("simKthin")}
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

## Removal of all points located at less than 5 km from another point
## Calculation of distance between points
dp <- dist(xy)  |>  as.matrix()
## identification, for eachpoint, of the index of the first point
## for which the distance is lower than 5 km (when the distance between
## i,j is lower than 5 km, if i<j, then the first point on the i-th row
## will be the point i, and the first point on the j-th row will also be
## the point i).
firstL5 <- sapply(1:nrow(dp),function(x) min(which(dp[x,]<5000)))
## Remove duplicated points
dup <- duplicated(firstL5)
xy_thin <- xy[!dup,]

## contour of the map of France
data(sfFrance)
oo <- sfFrance[[1]][[1]]

## The set of distances for which the function K is calculated
## (between 1 to 250 km)
vs <- seq(1000, 250000, length=100)

## Calculation of khat (division by 1000, since
## it is easier to interpret the function K on a
## scale expressed in kilometers instead of meters)
kh_thin <- splancs::khat(xy_thin/1000, oo/1000, vs/1000)

## Function L
kh_thin <- sqrt(kh_thin/pi)-vs/1000


## Fit the Thomas process (similarly, divide by 1000
## to work in kilometers instead of meters)
## Define a point pattern for the package spatstat
ppo_thin <- spatstat.geom::ppp(xy_thin[,1]/1000, xy_thin[,2]/1000,
                               window =
                               spatstat.geom::owin(poly=oo[nrow(oo):1,]/1000)) 

## Fit the Thomas process
ppro_thin <- spatstat.model::kppm(ppo_thin, ~1, "Thomas")

## Simulate this Thomas process 999 times,
## and recalculate the K function for each simulation.
## WARNING: Very slow!
\dontrun{
simKthin <- sapply(1:999, function(r) {
    cat(r,"\r")
    s <- simulate(ppro_thin,1)[[1]]
    khs <- splancs::khat(cbind(s$x,s$y), oo/1000, vs/1000)
    khs <- sqrt(khs/pi)-vs/1000
    return(khs)
})
}

## You can get the result with
data(simKthin)

## And show the fit with:
rk <- t(apply(simKthin,1,function(x) range(x)))
vsb <- vs/1000

## prepare the plot
plot(vsb,vsb, ty="n", ylim=range(c(c(rk),kh_thin)),
      ylab="Function L", xlab="Distance")

## show the envelope
polygon(c(vsb, rev(vsb)), c(rk[,1], rev(rk[,2])),
        col="grey", border=NA)
## and the function L
lines(vsb, kh_thin)

}
\keyword{datasets}
