\name{simEnv}
\alias{simEnv}
\docType{data}
\title{
  Monte Carlo Test of the Effect of Humid Areas and Human Population
  Density on the Occurrence of Dead Birds with Positive Usutu Test
}
\description{
  This dataset contains the mean of Humid areas and Human population
  density estimated for 999 simulations of 
  the Thomas process fitted to the Usutu location data after removal of
  all dead birds located at less than 5 km from another bird (see
  examples for the code used to obtain this object).
}
\usage{data("simEnv")}
\format{
  A matrix containing the value of the mean of humid area (first column)
  and human population density (second column) for each one of the 999
  simulations of the Thomas process (see also the vignette for more
  context). 
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

## Map of humid areas and human density
data(usutuEnvir)

## Fit the Thomas process (similarly, divide by 1000
## to work in kilometers instead of meters)
## Define a point pattern for the package spatstat
ppo_thin <- spatstat.geom::ppp(xy_thin[,1]/1000, xy_thin[,2]/1000,
                               window =
                               spatstat.geom::owin(poly=oo[nrow(oo):1,]/1000)) 

## Fit the Thomas process
ppro_thin <- spatstat.model::kppm(ppo_thin, ~1, "Thomas")

## observed mean values of humid areas and human population density 
sppo <- sp::SpatialPoints(xy_thin, proj4string=sp::CRS(sp::proj4string(usutuEnvir)))
obs_means <- colMeans(sp::over(sppo, usutuEnvir))

## Simulate the Thomas Process 999 times and calculate the mean values
## of humid areas and human density for each simulation
\dontrun{
simEnv <- t(sapply(1:999, function(r) {
    cat(r,"\r")
    s <- simulate(ppro_thin,1)[[1]]
    sppob <- sp::SpatialPoints(cbind(s$x,s$y)*1000,
                               proj4string=sp::CRS(sp::proj4string(usutuEnvir)))
    return(colMeans(sp::over(sppob, usutuEnvir), na.rm=TRUE))
}))
}

## To avoid this long calculation:
data(simEnv)

## Bivariate test
adehabitatHS::biv.test(as.data.frame(simEnv), obs_means)



}
\keyword{datasets}
