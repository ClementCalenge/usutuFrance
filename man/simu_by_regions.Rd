\name{simu_by_regions}
\alias{simu_by_regions}
\docType{data}
\title{
  Results of the Simulations of the Generalized Additive Model to
  Predict Blackbird Population Trends in Three Regions
}
\description{
  This dataset contains simulations of the Generalized Additive Model
  fitted to the bird counts carried out between 1996 and 2022 in
  continental France as a function of year, julian date, time since
  sunrise, sampling route, and the level of infection by the Usutu
  virus.
}
\usage{data("simu_by_regions")}
\format{
  The format is:
List of 3
 $ low   : num [1:1000, 1:27] 6.15 5.74 6.23 6.2 5.93 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:27] "1996" "1997" "1998" "1999" ...
 $ medium: num [1:1000, 1:27] 7.14 6.7 7.09 7.5 7.47 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:27] "1996" "1997" "1998" "1999" ...
 $ high  : num [1:1000, 1:27] 7.02 6.58 6.85 6.58 6.87 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:27] "1996" "1997" "1998" "1999" ...
}
\details{

  One thousand of coefficient vectors were simulated from a multivariate
  distribution from the results of the GAM (parameterized by the mean
  coefficient vector and the covariance matrix of the coefficient), and
  each simulated coefficient vector was used to predict the response
  variable (total number of male blackbirds counted on an average
  route). We carried out this prediction in three French regions defined
  by the level of infection by the usutu virus (low, medium and
  high). The result is a list with three components (one for each level
  of infection), each component being a data.frame with 1000 rows (the
  1000 simulations) and 27 columns (the 27 years).

  
}
\references{

  Bouchez-Zacria, M., Calenge, C., Villers, A., Lecollinet, S.,
  Gonzalez, G., Quintard, B., Leclerc, A., Baurier, F., Paty, M.-C.,
  Eraud, C. and Decors, A. 2024. Analysis of the Usutu episode of summer
  2018 in birds in France. In prep.
  
}
\examples{
data(simu_by_regions)

## Mean values:
mv <- lapply(simu_by_regions, colMeans)

plot(1996:2022, mv[[1]], ylim=range(unlist(mv)),
     xlab="Year", ylab="Mean count", ty="n")
points(1996:2022, mv[[1]], col="green", ty="b")
points(1996:2022, mv[[2]], col="orange", ty="b")
points(1996:2022, mv[[3]], col="red", ty="b")


}
\keyword{datasets}
