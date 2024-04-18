\name{compute_distri}
\alias{compute_distri}
\title{
  Simulate Predictions of the Generalized Additive Model Used to Predict
  the Number of Blackbirds on Sampling Routes in the ACT Survey
}
\description{
  This functions simulates the Generalized Additive Model
  fitted to the bird counts carried out between 1996 and 2022 in
  continental France as a function of year, julian date, time since
  sunrise, sampling route, and the level of infection by the Usutu
  virus.
}
\usage{
compute_distri(model, dfpred, time_var = "year", factor_name = NA, nreplicates = 1000, ...)
}
\arguments{
  \item{model}{
    the fitted model.
  }
  \item{dfpred}{
    a data frame in which to look for variables with
    which to predict.
  }
  \item{time_var}{
    a character string indicating the name of the year variable in
    \code{dfpred}.
}
  \item{factor_name}{
    optionaly, a character string indicating the name of a factor
    variable used to split the predictions (typically, the factor
    indicating the level of infection by the usutu virus).
  }
  \item{nreplicates}{
    the number of simulations.
  }
  \item{\dots}{
    additional arguments to be passed to and from other functions.
  }
}
\details{

  \code{nreplicate} coefficient vectors are simulated from a
  multivariate Gaussian distribution from the results of the GAM (parameterized
  by the mean coefficient vector and the covariance matrix of the
  coefficient, extracted from the model), and each simulated coefficient
  vector is used to predict the response variable (total number of male
  blackbirds counted on an average route). This prediction can be
  carried out for each level of \code{factor_name}. 

}
\value{
  A list with L components (one for each level of \code{factor_name}), each component
  being a data.frame with 1000 rows (the 1000 simulations) and Y
  columns (the Y years, in the variable \code{time_var}).
}
\references{

  Bouchez-Zacria, M., Calenge, C., Villers, A., Lecollinet, S.,
  Gonzalez, G., Quintard, B., Leclerc, A., Baurier, F., Paty, M.-C.,
  Eraud, C. and Decors, A. 2024. Analysis of the Usutu episode of summer
  2018 in birds in France. In prep.
}
\author{
  Alexandre Villers \email{alexandre.villers@ofb.gouv.fr}
}
\seealso{
  \code{\link{predict.gam}}
}
\examples{
## No examples. See
## vignette("usutuFrance")
## for an example of application
}
\keyword{ models }

