\name{KstFun}
\alias{KstFun}
\docType{data}
\title{
  Probability That the Space-Time K Function of Dead Bird Locations with
  Identification of Usutu is Greater than Expected by Chance
}
\description{
  
  This dataset contains the probability that the Space-Time K function
  of locations of dead birds diagnosed with usutu is greater than
  expected by chance, as simulated by Monte Carlo (randomization of the
  dates of collection, see examples).

}
\usage{data("KstFun")}
\format{

  A matrix containing the probability for each one of the 100 values of
  distance and each one of the 10 weeks of the study period (see
  examples; see also the vignette for more context).

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

## get define space-time points
po <- stpp::as.3dpoints(usutup$x/1000,usutup$y/1000,
                        as.numeric(usutup$date-min(usutup$date)))

## contour of the map of France
data(sfFrance)
oo <- sfFrance[[1]][[1]]

## The set of distances for which the function K is calculated
## (between 1 to 220 km)
vs <- seq(1000, 220000, length=100)
## The weeks
vt <- 1:10

## the dates
dates <- as.numeric(usutup$date-min(usutup$date))

## Space-time function
ch <- stpp::STIKhat(po, s.region = oo/1000,
                    t.region = c(min(dates-1), max(dates+1)),
                    times=vt, dist=vs/1000, infectious=TRUE)

## Randomization test: WARNING! THIS PART IS VERY SLOW
\dontrun{
si <- list()

## 999 randomization (slow calculation !!)
for (i in 1:999) {
    cat(i,"\r")
    pos <- po
    pos <- stpp::as.3dpoints(usutup$x/1000, usutup$y/1000,sample(dates))
    chs <- stpp::STIKhat(pos, s.region = oo/1000,
                         t.region = c(min(dates-1), max(dates+1)),
                         times=vt, dist=vs/1000, infectious=TRUE)
    si[[i]] <- chs$Khat
}


## Calculation of the proportion of simulations greater than observation

KstFun <- ch$Khat ## Initialization
for (i in 1:nrow(KstFun)) {
    cat(i,"\r")
    for (j in 1:ncol(KstFun)) {
        ot <- sapply(1:length(si), function(k) si[[k]][i,j])
        KstFun[i,j] <- mean(ot>ch$Khat[i,j])
    }
}

}

## If you want to avoid this long calculation step:
data(KstFun)

## Show the results:
dok <- data.frame(duration=rep(vt,length(vs)),
                  distance=rep(vs/1000,each=length(vt)),
                  Proba=as.vector(t(KstFun)))

ggplot2::ggplot(dok, ggplot2::aes(x=duration, y=distance))+
       ggplot2::geom_tile(ggplot2::aes(fill=Proba))+
       viridis::scale_fill_viridis(option="magma", limits=c(0,1))+
       ggplot2::geom_contour(ggplot2::aes(z=Proba), breaks = 0.05, col="white",
                    alpha=0.5, linewidth=1)+ggplot2::theme_bw()

}
\keyword{datasets}
