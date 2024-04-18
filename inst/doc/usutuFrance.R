## ----setup, include=FALSE, cache=FALSE--------------------
# set global chunk options
library('knitr')
opts_chunk$set(fig.path="usutuFrance-",
               fig.align="center",
               fig.show="hold",
               echo=TRUE,
               results="markup",
               fig.width=10,
               fig.height=10, out.width='\\linewidth',
               out.height='\\linewidth',
               cache=FALSE,
               dev='png',
               concordance=TRUE,
               error=FALSE)
opts_knit$set(aliases = c(h = 'fig.height',
              w = 'fig.width',
              wo='out.width',
              ho='out.height'))
#Sys.setlocale("LC_TIME", "C")
options(replace.assign=TRUE,width=60)
set.seed(9567)


## ----eval=FALSE-------------------------------------------
## ## If devtools is not yet installed, type
## install.packages("devtools")
## 
## ## Install the package usutuFrance
## devtools::install_github("ClementCalenge/usutuFrance", ref="main")


## ----load-package-----------------------------------------
library(usutuFrance)


## ----load-usutu-data--------------------------------------
data(usutu)
str(usutu)


## ----focus-positive---------------------------------------
usutup <- usutu[usutu$result=="Positive",]


## ----first-plot-of-the-data-------------------------------
## Required libraries
library(sf)
library(ggplot2)
library(ggspatial)

data(sfFrance)

theme_set(theme_bw())
ggplot()+
    geom_sf(data=sfFrance)+
    geom_point(aes(x=x,y=y),data=usutup,size=3)+
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed",
                                          linewidth = 0.5),
          panel.background = element_rect(fill = "aliceblue"))+
    theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none")+
    annotation_scale(location = "bl", width_hint = 0.5, plot_unit="m")


## ----estimate-K-function----------------------------------
## Distances for which the K function is desired (from 1 to 250 km)
vs <- seq(1000, 250000, length=100)

## The boundary of the study area
oo <- sfFrance[[1]][[1]]

## Coordinates of the points
xy <- as.matrix(usutup[,c("x","y")])

## Calculation of the K function (divide all coordinates and distance
## by 1000: representation in km is clearer than representation in
## meters)
kh <- splancs::khat(xy/1000, oo/1000, vs/1000)

## Calculation function L
kh <- sqrt(kh/pi)-vs/1000


## ----simulate-CSR-envelope, eval=FALSE--------------------
## ## Calculation enveloppe by simulation of complete spatial randomness
## ## (min and max for different values of distances)
## simKcsr <- splancs::Kenv.csr(nrow(xy),oo/1000,999,vs/1000,quiet=TRUE)


## ----load-result-simus------------------------------------
data(simKcsr)


## ----plot-L-function--------------------------------------
## transformation envelope of K into envelope of L
simKcsr$lower <- sqrt(simKcsr$lower/pi)-vs/1000
simKcsr$upper <- sqrt(simKcsr$upper/pi)-vs/1000

## Display the function L with envelope expected under CSR
plot(vs/1000, kh, ty="l",
     ylim=range(c(kh,unlist(simKcsr))), xlab="distance t (km)", ylab="L(t)")
polygon(cbind(c(vs/1000, rev(vs/1000)),
              c(simKcsr$lower,rev(simKcsr$upper))), col="grey", border = NA)
abline(h=0)


## ----fit-thomas-process-raw-------------------------------
## Define a point pattern (similarly, divide by 1000: it is clearer to
## work in kilometers).
ppo <- spatstat.geom::ppp(xy[,1]/1000, xy[,2]/1000,
                          window =
                              spatstat.geom::owin(poly=oo[nrow(oo):1,]/1000)) 

## Fit the Thomas process (no trend)
ppro <- spatstat.model::kppm(ppo, ~1, "Thomas")


## ----simulate-Thomas-1-envelope, eval=FALSE---------------
## simKraw <- sapply(1:999, function(r) {
##     cat(r,"\r")
##     s <- simulate(ppro,1)[[1]]
##     khs <- splancs::khat(cbind(s$x,s$y), oo/1000, vs/1000)
##     khs <- sqrt(khs/pi)-vs/1000
##     return(khs)
## })


## ----load-result-simus-1----------------------------------
data(simKraw)


## ----fit-Thomas-1-----------------------------------------
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


## ---------------------------------------------------------
## Removal of all points located at less than 5 km from another point
## Calculation of distance between points
dp <- as.matrix(dist(xy))
## identification, for eachpoint, of the index of the first point
## for which the distance is lower than 5 km (when the distance between
## i,j is lower than 5 km, if i<j, then the first point on the i-th row
## will be the point i, and the first point on the j-th row will also be
## the point i).
firstL5 <- sapply(1:nrow(dp),function(x) min(which(dp[x,]<5000)))
## Remove duplicated points
dup <- duplicated(firstL5)

## Number of pairs:
sum(dup)

## Number of pairs of the same species:
w <- which(dup)
sum(usutup$species[w]==usutup$species[w-1])

## Which corresponds to two pairs blackbirds and two pairs of Great
## grey owl
usutup$species[w][usutup$species[w]==usutup$species[w-1]]


## ----test-different-species-local-clustering--------------
set.seed(777)

ra <- sapply(1:999, function(r) {
    ## Random permutation
    w <- which(c(FALSE,sample(dup[-1])))
    ## Number of pairs of the same species
    sum(usutup$species[w]==usutup$species[w-1])
})

## The P-value of this test is
mean(ra<=4)


## ----thinning---------------------------------------------
xy_thin <- xy[!dup,]


## ----K-thin-----------------------------------------------
kh_thin <- splancs::khat(xy_thin/1000, oo/1000, vs/1000)

## Function L
kh_thin <- sqrt(kh_thin/pi)-vs/1000


## ----fit-thomas-2-----------------------------------------
## Define a point pattern for the package spatstat (divide by 1000 to work on km)
ppo_thin <- spatstat.geom::ppp(xy_thin[,1]/1000, xy_thin[,2]/1000,
                               window =
                               spatstat.geom::owin(poly=oo[nrow(oo):1,]/1000)) 

## Fit the Thomas process
ppro_thin <- spatstat.model::kppm(ppo_thin, ~1, "Thomas")


## ----simulate-Thomas-2-envelope, eval=FALSE---------------
## simKthin <- sapply(1:999, function(r) {
##     cat(r,"\r")
##     s <- simulate(ppro_thin,1)[[1]]
##     khs <- splancs::khat(cbind(s$x,s$y), oo/1000, vs/1000)
##     khs <- sqrt(khs/pi)-vs/1000
##     return(khs)
## })


## ----load-result-simus-2----------------------------------
data(simKthin)


## ----envelope-Thomas-2------------------------------------
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


## ----parameters-Thomas-Process----------------------------
## The log-density of clusters 
co <- coef(summary(ppro_thin))

## get the estimate and confidence limits, calculate the exponential
## to estimate the density, and multiply by 10000 to have the mean
## number of cases per 10000 squared km (mean number of cases in
## a square of 100 x 100 km)
exp(co[c(1,3,4)])*10000

## All the other parameters of the process
ppro_thin


## ----sd-distance-relationship-----------------------------
## Simulation of 1000 points drawn from a bivariate Gaussian
## distribution with standard deviation equal to 77 km, and
## calculation of the mean distance between pairs of points.
mean(dist(cbind(rnorm(1000,0,77),rnorm(1000,0,77))))


## ----raster-maps,fig.width=10, fig.height=5, out.width='\\linewidth',out.height='0.5\\linewidth'----
library(sp)
data(usutuEnvir)

par(mfrow = c(1,2),mar=c(0,0,2,0))
image(usutuEnvir[1], col=viridis::magma(10))
plot(sfFrance, add=TRUE)
kd <- MASS::kde2d(xy[,1], xy[,2], 125000, n=400, lims=c(range(oo[,1]),range(oo[,2])))
contour(kd$x, kd$y,kd$z, add=TRUE, lwd=2, col="lightgrey", nlevels=4,
        levels=seq(0,1.4e-11, length=10)[c(1,4,7,10)],drawlabels = FALSE)
title("Potential Wetlands")
image(usutuEnvir[2], col=viridis::magma(10))
plot(sfFrance, add=TRUE)
contour(kd$x, kd$y,kd$z, add=TRUE, lwd=2, col="lightgrey", nlevels=4,
        levels=seq(0,1.4e-11, length=10)[c(1,4,7,10)],drawlabels = FALSE)
title("Log-Human population density")


## ----mean-env-variables-----------------------------------
sppo <- sp::SpatialPoints(xy_thin, proj4string=sp::CRS(sp::proj4string(usutuEnvir)))
(obs_means <- colMeans(sp::over(sppo, usutuEnvir)))


## ----test-env-var, eval=FALSE-----------------------------
## simEnv <- t(sapply(1:999, function(r) {
##     cat(r,"\r")
##     s <- simulate(ppro_thin,1)[[1]]
##     sppob <- sp::SpatialPoints(cbind(s$x,s$y)*1000,
##                                proj4string=sp::CRS(sp::proj4string(usutuEnvir)))
##     return(colMeans(sp::over(sppob, usutuEnvir), na.rm=TRUE))
## }))


## ----load-simEnv------------------------------------------
data(simEnv)


## ----comparison-obs-sim-----------------------------------
## Mean and standard error for the index of probability of wetlands
c(mean(simEnv[,1]), sd(simEnv[,1]))
## Recall the observed mean:
obs_means[1]

## Mean and standard error for the log-human population density
c(mean(simEnv[,2]), sd(simEnv[,2]))
## Recall the observed mean:
obs_means[2]


## ----biv-test---------------------------------------------
adehabitatHS::biv.test(as.data.frame(simEnv), obs_means)


## ----time-series-plot-cases-------------------------------
library(lubridate)
seqPeriod <- seq(ymd('2018-07-15'),ymd('2018-08-31'), by = '1 day')
movAvg <- sapply(1:length(seqPeriod),  function(i) {
    x <- seqPeriod[i]
    ## days of the week centred on the day i
    d <- seq(x-3, x+3, by="1 day")
    ## number of days of the week in the study period
    ndays <- mean(d>=seqPeriod[1]&d<=seqPeriod[length(seqPeriod)])
    ## number of cases reported in this period
    ncases <- sum(usutup$date>=(seqPeriod[i]-3)&usutup$date<=(seqPeriod[i]+3))
    ## number of cases per day:
    return(ncases/ndays)
})

da <- data.frame(date=seqPeriod,
                 NumberCases=movAvg)

ggplot(da, aes(x=date,y=NumberCases))+geom_line()+xlab("Date")+
    ylab("Mean number of cases reported per day")


## ----map-spatio-temporal----------------------------------
theme_set(theme_bw())
labd <- pretty(usutup$date)
ggplot()+
    geom_sf(data=sfFrance)+
    geom_point(aes(x=x,y=y,col=as.numeric(date),group=as.factor(date)),data=usutup,size=3)+
    theme(panel.grid.major = element_line(color = gray(0.5),
                                          linetype = "dashed",
                                          linewidth = 0.5),
          panel.background = element_rect(fill = "aliceblue"))+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none")+
    viridis::scale_color_viridis(breaks = as.numeric(labd),
                                 labels = labd)+
    annotation_scale(location = "bl", width_hint = 0.5, plot_unit="m")


## ----stikhat-function-------------------------------------
## get define space-time points
po <- stpp::as.3dpoints(usutup$x/1000,usutup$y/1000,
                        as.numeric(usutup$date-min(usutup$date)))

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


## ----simus-stikhat, eval=FALSE----------------------------
## ## Randomization test: WARNING! THIS PART IS VERY SLOW
## si <- list()
## 
## ## 999 randomization (slow calculation !!)
## for (i in 1:999) {
##     cat(i,"\r")
##     pos <- po
##     pos <- stpp::as.3dpoints(usutup$x/1000, usutup$y/1000,sample(dates))
##     chs <- stpp::STIKhat(pos, s.region = oo/1000,
##                          t.region = c(min(dates-1), max(dates+1)),
##                          times=vt, dist=vs/1000, infectious=TRUE)
##     si[[i]] <- chs$Khat
## }
## 
## 
## ## Calculation of the proportion of simulations greater than observation
## KstFun <- ch$Khat ## Initialization
## for (i in 1:nrow(KstFun)) {
##     cat(i,"\r")
##     for (j in 1:ncol(KstFun)) {
##         ot <- sapply(1:length(si), function(k) si[[k]][i,j])
##         KstFun[i,j] <- mean(ot>ch$Khat[i,j])
##     }
## }


## ----load-dataset-----------------------------------------
data(KstFun)


## ----plot-stikhat-----------------------------------------
## Show the results:
dok <- data.frame(duration=rep(vt,length(vs)),
                  distance=rep(vs/1000,each=length(vt)),
                  Prob. =as.vector(t(KstFun)))

ggplot2::ggplot(dok, ggplot2::aes(x=duration, y=distance))+
    ggplot2::geom_tile(ggplot2::aes(fill=Prob.))+
    viridis::scale_fill_viridis(option="magma", limits=c(0,1))+
    ggplot2::geom_contour(ggplot2::aes(z=Prob.), breaks = 0.05, col="white",
                             alpha=0.5, linewidth=1)+ggplot2::theme_bw()+
    xlab("Duration (days)")+ylab("Distance (km)")



## ----synthesis-figure-------------------------------------
##png(filename="SynthesisFigure.png", width=1000, height=1000, pointsize=20)
library(gridBase)
library(grid)

par(mfrow=c(2, 2))
plot(vs/1000, kh, ty="l",
     ylim=range(c(kh,unlist(simKcsr))), xlab="Distance t (km)", ylab="L(t)",
     main="(A)")
polygon(cbind(c(vs/1000, rev(vs/1000)),
              c(simKcsr$lower,rev(simKcsr$upper))), col="grey", border = NA)
abline(h=0)

opar <- par(mar=c(1,1,4,1))
kd <- MASS::kde2d(xy[,1], xy[,2], 125000, n=400, lims=c(range(oo[,1]),range(oo[,2])))
plot(sfFrance, main="(B)", col=grey(0.8), border=NA)
contour(kd$x, kd$y,kd$z, add=TRUE, lwd=2, col=grey(0.6), nlevels=4,
        levels=seq(0,1.4e-11, length=10)[c(1,4,7,10)],drawlabels = FALSE)
usud <- unclass(usutup$date-min(usutup$date))
points(usutup[,c("x","y")], pch=3, col="black", cex=1.5)

da <- data.frame(date=seqPeriod,
                 NumberCases=movAvg)
par(opar)
plot(da$date, da$NumberCases, ty="l",xlab = "Date", ylab="Mean number of cases reported per day", main="(C)")
plot.new()              ## suggested by @Josh
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(1.8,1,3,1))
p <- ggplot2::ggplot(dok, ggplot2::aes(x=duration, y=distance))+
    ggplot2::geom_tile(ggplot2::aes(fill=Prob.))+
    viridis::scale_fill_viridis(option="magma", limits=c(0,1))+
    ggplot2::geom_contour(ggplot2::aes(z=Prob.), breaks = 0.05, col="white",
                          alpha=0.5, linewidth=0.8)+ggplot2::theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=13))+
    xlab("Duration (days)")+ylab("Distance (km)")+title("(D)")
print(p,vp = vp1) 
##dev.off()


## ----convert-to-spatRaster--------------------------------
## Load the package terra
library(terra)

## Convert to spatRaster
m <- kd
resx <- diff(kd$x)[1]
resy <- diff(kd$y)[1]
xmn <- min(m$x) - 0.5 * resx
xmx <- max(m$x) + 0.5 * resx
ymn <- min(m$y) - 0.5 * resy
ymx <- max(m$y) + 0.5 * resy
z <- m$z[,ncol(m$z):1]
r1 <- rast(ext(xmn, xmx, ymn, ymx), resolution=c(resx, resy))
r1[] <- as.numeric((z))

## Conversion of sfFrance to spatVector
sfv <- vect(sfFrance)

## Set the map to NA outside the limits of continental France
rr <- rasterize(sfv, r1)
r1 <- rr*r1

image(r1)


## ----location-of-points-and-routes-ACT--------------------
data(resop)
head(resop)

plot(resop[,c("x","y")], asp=1)


## ----calculate-usutu-density------------------------------
resop$usutu<-terra::extract(r1, as.matrix(resop[,c("x", "y")]))[,1]


## ----average-usutu-values-per-route-----------------------
val_agg <- resop  |>  dplyr::group_by(cd_route) |> 
dplyr::summarise(usutu_agg=mean(usutu, na.rm=T),
                 x=mean(x), y=mean(y))
head(val_agg)


## ----limits-density-usuv----------------------------------

## Median
usutu_med <- median(val_agg$usutu_agg)

## Third quartile
usutu_75 <- quantile(val_agg$usutu_agg,  probs=0.75)

## new variable in val_agg, defining the three levels of infection
val_agg$usutuF <- ifelse(val_agg$usutu_agg<=usutu_med, "low",
                  ifelse(val_agg$usutu_agg>usutu_75, "high", "medium"))
val_agg$usutuF <- factor(val_agg$usutuF, levels=c("low", "medium", "high"))

head(val_agg)


## ----load-merula------------------------------------------
data(merula)
head(merula)


## ----merula-join-usutu------------------------------------
## Join with the density map
ex2 <- terra::extract(r1, as.matrix(merula[,c("x","y")]))

## Store the results in the dataset:
merula$usutu_agg <- ex2[,1]

## median and third quartile limits:
limi <- c(usutu_med, usutu_75)

## Cut the density in three classes and store the result in merula:
cu <- cut(merula$usutu_agg, c(-10,limi,10))
cu <- factor(cu, labels=c("low","medium","high"))
merula$usutuF <- cu



## ----GAMM-fit, eval=FALSE---------------------------------
## ## Load the package mgcv
## library(mgcv)
## 
## ## Model fit with bam (large dataset)
## m_per <- bam(tot_MN ~  usutuF + s(year, k=10) +  s(year, by=usutuF, k=10) +
##                  s(cd_routeF, bs="re") + s(yearF, bs="re") +
##                  s(jd, k=6) + s(TsSR, k=6) , family=tw(),
##              data=merula, method="REML")
## 
## ## Check the model fit
## appraise(m_per)
## ## Nothing worrying there
## 
## ## Summary of the model
## summary(m_per)
## 
## ## Plot the smooth functions of explanatory variables:
## plot(m_per)
## 
## ## A data.frame for predictions, with unknown levels for the 2 random effects factors,
## ## to compute changes independently of these effects
## df_pred<-data.frame(expand.grid(usutuF=levels(merula$usutuF), cd_routeF="1",
##                                 year=1996:2022, jd=mean(merula$jd),
##                                 TsSR=mean(merula$TsSR), yearF=as.factor(2030)))
## 
## #Load the functions to compute simulations and estimate change
## simu_by_regions<-compute_distri(m_per, df_pred, time_var="year",
##                                 factor_name="usutuF", nreplicates=1000)


## ----load-simus-prediction--------------------------------
data(simu_by_regions)
str(simu_by_regions)


## ----simus-better-layout----------------------------------
#convert to a data frame from which to average values per relevant time periods
df_4analyses<- do.call("rbind", simu_by_regions) |> as.data.frame() |> 
dplyr::mutate(usutuF=rep(levels(merula$usutuF), each=1000),
       id_sim=rep(1:1000, 3))  |> 
tidyr::pivot_longer(!c(usutuF, id_sim), names_to="year",
                    values_to="abundance")  |> 
dplyr::mutate_at( 'year', as.numeric)

head(df_4analyses)


## ----calculate-trends-------------------------------------
## A data frame where abundance is averaged by simulation,
## by period (before after the 2018 Usutu episode) and
## by areas of Usutu "prevalence" (usutuF)
df_4trends<- df_4analyses  |> 
dplyr::filter(year>2014)  |> 
dplyr::mutate(usutu_BA=factor(ifelse(year<2019, "2015-2018", "2019-2022"),
                              levels=c("2015-2018", "2019-2022")))  |> 
dplyr::group_by(id_sim, usutu_BA, usutuF)  |> 
dplyr::summarise(mean_abund=mean(abundance))

## Define the order of the levels of infection
df_4trends$usutuF<-factor(df_4trends$usutuF, levels=c("low", "medium", "high"))

#Compute the trends for the 3 areas according to the previous values
df_trends<- df_4trends  |> 
dplyr::group_by(id_sim, usutu_BA, usutuF) |> 
tidyr::pivot_wider(names_from=usutu_BA, values_from=mean_abund) |> 
as.data.frame()

head(df_trends)


## ----mean-trend-with-ci-----------------------------------
## Low level of infection:
dl <- df_trends[df_trends$usutuF=="low",]
trl <- 100*(dl[,4]-dl[,3])/dl[,3]
tr_low <- c(mean(trl), quantile(trl, c(0.025,0.975)))
    
## Medium level of infection:
dm <- df_trends[df_trends$usutuF=="medium",]
trm <- 100*(dm[,4]-dm[,3])/dm[,3]
tr_medium <- c(mean(trm), quantile(trm, c(0.025,0.975)))

## High level of infection:
dh <- df_trends[df_trends$usutuF=="high",]
trh <- 100*(dh[,4]-dh[,3])/dh[,3]
tr_high <- c(mean(trh), quantile(trh, c(0.025,0.975)))

#the final data.frame which can be turned into a table for the MS
res_trends<-data.frame(UsutuArea=c("low", "medium", "high"),
                       as.data.frame(rbind(tr_low, tr_medium, tr_high)))
colnames(res_trends)<-c("Usutu area", "Mean trend", "2.5%CI", "97.5%CI")
rownames(res_trends)<-NULL

## Results:
res_trends


## ----final-summary-plot-----------------------------------

## format the data.frame df_4trends for an easier handling
df_stat_abund<-df_4trends  |> 
dplyr::group_by(usutu_BA, usutuF)  |> 
dplyr::summarise(averag_abund=mean(mean_abund),
                 CI_low=quantile(mean_abund, probs=0.025),
                 CI_high=quantile(mean_abund, probs=0.975))  |> 
dplyr::mutate_at(c('CI_low', 'CI_high'), as.numeric)


## correspondance between Usutu areas and colours
## for the following plots (scale_colour_manual)
col_low<-"palegreen4"
col_medium<-"darkorange2"
col_high<-"deeppink4"
val_colours<-c(col_low, col_medium, col_high)


## MAIN PLOT
## A plot of each simulated average values of abundance before and
## after, and the mean and 95%CI for each area
p_trends <- ggplot2::ggplot(df_4trends) +
    ggplot2::geom_line(ggplot2::aes(x=usutu_BA,  y=mean_abund,
                  group=id_sim, colour=usutuF), alpha=0.075) +
    ggplot2::geom_point(data=df_stat_abund,
                        ggplot2::aes(x=usutu_BA, y=averag_abund,
                                     colour=usutuF), size=4) +
    ggplot2::geom_errorbar(data=df_stat_abund,
                           ggplot2::aes(x=usutu_BA, ymin=CI_low, ymax=CI_high,
                                        colour=usutuF),  width=0.1, linewidth=1) +
    ggplot2::scale_colour_manual(name= "Usutu area", values=val_colours)+
    ggplot2::xlab("Period") + ggplot2::ylab("Abundance index") +
    ggplot2::facet_wrap(~usutuF) +
    ggplot2::theme( strip.text.x = ggplot2::element_blank() ) 




## Histogram of Usutu densities values for the centroid of each route
## and categorisation in 3 groups
usutu_histo <- ggplot2::ggplot(val_agg) +
    ggplot2::geom_histogram(ggplot2::aes(x=usutu_agg, fill=usutuF)) +
    ggplot2::geom_vline(xintercept=c(usutu_med, usutu_75),
                        colour="black", linewidth=1,
                        linetype=c("solid", "dashed")) +
    ggplot2::ylab("Count") + ggplot2::xlab("Density of Usutu cases") +
    ggplot2::scale_fill_manual(name="Usutu area", values=val_colours)



## Spatial plot of ACT and density of Usutu cases
library(tidyterra)
usutu_sp_act<-ggplot2::ggplot() + tidyterra::geom_spatraster(data=r1) +
    ggplot2::geom_point(data=val_agg, ggplot2::aes(x=x, y=y, colour=usutuF)) +
    ggplot2::scale_colour_manual(name="Usutu area", values=val_colours) +
    tidyterra::scale_fill_whitebox_c(name="Density of Usutu cases", palette="arid")+
    ggplot2::theme(aspect.ratio=1,axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),)



#the figure for the blackbird part
#
##png("Fig2.png", width=3000, height=3000, pointsize=3,res=300)
library(ggpubr)
ggpubr::ggarrange(ggpubr::ggarrange(usutu_histo,
                                    usutu_sp_act, ncol=2,
                                    labels=c("A", "B"), widths = c(0.6, 1)),
                  p_trends, nrow=2, labels=c("", "C"))
##dev.off()


