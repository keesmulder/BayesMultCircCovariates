rm(list=ls())

source("Code/generateCircularGLMData.R")
source("Code/describeCirc.R")
source("Code/CircularVonMisesRegressionTest.R")
library(circular)


# Link functions
linkfun    <- function(x, c=2) c * atan(x)
invlinkfun <- function(x, c=2) tan(x/c)


runCircGLM <- function(d, method, ...) {
  plot(circular(d[,1]))

  attd <- attributes(d)

  truepars <- c(attd$truebeta0 %% (2*pi), attd$residkappa, attd$truebeta)

  r <- method(th=d[,1], X=d[,-1, drop=FALSE], ...)

  est <- c(beta_0=meanDir(r[, 1]), colMeans(r[, -1, drop=FALSE]))

  names(truepars) <- names(est)

  rbind(true=truepars, estimate=est)
}





d <- generateCircGLMData(truebeta0="random", truebeta=1,
                         nconpred=1, ncatpred=0, n=100,
                         linkfun=function(x) 2 * atan(x),
                         invlinkfun=function(x) tan(x/2))


runCircGLM(d, mcmcGCM,
           linkfun=function(x)2*atan(x), invlinkfun=function(x)tan(x/2),
           Q=1000)


