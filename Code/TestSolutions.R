rm(list=ls())

source("Code/generateCircularGLMData.R")
source("Code/describeCirc.R")
source("Code/CircularVonMisesRegressionTest.R")
library(circular)


# Link functions
linkfun    <- function(x, c=2) c * atan(x)
invlinkfun <- function(x, c=2) tan(x/c)

# method: The method of analysis
# dpar: Arguments to be passed to the data generation method.
# mpar: Arguments to be passed to the data analysis method.
runCircGLM <- function(method=mcmcGCM, dpar, mpar) {

  # Generate data
  d <- do.call(generateCircGLMData, dpar)

  # Save the true parameters
  attd <- attributes(d)
  truepars <- c(attd$truebeta0 %% (2*pi), attd$residkappa, attd$truebeta)

  # Analyze the dataset
  r <- do.call(method, args=c(list(th=d[,1], X=d[,-1, drop=FALSE]), mpar))

  # Obtain estimates
  est <- c(beta_0=meanDir(r[, 1]), colMeans(r[, -1, drop=FALSE]))

  # Return results.
  rbind(true=truepars, estimate=est)
}

dpar <- list(truebeta0="random", truebeta=1, nconpred=1, ncatpred=0, n=100,
             linkfun=function(x) 2 * atan(x), invlinkfun=function(x) tan(x/2))

mpar <- list(linkfun=function(x)2*atan(x), invlinkfun=function(x)tan(x/2),
             Q=1000)

d <- generateCircGLMData()
th <- d[,1]
X <- d[,-1]


runCircGLM(mcmcGCM, dpar, mpar)


