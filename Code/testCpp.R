
source("Code/generateCircularGLMData.R")
# source("Code/describeCirc.R")
# source("Code/CircularVonMisesRegressionTest.R")
sourceCpp('Code/circGLM.cpp')

library(circular)

set.seed(1337)
d  <- generateCircGLMData()
th <- d[,  1]
X  <- d[, -1]

attributes(d)
#
b0 <- attr(d, "truebeta0")
bt <- attr(d, "truebeta")
kp <- attr(d, "residkappa")
#
# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(b0, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", bt))))
  }
}
# # Link functions
linkfun    <- function(x, c=2) c * atan(x)
invlinkfun <- function(x, c=2) tan(x/c)

# conj_prior

K <- ncol(X)

circGLM(th, X,
        conj_prior = rep(0, 3),
        bt_prior = matrix(c(0,1), nc=2, nr=K, byrow = TRUE),
        starting_values = c(0, 1, rep(0, K)),
        burnin = 0,
        lag = 1,
        bwb=.05,
        kappaModeEstBandwith=.05,
        CIsize=.95,
        Q=10000,
        c=2,
        returnPostSample=TRUE)

Rll <- buildLogLikfun(th, X, linkfun)



Rll(b0, kp, bt)

btlik <- function(bt1) computeLogBtPost(b0+pi, kp, c(bt1, 0.1, 0.1, 0.1), th, X, 2, rep(0, 4), rep(1, 4))

plot(Vectorize(btlik), xlim=c(-10,10))



