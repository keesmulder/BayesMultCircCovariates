
source("Code/generateCircularGLMData.R")
source("Code/describeCirc.R")
source("Code/CircularVonMisesRegressionTest.R")
sourceCpp('Code/circGLM.cpp')

library(circular)

set.seed(1337)
d  <- generateCircGLMData()
th <- d[,  1]
X  <- d[, -1]

attributes(d)

b0 <- attr(d, "truebeta0")
bt <- attr(d, "truebeta")
kp <- attr(d, "residkappa")

# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(b0, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", bt))))
  }
}
# Link functions
linkfun    <- function(x, c=2) c * atan(x)
invlinkfun <- function(x, c=2) tan(x/c)



llr <- buildLogLikfun(th, X, linkfun)


llr(b0, kp, bt)
ll(b0 = b0, kp = kp, bt = bt, X = X, th = th, c = 2)
llincl(b0 = b0, kp = kp, bt = bt, X = X, th = th, c = 2)

Q <- 1000000

a2 <- system.time(replicate(Q, rhsll(b0 = b0, kp = kp, bt = bt, X = X, th = th, c = 2)))
a3 <- system.time(replicate(Q, ll(b0 = b0, kp = kp, bt = bt, X = X, th = th, c = 2)))
a4 <- system.time(replicate(Q, llincl(b0 = b0, kp = kp, bt = bt, X = X, th = th, c = 2)))

a2
a3
a4





