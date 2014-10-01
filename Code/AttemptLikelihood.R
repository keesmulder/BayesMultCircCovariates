library(Rcpp)

rm(list=ls())

sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

n <- 10
mu <- 2
kp <- 3

linkfun <- function(x) 2 * atan(x)
plot(linkfun, xlim=c(-10, 10))

th <- rvmc(n, mu, kp)
X  <- matrix(2*tan(th) + rnorm(n))


plot(X, th)
plot(linkfun(X), th)

# Joint Log-likelihood function, with data already built in.
buildLikfunGivenData <- function(th, X, linkfun){
  function(mu, kp, bt){
    n <- length(th)
    - n * log (besselI(kp, 0)) + kp * sum(cos(th - mu - linkfun(apply(X, 1, "%*%", bt))))
  }
}

Q <- 100

lfun <- buildLikfunGivenData(th, X, linkfun)

sq <- seq(-1, 8, 0.01)

R <- sqrt(sum(cos(th))^2 + sum(sin(th))^2)
mu_n <- meanDir(th)


bt <- 0.1
kp <- 3
# These plots are the same.
plot(sq, sapply(sq, function(x) exp(lfun(x, kp, bt))), type="l", xlim=c(0, 2*pi))
plot(sq, sapply(sq, function(x) dvm(x, mu = meanDir(th) + linkfun(bt %*% x)), type="l", xlim=c(0, 2*pi))
plot(sq, sapply(sq, function(x) besselI(kp, 0)^(-n) * exp(R * kp * cos(x + linkfun(bt %*% x)- mu_n))), type="l", xlim=c(0, 2*pi))

abline(v=c(0, 2*pi), col="red")


