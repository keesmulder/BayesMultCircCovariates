# Re-find Gill & Hangartner's simulation
# rm(list=ls())

library(Rcpp)
library(circular)
library(mvtnorm)
library(ggplot2)
library(rootSolve)

sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

par(cex=14)

# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

# Data generation
n       <- 1000
true_b0 <- pi
true_bt <- 1
true_kp <- 20
X       <- matrix(rnorm(n, sd=.5))
err     <- rvmc(n, mu = 0, kp = true_kp)
th      <- (true_b0 + linkfun(true_bt * X) + err) %% (2*pi)

# plot.new()
par(mfrow=c(2, 2))

# Show generated data
plotCircular(th)
plot(X, th)

# Limits
xl <- c(-12, 12)
xl <- c(-100, 100)
xlC <- c(-pi, pi)

b0 <- pi/2

# Main part of the loglikelihood function.
Rllfun <- function(bt) sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", bt))))
RllfunZeta <- function(zeta) sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", tan(zeta/2)))))

# Plot likelihood of beta at true mean and true kp, and add the extrema.
plot(Vectorize(function(x) Rllfun(x)), xlim=xl, main="LogLik of Beta", xlab=expression(beta), ylab="Log-Likelihood")
plot(Vectorize(function(x) RllfunZeta(x)), xlim=xlC, main="LogLik of Zeta", xlab=expression(zeta), ylab="Log-Likelihood")
abline(v=bts, col=2:(length(bts)+1), lwd=3.5)

# Plot predictions for different beta's
plot(X, th, xlim=c(min(X), max(X)), ylim = c(0, 2*pi),
     main=paste("Given mu=", round(b0)), legend=list())

for (i in seq(mus)) {
  estth <- function(x) (mus[i] + linkfun(bts[i]*x)) %% (2*pi)
  curve(estth, add=TRUE, col=i+1, lwd=3.5)
}

abline(h=b0, col="blue", lty = 2)

par(mfrow=c(1, 1))


logit <- function(x)   x / (1-x)
logistic <- function(x)   exp(x) / (1+exp(x))
plot(logistic, xlim=xl)


