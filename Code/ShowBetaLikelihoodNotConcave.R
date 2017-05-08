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
true_mu <- pi
true_bt <- 1
true_kp <- 20
X       <- matrix(rnorm(n, sd=4))
err     <- rvmc(n, mu = 0, kp = true_kp)
th      <- (true_mu + linkfun(true_bt * X) + err) %% (2*pi)

# plot.new()
par(mfrow=c(2, 2))

# Show generated data
plotCircular(th)
plot(X, th)

# Limits
xl <- c(-12, 12)

test_mu <- pi/2
test_mu <- pi
test_mu <- true_mu


# Main part of the loglikelihood function.
Rllfun <- function(bt) sum(cos(th - test_mu - linkfun(apply(X, 1, "%*%", bt))))

llmax <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = TRUE)$maximum
llmin <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = FALSE)$minimum

mus <- c(test_mu, test_mu, test_mu, test_mu)
bts <- c(llmin, llmax, -5, 5)

# Plot likelihood of beta at true mean and true kp, and add the extrema.
plot(Vectorize(function(x) Rllfun(x)), xlim=xl, main="LogLik of Beta", xlab=expression(beta), ylab="Log-Likelihood")
abline(v=bts, col=2:(length(bts)+1), lwd=3.5)

# Plot predictions for different beta's
plot(X, th, xlim=c(min(X), max(X)), ylim = c(0, 2*pi),
     main=paste("Given mu=", round(test_mu)), legend=list())

for (i in seq(mus)) {
  estth <- function(x) (mus[i] + linkfun(bts[i]*x)) %% (2*pi)
  curve(estth, add=TRUE, col=i+1, lwd=3.5)
}
abline(h=test_mu, col="blue", lty = 2)

par(mfrow=c(1, 1))

