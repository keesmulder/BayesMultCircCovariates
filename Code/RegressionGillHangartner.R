rm(list=ls())

library(Rcpp)
library(circular)
library(MASS)
sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

# Link function
linkfun <- function(x) 2 * atan(x)

# Cov-matrix
estCovMat <- function(X, bt, kp, linkfun) {
  g      <- linkfun(apply(X, 1, "%*%", bt))
  G      <- diag(g)
  Gsq    <- G %*% G
  nc     <- 1 / (kp * A1(kp))
  pt1    <- solve(t(X) %*% Gsq %*% X)
  pt2num <- pt1 %*% t(X) %*% g %*% t(g) %*% X %*% pt1
  pt2den <- n - (t(g) %*% X  %*% pt1  %*% t(X) %*% g)

  nc * (pt1 + pt2num / as.numeric(pt2den))
}


# Data generation
n <- 30
th <- rvmc(n, mu = 2, kp = 5)
X  <- matrix(2*tan(th) + rnorm(n))

# Data summary statistics
R <- sqrt(sum(cos(th))^2 + sum(sin(th))^2)
mu_hat <- meanDir(th)

# Settings for sampler
mu_start <- 2
kp_start <- 1
bt_start <- 0.1
tau_mu <- 1
tau_bt <- 1
Q <- 20
p <- length(bt_start)

# Empty matrices and initialization.
mu <- c(mu_start, rep(NA, Q-1))
kp <- c(kp_start, rep(NA, Q-1))
bt <- rbind(bt_start, matrix(nc=p, nr=Q-1))
mu_cur <- mu[1]
kp_cur <- kp[1]
bt_cur <- bt[1, ]

# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(mu, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) + kp * sum(cos(th - mu - linkfun(apply(X, 1, "%*%", bt))))
  }
}

# Log-likelihood for this dataset and link function.
llfun <- buildLogLikfun(th, X, linkfun)


for(it in 2:Q) {

  # linkxbeta is what the set of covariates adds in terms of predicted angle.
  linkxbeta <- linkfun(apply(X, 1, "%*%", bt[itm ]))

  # New properties
  S_fromBeta <- sum(sin(th - linkxbeta))
  C_fromBeta <- sum(cos(th - linkxbeta))
  R_fromBeta <- sqrt(C_fromBeta^2 + S_fromBeta^2)

  # Generate a new value for mu.
  mu_can <- rvmc(1, mu_cur, R_fromBeta*kp_cur)


  estCov_cur   <- estCovMat(X = X, bt = bt_cur, kp = kp_cur, linkfun = linkfun)
  bt_can     <- mvrnorm(n = 1, mu = bt_cur, Sigma = estCov_cur)


#   kp_can     <-

  rho <- min(1, llfun(mu_can, kp_can, bt_can) / llfun(mu[it-1], kp[it-1], bt[it-1, ]))
}

