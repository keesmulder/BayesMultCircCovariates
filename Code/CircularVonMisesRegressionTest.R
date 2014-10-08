# rm(list=ls())

library(Rcpp)
library(circular)
library(mvtnorm)
sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

# Link function
linkfun <- function(x) 2 * atan(x)

# Cov-matrix
estCovMat <- function(X, bt, kp, linkfun, linkxbeta = apply(X, 1, "%*%", bt)) {
  g      <- linkfun(linkxbeta)
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
mu_start <- mu_hat
kp_start <- 100
bt_start <- 10
tau_mu <- 1
tau_bt <- 1
Q <- 1000
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
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - mu - linkfun(apply(X, 1, "%*%", bt))))
  }
}

# Log-likelihood for this dataset and link function.
llfun <- buildLogLikfun(th, X, linkfun)

it <- 2
for(it in 2:Q) {

  # linkxbeta is what the set of covariates adds in terms of predicted angle.
  linkxbeta <- linkfun(apply(X, 1, "%*%", bt_cur))

  # New properties
  S_fromBeta  <- sum(sin(th - linkxbeta))
  C_fromBeta  <- sum(cos(th - linkxbeta))
  R_fromBeta  <- sqrt(C_fromBeta^2 + S_fromBeta^2)
  mu_fromBeta <- atan2(S_fromBeta, C_fromBeta)

  # Generate a new value for mu.
  ### TODO: Mathematically prove this to be equal to the ll:
#   plot(circSeq, sapply(circSeq, function(x)
       # exp(llfun(x, kp=kp_cur, bt=bt[it-1, , drop=FALSE]))), xlim=c(-0.5, 7), type="l")
  mu_cur <- rvmc(1, mu_fromBeta, R_fromBeta*kp_cur)

#   plot(xlim=c(-0.5, 7),
#        function(x) dvm(x, mu=mu_fromBeta, kappa = R_fromBeta*kp_cur))
#
#   plot(xlim=c(-0.5, 7), type="l",
#        circSeq,
#        sapply(circSeq, function(x) exp(llfun(x, kp=kp_cur, bt=bt[it-1, , drop=FALSE]))))
#

  estCov_cur <- estCovMat(X = X, bt = bt_cur, kp = kp_cur,
                          linkfun = linkfun, linkxbeta = linkxbeta)
  bt_can     <- mvtnorm::rmvnorm(n = 1, mean = bt_cur, sigma = estCov_cur)

  bt_lratio <- llfun(mu = mu_cur, kp = kp_cur, bt = bt_can) +
                dmvnorm(bt_can, log = TRUE, mean = mu_cur, sigma = estCov_cur) -
                llfun(mu = mu_cur, kp = kp_cur, bt = bt_cur) -
                dmvnorm(bt_cur, log = TRUE, mean = mu_cur, sigma = estCov_cur)

  if (bt_lratio > log(runif(1))) {
    bt_cur <- bt_can
  }



  ### KAPPA
  kp_can     <- rchisq(1, df = kp_cur)

  kp_lratio <- llfun(mu = mu_cur, kp = kp_can, bt = bt_cur) +
                dchisq(kp_can, df = kp_cur, log = TRUE) -
                llfun(mu = mu_cur, kp = kp_cur, bt = bt_cur) -
                dchisq(kp_cur, df = kp_can, log = TRUE)


  if (kp_lratio > log(runif(1))) {
    kp_cur <- kp_can
  }


#   lim <- c(-0.05, 0.4)
#   sq  <- seq(lim[1], lim[2], length.out = 2000)
#
#   plot(xlim=lim, type="l",
#         sq,
#         sapply(sq, function(x) exp(llfun(mu_cur, kp=kp_cur, bt=x))) )
#  plot(xlim=lim, type="l",
#         sq,
#         sapply(sq, function(x) dnorm(x, mean = bt_cur, sd = estCov_cur)))



  mu[it]   <- mu_cur
  kp[it]   <- kp_cur
  bt[it, ] <- bt_cur

#   rho <- min(1, llfun(mu_can, kp_can, bt_can) / llfun(mu[it-1], kp[it-1], bt[it-1, ]))
it <- it+1

}




plot(mu, type='l')
plot(kp, type='l')
plot(bt, type='l')





