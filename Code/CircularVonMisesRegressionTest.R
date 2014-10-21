# rm(list=ls())

library(Rcpp)
library(circular)
library(mvtnorm)
library(ggplot2)
library(rootSolve)

sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

printcount <- FALSE
verbose <- FALSE

# Link function
linkfun    <- function(x) 4 * atan(x)
invlinkfun <- function(x) tan(x/4)

plot(linkfun, xlim=c(-10, 10))


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
n       <- 100
X       <- matrix(rnorm(n))
true_mu <- 1
true_bt <- -0.05
err     <- rnorm(n, sd = 0.1)
th      <- true_mu + linkfun(true_bt * X) + err


plot(circular(th))
hist(X, breaks=30)
plot(X, atan(th))

# Data summary statistics
R <- sqrt(sum(cos(th))^2 + sum(sin(th))^2)
mu_hat <- meanDir(th)

# Settings for sampler
mu_start <- mu_hat
kp_start <- 1
bt_start <- 2
Q <- 1000
p <- length(bt_start)


# Empty matrices and initialization.
mu <- c(mu_start, rep(NA, Q-1))
kp <- c(kp_start, rep(NA, Q-1))
bt <- rbind(bt_start, matrix(nc=p, nr=Q-1))
row.names(bt) <- NULL
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

# Bandwith for the proposal for beta
bwb <- .05

it <- 2
for(it in 2:Q) {

  if (printcount) print(paste("Start of iteration", it))

  # linkxbeta is what the set of covariates adds in terms of predicted angle.
  linkxbeta <- linkfun(apply(X, 1, "%*%", bt_cur))

  # New properties
  S_fromBeta  <- sum(sin(th - linkxbeta))
  C_fromBeta  <- sum(cos(th - linkxbeta))
  R_fromBeta  <- sqrt(C_fromBeta^2 + S_fromBeta^2)
  mu_fromBeta <- atan2(S_fromBeta, C_fromBeta)

  # Generate a new value for mu.
  ### TODO: Mathematically prove this to be equal to the ll:
  # Plots showing that this is true:
  #   plot(xlim=c(-0.5, 7),
  #        function(x) dvm(x, mu=mu_fromBeta, kappa = R_fromBeta*kp_cur))
  #   plot(xlim=c(-0.5, 7), type="l",
  #        circSeq,
  #        sapply(circSeq, function(x) exp(llfun(x, kp=kp_cur, bt=bt_cur))))

  if (is.na(mu_fromBeta) | kp_cur<1e-20) stop(paste("NA mu at it", it))
  mu_cur <- rvmc(1, mu_fromBeta, R_fromBeta*kp_cur)

  if (verbose) {
    cat("\014")
    print(paste("Start of iteration", it))
    printVars(bt_cur, kp_cur, mu_cur)
  }


  #### BETA
  #
  #   estCov_cur <- estCovMat(X = X, bt = bt_cur, kp = kp_cur,
  #                           linkfun = invlinkfun, linkxbeta = linkxbeta)
  #   bt_can     <- mvtnorm::rmvnorm(n = 1, mean = bt_cur, sigma = estCov_cur)
  #
  #   bt_lratio <- llfun(mu = mu_cur, kp = kp_cur, bt = bt_can) +
  #     dmvnorm(bt_cur, log = TRUE, mean = bt_can, sigma = estCov_cur) -
  #     llfun(mu = mu_cur, kp = kp_cur, bt = bt_cur) -
  #     dmvnorm(bt_can, log = TRUE, mean = bt_cur, sigma = estCov_cur)
  bt_can     <- bt_cur + runif(1, -bwb, bwb)

  bt_lratio <- llfun(mu = mu_cur, kp = kp_cur, bt = bt_can)  -
    llfun(mu = mu_cur, kp = kp_cur, bt = bt_cur)

  if (bt_lratio > log(runif(1))) {
    bt_cur <- bt_can
  }

  if (verbose) {
    print("After choosing new beta")
    printVars(bt_cur, bt_can, estCov_cur)
  }

  ### KAPPA
  kp_can     <- rchisq(1, df = kp_cur)

  kp_lratio <- llfun(mu = mu_cur, kp = kp_can, bt = bt_cur) +
    dchisq(kp_cur, df = kp_can, log = TRUE) -
    llfun(mu = mu_cur, kp = kp_cur, bt = bt_cur) -
    dchisq(kp_can, df = kp_cur, log = TRUE)


  if (kp_lratio > log(runif(1))) {
    kp_cur <- kp_can
  }

  if (verbose) {
    print("After choosing Kappa")
    printVars(kp_cur, kp_can)
  }




  #### PRINTING ####
  if(FALSE){


    # Likelihood
    bandwith <- .2
    lim <- bt_cur + c(-1, 1) * bandwith

    # # # # Beta
    plot(Vectorize(function(x) exp(llfun(mu_cur, kp=kp_cur, bt=x))), xlim=lim)
    # Proposal
    plot(Vectorize(function(x) dunif(x, bt_cur-bwb, bt_cur+bwb)), xlim=lim)

    #   #
    lim <- c(-0.05, 200)

    # # Kappa
    #     Likelihood
    plot(Vectorize(function(x) exp(llfun(mu_cur, kp=x, bt=bt_cur))), xlim=lim)
    plot(Vectorize(function(x) llfun(mu_cur, kp=x, bt=bt_cur)), xlim=lim)
    # Proposal
    plot(Vectorize(function(x) dchisq(x, df=kp_cur)), xlim=lim)
    plot(Vectorize(function(x) dchisq(x, df=kp_can)), xlim=lim)
  }

  mu[it]   <- mu_cur
  kp[it]   <- kp_cur
  bt[it, ] <- bt_cur

  #   lines(x = c(it-1, it), y=c(bt[it-1], bt[it]))

  #   rho <- min(1, llfun(mu_can, kp_can, bt_can) / llfun(mu[it-1], kp[it-1], bt[it-1, ]))

  # it <- it+1

}



mu_cur
kp_cur
bt_cur

plot.ts(mu)
plot.ts(mu[-(1:100)])
plot.ts(kp)
plot.ts(bt)
# plot(circular(th))

# mu_cur <- mean(circular(th))


mlebeta <- function(b) sum( (X * sin(th - mu_cur - linkfun(b * X))) / (1+linkfun(b * X)^2))

xl <- c(-3, 4)
plot(Vectorize(mlebeta), xlim=xl) ; abline(h=0)
plot(Vectorize(function(x) llfun(mu_cur, kp=kp_cur, bt=x)), xlim=xl)

uniroot.all(Vectorize(mlebeta), interval = xl)

optimize(f = function(x) llfun(mu_cur, kp=kp_cur, bt=x), interval = xl, maximum = TRUE)




