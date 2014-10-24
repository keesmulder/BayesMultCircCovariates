rm(list=ls())

library(Rcpp)
library(circular)
library(mvtnorm)
library(ggplot2)
library(rootSolve)

sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")

# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(mu, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - mu - linkfun(apply(X, 1, "%*%", bt))))
  }
}

# bwb: Bandwith for the proposal for beta
mcmcGCM <- function(th, X, linkfun, invlinkfun,
                    mu_start=0, kp_start=1, bt_start=2,
                    bwb=.05, Q=1000,
                    printcount=FALSE, verbose=FALSE) {

  # Empty matrices and initialization.
  p <- length(bt_start)
  mu <- c(mu_start, rep(NA, Q-1))
  kp <- c(kp_start, rep(NA, Q-1))
  bt <- rbind(bt_start, matrix(nc=p, nr=Q-1))
  rownames(bt) <- NULL
  colnames(bt) <- paste0("beta_", 1:p)
  mu_cur <- mu[1]
  kp_cur <- kp[1]
  bt_cur <- bt[1, ]

  # Log-likelihood for this dataset and link function.
  llfun <- buildLogLikfun(th, X, linkfun)

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
      plot(Vectorize(function(x) exp(llfun(0.1, kp=x, bt=bt_cur))), xlim=lim)
      plot(Vectorize(function(x) exp(llfun(0.1, kp=x, bt=0))), xlim=lim)
      plot(Vectorize(function(x) exp(llfun(0.4, kp=x, bt=bt_cur))), xlim=lim)
      plot(Vectorize(function(x) llfun(mu_cur, kp=x, bt=bt_cur)), xlim=lim)
      # Proposal
      plot(Vectorize(function(x) dchisq(x, df=kp_cur)), xlim=lim)
      plot(Vectorize(function(x) dchisq(x, df=kp_can)), xlim=lim)
    }

    mu[it]   <- mu_cur
    kp[it]   <- kp_cur
    bt[it, ] <- bt_cur


  }

  sapply(ls(), function(obj) assign(obj, get(obj), envir = .GlobalEnv))
  return(data.frame(mu, kp, bt))
}



# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

# Data generation
n       <- 400

# Re-find Gill & Hangartner's simulation
true_mu <- 0
true_bt <- 2
true_kp <- 2

X       <- matrix(runif(n, -1, 1))
err     <- rvmc(n, mu = 0, kp = true_kp)

th      <- (true_mu + linkfun(true_bt * X) + err) %% (2*pi)

# Data summary statistics
R      <- sqrt(sum(cos(th))^2 + sum(sin(th))^2)
R_bar  <- R/n
mu_hat <- meanDir(th)

plot.new()
par(mfrow=c(3, 3))


# Show generated data
plot(circular(th))
# plot(linkfun, xlim=c(-10, 10))
hist(X, breaks=30)
plot(X, th)

# Run sampler
 s <- mcmcGCM(th = th, X = X, linkfun = linkfun, invlinkfun = invlinkfun)

# Here, we can see that the estimate of kappa coincides with the residual variance
# plot(circular(th))
# plot(circular(th - linkfun(true_bt * X)))

# Approx kappa of the full sample
approxKappaML(th)
A1inv(R_bar)

# Approx residual kappa
approxKappaML(circular(th - linkfun(true_bt * X)))

# Plot sample results
plot.ts(s$mu)
plot.ts(s$kp)
plot.ts(s$beta_1)


# MLE BETA EXPLORATION
mlebeta <- function(b) sum( (X * sin(th - mu_cur - linkfun(b * X))) / (1 + linkfun(b * X)^2))

xl <- c(-150, 150)
xl <- c(-15, 15)

mu_cur <- true_mu
kp_cur <- true_kp
# mu_cur <- pi/2
# kp_cur <- 300


Rllfun <- function(bt) sum(cos(th - mu_cur - linkfun(apply(X, 1, "%*%", bt))))

mleroot <- uniroot.all(Vectorize(mlebeta), interval = xl)


llmax <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = TRUE)$maximum
llmin <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = FALSE)$minimum

mleroot
c(llmin, llmax)


plot(Vectorize(mlebeta), xlim=xl, main="MLEBETA function")
abline(h=0, col="green")
abline(v=mleroot, col="blue")
abline(v=c(llmin, llmax), col="red")
plot(Vectorize(function(x) Rllfun(x)), xlim=xl, main="LogLik of Beta")
abline(v=mleroot, col="blue")
abline(v=c(llmin, llmax), col="red")


par(mfrow=c(1, 1))


mus <- c(true_mu, true_mu, true_mu)
bts <- c(llmin, llmax, -1500)

plot(X, th)
i <- 3
for (i in seq(mus)) {
  estth <- function(x) (mus[i] + linkfun(bts[i]*x)) %% (2*pi)
  curve(estth, xlim=c(-1, 1), ylim = c(-pi, pi), add=TRUE, col=i)
}






