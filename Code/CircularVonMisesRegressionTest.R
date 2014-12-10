# rm(list=ls())

library(Rcpp)
library(circular)
library(mvtnorm)
library(ggplot2)

sourceCpp("Code/rvmc.cpp")

# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(b0, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", bt))))
  }
}



mlebeta <- function(X, k=1, xl = c(-20, 20)) {
  require(rootSolve)

  # Amount of parameters.
  K <- ncol(X)

  ddbtloglik <- function(b) {
    # Function of the derivative of the log likelihood with respect to beta.
    # Provides MLE's at 0, for beta_k.
    sum( (X[, k] * sin(th - th_bar - linkfun(b * X[, k]))) / (1+(linkfun(b * X[, k])^2)))
  }

  # Log-likelihood for this dataset and link function.
  llfun <- buildLogLikfun(th, X, linkfun)

  # The optima found in the log-likelihood of beta of the current predictor.
  lloptima <- uniroot.all(Vectorize(ddbtloglik), interval = xl)

  # Obtain the likelihood if we set beta_k to a specific value, and the rest to
  # zero.
  getllOneBt <- function(x) {
    bts    <- integer(K)
    bts[k] <- x
    llfun(th_bar, 1, bts)
  }

  optll    <- sapply(lloptima, getllOneBt)
  whichmle <- which.max(optll)

  lloptima[whichmle]
}


# bwb: Bandwith for the proposal for beta
circGLMR <- function(th, X, linkfun, invlinkfun,
                    b0_start=0, kp_start=1, bt_start=rep(0, ncol(X)),
                    bwb=rep(.5, ncol(X)), Q=1000) {

  # Log-likelihood for this dataset and link function.
  llfun <- buildLogLikfun(th, X, linkfun)

  # Empty matrices and initialization.
  K  <- length(bt_start)
  b0 <- c(b0_start, rep(NA, Q-1))
  kp <- c(kp_start, rep(NA, Q-1))
  bt <- rbind(bt_start, matrix(nc=K, nr=Q-1))
  rownames(bt) <- NULL
  colnames(bt) <- paste0("beta_", 1:K)
  b0_cur <- b0[1]
  kp_cur <- kp[1]
  bt_cur <- bt[1, ]


  for(it in 2:Q) {

    # linkxbeta is what the set of covariates adds in terms of predicted angle.
    linkxbeta <- linkfun(apply(X, 1, "%*%", bt_cur))

    # New properties
    psi    <- th - linkxbeta
    S_psi  <- sum(sin(psi))
    C_psi  <- sum(cos(psi))
    R_psi  <- sqrt(C_psi^2 + S_psi^2)
    b0_psi <- atan2(S_psi, C_psi)

    # Generate a new value for b0.
    if (is.na(b0_psi) | kp_cur<1e-20) stop(paste("NA b0 at it", it))
    b0_cur <- rvmc(1, b0_psi, R_psi*kp_cur)
#     plot(function(x) dvonmises(x, b0_psi, R_psi*kp_cur), xlim=c(0, 6))

    #### BETA
    for (k in 1:K) {

      # Proposal for the beta we are currently working on.
      btk_can     <- bt_cur[k] + runif(1, -bwb[k], bwb[k])

      # Resulting vector of beta's including our new vector.
      bt_can      <- bt_cur
      bt_can[k]   <- btk_can

      # Log-MH-ratio for
      bt_lratio   <- llfun(b0 = b0_cur, kp = kp_cur, bt = bt_can) -
                     llfun(b0 = b0_cur, kp = kp_cur, bt = bt_cur)

#       xl <- c(-100, 100)
# plot(Vectorize(function(x) llfun(b0=b0_cur, kp=kp_cur, bt=c(1, x))), xlim=xl)

      # MH acceptation
      if (bt_lratio > log(runif(1))) {
        bt_cur[k] <- btk_can
      }
    }

    ### KAPPA
    kp_can     <- rchisq(1, df = kp_cur)

    kp_lratio <- llfun(b0 = b0_cur, kp = kp_can, bt = bt_cur) +
                 dchisq(kp_cur, df = kp_can, log = TRUE) -
                 llfun(b0 = b0_cur, kp = kp_cur, bt = bt_cur) -
                 dchisq(kp_can, df = kp_cur, log = TRUE)

    if (kp_lratio > log(runif(1))) {
      kp_cur <- kp_can
    }

    b0[it]   <- b0_cur
    kp[it]   <- kp_cur
    bt[it, ] <- bt_cur

  }


#   sapply(ls(), function(obj) assign(obj, get(obj), envir = .GlobalEnv))]
  postsample <- data.frame(b0, kp, bt)

  estimates  <- c(beta_0=meanDir(postsample[, 1]), colMeans(postsample[, -1, drop=FALSE]))

  result     <- structure(estimates, postsample=postsample)
}


# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

# # Data generation
# n       <- 100
#
# true_b0 <- pi
# true_bt <- c(3, 6)
# true_kp <- 20
# X       <- cbind(rnorm(n, sd=.2), rnorm(n, sd=.2))
# err     <- rvmc(n, mu = 0, kp = true_kp)
# th      <- (true_b0 + linkfun(apply(X, 1, "%*%", true_bt)) + err) %% (2*pi)
#
# th_bar <- as.vector(mean(circular(th)))
# K       <- ncol(X)
#
# beta_hat <- sapply(1:K, function(k) mlebeta(X=X, k=k, xl = c(-20, 20)))
#
# # Run sampler
# se <- mcmcGCM(th = th, X = X, linkfun = linkfun, invlinkfun = invlinkfun,
#               b0_start=true_b0, bt_start=c(0, 0))

# # Plot sample results
# plot.ts(se$b0)
# plot.ts(se$kp)
# plot.ts(se$beta_1)
# plot.ts(se$beta_2)

#
# estb0 <- mean(s$b0[200:1000])
# estb0
# th_bar %% (2*pi)
# estbt <- mean(s$beta_1[200:1000])
# btind <- 2
#
# predictth <- function(x) (estb0 + linkfun(x * estbt))  %% (2*pi)
#
# plot(X[, 2], th)
# curve(predictth, add = TRUE, col=3, lwd=3)
#
#
# c(1, 2) %*% c(3, 5)
# c(3, 5) %*% c(1, 2)







