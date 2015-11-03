library(Rcpp)
library(dplyr)
sourceCpp('DataAnalysis/circGLM.cpp')

logBesselI <- function(nu, x) log(besselI(x, nu, expon.scaled = TRUE)) + x

# Finding nice dimensions for printing a grid of plots.
# n is the number of plots to be printed, tasp is the target aspect ratio.
findPlotGridForm <- function(n, tasp) {
  div <- seq_len(n)
  ftr <- div[n %% div == 0]

  if (n==1) return(c(1, 1))
  if (n==2) return(c(1, 2))

  # If prime, try again with n+1
  if (length(ftr) < 3) {
    return(findPlotGridForm(n+1, tasp))

    # Otherwise figure out which factor provides an aspect ratio closest to the
    # target.
  } else {
    asps <- ftr^2/n
    chs <- ftr[which.min(abs(log(asps) - log(tasp)))]
    return(c(n/chs, chs))
  }
}

# Check if a predictor is dichotomous.
is.dichotomous <- function(x) {
  if (length(unique(x)) == 2) {
    if (all(x == 0 | x == 1)) {
      return(TRUE)
    } else {
      warning("A predictor might be dichotomous but not 0|1.")
    }
  }
  FALSE
}

# Coef is beta or zeta depending on which should be used.
plot.circGLM <- function(m, tasp=1, coef="Beta") {

  if (is.null(m$b0_chain)) {
    warning("No posterior sample saved for this result.")
  } else {


    if (coef == "Beta") {
      coef_chain <- m$bt_chain
    } else if (coef == "Zeta") {
      coef_chain <- m$zt_chain
    } else {
      stop("Coef type not found")
    }

    # Find nice dimensions of the grid of plots.
    npanel <- min(ncol(m$bt_chain)+ncol(m$dt_chain)+2, 16)
    mfr <- findPlotGridForm(npanel, tasp)
    old.par <- par(mfrow = mfr, mar=c(2.5,4,2,1))

    # If the chain is reasonably small, plot the whole chain.
    if (m$SavedIts < 5000) {
      plot.ts(m$b0_chain, xlab="Iteration", ylab="Beta_0", main="Beta_0 chain")
      plot.ts(m$kp_chain, xlab="Iteration", ylab="Kappa", main="Kappa chain")
      apply(coef_chain, 2, plot.ts,
            xlab="Iteration", ylab=coef, main=paste(coef, "chain"))

      apply(m$dt_chain, 2, plot.ts,
            xlab="Iteration", ylab="Delta", main="Delta chain")

    # Otherwise plot only 5000 values, otherwise the plotting will be too slow.
    } else {
      idx <- round(seq(1, m$SavedIts, m$SavedIts/5000))
      plot.ts(x=idx, y=m$b0_chain[idx], xy.lines=TRUE,
              xlab="Iteration", ylab="Beta_0", pch=NA, main="Beta_0 chain")
      plot.ts(x=idx, y=m$kp_chain[idx], xy.lines=TRUE,
              xlab="Iteration", ylab="Kappa", pch=NA, main="Kappa chain")
      apply(coef_chain[idx, , drop=FALSE], 2,
            function(cfchn) plot.ts(x=idx, y=cfchn, xy.lines=TRUE,
                                      xlab="Iteration", ylab=coef,
                                      main=paste(coef, "chain"), pch=NA))
      apply(m$dt_chain[idx, , drop=FALSE], 2,
            function(cfchn) plot.ts(x=idx, y=cfchn, xy.lines=TRUE,
                                      xlab="Iteration", ylab="Delta",
                                      main="Delta chain", pch=NA))
    }
    par(old.par)
  }

}


# Print the results, but never the elements with the full posteriors.
print.circGLM <- function(m) {
  print(m[-grep("chain|Call|data_", names(m))])
  print(lapply(m[grep("chain", names(m))], head))
  invisible(return(NULL))
}


fixResultNames <- function(nms){

  nbts <- length(grep("bt_mean", nms))

  nms[grep("b0_CCI", nms)]     <- c("b0_CCI_LB", "b0_CCI_UB")
  nms[grep("kp_HDI", nms)]     <- c("kp_HDI_LB", "kp_HDI_UB")
  nms[grep("bt_mean", nms)]    <- paste0("bt_", 1:nbts, "_mean")
  nms[grep("zt_mean", nms)]    <- paste0("zt_", 1:nbts, "_mean")
  nms[grep("zt_mdir", nms)]    <- paste0("zt_", 1:nbts, "_mdir")
  nms[grep("bt_propacc", nms)] <- paste0("bt_", 1:nbts, "_propacc")
  nms[grep("bt_CCI", nms)] <- paste0("bt_",
                                     rep(1:nbts, each=2),
                                     c("_LB", "_UB"))
  nms[grep("zt_CCI", nms)] <- paste0("zt_",
                                     rep(1:nbts, each=2),
                                     c("_LB", "_UB"))
  nms
}

# Plot the plot with the first predictor and the first grouping.
predict.plot.circGLM <- function(m) {
  plot(m$data_stX[, 1], m$data_th, pch = 16, col = rgb(m$data_d, 0, 0, 0.6))

  xmin <- min(m$data_stX)
  xmax <- max(m$data_stX)

  sq <- seq(xmin, xmax, length.out = 101)

  pred.grp1 <- m$b0_meandir + atanLF(m$bt_mean %*% t(sq), 2)
  pred.grp2 <- m$b0_meandir + m$dt_meandir[1] + atanLF(m$bt_mean %*% t(sq), 2)

  lines(sq, pred.grp1, col = 1, lwd = 2)
  lines(sq, pred.grp1 + 2*pi, col = 1, lwd = 2)
  lines(sq, pred.grp2, col = 2, lwd = 2)
  lines(sq, pred.grp2 + 2*pi, col = 2, lwd = 2)


  abline(h = c(0, 2*pi), col = "gray80")
  polygon(c(xmin - 1, xmax + 1, xmax + 1, xmin + 1), c(2*pi, 2*pi, 8, 8), col = "gray80")
  polygon(c(xmin - 1, xmax + 1, xmax + 1, xmin + 1), c(0, 0, -1, -1), col = "gray80")
}



circGLM <- function(th, X,
                    conj_prior = rep(0, 3),
                    bt_prior = matrix(0:1, nrow=ncol(X), ncol=2, byrow=TRUE),
                    starting_values = c(0, 1, rep(0, ncol(X))),
                    burnin = 1000,
                    lag = 1,
                    bwb = rep(.05, ncol(X)),
                    kappaModeEstBandwith = .1,
                    CIsize = .95,
                    Q = 10000,
                    r = 2,
                    returnPostSample = FALSE,
                    bt_prior_type=1,
                    output = "list",
                    reparametrize = FALSE,
                    debug = FALSE,
                    loopDebug = FALSE) {

  # Check if the inputs are matrices.
  if (!is.matrix(th)) th <- as.matrix(th)
  if (!is.matrix(X))  X <- as.matrix(X)

  # Check if theta is in radians
  if (any(th > 7)) {
    cat("Setting theta to radians instead of degrees.")
    th <- th * pi / 180
  }

  # Indices of dichotomous predictors.
  dichInd <- apply(X, 2, is.dichotomous)

  # Standardize continuous predictors.
  stX <- scale(X[, !dichInd, drop=FALSE])

  d <- X[, dichInd, drop=FALSE]

  bt_prior <- matrix(0:1, nrow=ncol(stX), ncol=2, byrow=TRUE)

  res <- circGLMC(th=th, X=stX, D=d,
                  conj_prior=conj_prior, bt_prior=bt_prior, #btd_prior=btd_prior,
                  starting_values=starting_values,
                  burnin=burnin, lag=lag, bwb=bwb,
                  kappaModeEstBandwith=kappaModeEstBandwith,
                  CIsize=CIsize,
                  Q=Q, r=r,
                  returnPostSample=returnPostSample,
                  bt_prior_type=bt_prior_type,
                  reparametrize=reparametrize, debug=debug, loopDebug=loopDebug)

  # Set some names for clarity in the output.
  names(res$b0_CCI)     <- c("LB", "UB")
  names(res$kp_HDI)     <- c("LB", "UB")
  names(res$bt_propacc) <- paste0("bt_", 1:length(res$bt_propacc))
  colnames(res$bt_CCI)  <- paste0("bt_", 1:ncol(res$bt_CCI))
  rownames(res$bt_CCI)  <- c("LB", "UB")
  colnames(res$zt_CCI)  <- paste0("zt_", 1:ncol(res$zt_CCI))
  rownames(res$zt_CCI)  <- c("LB", "UB")

  # Add a class 'circGLM', which will make the defined plot and print methods
  # for this class work.
  class(res) <- c("circGLM", class(res))

  # Choose how to return the output.
  if (output == "list") {
    res$Call <- match.call()
    res$data_th   <- th
    res$data_X    <- X
    res$data_d    <- d
    res$data_stX  <- stX

    return(res)

  } else if (output == "vector") {
    if (returnPostSample == TRUE) {message("Vector output with full chains.")}

    out <- unlist(res)
    names(out) <- fixResultNames(names(out))
    return(out)

  } else {
    stop(paste("Output type", output, "not found"))
  }
}
