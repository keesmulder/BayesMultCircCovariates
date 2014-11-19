library(Rcpp)

sourceCpp("Code/rvmc.cpp")

generateCircGLMData <- function(n=30, residkappa=5, nconpred=2, ncatpred=2,
                                truebeta0=1,
                                truebeta=rep(.1, nconpred+ncatpred),
                                linkfun    = function(x) 2 * atan(x),
                                invlinkfun = function(x) tan(x/2)) {
  if (!is.numeric(truebeta0)  && truebeta0 == "random") {
    truebeta0 <- runif(1, -2, 8)
  }
  if (!is.numeric(truebeta)   && truebeta == "random") {
    truebeta  <- rnorm(nconpred+ncatpred, sd = 0.5)
  }
  if (!is.numeric(residkappa) && residkappa == "random") {
    residkappa  <- rchisq(1, 10)
  }


  if (nconpred > 0) {
    Xcon <- sapply(1:nconpred, function(x) scale(rnorm(n)))
    colnames(Xcon)  <- paste0("l", 1:nconpred)
  }

  if (ncatpred > 0) {
    Xcat  <- sapply(1:ncatpred, function(x) sample(0:1, size=n, replace=TRUE))
    colnames(Xcat) <- paste0("c", 1:ncatpred)
  }

  if (nconpred < 1 & ncatpred < 1) {
    stop("No predictors.")
  } else if (nconpred > 0 & ncatpred < 1) {
    X <- Xcon
  } else if (nconpred < 1 & ncatpred > 0) {
    X <- Xcat
  } else {
    X <- cbind(Xcon, Xcat)
  }

  thpred <- truebeta0 + linkfun(apply(X, 1, "%*%", truebeta))
  therr  <- rvmc(n, 0, residkappa)
  th     <- thpred + therr

  dmat   <- cbind(th, X)

  attr(dmat, "truebeta0")  <- truebeta0
  attr(dmat, "truebeta")   <- truebeta
  attr(dmat, "linkfun")    <- linkfun
  attr(dmat, "residkappa") <- residkappa

  return(dmat)
}
