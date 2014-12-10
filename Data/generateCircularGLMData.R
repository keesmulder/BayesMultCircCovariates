library(Rcpp)

sourceCpp("Code/rvmc.cpp")

generateCircGLMData <- function(n=30, residkappa=5, nconpred=2, ncatpred=2,
                                truebeta0 = pi/2,
                                truebeta  = rep(.1, nconpred+ncatpred),
                                linkfun   = function(x) 2 * atan(x),
                                seed ="random") {
  if (!is.numeric(seed) && seed == "random") {
    seed  <- sample(1:1000000, 1)
  }
  set.seed(seed)

  # Check whether true parameter values must be drawn.
  if (!is.numeric(truebeta0)  && truebeta0 == "random") {
    truebeta0 <- runif(1, -2, 8)
  }
  if (!is.numeric(truebeta)   && truebeta == "random") {
    truebeta  <- rnorm(nconpred+ncatpred, sd = 0.5)
  }
  if (!is.numeric(residkappa) && residkappa == "random") {
    residkappa  <- rchisq(1, 10)
  }

  # Generate predictors
  if (nconpred > 0) {
    Xcon <- sapply(1:nconpred, function(x) scale(rnorm(n)))
    colnames(Xcon)  <- paste0("l", 1:nconpred)
  }

  if (ncatpred > 0) {
    Xcat  <- sapply(1:ncatpred, function(x) sample(0:1, size=n, replace=TRUE))
    colnames(Xcat) <- paste0("c", 1:ncatpred)
  }

  # Combine continuous and categorical predictors.
  if (nconpred < 1 & ncatpred < 1) {
    stop("No predictors.")
  } else if (nconpred > 0 & ncatpred < 1) {
    X <- Xcon
  } else if (nconpred < 1 & ncatpred > 0) {
    X <- Xcat
  } else {
    X <- cbind(Xcon, Xcat)
  }

  # Generate values for the circular outcome.
  thpred <- truebeta0 + linkfun(apply(X, 1, "%*%", truebeta))
  therr  <- rvmc(n, 0, residkappa)
  th     <- thpred + therr

  dmat   <- cbind(th, X)

  # Add the true values as attributes.
  attr(dmat, "truebeta0")  <- truebeta0
  attr(dmat, "truebeta")   <- truebeta
  attr(dmat, "residkappa") <- residkappa
  attr(dmat, "linkfun")    <- linkfun

  return(dmat)
}


# Generate circular outcome GLM datasets, and save them to disk.
saveCircGLMDatasets <- function (truens, truekps, truebts,
                                 truebeta0 = pi/2, nsim = 100) {

  # All designs we want to generate data for.
  designs <- expand.grid(n=truens, kp=truekps, bt=truebts,
                         stringsAsFactors=FALSE)

  # Go through all designs.
  for (ides in 1:nrow(designs)) {

    curDesign    <- designs[ides, ]

    # The amount of categorical and continuous predictors is passed as the names
    # of each set of betas.
    curbttype <- strsplit(names(curDesign$bt), split = "")[[1]]
    ncatpred <- sum(curbttype == "c")
    nconpred <- sum(curbttype == "l")

    # Directory for placing the datasets.
    saveDirName       <- paste0(getwd(),
                            "/Data/Datasets/",
                            "n=",         curDesign$n,
                            "_kp=",       curDesign$kp,
                            "_bt=",       curDesign$bt,
                            "_bttype=",   paste0(curbttype))

    # Create a folder for the datasets
    dir.create(saveDirName, showWarnings=FALSE)

    # Repeat generation nsim times.
    for (isim in 1:nsim) {

      writefilename <- paste0(saveDirName, "/nr", isim, ".csv")


      # Actually generate the data.
      d <- generateCircGLMData(n = curDesign$n,
                               residkappa = curDesign$kp,
                               truebeta   = curDesign$bt[[1]],
                               truebeta0  = pi/2,
                               nconpred   = nconpred,
                               ncatpred   = ncatpred)

      write.table(d, writefilename, sep = ",", row.names=FALSE, col.names=FALSE)
    }
  }
}







