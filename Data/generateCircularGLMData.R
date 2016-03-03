library(Rcpp)

sourceCpp("Code/rvmc.cpp")

generateCircGLMData <- function(n=30, residkappa=5, nconpred=2, ncatpred=2,
                                truebeta0 = pi/2,
                                truebeta  = rep(.1, nconpred),
                                truedelta = rep(.1, ncatpred),
                                linkfun   = function(x) 2 * atan(x)) {

  dtpart <- btpart <- 0

  Xcon <- matrix(nr = n, nc = nconpred)
  Xcat <- matrix(nr = n, nc = ncatpred)

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
    btpart <- linkfun(apply(Xcon, 1, "%*%", truebeta))
  }

  if (ncatpred > 0) {
    Xcat  <- sapply(1:ncatpred, function(x) sample(0:1, size=n, replace=TRUE))
    colnames(Xcat) <- paste0("c", 1:ncatpred)
    dtpart <- apply(Xcat, 1, "%*%", truedelta)
  }



  # Generate values for the circular outcome.
  thpred <- truebeta0 + dtpart + btpart

  therr  <- rvmc(n, 0, residkappa)
  th     <- (thpred + therr) %% (2*pi)

  dmat   <- cbind(th, Xcon, Xcat)

  # Save the percentage of data that is found around the true beta_0.
  uLB <- truebeta0-pi/2 %% (2*pi)
  uUB <- truebeta0+pi/2 %% (2*pi)
  if (uLB < uUB) {
    u <- mean(uLB < th & th < uUB)
  } else {
    u <- mean(uLB < th | th < uUB)
  }

  # Add the true values as attributes.
  attr(dmat, "truebeta0")  <- truebeta0
  attr(dmat, "truebeta")   <- truebeta
  attr(dmat, "truezeta")   <- (2/pi)*atan(truebeta)
  attr(dmat, "residkappa") <- residkappa
  attr(dmat, "linkfun")    <- linkfun
  attr(dmat, "u")          <- u

  return(dmat)
}


# Generate circular outcome GLM datasets, and save them to disk.
saveCircGLMDatasets <- function (truens, truekps, truebts,
                                 truebeta0 = pi/2, nsim = 100, seed = 139738) {

  set.seed(seed)

  # All designs we want to generate data for.
  designs <- expand.grid(n=truens, kp=truekps, bt=truebts,
                         stringsAsFactors=FALSE)

  alreadyExistingCounter <- 0

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
                                "_bttype=",   paste0(names(curDesign$bt)))

    # Create a folder for the datasets
    dir.create(saveDirName, showWarnings=FALSE)


    # Repeat generation nsim times.
    for (isim in 1:nsim) {

      writefilename <- paste0(saveDirName, "/nr", isim, ".csv")

      if (!file.exists(writefilename)) {

        # Actually generate the data.
        d <- generateCircGLMData(n = curDesign$n,
                                 residkappa = curDesign$kp,
                                 truebeta   = curDesign$bt[[1]],
                                 truebeta0  = pi/2,
                                 nconpred   = nconpred,
                                 ncatpred   = ncatpred)

        # Write data to chosen filename.
        write.table(d, writefilename, sep = ",", row.names=FALSE, col.names=FALSE)
      } else {
        alreadyExistingCounter <- alreadyExistingCounter + 1
      }

    }


  }
  if (alreadyExistingCounter > 0) {
    cat("\n[Data generation: ", alreadyExistingCounter, "/", nsim*nrow(designs),
        " datasets already existed.]\n")
  }
}







