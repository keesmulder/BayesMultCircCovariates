
# source("Code/describeCirc.R")
source("Data/generateCircularGLMData.R")
source("DataAnalysis/circGLM.R")

library(circular)
library(magrittr)
library(dplyr)

makeEmptyNsimArray <- function(obj, nsim){
  olddim <- dim(as.array(obj))
  array(NA, c(nsim, olddim))
}

fixResultNames <- function(nms){
  nms[grep("b0_CCI", nms)] <- c("b0_CCI_LB", "b0_CCI_UB")
  nms[grep("kp_HDI", nms)] <- c("kp_HDI_LB", "kp_HDI_UB")
  nbts                     <- length(grep("bt_CCI", nms))/2
  nms[grep("bt_CCI", nms)] <- paste0("bt_CCI_",
                                     rep(1:nbts, each=2),
                                     c("_LB", "_UB"))
  nms
}


# Returns whether an angle th is in circular interval CI. The interval must be
# supplied as {lower bound, upper bound}.
angleInCircInterval <- function (th, CI) {
  if(CI[1] < CI[2]) return(th > CI[1] & th < CI[2])
  if(CI[1] > CI[2]) return(th > CI[1] | th < CI[2])
}


# output options: 'array', 'df', and 'full'.
simStudyCircGLM <- function(truens, truekps, betaDesigns,
                            nsim = 100, output="array", seed=1,
                            mcmcpar=list(conj_prior = rep(0, 3),
                                         Q=10000,
                                         burnin = 0,
                                         lag = 1,
                                         kappaModeEstBandwith=.05,
                                         CIsize=.95,
                                         r=2,
                                         bt_prior_type=1)) {

  set.seed(seed)

  # Set those values to list that may have been provided by themselves.
  if(!is.list(betaDesigns[[1]])) betaDesigns %<>% list

  # The directory containing the datasets.
  datasetsDir <- paste0(getwd(), "/Data/Datasets/")

  # Final resulting vector.
  simStudyResults <- list()

  # The outermost loop loops through different values for K. If we have a
  # different number of predictors, the output will look different.
  for(btdes in betaDesigns) {

    # Current beta design.
    curTrueBts <- btdes[1]
    nbts       <- length(curTrueBts[[1]])

    designs    <- expand.grid(n=truens,
                              kp=truekps,
                              bt=curTrueBts,
                              starting_values=btdes[2],
                              bt_prior=btdes[3],
                              bwb=btdes[4])
    ndesigns   <- nrow(designs)

    curPar     <- c(mcmcpar, btdes[-1], returnPostSample=FALSE)

    # Perform one test-run to obtain the dimensions of all outputs of circGLM.
    testDesign     <- designs[1, ]
    testrunpar     <- curPar
    testrunpar$lag <- 1; testrunpar$Q <- 5; testrunpar$burnin <- 0
    testrunth      <- runif(testDesign$n)
    testrunX       <- matrix(runif(testDesign$n*length(testDesign$bt[[1]])),
                             nr=testDesign$n)
    testrunres     <- do.call(circGLM, c(list(th = testrunth,
                                              X  = testrunX), testrunpar))

    # Make a template in which we can store output for this betadesign.
    outputTemplate <- c(lapply(testrunres, makeEmptyNsimArray, nsim=nsim),
                        list("biasCvgRes"=matrix(NA, nr=nsim, nc=7+nbts*2)))

    # Number of output arrays to fill.
    nout <- length(outputTemplate)

    # Additional results, to be calculated for each simulation in relation to
    # the true scores.
    biasCvgResNames <- c("b0_bias",
                        "b0_coverage",
                        "kp_mean_bias",
                        "kp_mode_bias",
                        "kp_coverage",
                        paste0("bt_bias_", 1:nbts),
                        paste0("bt_coverage_", 1:nbts))

    totalResultNames <- fixResultNames(c(names(unlist(testrunres)),
                                         biasCvgResNames))



    # Make an empty array for the output if asked for, otherwise empty matrix.
    if (output == "array") {
      thisBtdesArray <- array(NA, c(length(truens),
                                    length(truekps),
                                    length(totalResultNames)),
                              dimnames=list(truens,
                                            truekps,
                                            totalResultNames))
    } else if (output == "df" | output == "dataframe") {
      thisBtdesRight <- data.frame(matrix(NA,
                                          nr=ndesigns,
                                          nc=length(totalResultNames)))
      colnames(thisBtdesRight) <- totalResultNames
      thisBtdesDf    <- cbind(designs, thisBtdesRight)

    } else if (output == "full"){

    } else {
      stop("Unrecognized output type")
    }


    # Now go through all the other elements of the design.
    for (idesign in 1:ndesigns) {

      curDesignResults <- outputTemplate

      curDesign <- designs[idesign, ]

      curReadDir <- paste0("n=",         curDesign$n,
                           "_kp=",       curDesign$kp,
                           "_bt=",       curDesign$bt,
                           "_bttype=", names(curDesign$bt))

      # Progress bar
      pb <- winProgressBar(title = paste0("SimStudyCircGLM with nsim  =  ",
                                          nsim,         "  |  n = ",
                                          curDesign[1], "  |  kp = " ,
                                          curDesign[2], "  |  bt = " ,
                                          curDesign[3]),
                           label=paste0("Design ", idesign, "/", ndesigns, ":",
                                        "0% done. Starting now."),
                           min=0, max=nsim, initial=0,
                           width=1000)

      # Save time to analyze each dataset in order to give an estimate of
      # duration.
      simTimes <- rep(NA, nsim)

      for (isim in 1:nsim) {

        # Filename to read dataset to analyze from.
        readFileName <- paste0(datasetsDir, curReadDir, "/nr", isim, ".csv")

        d  <- read.table(readFileName, sep=",")
        th <- d[, 1]
        X  <- as.matrix(d[, -1])
        K  <- ncol(X)

        curSimRes <- do.call(circGLM, c(list(th=th, X=X, output="vector"), curPar))


        bt_bias        <- drop(curSimRes$bt_mean - curDesign$bt[[1]])
        names(bt_bias) <- paste0("bt_bias_", 1:nbts)
        bt_in_CCI      <- curDesign$bt[[1]] > curSimRes$bt_CCI[1] &
                          curDesign$bt[[1]] < curSimRes$bt_CCI[2]
        names(bt_in_CCI) <- paste0("bt_in_CCI_", 1:nbts)


        # True beta_0 should always be generated as pi/2.
        biasCvgRes <- c("b0_bias"      = curSimRes$b0_meandir - pi/2,
                        "b0_in_CCI"    = angleInCircInterval(pi/2,
                                                             curSimRes$b0_CCI),
                        "kp_mean_bias" = curSimRes$kp_mean - curDesign$kp,
                        "kp_mode_bias" = curSimRes$kp_mode - curDesign$kp,
                        "kp_in_HDI"    = curDesign$kp > curSimRes$kp_HDI[1] &
                                         curDesign$kp < curSimRes$kp_HDI[2],
                        bt_bias,
                        bt_bias_meanOverBt = mean(bt_bias),
                        bt_in_CCI,
                        bt_in_CCI_meanOverBt = mean(bt_in_CCI))

        curCombinedRes <- c(curSimRes, list(biasCvgRes))


        # Place the outputs in the correct locations.
        for(iout in 1:nout) {
          if(length(dim(curDesignResults[[iout]])) == 2)
          {
            curDesignResults[[iout]][isim, ]   <- curCombinedRes[[iout]]
          } else if(length(dim(curDesignResults[[iout]])) == 3)
          {
            curDesignResults[[iout]][isim, , ] <- curCombinedRes[[iout]]
          }
        }


        # Progress bar
        simTimes[isim]    <- curSimRes$TimeTaken
        avgTimeTaken      <- mean(simTimes, na.rm = TRUE)
        timeRemaining     <- (nsim - isim) * avgTimeTaken
        dispTimeRemaining <- ifelse(timeRemaining > 120,
                                    paste(round(timeRemaining/60, 0),"minutes"),
                                    paste(round(timeRemaining, 1),"seconds"))

        info <- paste0("\t Design ", idesign, "/", ndesigns, ":  ",
                       round((isim/nsim)*100), "% done. \t Mean time",
                       " per dataset ", round(avgTimeTaken, 2),
                       " sec. \t Est. to finish this design in ",
                       dispTimeRemaining, ", at ",
                       Sys.time() + timeRemaining, ".")
        setWinProgressBar(pb, isim, label=info)

      } # End of 1:nsim loop.

      colnames(curDesignResults$biasCvgRes) <- names(biasCvgRes)

      # AGGREGATE RESULTS INTO ONE OBJECT (ARRAY/MATRIX)

      btdesOutName <- paste0(names(curTrueBts), "=", curTrueBts)

      if (output == "full") {
        simStudyResults[[btdesOutName]] <- curDesignResults
      }

      close(pb)

    } # End of 1:ndesigns loop.

  } # End of betaDesigns loop.

}






