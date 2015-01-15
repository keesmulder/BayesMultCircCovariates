source("Data/generateCircularGLMData.R")
source("DataAnalysis/circGLM.R")

library(circular)
library(magrittr)
library(dplyr)


# Method for printing simulation results.
print.cGLMSim <- function(res,
                          selection=c('n', 'kp',
                                      'b0_meandir', 'b0_in_CCI'
                                      'kp_mean', 'kp_mode', 'kp_in_HDI',
                                      'bt_1_mean', 'bt_1_in_CCI', 'crashed')) {
  for (i in 1:length(res)) {
    cat("\n\n(", i, ") Betadesign: ", names(res)[i], "\n", sep = "")
    print(res[[i]] %>%
#             tbl_df() %>%
            select(n, kp,
                   kp_mean, kp_mode, kp_in_HDI, b0_bias, bt_1_mean, crashed))
  }
}

# Returns whether an angle th is in circular interval CI. The interval must be
# supplied as {lower bound, upper bound}.
angleInCircInterval <- function (th, CI) {
  if(CI[1] < CI[2]) return(unname(th > CI[1] & th < CI[2]))
  if(CI[1] > CI[2]) return(unname(th > CI[1] | th < CI[2]))
}


# output options: 'array', 'df', and 'full'.
simStudyCircGLM <- function(truens, truekps, betaDesigns,
                            nsim = 100,
                            saveResults=TRUE, overwrite = TRUE,
                            seed=489734,
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

  # Final result.
  simStudyResults  <- list()


  btDesName <-  paste0("(bt)",
                       paste0(sapply(betaDesigns,
                                     function(x) paste0(names(x)[1], "=", paste(x[[1]], collapse=","))),
                              collapse="_"))
  simFileName <- paste("[simStudCircGLM",
                        paste0("(nsim)", nsim),
                        paste0("(n)", paste(truens,  collapse=",")),
                        paste0("(kp)", paste(truekps, collapse=",")),
                        paste0(btDesName, "].rda"), sep = "]__[")
  simFilePath <- paste0(getwd(), "/Simulation/Results/", simFileName)

  if(file.exists(simFilePath)&saveResults) {
    cat("Simulation file:\n\n", simFileName, "\n\nalready existed. ")
    if (!overwrite) {
      cat("Returning previous file.")
      load(simFilePath)
      return(simStudyResults)
    } else {
      cat("Overwriting.")
    }
  }

  # Total number of designs and a counter, to be used for progressbar.
  totalNDesigns      <- length(truens) * length(truekps) * length(betaDesigns)
  totalDesignCounter <- 0

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



    # Amount of outcomes we get for each beta: mean, Lower Bound, Upper Bound,
    # and proportion accepted.
    nBtOutcomes    <- 6

    # Amount of outcomes that do not involve beta.
    nNonBtOutcomes <- 19

    # Total number outcomes.
    nOutcomes      <- nNonBtOutcomes + nbts*nBtOutcomes

    # Make a template in which we can store output for this betadesign.
    outputTemplate <- matrix(NA, nr=nsim, nc=nOutcomes)

    thisBtdesDf    <- cbind(designs, data.frame(matrix(NA,
                                                       nr=ndesigns,
                                                       nc=nOutcomes)))
    thisBtdesFullArray <- array(NA, c(nsim,
                                      nOutcomes,
                                      length(truens),
                                      length(truekps)),
                                dimnames=list(1:nsim,
                                              NULL,
                                              truens,
                                              truekps))

    # Now go through all the other elements of the design.
    for (idesign in 1:ndesigns) {

      curDgnRes <- outputTemplate

      totalDesignCounter <- totalDesignCounter + 1

      curDesign <- designs[idesign, ]

      curReadDir <- paste0("n=",       curDesign$n,
                           "_kp=",     curDesign$kp,
                           "_bt=",     curDesign$bt,
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


        tryCatch({

          curSimRes <- do.call(circGLM,
                               c(list(th=th, X=X, output="vector"), curPar))

          bt_means <- curSimRes[grep("bt.*mean", names(curSimRes))]
          bt_LBs   <- curSimRes[grep("bt.*LB",   names(curSimRes))]
          bt_UBs   <- curSimRes[grep("bt.*UB",   names(curSimRes))]

          bt_bias        <- bt_means - curDesign$bt[[1]]
          names(bt_bias) <- paste0("bt_", 1:nbts, "_bias")

          bt_in_CCI <- curDesign$bt[[1]] > bt_LBs & curDesign$bt[[1]] < bt_UBs
          names(bt_in_CCI) <- paste0("bt_", 1:nbts, "_in_CCI")

          b0_CCI <- curSimRes[grep("^b0_CCI", names(curSimRes))]

          # True beta_0 should always be generated as pi/2.
          biasCvgRes <- c("b0_bias"      = curSimRes[['b0_meandir']] - pi/2,
                          "b0_in_CCI"    = angleInCircInterval(pi/2, b0_CCI),
                          "kp_mean_bias" = curSimRes[['kp_mean']] - curDesign$kp,
                          "kp_mode_bias" = curSimRes[['kp_mode']] - curDesign$kp,
                          "kp_in_HDI"    = curDesign$kp > curSimRes[['kp_HDI_LB']] &
                            curDesign$kp < curSimRes[['kp_HDI_UB']],
                          bt_bias,
                          bt_bias_meanOverBt = mean(bt_bias),
                          bt_in_CCI,
                          bt_in_CCI_meanOverBt = mean(bt_in_CCI),
                          crashed = FALSE)

          curCombinedRes <- c(curSimRes, biasCvgRes)

        }, error = function (e) {

          cat("\n Error in ", curReadDir, "/nr", isim, ".csv", ":\n", paste(e))

          curCombinedRes <<- rep(NA, nOutcomes)
          curCombinedRes[nOutcomes] <<- TRUE
        })


        curDgnRes[isim, ] <- curCombinedRes

        # Progress bar
        simTimes[isim]     <- curSimRes['TimeTaken']
        avgTimeTaken       <- mean(simTimes, na.rm = TRUE)
        timeRemaining      <- (nsim - isim) * avgTimeTaken
        dispTimeRemaining  <- ifelse(timeRemaining > 120,
                                     paste(round(timeRemaining/60, 0),"minutes"),
                                     paste(round(timeRemaining, 1),"seconds"))

        info <- paste0("\t Design ", totalDesignCounter, "/", totalNDesigns, ":  ",
                       round((isim/nsim)*100), "% done. \t Mean time",
                       " per dataset ", round(avgTimeTaken, 2),
                       " sec. \t Est. to finish this design in ",
                       dispTimeRemaining, ", at ",
                       Sys.time() + timeRemaining, ".")
        setWinProgressBar(pb, isim, label=info)

      } # End of 1:nsim loop.

      # Set column names for results
      colnames(curDgnRes) <- names(curCombinedRes)

      close(pb)

      # Unless we are saving information from every simulation, we summarize the
      # results over 1:nsim results which are in curDgnRes.
      summaryThisDesignRes <- c("b0_meandir" =
                                  computeMeanDirection(na.omit(curDgnRes[, 1])),
                                apply(curDgnRes[, -1], 2, mean, na.rm=TRUE))

      thisBtdesDf[idesign, -(1:6)] <- summaryThisDesignRes

      # For the full output, we place the full results into the large array.
      thisBtdesFullArray[, , as.character(curDesign$n),
                         as.character(curDesign$kp)] <- curDgnRes

    } # End of 1:ndesigns loop.

    btdesName <- paste(names(curTrueBts), "=", curTrueBts)

    colnames(thisBtdesDf)[-(1:6)] <- names(summaryThisDesignRes)
    simStudyResults[[btdesName]] <- thisBtdesDf

    dimnames(thisBtdesFullArray)[[2]] <- names(summaryThisDesignRes)
    attr(simStudyResults[[btdesName]], "full") <- thisBtdesFullArray

  } # End of betaDesigns loop.

  class(simStudyResults) <- c("cGLMSim", class(simStudyResults))

  if (saveResults) {
    save(simStudyResults, file = simFilePath)
  }

  return(simStudyResults)
}



