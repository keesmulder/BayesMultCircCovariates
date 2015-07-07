source("Data/generateCircularGLMData.R")
source("DataAnalysis/circGLM.R")

library(circular)
library(dplyr)


# Method for printing simulation results.
print.cGLMSim <- function(res, header=TRUE, digits=2,
                          selection=list('n', 'kp',
                                         'b0_meandir', 'b0_in_CCI',
                                         'kp_mean', 'kp_mode', 'kp_in_HDI',
                                         'bt_1_mean', 'bt_1_in_CCI',
                                         #     'zt_1_mean', 'zt_1_mdir', 'zt_1_in_CCI',
                                         'crashed')) {

  oldDgt <- options()$digits
  options(digits=digits)

  resdes <- attr(res, 'args')
  if (header) {
    cat(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",
        "Simulation study results with number of simulations: ", resdes$nsim,
        ". Iterations per dataset Q: ", resdes$mcmcpar$Q, ".",
        "\nRange of the link function is r*pi, here r = ", resdes$mcmcpar$r,
        ". Using a ", ifelse(resdes$mcmcpar$bt_prior_type,
                             paste0("N(", resdes$betaDesigns[[1]]$bt_prior[1], ", ",
                                    resdes$betaDesigns[[1]]$bt_prior[2], ")"),
                             "constant"),
        " prior for Beta.\n", "Reparametrization was ",
        ifelse(resdes$mcmcpar$reparametrize, "", "not "), "performed.\n",
        " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",
        sep="")
  }

  for (i in 1:length(res)) {
    cat("\n(", i, ") Betadesign: ",
        gsub('c\\(', "(", names(res)[i]), "\n", sep = "")
    if (selection == "all"){
      print(res[[i]])
    } else {
      print(do.call(select_, c(list(res[[i]]), .dots = list(selection))))
    }
  }
  options(digits=oldDgt)
}

# Method for obtaining a single vector of length nsim, containing results from a
# single statistic for a design of the simulation study.
slice.cGLMSim <- function(res, btDesNumber, n, kp, stat) {
  attr(res[[btDesNumber]], "full")[,
                                   stat,
                                   as.character(n),
                                   as.character(kp)]
}

# Method for plotting simulation results.
plot.cGLMSim <- function(res, btDesNumber, n, kp, stat, ...) {
  mn <- paste0("Histogram of ", stat, " from Betadesign (", btDesNumber,
               "), n=", n, ", ", "kappa=", kp, ". ")

  hist(slice.cGLMSim(res, btDesNumber, n, kp, stat),
       main=mn, xlab=stat, ...)
}

# Returns whether an angle th is in circular interval CI. The interval must be
# supplied as {lower bound, upper bound}.
angleInCircInterval <- function (th, CI) {
  if(CI[1] < CI[2]) return(unname(th > CI[1] & th < CI[2]))
  if(CI[1] > CI[2]) return(unname(th > CI[1] | th < CI[2]))
}

# returns is a character vector of regular expressions which will match
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
                                         reparametrize   = TRUE,
                                         bt_prior_type=1),
                            returns = c()) {

  set.seed(seed)

  # Set those values to list that may have been provided by themselves.
  if(!is.list(betaDesigns[[1]])) betaDesigns %<>% list

  # The directory containing the datasets.
  datasetsDir <- paste0(getwd(), "/Data/Datasets/")

  # Final result.
  simStudyResults  <- list()


  btDesName <-  paste0("bt",
                       paste0(sapply(betaDesigns,
                                     function(x) {
                                       if (length(unique(x[[1]])) == 1) {
                                         return(paste0(names(x)[1],
                                                       "=", x[[1]][1], collapse=","))
                                       } else {
                                         return(paste0(names(x)[1], "=",
                                                       paste(x[[1]], collapse=",")))
                                       }
                                     }
                       ),
                       collapse=","))
  btDesName <-  paste0("bt",
                       paste0(unique(unlist(sapply(betaDesigns,
                                                   function(x) names(x)[1]))),
                              collapse=","), ",",
                       paste0(unique(unlist(sapply(betaDesigns,
                                                   function(x) x[[1]]))),
                              collapse=","),
                       collapse=",")
  simFileName <- paste("[simStudCircGLM",
                       paste0("nsim", nsim),
                       paste0("Q",        mcmcpar$Q),
                       paste0("burnin",   mcmcpar$burnin),
                       paste0("r",        mcmcpar$r),
                       paste0("bt_prior", mcmcpar$bt_prior_type),
                       paste0("seed",     seed),
                       paste0("n",        paste(truens,  collapse=",")),
                       paste0("kp",       paste(truekps, collapse=",")),
                       paste0(btDesName, "].rda"), sep = "_")
  simFilePath <- paste0(getwd(), "/Simulation/Results/", simFileName)

  print(simFilePath)

  if(file.exists(simFilePath)&saveResults) {
    cat("\n ------ \nSimulation file:\n\n", simFileName, "\n\nalready existed. ")
    if (!overwrite) {
      cat("Returning previous file.\n ------ ")
      load(simFilePath)
      return(simStudyResults)
    } else {
      cat("Overwriting.\n ------ ")
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


    # Now go through all the other elements of the design.
    for (idesign in 1:ndesigns) {

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
        th <- as.matrix(d[, 1,  drop=FALSE])
        X  <- as.matrix(d[, -1, drop=FALSE])
        K  <- ncol(X)



        tryCatch({

          # Run the MCMC sampler.
          curSimRes <- do.call(circGLM,
                               c(list(th=th, X=X, output="vector"), curPar))

          # Calculate bias and coverage results.
          #### BETA0 ####
          b0_CCI <- curSimRes[grep("^b0_CCI", names(curSimRes))]
          b0_res <- c("b0_bias"      = curSimRes[['b0_meandir']] - pi/2,
                      "b0_in_CCI"    = angleInCircInterval(pi/2, b0_CCI))
          ###############



          #### KAPPA ####
          kp_res <- c("kp_mean_bias" = curSimRes[['kp_mean']] - curDesign$kp,
                      "kp_mode_bias" = curSimRes[['kp_mode']] - curDesign$kp,
                      "kp_in_HDI"    = curDesign$kp > curSimRes[['kp_HDI_LB']] &
                        curDesign$kp < curSimRes[['kp_HDI_UB']])
          ###############



          #### BETA  ####
          bt_means <- curSimRes[grep("bt.*mean", names(curSimRes))]
          bt_LBs   <- curSimRes[grep("bt.*LB",   names(curSimRes))]
          bt_UBs   <- curSimRes[grep("bt.*UB",   names(curSimRes))]

          bt_bias        <- bt_means - curDesign$bt[[1]]
          names(bt_bias) <- paste0("bt_", 1:nbts, "_bias")

          bt_in_CCI <- curDesign$bt[[1]] > bt_LBs & curDesign$bt[[1]] < bt_UBs
          names(bt_in_CCI) <- paste0("bt_", 1:nbts, "_in_CCI")

          bt_res <- c(bt_bias,
                      bt_bias_meanOverBt = mean(bt_bias),
                      bt_in_CCI,
                      bt_in_CCI_meanOverBt = mean(bt_in_CCI))
          ###############



          #### ZETA  ####
          if (any(grepl("zt.*mean", names(curSimRes)))) {

            true_zt        <- atanLF(curDesign$bt[[1]], 2/pi)
            zt_means       <- curSimRes[grep("zt.*mean", names(curSimRes))]
            zt_meandirs    <- curSimRes[grep("zt.*mdir", names(curSimRes))]
            zt_LBs         <- curSimRes[grep("zt.*LB",   names(curSimRes))]
            zt_UBs         <- curSimRes[grep("zt.*UB",   names(curSimRes))]

            zt_meanbias        <- zt_means    - true_zt
            zt_mdirbias        <- zt_meandirs - true_zt
            names(zt_meanbias) <- paste0("zt_", 1:nbts, "_meanbias")
            names(zt_mdirbias) <- paste0("zt_", 1:nbts, "_mdirbias")

            zt_in_CCI <- true_zt > zt_LBs & true_zt < zt_UBs
            names(zt_in_CCI) <- paste0("zt_", 1:nbts, "_in_CCI")

            zt_res <- c(zt_meanbias,
                        zt_mdirbias,
                        zt_bias_meanOverBt = mean(zt_meanbias),
                        zt_in_CCI,
                        zt_in_CCI_meanOverBt = mean(zt_in_CCI))
          } else {
            zt_res <- NULL
          }
          ###############



          ###############
          # Gather results.
          # True beta_0 should always be generated as pi/2.
          biasCvgRes <- c(b0_res,
                          kp_res,
                          bt_res,
                          zt_res,
                          crashed = FALSE)

          curCombinedRes <- c(curSimRes, biasCvgRes)

        }, error = function (e) {

          cat("\n[[!]] Error in ", curReadDir,
              "/nr", isim, ".csv", ":\n", paste(e), sep = "")

          curCombinedRes <<- rep(NA, ncol(curDgnRes))
          curCombinedRes["crashed"] <<- TRUE
        })

        # At the first dataset, make a current design result object that can
        # hold the length of the output currently given.
        if (isim == 1) {
          curDgnRes <- matrix(nrow = nsim, ncol = length(curCombinedRes))

          # Set column names for results
          colnames(curDgnRes) <- names(curCombinedRes)
        }

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
      close(pb)

      # With the first design, create output objects for this Betadesign.
      if (idesign == 1) {

        thisBtdesDf    <- cbind(designs,
                                data.frame(matrix(NA,
                                                  nrow=ndesigns,
                                                  ncol=ncol(curDgnRes))))
        thisBtdesFullArray <- array(NA, c(nsim,
                                          ncol(curDgnRes),
                                          length(truens),
                                          length(truekps)),
                                    dimnames=list(1:nsim,
                                                  NULL,
                                                  truens,
                                                  truekps))
      }

      # Indices of the columns for which we must take the circular mean.
      circIdx <- which(colnames(curDgnRes) %in% c("b0_meandir", "b0_bias"))

      # Unless we are saving information from every simulation, we summarize the
      # results over 1:nsim results which are in curDgnRes.
      summaryThisDesignRes <- c("b0_meandir" =
                                  computeMeanDirection(na.omit(curDgnRes[, "b0_meandir"])),
                                "b0_bias" =
                                  computeMeanDirection(na.omit(curDgnRes[, "b0_bias"])),
                                apply(curDgnRes[, -circIdx], 2, mean, na.rm=TRUE))

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
  attr(simStudyResults, 'call') <- match.call()
  attr(simStudyResults, 'args') <- lapply(match.call()[-1], eval)


  if (saveResults) {
    save(simStudyResults, file = simFilePath)
  }

  return(simStudyResults)
}



