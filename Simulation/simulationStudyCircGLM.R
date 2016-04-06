source("Data/generateCircularGLMData.R")
source("DataAnalysis/circGLM.R")

library(circular)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)


# Method for printing simulation results.
print.cGLMSim <- function(res, header=TRUE, digits=2,
                          selection='automatic') {

  oldDgt <- options()$digits
  options(digits=digits)

  if (selection == 'default') {
    selcols <- list('n', 'kp',
                    'b0_meandir', 'b0_in_CCI',
                    'kp_mean', 'kp_mode', 'kp_in_HDI',
                    'crashed')
  }

  resdes <- attr(res, 'args')
  if (header) {
    cat(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",
        "Simulation study results with number of simulations: ", resdes$nsim,
        ". Iterations per dataset Q: ", resdes$mcmcpar$Q, ".",
        "\nRange of the link function is r*pi, here r = ", resdes$mcmcpar$r,
        ". Using an ", ifelse(resdes$predDesigns[[2]]$bt_prior_type,
                              paste0("N(", resdes$predDesigns[[1]]$bt_prior_musd[1], ", ",
                                     resdes$predDesigns[[1]]$bt_prior_musd[2], ")"),
                              "constant"),
        " prior for Beta.\n", "Reparametrization was ",
        ifelse(resdes$mcmcpar$reparametrize, "", "not "), "performed.\n",
        " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n",
        sep="")
  }

  for (i in 1:length(res)) {
    cat("\n(", i, ") predDesign: ",
        gsub('c\\(', "(", names(res)[i]), "\n", sep = "")
    if (selection == "all"){
      print(res[[i]])
    } else if (selection == 'automatic') {
      btcols <- dtcols <- NULL
      if(any(grepl("bt.*mean", names(res[[i]])))) {
        btcols <- c("bt_1_mean", "bt_1_in_CCI")
      }
      if(any(grepl("dt.*mdir", names(res[[i]])))) {
        dtcols <- c("dt_1_mdir", "dt_1_in_CCI")
      }
      selcols <-  c('n', 'kp',
                    'b0_meandir', 'b0_in_CCI',
                    'kp_mean', 'kp_mode', 'kp_in_HDI',
                    btcols, dtcols,
                    'crashed')


      print(do.call(select_, c(list(res[[i]]), .dots = list(selcols))))
    } else if (selection == 'default') {
      print(do.call(select_, c(list(res[[i]]), .dots = list(selcols))))
    }
  }
  options(digits=oldDgt)
}

# Method for obtaining a single vector of length nsim, containing results from a
# single statistic for a design of the simulation study.
slice.cGLMSim <- function(res, pdDesNumber, n, kp, stat) {
  attr(res[[pdDesNumber]], "full")[,
                                   stat,
                                   as.character(n),
                                   as.character(kp)]
}

# Method for plotting simulation results.
plot.cGLMSim <- function(res, pdDesNumber, n, kp, stat, ...) {
  mn <- paste0("Histogram of ", stat, " from predDesign (", pdDesNumber,
               "), n=", n, ", ", "kappa=", kp, ". ")

  hist(slice.cGLMSim(res, pdDesNumber, n, kp, stat),
       main=mn, xlab=stat, ...)
}

# Returns whether an angle th is in circular interval CI. The interval must be
# supplied as {lower bound, upper bound}.
angleInCircInterval <- function (th, CI) {
  if(CI[1] < CI[2]) return(unname(th > CI[1] & th < CI[2]))
  if(CI[1] > CI[2]) return(unname(th > CI[1] | th < CI[2]))
}


# dataFilename is a string denoting where on disk the dataset to be analyzed is.
# params are the parameters to pass to the circGLM analysis function.
runCircGLMSamplerOnce <- function(dataFilename, params,
                                  curTrueBts, curTrueDts) {

  d  <- read.table(dataFilename, sep=",")
  th <- as.matrix(d[, 1,  drop=FALSE])
  X  <- as.matrix(d[, -1, drop=FALSE])
  K  <- ncol(X)
  nbts <- length(curTrueBts)
  ndts <- length(curTrueDts)

  curCombinedRes <- NA

  tryCatch({

    # Run the MCMC sampler.
    curSimRes <- do.call(circGLM,
                         c(list(th=th, X=X, output="vector"), params))

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
    if (any(grepl("bt.*mean", names(curSimRes)))) {
      bt_means <- curSimRes[grep("bt.*mean", names(curSimRes))]
      bt_LBs   <- curSimRes[grep("bt.*LB",   names(curSimRes))]
      bt_UBs   <- curSimRes[grep("bt.*UB",   names(curSimRes))]

      bt_bias        <- bt_means - curTrueBts
      names(bt_bias) <- paste0("bt_", 1:nbts, "_bias")

      bt_in_CCI <- curTrueBts > bt_LBs & curTrueBts < bt_UBs
      names(bt_in_CCI) <- paste0("bt_", 1:nbts, "_in_CCI")

      bt_res <- c(bt_bias,
                  bt_bias_meanOverBt = mean(bt_bias),
                  bt_in_CCI,
                  bt_in_CCI_meanOverBt = mean(bt_in_CCI))
    } else {
      bt_res <- NULL
    }
    ###############





    #### ZETA  ####
    if (any(grepl("zt.*mean", names(curSimRes)))) {

      true_zt        <- atanLF(curTrueBts, 2/pi)
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






    #### DELTA  ####
    if (any(grepl("dt.*mdir", names(curSimRes)))) {
      dt_meandirs <- curSimRes[grep("dt.*mdir", names(curSimRes))]
      dt_LBs      <- curSimRes[grep("dt.*LB",   names(curSimRes))]
      dt_UBs      <- curSimRes[grep("dt.*UB",   names(curSimRes))]

      dt_bias        <- dt_meandirs - curTrueDts
      names(dt_bias) <- paste0("dt_", 1:ndts, "_bias")

      dt_in_CCI <- numeric(ndts)

      for (dti in 1:ndts) {
        dt_in_CCI[dti] <- angleInCircInterval(curTrueDts[dti],
                                              c(dt_LBs[dti], dt_UBs[dti]))
      }

      names(dt_in_CCI) <- paste0("dt_", 1:ndts, "_in_CCI")

      dt_res <- c(dt_bias,
                  dt_bias_meanOverDt = mean(dt_bias),
                  dt_in_CCI,
                  dt_in_CCI_meanOverDt = mean(dt_in_CCI))
    } else {
      dt_res <- NULL
    }
    ###############





    ###############
    # Gather results.
    # True beta_0 should always be generated as pi/2.
    biasCvgRes <- c(b0_res,
                    kp_res,
                    bt_res,
                    zt_res,
                    dt_res,
                    crashed = FALSE)

    curCombinedRes <- c(curSimRes, biasCvgRes)

  }, error = function (e) {
    cat("\n[[!]] Error in ", curReadDir,
        "/nr", isim, ".csv", ":\n", paste(e), sep = "")
  })

  curCombinedRes

}


# returns is a character vector of regular expressions which will match
simStudyCircGLM <- function(truens, truekps, predDesigns,
                            nsim = 100,
                            saveResults=TRUE, overwrite = TRUE,
                            saveLoadEachDgn = TRUE,
                            seed=489734,
                            mcmcpar=list(conj_prior = rep(0, 3),
                                         Q=10000,
                                         burnin = 0,
                                         lag = 1,
                                         kappaModeEstBandwith=.05,
                                         CIsize=.95,
                                         r=2,
                                         reparametrize = TRUE),
                            returns = c()) {


  # Register a parallel computing cluster.
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)


  set.seed(seed)

  # Set those values to list that may have been provided by themselves.
  if(!is.list(predDesigns[[1]])) predDesigns %<>% list

  # The directory containing the datasets.
  datasetsDir <- paste0(getwd(), "/Data/Datasets/")

  # Final result.
  simStudyResults  <- list()


  pdDesName <-  paste0("bt",
                       paste0(unique(unlist(sapply(predDesigns,
                                                   function(x) names(x)[1]))),
                              collapse=","), ",",
                       paste0(unique(unlist(sapply(predDesigns,
                                                   function(x) x[[1]]))),
                              collapse=","),
                       collapse=",")
  baseName <- paste("[simStudCircGLM",
                    paste0("nsim", nsim),
                    paste0("Q",        mcmcpar$Q),
                    paste0("burnin",   mcmcpar$burnin),
                    paste0("r",        mcmcpar$r),
                    paste0("seed",     seed),sep = "_")
  designInFileName <- paste(paste0("n",        paste(truens,  collapse=",")),
                            paste0("kp",       paste(truekps, collapse=",")),
                            paste0(pdDesName, "].rda"), sep = "_")

  simFilePath <- paste0(getwd(), "/SimulationResults/",
                        baseName, designInFileName)

  print(simFilePath)

  if(file.exists(simFilePath)&saveResults) {
    cat("\n ------ \nSimulation file:\n\n", simFileName, "\n\nalready existed. ")
    if (!overwrite) {
      cat("Returning previous file.\n ------ ")
      load(simFilePath)
      return(simStudyResults)
    } else {
      cat("Overwriting.\n ------ \n ")
    }
  }


  # Total number of designs and a counter, to be used for progressbar.
  totalNDesigns      <- length(truens) * length(truekps) * length(predDesigns)
  totalDesignCounter <- 0




  # The outermost loop loops through different values for K. If we have a
  # different number of predictors, the output will look different.
  for(pddes in predDesigns) {

    # pddes <- predDesigns[[1]]

    # Current predictor design.
    curTruePds     <- pddes[1]
    curTruePdTypes <- strsplit(names(curTruePds), "")[[1]]
    conidx         <- curTruePdTypes == "l"
    catidx         <- curTruePdTypes == "c"
    npds           <- length(curTruePds[[1]])
    nbts           <- sum(conidx)
    ndts           <- sum(catidx)

    # Current true values for predictors.
    curTrueBts     <- curTruePds[[1]][conidx]
    curTrueDts     <- curTruePds[[1]][catidx]

    designs    <- expand.grid(n=truens,
                              kp=truekps,
                              pred=curTruePds,
                              starting_values=pddes[2],
                              bt_prior_musd=pddes[3],
                              bt_prior_type=pddes[4],
                              bwb=pddes[5])
    ndesigns   <- nrow(designs)


    curPar     <- c(mcmcpar, pddes[-1], returnPostSample=FALSE)


    # Now go through all the other elements of the design.
    for(idesign in 1:ndesigns) {


      totalDesignCounter <- totalDesignCounter + 1

      curDesign <- designs[idesign, ]

      curReadDir <- paste0("n=",       curDesign$n,
                           "_kp=",     curDesign$kp,
                           "_pred=",   curDesign$pred,
                           "_predtype=", names(curDesign$pred))

      # Progress bar
      pb <- winProgressBar(title = paste0("SimStudyCircGLM with nsim  =  ",
                                          nsim,         "  |  n = ",
                                          curDesign$n,  "  |  kp = " ,
                                          curDesign$kp, "  |  predtype = " ,
                                          names(curDesign$pred),"  |  pred = " ,
                                          curDesign$pred),
                           label=paste0("Design ", idesign, "/", ndesigns, ":",
                                        "0% done. Starting now."),
                           min=0, max=nsim, initial=0,
                           width=1000)

      # Save time to analyze each dataset in order to give an estimate of
      # duration.
      simTimes <- rep(NA, nsim)

      # A string name for the current design, for the current temp file
      desString <- paste0("SimulationResults/Temp/TEMPFILE_", baseName, "_", curReadDir, "].rda")

      if (file.exists(desString) && saveLoadEachDgn) {
        load(file = desString)
        cat("[Tempfile", totalDesignCounter, "loaded]\n")
      } else {

        curDgnResList <- list()

        for (isim in 1:nsim) {



          # Filename to read dataset to analyze from.
          readFileName <- paste0(datasetsDir, curReadDir, "/nr", isim, ".csv")


          # Run the sampler, obtain results on average and coverage
          curSimRes <- runCircGLMSamplerOnce(dataFilename = readFileName,
                                             params = curPar,
                                             curTrueBts, curTrueDts)


          curDgnResList[[isim]] <- curSimRes


          # Progress bar
          simTimes[isim]     <- curSimRes['TimeTaken4']
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




        # Identify errors because they will be length 1 (length(NA) == 1)
        curDgnResLens <- sapply(curDgnResList, length)

        # If all are correct, we can just combine the list to a dataframe.
        if (all(curDgnResLens > 1)) {
          curDgnRes <- t(as.data.frame(curDgnResList))
          rownames(curDgnRes) <- 1:nsim

        } else if (any(curDgnResLens == 1) & !all(curDgnResLens == 1)) {

          curDgnRes <- matrix(NA, nrow=nsim, ncol=max(curDgnResLens))

          # Expand the list if there are errors.
          for (i in 1:nsim) {
            if (curDgnResLens == 1) {
              curDgnRes[i, ] <- rep(NA, max(curDgnResLens))
              curDgnRes[i, "crashed"] <- TRUE
            } else {
              curDgnRes[i, ] <- curDgnResList[[i]]
            }
          }


          # Otherwise, we stop because this whole thing has gone wrong.
        } else {
          stop(paste("In", desString, "a fatal error exists. (No valid runs.)"))
        }

      }


      close(pb)



      if (!file.exists(desString) & saveLoadEachDgn) {
        save(file = desString, list = "curDgnRes")
      }

      # With the first design, create output objects for this predDesign.
      if (idesign == 1) {

        thisPredDesDf    <- cbind(designs,
                                  data.frame(matrix(NA,
                                                    nrow=ndesigns,
                                                    ncol=ncol(curDgnRes))))
        thisPredDesFullArray <- array(NA, c(nsim,
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

      # For the full output, we place the full results into the large array.
      thisPredDesFullArray[, , as.character(curDesign$n),
                           as.character(curDesign$kp)] <- as.matrix(curDgnRes)

      # Return the summary, as this object will be combined by foreach()
      thisPredDesDf[idesign, ] <- c(curDesign, summaryThisDesignRes)

    } # End of 1:ndesigns loop.
#
#     # Combine the results with information about the what design generated the
#     # results.
#     thisPredDesDf <- cbind(designs, thisPredDesResDf)

    # The name of the predictor design.
    pddesName <- paste(names(curTruePds), "=", curTruePds)

    colnames(thisPredDesDf)[-(1:7)] <- names(summaryThisDesignRes)
    simStudyResults[[pddesName]] <- thisPredDesDf

    dimnames(thisPredDesFullArray)[[2]] <- names(summaryThisDesignRes)
    attr(simStudyResults[[pddesName]], "full") <- thisPredDesFullArray

  } # End of predDesigns loop.




  # Stop the parallel cluster
  stopCluster(cl)

  class(simStudyResults) <- c("cGLMSim", class(simStudyResults))
  attr(simStudyResults, 'call') <- match.call()
  attr(simStudyResults, 'args') <- lapply(match.call()[-1], eval)


  if (saveResults) {
    save(simStudyResults, file = simFilePath)
  }

  return(simStudyResults)
}



