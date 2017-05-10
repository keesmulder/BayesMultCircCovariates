library(ggplot2)


obtainBFResults <- function(SSRlist, selNum = 1) {

  allBFDfs <- list()

  nms <- names(SSRlist)

  for (SSRi in 1:length(SSRlist)) {

    thisSSR <- SSRlist[[SSRi]]

    for (pdi in 1:length(thisSSR)) {

      thisPd <- thisSSR[[pdi]]

      thisPdName <- thisSSR[[pdi]]$pred[[1]][1]

      thisFull <- attr(thisPd, "full")


      for (j in 1:nrow(thisPd)) {
        thisN  <- as.character(thisPd[j, "n"])
        thisKp <- as.character(thisPd[j, "kp"])

        thisRes <- data.frame(thisFull[, , thisN, thisKp])

        thisFirstIneqBF <- grep("IneqBayesFactors", colnames(thisRes))[1]
        thisFirstEqBF   <- grep("SDDBayesFactors", colnames(thisRes))[1]

        thisIneqBFs  <- thisRes[, thisFirstIneqBF]
        thisEqBFs    <- thisRes[, thisFirstEqBF]


        thisBFDf <- data.frame(Model = nms[SSRi],
                               IneqBFName = colnames(thisRes)[thisFirstIneqBF],
                               EqBFName = colnames(thisRes)[thisFirstEqBF],
                               n = thisN,
                               Kappa = thisKp,
                               Coefficient = thisPdName,
                               Ineq = thisIneqBFs,
                               Eq = thisEqBFs,
                               stringsAsFactors = FALSE)

        thisBFDfName <- paste0(nms[SSRi], thisN, thisKp, thisPdName)

        allBFDfs[[thisBFDfName]] <- thisBFDf
      }
    }
  }


  NDfs      <- length(allBFDfs)
  eachNrow  <- nrow(allBFDfs[[1]])
  superNcol <- 8
  superNrow <- NDfs * eachNrow

  superDF <- data.frame(matrix(NA, nrow = superNrow, ncol = superNcol))
  colnames(superDF) <- colnames(allBFDfs[[1]])

  for (si in 1:NDfs) {
    start <- ((si - 1) * eachNrow) + 1

    superDF[start:(start+eachNrow-1), ] <- allBFDfs[[si]]
  }

  superDF$n     <- factor(superDF$n,     levels = rev(sort(unique(superDF$n))))
  superDF$Model <- factor(superDF$Model, levels = rev(sort(unique(superDF$Model))))

  superDF
}




