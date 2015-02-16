source('Simulation/simulationStudyCircGLM.R')

print.comparecGLMSim <- function (obj) {
  cat("Comparison of ",
      paste(attr(obj, 'objnames'), collapse=", "), ".\n", sep = "")

  cat("Simulation studies have number of simulations: ",
      paste0("(", paste0(attr(obj, 'nsims'), collapse="|"), ")"),
      ".\nIterations per dataset Q: ",
      paste0("(", paste0(attr(obj, 'Qs'), collapse="|"), ")"), ".",
      "\nRange of the link function is r*pi, here r = ",
      paste0("(", paste0(attr(obj, 'rs'), collapse="|"), ")"),
      ".\nUsing a ",
      paste0("(", paste0(attr(obj, 'priors'), collapse="|"), ")"),
      " prior for Beta.\n", sep="")
  print.cGLMSim(obj, header=FALSE)
}

compareSimRes <- function(..., type = "meansd", digits = 2) {
  # Compare two simulation study files in various ways.


  # A list of all the sim study result objects.
  rs <- list(...)

  # Outputfile will have the same structure as the first result object.
  out <- rs[[1]]

  # The amount of results to be compared.
  nrs <- length(rs)

  # Check if type is applicable.
  if (type == "diff" & nrs > 2) {
    stop("Use type = meansd with more than two objects.")
  }

  # Number of beta designs.
  nbts <- length(attr(out, "args")$betaDesigns)



  for(i in 1:nbts) {

    # Indices of the columns that are numeric.
    numericIdx <- which(sapply(out[[i]], class) == 'numeric')

    # Iterate over each numeric cell in the data.frame.
    for(icol in numericIdx) {
      for(irow in 1:nrow(out[[i]])) {

        # The values for this particular cell in the given result objects.
        vals <- sapply(rs, function(r) r[[i]][irow, icol])

        if (type == "meansd") {

          # Calculate mean and sd.
          mn    <- round(mean(vals), digits)
          stdev <- round(sd(vals), digits)

          # Save the standard deviation of the objects.
          out[[i]][irow, icol] <- paste0(mn, " (", stdev, ")")

        } else if (type == "diff") {
          out[[i]][irow, icol] <- abs(vals[1] - vals[2])

        } else if (type == "sidebyside") {
          roundvals <- round(vals, digits)
          out[[i]][irow, icol] <- paste0("(", paste0(roundvals, collapse="|"), ")")
        } else {stop("Unknown type.")}

      }
    }
  }
  class(out) <- c("comparecGLMSim", class(out))
  cl <- match.call()
  attr(out, "objnames")  <- as.character(cl[-c(1, length(cl)-c(1, 0))])
  attr(out, 'nsims')   <- sapply(rs, function(r) attr(r, 'args')$nsim)
  attr(out, "Qs")     <- sapply(rs, function(r) attr(r, 'args')$mcmcpar$Q)
  attr(out, "rs")     <- sapply(rs, function(r) attr(r, 'args')$mcmcpar$r)
  attr(out, "priors") <- sapply(rs, function(r) ifelse(attr(r, 'args')$mcmcpar$bt_prior_type,
                                                       paste0("N(", attr(r, 'args')$betaDesigns[[1]]$bt_prior[1], ", ",
                                                              attr(r, 'args')$betaDesigns[[1]]$bt_prior[2], ")"),
                                                       "constant"))
  out
}
