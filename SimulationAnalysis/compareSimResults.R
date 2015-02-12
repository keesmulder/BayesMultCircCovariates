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
  out
}
