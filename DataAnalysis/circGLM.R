library(Rcpp)
library(dplyr)
sourceCpp('DataAnalysis/circGLM.cpp')


plot.circGLM <- function(m) {
  par(mfrow = c(2,3))
  plot.ts(m$b0_chain, main="beta_0")
  plot.ts(m$kp_chain, main="kappa")
  apply(m$bt_chain, 2, plot.ts, main="beta")
  par(mfrow = c(1,1))
}



# Print the results, but never the elements with the full posteriors.
print.circGLM <- function(m) {
  print(m[-grep("chain|Call", names(m))])
  print(lapply(m[grep("chain", names(m))], head))
  return(NULL)
}


fixResultNames <- function(nms){

  nbts <- length(grep("bt_mean", nms))

  nms[grep("b0_CCI", nms)]     <- c("b0_CCI_LB", "b0_CCI_UB")
  nms[grep("kp_HDI", nms)]     <- c("kp_HDI_LB", "kp_HDI_UB")
  nms[grep("bt_mean", nms)]    <- paste0("bt_", 1:nbts, "_mean")
  nms[grep("bt_propacc", nms)] <- paste0("bt_", 1:nbts, "_propacc")
  nms[grep("bt_CCI", nms)] <- paste0("bt_",
                                     rep(1:nbts, each=2),
                                     c("_LB", "_UB"))
  nms
}

circGLM <- function(th, X, conj_prior, bt_prior, starting_values,
                    burnin, lag, bwb, kappaModeEstBandwith, CIsize,
                    Q, r, returnPostSample, bt_prior_type,
                    output = "list") {

  res <- circGLMC(th, X, conj_prior, bt_prior, starting_values,
                  burnin, lag, bwb, kappaModeEstBandwith, CIsize,
                  Q, r, returnPostSample, bt_prior_type)

  # Drop each superfluous dimension, which are a result of Rcpp's treatment of
  # Armadillo's arma::mat and arma::vec objects.
  res %<>% lapply(drop)

  # Set some names for clarity in the output.
  names(res$b0_CCI)     <- c("LB", "UB")
  names(res$kp_HDI)     <- c("LB", "UB")
  names(res$bt_propacc) <- paste0("bt_", 1:length(res$bt_propacc))
  colnames(res$bt_CCI)  <- paste0("bt_", 1:ncol(res$bt_CCI))
  rownames(res$bt_CCI)  <- c("LB", "UB")

  # Add a class 'circGLM', which will make the defined plot and print methods
  # for this class work.
  class(res) <- c("circGLM", class(res))

  # Choose how to return the output.
  if (output == "list") {
    res$Call <- match.call()
    return(res)

  } else if (output == "vector") {
    out <- unlist(res)
    names(out) <- fixResultNames(names(out))
    return(out)

  } else {
    stop(paste("Output type", output, "not found"))
  }
}
