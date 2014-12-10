library(Rcpp)
library(dplyr)
sourceCpp('DataAnalysis/circGLM.cpp')



fixResultNames <- function(nms){
  nms[grep("b0_CCI", nms)] <- c("b0_CCI_LB", "b0_CCI_UB")
  nms[grep("kp_HDI", nms)] <- c("kp_HDI_LB", "kp_HDI_UB")
  nbts                     <- length(grep("bt_CCI", nms))/2
  nms[grep("bt_CCI", nms)] <- paste0("bt_CCI_",
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
  colnames(res$bt_CCI)  <- paste0("bt_", 1:nrow(res$bt_CCI))
  rownames(res$bt_CCI)  <- c("LB", "UB")

  # Choose how to return the output.
  if (output == "list") {
    return(res)

  } else if (output == "vector") {
    out <- unlist(res)
    names(out) <- fixResultNames(names(out))
    return(out)

  } else {
    stop(paste("Output type", output, "not found"))
  }
}
