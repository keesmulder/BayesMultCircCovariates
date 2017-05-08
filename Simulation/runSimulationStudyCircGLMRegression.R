source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')
source('Simulation/parallelizeSimulationStudyCircGLM.R')


library(dplyr)

makePredList <- function(pg) {
  pl <- list()
  for (i in 1:nrow(pg)) {
    pl[[i]] <- rep(pg$pred[i], pg$n[i])

    if (pg[i, "type"] == "anc") {
      names(pl)[i] <- paste0(c("c", rep("l", pg$n[i]-1)), collapse="")
    } else {
      names(pl)[i] <- paste0(rep(pg$type[i], pg$n[i]), collapse="")
    }
  }
  pl
}




# Properties of the simulation study.
truens  <- c(20, 100)
truekps <- c(2, 20)
# truebeta0 is always pi/2 because it should be irrelevant.

npreds   <- c(1)
truepred <- c(.05, .8)
type     <- c("l")

pg <- expand.grid(n=npreds, pred=truepred, type=type)

pl <- makePredList(pg)

nsim <- 5000

# Save the datasets as .csv files.
saveCircGLMDatasets(truens = truens, truekps = truekps, truepreds = pl,
                    truebeta0 = truebeta0, nsim = nsim)

# General MCMC parameters.
mcmcpar=list(conj_prior = rep(0, 3),
             Q=20000, burnin = 1000, lag = 1,
             kappaModeEstBandwith=.05, CIsize=.95,
             r=2, reparametrize=TRUE)

# Generate the designs for predictors
predDesigns <- lapply(1:nrow(pg), function(i){
  c(pl[i],
    list(starting_values=c(0, 1, rep(0, pg$n[i])),
         bt_prior_musd=c("mu"=0, "sd"=1),
         bt_prior_type=1,
         bwb=rep(.05, pg$n[i])))
})

# Run the simulation study.
simres <- simStudyCircGLM(nsim = nsim, saveLoadEachDgn = TRUE,
                       truens = truens, truekps = truekps,
                       predDesigns = predDesigns, overwrite=FALSE,
                       seed = 38944, mcmcpar = mcmcpar)













