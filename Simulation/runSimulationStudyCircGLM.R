source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

# Properties of the simulation study.
truens  <- c(20, 50, 100)
truekps <- c(1, 4, 32)
# truebeta0 is always pi/2 because it should be irrelevant.

nbts       <- rep(c(1, 3, 6), each=4)
truebtvals <- rep(c(0.05, 0.3), times=6)
type       <- rep(c("l","c"), times=3, each=2)

truebts <- sapply(1:length(nbts), function(i) {
  out        <- list(rep(truebtvals[i], nbts[i]))
  names(out) <- paste0(rep(type[i], nbts[i]), collapse="")
  out
})

nsim <- 1000

# Save the datasets as .csv files.
saveCircGLMDatasets(truens = truens, truekps = truekps, truebts = truebts,
                    truebeta0 = truebeta0, nsim = nsim)

# General MCMC parameters.
mcmcpar=list(conj_prior = rep(0, 3), bt_prior_type=1,
             Q=20000, burnin = 1000, lag = 1,
             kappaModeEstBandwith=.05, CIsize=.95,
             r=2, reparametrize=TRUE)

# Generate the designs for Beta
betaDesigns <- lapply(1:length(nbts), function(i){
  c(truebts[i],
    list(starting_values=c(0, 1, rep(0, nbts[i])),
         bt_prior=matrix(c(0,2), nc=2, nr=nbts[i], byrow = TRUE),
         bwb=rep(.05, nbts[i])))
})

# Run the simulation study.
simres <- simStudyCircGLM(nsim = nsim,
                       truens = truens, truekps = truekps,
                       betaDesigns = betaDesigns, overwrite=FALSE,
                       seed = 389238, mcmcpar = mcmcpar)




