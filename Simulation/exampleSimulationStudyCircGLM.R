source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

# Properties of the simulation study.
truens  <- c(12, 100)
truekps <- c(0.5, 30)
# truebeta0 is pi/2
nbts       <- c(1, 1, 4, 3, 6, 6)
truebtvals <- c(0.1, -1, 0.1, -1, 0.1, -1)
type       <- c("c", "c", "c", "l", "c", "c")

truebts <- sapply(1:length(nbts), function(i) {
  out        <- list(rep(truebtvals[i], nbts[i]))
  names(out) <- paste0(rep(type[i], nbts[i]), collapse="")
  out
})

nsim <- 10

# Save the datasets as .csv files.
saveCircGLMDatasets(truens = truens, truekps = truekps, truebts = truebts,
                    truebeta0 = truebeta0, nsim = nsim)

# General MCMC parameters.
mcmcpar=list(conj_prior = rep(0, 3), bt_prior_type=0,
             Q=10, burnin = 100, lag = 10,
             kappaModeEstBandwith=.05, CIsize=.95,
             r=2, reparametrize=FALSE)

# Generate the designs for Beta
betaDesigns <- lapply(1:length(nbts), function(i){
  c(truebts[i],
    list(starting_values=c(0, 1, rep(0, nbts[i])),
         bt_prior=matrix(c(0,1), ncol=2, nrow=nbts[i], byrow = TRUE),
         bwb=rep(.05, nbts[i])))
})

# Run the simulation study.
simres <- simStudyCircGLM(nsim = nsim,
                          truens = truens, truekps = truekps,
                          betaDesigns = betaDesigns, overwrite=TRUE,
                          seed = 389243, mcmcpar = mcmcpar)
"F:/Dropbox/Research/BayesMultCircCovariates/Simulation/Results/[simStudCircGLM][nsim10][Q10][burnin100][r2][bt_prior0][seed389243][n12,100][kp0.5,30][btc=0.1,c=-1,cccc=0.1,lll=-1,cccccc=0.1,cccccc=-1].rda"
"F:/Dropbox/Research/BayesMultCircCovariates/Simulation/Results/[simStudCircGLM][nsim2][Q20000][burnin100][r2][bt_prior1][seed389238][n20,50,100][kp1,4,32][btl=0.05,l=0.3,c=0.05,c=0.3,lll=0.05,lll=0.3,ccc=0.05,ccc=0.3,llllll=0.05,llllll=0.3,cccccc=0.05,cccccc=0.3].rda"
simres
