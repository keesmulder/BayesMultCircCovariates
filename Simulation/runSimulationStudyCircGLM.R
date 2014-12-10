source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')


truens  <- c(10, 30, 100)
truekps <- c(0.5, 5, 30)

# The true Betas MUST be passed with the beta's for continuous predictors first.
truebts        <- list(rep(.1, 2), rep(1, 4), rep(5, 10))
names(truebts) <- c("ll", "llll", "lllll lllll")

truebeta0 <- pi/2
nsim      <- 100
Q         <- 10000

# #
# saveCircGLMDatasets(truens = truens,
#                     truekps = truekps,
#                     truebts = truebts,
#                     truebeta0 = truebeta0,
#                     nsim = nsim)



mcmcpar=list(conj_prior = rep(0, 3),
             burnin = 0,
             lag = 1,
             kappaModeEstBandwith=.05,
             CIsize=.95,
             Q=Q,
             r=2,
             bt_prior_type=1)


bt_prior = list(matrix(c(0,1),
                       nc=2, nr=2,
                       byrow = TRUE),
                matrix(c(0,1),
                       nc=2, nr=4,
                       byrow = TRUE),
                matrix(c(0,1),
                       nc=2, nr=10,
                       byrow = TRUE))

starting_values = list(c(pi+1, 1, rep(0, 2)),
                       c(pi+1, 1, rep(0, 4)),
                       c(pi+1, 1, rep(0, 10)))
bwb=list(rep(.05, 2),
         rep(.05, 4),
         rep(.05, 10))



nbt1 <- 2
betades1 <- list(ll=rep(0.1, nbt1),
                 starting_values=c(0, 1, rep(0, nbt1)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt1, byrow = TRUE),
                 bwb=rep(.05, 2))

nbt2 <- 4
betades2 <- list(llll=rep(0.1, nbt2),
                 starting_values=c(0, 1, rep(0, nbt2)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt2, byrow = TRUE),
                 bwb=rep(.05, nbt2))

nbt3 <- 10
betades3 <- list('lllll lllll'= rep(0.1, nbt3),
                 starting_values=c(0, 1, rep(0, nbt3)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt3, byrow = TRUE),
                 bwb=rep(.05, nbt3))

betaDesigns <- list(betades1, betades2, betades3)






