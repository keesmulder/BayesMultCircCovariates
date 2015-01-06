source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

truens  <- c(10, 30, 50)
truekps <- c(0.5, 5, 30)

# The true Betas MUST be passed with the beta's for continuous predictors first.
truebts        <- list(rep(.1, 2), rep(-2, 2), rep(.1, 10), rep(-2, 10))
names(truebts) <- c("ll", "ll", "llllllllll", "llllllllll")

truebeta0 <- pi/2
nsim      <- 10
Q         <- 10000

# #
saveCircGLMDatasets(truens = truens,
                    truekps = truekps,
                    truebts = truebts,
                    truebeta0 = truebeta0,
                    nsim = nsim)



mcmcpar=list(conj_prior = rep(0, 3),
             burnin = 0,
             lag = 1,
             kappaModeEstBandwith=.05,
             CIsize=.95,
             Q=Q,
             r=2,
             bt_prior_type=1)



nbt1 <- 2
betades1 <- list(ll = truebts[[1]],
                 starting_values=c(0, 1, rep(0, nbt1)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt1, byrow = TRUE),
                 bwb=rep(.05, nbt1))

nbt2 <- 2
betades2 <- list(ll = truebts[[2]],
                 starting_values=c(0, 1, rep(0, nbt2)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt2, byrow = TRUE),
                 bwb=rep(.05, nbt2))

nbt3 <- 10
betades3 <- list(llllllllll = truebts[[3]],
                 starting_values=c(0, 1, rep(0, nbt3)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt3, byrow = TRUE),
                 bwb=rep(.05, nbt3))

nbt4 <- 10
betades4 <- list(llllllllll = truebts[[4]],
                 starting_values=c(0, 1, rep(0, nbt4)),
                 bt_prior=matrix(c(0,1), nc=2, nr=nbt4, byrow = TRUE),
                 bwb=rep(.05, nbt4))

betaDesigns <- list(betades1, betades2, betades3, betades4)


res <- simStudyCircGLM(truens = truens, truekps = truekps,
                       betaDesigns = betaDesigns, nsim = nsim,
                       output = "df", seed = 389238, mcmcpar = mcmcpar)


names(res[[1]])
res[[1]] %>%
  tbl_df() %>%
  select(n, kp, kp_mean_bias, kp_mode_bias, kp_in_HDI, b0_bias, bt_1_bias)



############ RUN ONE DATASET
set.seed(1)
d <- generateCircGLMData()

mcmcpar=list(conj_prior = rep(0, 3),
             burnin = 0,
             lag = 1,
             kappaModeEstBandwith=.05,
             CIsize=.95,
             Q=1000000,
             r=2,
             bt_prior_type=1)

res <- do.call(circGLM, c(list(th=d[,1], X=d[,-1],
                               output="list",
                               starting_values=c(0, 1, rep(0, 4)),
                               bt_prior=matrix(c(0,1), nc=2, nr=4, byrow = TRUE),
                               bwb=rep(.05, 4),
                               returnPostSample=TRUE),
                          mcmcpar))

unlist(attributes(d)[-c(2, 3, 6)])
res
pi/2
plot(res)


