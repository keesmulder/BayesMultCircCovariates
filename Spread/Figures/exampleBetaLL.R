rm(list=ls())

source('Data/generateCircularGLMData.R')
source("DataAnalysis/circGLM.R")

# source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

############ RUN ONE DATASET
set.seed(1)

X  <- -3:3
bt <- 0.8
th <- pi/2 + atanlf(bt*X, 2)

# plot(circular(th))

plot(X, th)

betall <- Vectorize(function (b) ll(b0 = pi/2, kp = 1, bt = b,
                                    th = th, X = as.matrix(X), r = 2))

betall(10)
plot(betall, xlim = c(-10, 10))

mcmcpar <- list(conj_prior = rep(0, 3),
                burnin = 0,
                lag = 1,
                kappaModeEstBandwith=.05,
                CIsize=.95,
                Q=3000,
                r=2,
                bt_prior_type=0)
res <- do.call(circGLM, c(list(th=d[,1, drop=FALSE], X=d[,-1, drop=FALSE],
                               output="list",
                               starting_values=c(0, 1, rep(0, nb)),
                               bt_prior=matrix(c(0,1), nc=2, nr=nb, byrow=TRUE),
                               bwb=rep(.05, nb), returnPostSample=TRUE),
                          mcmcpar))

