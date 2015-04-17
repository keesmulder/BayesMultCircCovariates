source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

############ RUN ONE DATASET
set.seed(1)
ncon <- 1
ncat <- 0
nb <- ncon + ncat
n <- 30
d <- generateCircGLMData(n = n, residkappa = 30,
                         nconpred = ncon, ncatpred = ncat,
                         truebeta=seq(-10, 1, length.out = nb))

mcmcpar <- list(conj_prior = rep(0, 3),
             burnin = 0,
             lag = 1,
             kappaModeEstBandwith=.05,
             CIsize=.95,
             Q=30000,
             r=2,
             bt_prior_type=0)
res <- do.call(circGLM, c(list(th=d[,1, drop=FALSE], X=d[,-1, drop=FALSE],
                               output="list",
                               starting_values=c(0, 1, rep(0, nb)),
                               bt_prior=matrix(c(0,1), nc=2, nr=nb, byrow=TRUE),
                               bwb=rep(.05, nb),
                               returnPostSample=TRUE, reparametrize=TRUE),
                          mcmcpar))

# After a while, the chains run amok. This is not good.
plot(res, coef="Beta")
plot(res, coef="Zeta")




# The problem is solved by using the normal prior.
mcmcpar$bt_prior_type <- 1
res <- do.call(circGLM, c(list(th=d[,1, drop=FALSE], X=d[,-1, drop=FALSE],
                               output="list",
                               starting_values=c(0, 1, rep(0, nb)),
                               bt_prior=matrix(c(0,1), nc=2, nr=nb, byrow=TRUE),
                               bwb=rep(.05, nb), returnPostSample=TRUE),
                          mcmcpar))
plot(res)



