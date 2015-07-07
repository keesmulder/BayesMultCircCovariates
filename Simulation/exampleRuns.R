source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

############ RUN ONE DATASET
set.seed(2)
ncon <- 0
ncat <- 1
nb <- ncon + ncat
n <- 100
d <- generateCircGLMData(n = n, residkappa = 30,
                         nconpred = ncon, ncatpred = ncat,
                         truebeta=seq(0.8, 1, length.out = nb))

# plot(d[,2], d[,1])

mcmcpar <- list(conj_prior      = rep(0, 3),
                burnin          = 0,
                lag             = 1,
                CIsize          = .95,
                Q               = 100000,
                r               = 2,
                bt_prior_type   = 0,
                starting_values = c(pi/2, 1, rep(-10, nb)),
                bt_prior        = matrix(c(0,1), ncol=2, nrow=nb, byrow=TRUE),
                bwb             = rep(.05, nb),
                reparametrize   = TRUE,
                kappaModeEstBandwith = .05)

functionpar <- list(th     = d[,1,  drop=FALSE],
                    X      = d[,-1, drop=FALSE],
                    output = "list",
                    returnPostSample = TRUE)

set.seed(2)
tpar <- c(functionpar,mcmcpar)
res1 <- do.call(circGLM, c(functionpar, mcmcpar))



# RUN 2
mcmcpar2     <- mcmcpar
functionpar2 <- functionpar
mcmcpar2[["Q"]] <- 100
functionpar2[["output"]] <- "list"
functionpar2[["returnPostSample"]] <- TRUE

set.seed(2)
res2 <- do.call(circGLM, c(functionpar2,
                           mcmcpar2))
res2


# RUN 3
mcmcpar3     <- mcmcpar
functionpar3 <- functionpar
mcmcpar3[["bwb"]] <- rep(.5, nb)

set.seed(2)
res3 <- do.call(circGLM, c(functionpar3,
                           mcmcpar3))


# RUN 4
mcmcpar4     <- mcmcpar
functionpar4 <- functionpar
mcmcpar4[["bwb"]] <- rep(1, nb)

set.seed(2)
res4 <- do.call(circGLM, c(functionpar4,
                           mcmcpar4))


str(res1)

#
#
plot(res1, coef="Zeta")
plot(res2, coef = "Zeta")
plot(res3, coef = "Zeta")
plot(res4, coef = "Zeta")
#
# plot(res1, coef="Beta")
# plot(res2, coef = "Beta")
# plot(res3, coef = "Beta")
# plot(res4, coef = "Beta")
#
par(mfrow=c(2,2))
truezeta <- as.numeric(atanLF(attr(d, "truebeta"), 2/pi))
c(truezeta, attr(d, "truezeta"))
hist(res1$zt_chain, breaks=2000, xlim = c(0.4, 0.46))
abline(v=truezeta, col="green")
hist(res2$zt_chain, breaks=100)
abline(v=truezeta, col="green")
hist(res3$zt_chain, breaks=100)
abline(v=truezeta, col="green")
hist(res4$zt_chain, breaks=100)
abline(v=truezeta, col="green")



par(mfrow=c(2,2))
truezeta <- as.numeric(atanLF(attr(d, "truebeta"), 2/pi))
c(truezeta, attr(d, "truezeta"))
plot.ts(res1$zt_chain)
abline(h=truezeta, col="green")
plot.ts(res2$zt_chain)
abline(h=truezeta, col="green")
plot.ts(res3$zt_chain)
abline(h=truezeta, col="green")
plot.ts(res4$zt_chain)
abline(h=truezeta, col="green")



par(mfrow=c(1,1))
#

# hist(res1$bt_chain, breaks=100)
# abline(v=attr(d, "truebeta"), col="green")
#
# hist(res1$zt_chain, breaks=100)
# abline(v=truezeta, col="green")
#
#
#






functionpar <- list(th     = d[,1, drop=FALSE],
                    X      = d[,-1, drop=FALSE],
                    output = "vector",
                    returnPostSample = FALSE)

set.seed(2)
curSimRes <- do.call(circGLM, c(functionpar,
                           mcmcpar))
curSimRes





