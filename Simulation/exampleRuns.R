source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')

library(dplyr)

############ RUN ONE DATASET
set.seed(2)
ncon <- 1
ncat <- 0
nb <- ncon + ncat
n <- 100
d <- generateCircGLMData(n = n, residkappa = 30,
                         nconpred = ncon, ncatpred = ncat,
                         truebeta=seq(0.8, 1, length.out = nb))

plot(d[,2], d[,1])

mcmcpar <- list(conj_prior      = rep(0, 3),
                burnin          = 2000,
                lag             = 1,
                CIsize          = .95,
                Q               = 100000,
                r               = 2,
                bt_prior_type   = 0,
                starting_values = c(pi/2, 1, rep(-3, nb)),
                bt_prior        = matrix(c(0,1), ncol=2, nrow=nb, byrow=TRUE),
                bwb             = rep(.05, nb),
                reparametrize   = TRUE,
                kappaModeEstBandwith = .05)

functionpar <- list(th     = d[,1, drop=FALSE],
                    X      = d[,-1, drop=FALSE],
                    output = "list",
                    returnPostSample = TRUE)

set.seed(2)
res1 <- do.call(circGLM, c(functionpar,
                           mcmcpar))



# RUN 2
mcmcpar2     <- mcmcpar
functionpar2 <- functionpar
mcmcpar2[["reparametrize"]] <- FALSE

set.seed(2)
res2 <- do.call(circGLM, c(functionpar2,
                           mcmcpar2))

plot(res1, coef="Beta")
plot(res1, coef="Zeta")
plot(res2, coef = "Zeta")


hist(res1$bt_chain, breaks=100)
hist(res2$bt_chain, breaks=100)
hist(res2$zt_chain, breaks=100)


truezeta <- as.numeric(invAtanLF(attr(d, "truebeta"), pi/2))
attr(d, "truezeta")
hist(res1$bt_chain, breaks=100)
abline(v=attr(d, "truebeta"), col="green")

hist(res1$zt_chain, breaks=100)
abline(v=truezeta, col="green")



