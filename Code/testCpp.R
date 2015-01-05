
source("Data/generateCircularGLMData.R")
source("Code/describeCirc.R")
# source("Code/CircularVonMisesRegressionTest.R")
sourceCpp('DataAnalysis/circGLM.cpp')

library(circular)


n <- 30
K <- 2
set.seed(1335)
d  <- generateCircGLMData(n = n, nconpred = K, ncatpred = 0, truebeta = rep(.1, K),
                          residkappa=4)
th <- d[,  1]
X  <- d[, -1]



Q <- 10000
set.seed(202)
m1 <- circGLM(th, X,
              conj_prior = rep(0, 3),
              bt_prior = matrix(c(0,1), nc=2, nr=K, byrow = TRUE),
              starting_values = c(pi+1, 1, rep(0, K)),
              burnin = 0,
              lag = 1,
              bwb=rep(.05, K),
              kappaModeEstBandwith=.05,
              CIsize=.95,
              Q=Q,
              r=2,
              returnPostSample=FALSE,
              bt_prior_type=1,
              output="vector")

m1
par(mfrow=c(1,1))

plot(m1)
hist(m1$b0_chain)
abline(v=circularQuantile(m1$b0_chain, c(0.025, .975)))

p <- m1$ProportionAccepted


# plot.new()

# # plot.new()
# par(mfrow = c(2,3))
# hist(m1$b0_chain, main="beta_0")
# hist(m1$kp_chain, main="kappa")
# apply(m1$bt_chain, 2, hist, main="beta")

# acf(m1$b0_chain)

