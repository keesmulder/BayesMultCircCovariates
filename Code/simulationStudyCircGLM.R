
# source("Code/describeCirc.R")
source("Code/generateCircularGLMData.R")
sourceCpp('Code/circGLM.cpp')

library(circular)

truens  <- c(10, 30, 100)
truekps <- c(0.5, 5, 30)
truebts <- list(rep(.1, 4), rep(1, 4), rep(5, 4))



simulationStudyCircGLM <- function(truens, truekps, truebts, predtypes, nsim = 100,
                                   mcmcpar=list(conj_prior = rep(0, 3),
                                                bt_prior = matrix(c(0,1),
                                                                  nc=2, nr=K,
                                                                  byrow = TRUE),
                                                starting_values = c(pi+1,
                                                                    1,
                                                                    rep(0, K)),
                                                burnin = 0,
                                                lag = 1,
                                                bwb=rep(.05, K),
                                                kappaModeEstBandwith=.05,
                                                CIsize=.95,
                                                Q=Q,
                                                c=2,
                                                returnPostSample=TRUE,
                                                bt_prior_type=1)) {

  K       <- length(truebts[[1]])

  designs <- expand.grid(n=truens, kp=truekps, bt=truebts, predtype=predtypes)

#   set.seed(1338)
#   d  <- generateCircGLMData(n = n, nconpred = K, ncatpred = 0,
#                             truebeta = rep(.1, K), residkappa = 2)
#   th <- d[,  1]
#   X  <- d[, -1]
#
#   Q <- 10000
#   set.seed(204)

  for (isim in 1:nsim) {
    curDesign <- designs[i, ]

    readfilename <- paste0(wd, "/Data/Datasets/",
                           "n=",         curDesign$n,
                           "_kp=",       curDesign$kp,
                           "_bt=",       curDesign$bt,
                           "_predtype=", curDesign$predtype,
                           "/nr", i, ".csv")
    th <- read.table(readfilename, sep=",")
    d <- read.table()



  }

  results <- do.call(circGLM, c(list(th, x), mcmcpar))

  m1 <- circGLM(th, X,
                conj_prior           = conj_prior,
                bt_prior             = bt_prior,
                starting_values      = starting_values,
                burnin               = burnin,
                lag                  = lag,
                bwb                  = bwb,
                kappaModeEstBandwith = kappaModeEstBandwith,
                CIsize               = CIsize,
                Q                    = Q,
                r                    = r,
                returnPostSample     = FALSE,
                bt_prior_type        = bt_prior_type)
  class(m1) <- c("circGLM", class(m1))
  m1

}














docalltest <- function(a, b, d=3) {
  cat("a=", a, "\t, b=", b, "\t, d=", d, "\n")
}

l1 <- list(1, 2, 3)
l2 <- list(b=1, d=2, a=3)
l3 <- list(b=1, a=3)
l3 <- list(a=3)

do.call(docalltest, l1)
do.call(docalltest, l2)
do.call(docalltest, l3)








