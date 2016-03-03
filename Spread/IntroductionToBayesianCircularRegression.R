
## ----setup, include=FALSE, cache=FALSE, warning=FALSE--------------------
library(knitr)

# rm(list=ls())
library(Rcpp)
library(circular)
library(mvtnorm)
library(ggplot2)
library(rootSolve)

setwd("C:/Dropbox/Research/BayesMultCircCovariates")
sourceCpp("Code/rvmc.cpp")
source("Code/describeCirc.R")
source("Code/vonMises.R")
source('Code/CircularVonMisesRegressionTest.R')

# set global chunk options
opts_chunk$set(fig.path='figure/', fig.align='center', fig.show='hold', fig.height=10, fig.width=11, warning=FALSE, message=FALSE)
options(formatR.arrow=TRUE,width=90)

# ggplot theme to be used
plotTheme <- theme(
  panel.grid.major = element_line(
    colour = "white",
    linetype = "solid"),
  panel.grid.minor = element_blank(),
  axis.ticks = element_line(
    colour = "black"
  ),
  title = element_text(
    size = rel(1.3),
    face = "bold")
)



## ----plotfunc, include=FALSE---------------------------------------------

par(cex=14)

makeProp <-function(x) (x-min(x))/(max(x)-min(x))

# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)



generateBetaShapePlot <- function(b0_cur = pi, true_b0 = pi, true_bt = 1, true_kp = 20, Xsd=4, n = 100, xl = c(-8, 8), add=1.5) {

  # Data generation
  X       <- matrix(rnorm(n, sd=Xsd))
  err     <- rvmc(n, mu = 0, kp = true_kp)
  th      <- (true_b0 + linkfun(true_bt * X) + err) %% (2*pi)

  par(mfrow=c(2, 2))

  # Show generated data
  xh <- makeProp(X)

  plotCircular(th)
  plot(X, th, main = paste0("True: beta_0=", round(true_b0, 2), ", beta=", round(true_bt, 2), ", kappa=", round(true_kp, 2)), ylim=c(0, 2*pi))

  # Main part of the loglikelihood function.
  Rllfun <- function(bt) sum(cos(th - b0_cur - linkfun(apply(X, 1, "%*%", bt))))

  llmax <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = TRUE)$maximum
  llmin <- optimize(f = function(x) Rllfun(x), interval = xl, maximum = FALSE)$minimum

  b0s <- c(b0_cur, b0_cur, b0_cur, b0_cur)
  bts <- c(llmin, llmax, -5, 5)

  # Plot likelihood of beta at true mean and true kp, and add the extrema.
  plot(Vectorize(function(x) Rllfun(x)), xlim=xl, main=paste("Shape of Log-Likelihood of Beta when Beta_0 =", round(b0_cur, 2)), xlab=expression(beta), ylab="Log-Likelihood")
  abline(v=bts, col=2:(length(bts)+1), lwd=2.5)


  shift <- function(th, shift) ((th + shift) %% (2*pi)) - shift

  dif <- pi - b0_cur

  sth <- shift(th, dif)

  underth <- sth < (add - dif)
  overth  <- sth > (2*pi-(add+dif))

  dth <- c(sth[underth] + (2*pi), sth, sth[overth] - (2*pi))
  dX  <- c(X[underth, ], as.vector(X), X[overth, ])

  # Plot predictions for different beta's
  plot(dX, dth, xlim=c(min(X), max(X)), ylim = b0_cur + c(-(pi+add), pi+add),
       main=paste("Given current beta_0 =", round(b0_cur, 2)),
       ylab="Theta/Predicted Theta (Shifted)", xlab="X")
  abline(h=b0_cur + c(-pi, pi), col="grey60", lty = 3, lwd=3)

  for (i in seq(b0s)) {
    estth <- function(x) shift(b0s[i] + linkfun(bts[i]*x), dif)
    curve(estth, add=TRUE, col=i+1, lwd=2.5)
  }
  abline(h=b0_cur, col="purple", lty = 2)

  par(mfrow=c(1, 1))
}



plot3DLikelihood <- function (n, true_b0, true_bt, true_kp, loglik=TRUE,
                              bt0lim = c(0, 2*pi), bt1lim = c(-10, 10)) {

  # Link functions
  linkfun    <- function(x) 2 * atan(x)
  invlinkfun <- function(x) tan(x/2)

  X       <- cbind(runif(n, -1, 1))
  err     <- rvmc(n, mu = 0, kp = true_kp)
  th      <- (true_b0 + linkfun(apply(X, 1, "%*%", true_bt)) + err) %% (2*pi)

  # Log-likelihood for this dataset and link function.
  llfun <- buildLogLikfun(th, X, linkfun)

  b0s   <- seq(bt0lim[1], bt0lim[2], length = 80)
  bt1   <- seq(bt1lim[1], bt1lim[2], length = 80)
  if(loglik)   f <- Vectorize(function(x, y) llfun(x, true_kp, y))
  if(!loglik)  f <- Vectorize(function(x, y) exp(llfun(x, true_kp, y)))
  z <- outer(b0s, bt1, f)
  op <- par(bg = "white")

  # 3D-plot
  persp(b0s, bt1, z, theta = 50, phi = 40, expand = 0.4, col = "lightblue",
        ltheta = 120, shade = 0.75, ticktype = "detailed",
        xlab = "Beta_0", ylab = "Beta_1", zlab = "Density")

}




## ----Link, echo=FALSE, fig.height=4, fig.width=7, out.width="\\linewidth"----
# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

xl <- c(-20, 20)
yl <- c(-pi, pi)
ggplot(data.frame(x = xl), aes(x)) +
  stat_function(fun = linkfun, size=1) +
  scale_x_continuous(breaks=seq(xl[1], xl[2], length.out = 11)) +
  xlab("x") + ylab(expression(paste(g(x),"=", 2, tan^-1, (x)))) +
  scale_y_continuous(breaks=round(seq(yl[1], yl[2], length.out = 7))) + plotTheme


## ----examplerun, echo=FALSE----------------------------------------------

# Data generation
n       <- 100

true_b0 <- pi
true_bt <- c(3, 6)
true_kp <- 20
X       <- cbind(rnorm(n, sd=.2), rnorm(n, sd=.2))
err     <- rvmc(n, mu = 0, kp = true_kp)
th      <- (true_b0 + linkfun(apply(X, 1, "%*%", true_bt)) + err) %% (2*pi)

th_bar <- as.vector(mean(circular(th)))
K       <- ncol(X)


# Run sampler
se <- mcmcGCM(th = th, X = X, linkfun = linkfun, invlinkfun = invlinkfun,
              b0_start=true_b0, bt_start=c(0, 0))

par(mfrow=c(2, 2))

# # Plot sample results
plot.ts(se$b0 %% (2*pi), main="Beta_0")
plot.ts(se$kp, main="Kappa")
plot.ts(se$beta_1, main="Beta_1")
plot.ts(se$beta_2, main="Beta_2")
par(mfrow=c(1,1))


## ----plot3d, echo=FALSE--------------------------------------------------
# Properties of the data.
n       <-  100
true_b0 <-  pi
true_bt <-  0.2
true_kp <-  20

# Current value of beta_0
b0_cur  <- pi

set.seed(7)
plot3DLikelihood(n, true_b0 = true_b0, true_bt = true_bt, true_kp = true_kp,
                 loglik=FALSE, bt0lim=c(pi-0.5, pi+0.5), bt1lim=c(0.2-0.5, 0.2+0.5))


## ----plot3dlog, echo=FALSE-----------------------------------------------
set.seed(7)
plot3DLikelihood(100, true_b0 = true_b0, true_bt = true_bt, true_kp = true_kp)


## ----plot1, echo=FALSE---------------------------------------------------
set.seed(7)
generateBetaShapePlot(b0_cur = pi,
                      true_b0 = pi,
                      true_bt = 0.2,
                      true_kp = 20)



## ----plot2, echo=FALSE---------------------------------------------------
set.seed(7)
generateBetaShapePlot(b0_cur = true_b0+1, true_b0 = true_b0, true_bt = 0.1, true_kp = true_kp)


## ----plot3, echo=FALSE---------------------------------------------------
generateBetaShapePlot(b0_cur = true_b0+pi, true_b0 = true_b0, true_bt = true_bt, true_kp = true_kp)


## ----plot4, echo=FALSE---------------------------------------------------
generateBetaShapePlot(b0_cur = 3, true_b0 = 3, true_bt = .5, true_kp = 10)


## ----plot5, echo=FALSE---------------------------------------------------
generateBetaShapePlot(b0_cur = 3, true_b0 = 0, true_bt = .1, true_kp = 3)


