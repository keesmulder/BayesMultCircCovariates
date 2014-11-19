# rm(list=ls())

library(Rcpp)
library(circular)
library(ggplot2)
library(rgl)

sourceCpp("Code/rvmc.cpp")

# Joint Log-likelihood function, with data already built in.
buildLogLikfun <- function(th, X, linkfun){
  function(b0, kp, bt){
    n <- length(th)
    - n * log(besselI(kp, 0)) +
      kp * sum(cos(th - b0 - linkfun(apply(X, 1, "%*%", bt))))
  }
}

# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

# Data generation
n       <- 100
true_b0 <- 1
true_bt <- c(2, 3)
true_kp <- 2



# plot3DLikelihood <- function (n, true_kp, true_b0, true_bt) {

  # Link functions
  linkfun    <- function(x) 2 * atan(x)
  invlinkfun <- function(x) tan(x/2)

  X       <- cbind(runif(n, -1, 1), runif(n, -1, 1))
  err     <- rvmc(n, mu = 0, kp = true_kp)
  th      <- (true_b0 + linkfun(apply(X, 1, "%*%", true_bt)) + err) %% (2*pi)

  # Log-likelihood for this dataset and link function.
  llfun <- buildLogLikfun(th, X, linkfun)

  b0s   <- seq(0, 2*pi, length = 50)
  bt1   <- seq(-10, 10, length = 50)
  f <- Vectorize(function(x, y) llfun(x, true_kp, c(y, true_bt[2])))
  z <- outer(b0s, bt1, f)
  op <- par(bg = "white")

  # 3D-plot
  persp(b0s, bt1, z, theta = 50, phi = 40, expand = 0.4, col = "lightblue",
        ltheta = 120, shade = 0.75, ticktype = "detailed",
        xlab = "Beta_0", ylab = "Beta", zlab = "Density")

# }

# Interactive 3D-plot
persp3d(b0s, bt1, z, col='lightgreen',
        xlab="Beta_0", ylab="Beta", zlab="Density", smooth=FALSE, bg="grey")
lines3d(b0s, 3, outer(b0s, 3, f)+1, col = "blue", lwd=20)
lines3d(b0s, -10, outer(b0s, -10, f)+1, col = "blue", lwd=20)

lines3d(1, bt1, outer(1, bt1, f)+1, col = "red", lwd=20)
lines3d(4, bt1, outer(4, bt1, f)+1, col = "red", lwd=20)





