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
