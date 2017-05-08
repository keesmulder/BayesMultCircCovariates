# Link functions
linkfun    <- function(x) 2 * atan(x)
invlinkfun <- function(x) tan(x/2)

plot(dnorm, xlim=c(-4, 4))
plot(pnorm, xlim=c(-4, 4))
plot(invlinkfun, xlim=c(-4, 4))

x <- seq(-4, 4, .1)
px <- dnorm(x)

plot(x, px)

th  <- seq(-3, 3, .1)
xth <- invlinkfun(th)

pth <- dnorm(xth)

plot(th, pth)

# Middle band of size pi on the circle.
testbt <- function (bt = 1) {

  LBth <- -pi/2
  UBth <- pi/2

  LBx <- invlinkfun(LBth)
  UBx <- invlinkfun(UBth)

  LBbtx <- invlinkfun(LBth)/bt
  UBbtx <- invlinkfun(UBth)/bt

  return(abs(pnorm(UBbtx) - pnorm(LBbtx)))
}

plot(Vectorize(testbt), xlim=c(-2, 2), xlab="Beta")

abline(h=.5, col="red")

lower <- uniroot(function(x) testbt(x) - .5, interval = c(-2, 0))
upper <- uniroot(function(x) testbt(x) - .5, interval = c(0, 2))

abline(v=lower$root)
abline(v=upper$root)

abs(lower$root)


# Multivariate

require(mvtnorm)

# The same, so we can use pnorm(sum_of_xs).
pmvnorm(lower = c(-Inf, -Inf), upper = c(Inf, 2))
pnorm(2)





# The same, but simpler.

fun <- function(sz) pnorm(1 / sz) - pnorm(-1 / sz) - .5
uniroot(fun, interval = c(0, 10))

plot(fun, xlim=c(0,10))
abline(h=0)







# Testing if u = .5 with sigma^2_Z = 1.483, not only for kappa = 0 and kappa =
# inf, but for an arbitrary kappa.



getPropCircGLM <- function(kp, bt = 1.483, b0 = 0, n=10000, scaled = FALSE) {

  # Data generation
  X       <- matrix(rnorm(n, sd=1))
  if (scaled) X <- scale(X)
  err     <- rvmc(n, mu = 0, kp = kp)
  th      <- (b0 + linkfun(bt * X) + err) %% (2*pi)

  sum(c(th > (2*pi - pi/2), th < (pi/2)))/n
}

# getPropCircGLM(10, n = 100)

# plot(Vectorize(getPropCircGLM), xlim=c(0, 10))
plot(Vectorize(function(kp) getPropCircGLM(kp, scaled  = TRUE, bt = 1.5, n=100000)), xlim=c(0, 10))

















