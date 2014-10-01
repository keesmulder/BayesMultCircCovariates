require(ggplot2)
rm(list=ls())

set.seed(200)
n <- 1000
genx <- 4 + rnorm(n)
gend <- sample(rep(c(0, 1), n/2), n)
tBx  <- 5
tInt <- 100
tDif <- 15
dat <- data.frame(y = 5 + tBx * genx + tInt * genx * gend + tDif * gend + rnorm(n, sd = 8), 
                  x = genx, 
                  d = factor(gend))
dat

p <- ggplot(aes(x=x, y=y, group=d, color=d), data = dat) + geom_point()
p



cffull   <- lm(y ~ d + x, data=dat)$coefficients
cfglobal <- lm(y ~ x, data=dat)$coefficients
cfd0     <- lm(y ~ x, data=dat[dat$d==0, ])$coefficients
cfd1     <- lm(y ~ x, data=dat[dat$d==1, ])$coefficients

p + 
  geom_abline(intercept=cffull[1], slope=cffull[3], col='red') + 
  geom_abline(intercept=cffull[1]+cffull[2], slope=cffull[3], col='blue') +
  geom_abline(intercept=cfglobal[1], slope=cfglobal[2], col='black') +
  geom_abline(intercept=cfd0[1], slope=cfd0[2], col='pink') +
  geom_abline(intercept=cfd1[1], slope=cfd1[2], col='lightblue') 

# The regression lines in the full model are clearly not parallel to any of the 
# others, and so the ANCOVA is neither a change on the global or separate
# models.

cfavg <- (cfd0 + cfd1)/2

p + 
  geom_abline(intercept=cffull[1], slope=  cffull[3], col='red') + 
  geom_abline(intercept=cffull[1]+cffull[2], slope=cffull[3], col='blue') +
  geom_abline(intercept=cfavg[1], slope=cfavg[2], col='orange')


cffull


x_bar_j <- as.vector(by(dat$x, dat$d, mean))
y_bar_j <- as.vector(by(dat$y, dat$d, mean))

beta_try1 <- sum(outer(dat$x, x_bar_j, '-') * outer(dat$y, y_bar_j, '-')) / 
  sum(outer(dat$x, x_bar_j, '-')^2)

beta_try2 <- sum(outer(dat$x, x_bar_j, '-') * outer(dat$y, y_bar_j, '-')) / 
  sum(outer(dat$x, x_bar_j, '-')^2)

cffull
cfglobal
beta_try
cfavg
