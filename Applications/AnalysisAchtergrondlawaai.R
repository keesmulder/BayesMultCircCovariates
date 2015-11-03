require(Rcpp)
require(ggplot2)

source("DataAnalysis/circGLM.R")

wd <- read.table("C:/Dropbox/LiteratureCircular/AppliedDatasets/achtergrondlawaai.txt", header = TRUE)
wd$maand <- as.factor(wd$maand)

wdmaandnum <- as.numeric(wd$maand) - 1

#
# str(wd)
#
# plot(windhoeken ~ windsnelheid, data = wd, col = maand, pch=16)
#
# # Compare months only.
# res2 <- circGLM(wd$windhoeken, wdmaandnum, returnPostSample = TRUE, Q = 10000, burnin = 0)
#
# plot(res2)
#
#
# # Show difference between groups
# ggplot(wd, aes(x = windhoeken, group = maand, col = maand)) + geom_density()
#
#
# # Plot output with predictions.
# plot(wdmaandnum, wd$windhoeken * pi / 180, pch = 16, col = rgb(0, 0, 0, 0.1))
#
# sq <- seq(min(), 1, length.out = 101)
#
# lines(sq, res2$b0_meandir + atanLF(res2$bt_mean * sq, 2))
# lines(sq, res2$b0_meandir + atanLF(res2$bt_mean * sq, 2) + 2*pi)
# abline(h = c(0, 2*pi), col = "gray80")
#
#
#





# ANCOVA
res3 <- circGLM(wd$windhoeken, cbind(wdmaandnum, geluid = wd$achtergrondgeluid),
                returnPostSample = TRUE, Q = 10000, burnin = 0, debug = TRUE, loopDebug = TRUE)

plot(res3)

# Plot output with predictions.
plot(scale(wd$achtergrondgeluid), wd$windhoeken * pi / 180, pch = 16, col = rgb(wdmaandnum, 0, 0, 0.6))

sq <- seq(-3, 3, length.out = 101)

pred.mnd1 <- res3


lines(sq, res3$b0_meandir + atanLF(res3$bt_mean %*% t(cbind(1, sq)), 2), col = "red")
lines(sq, res3$b0_meandir + atanLF(res3$bt_mean %*% t(cbind(0, sq)), 2) + 2*pi, col = "black")

abline(h = c(0, 2*pi), col = "gray80")

res3$bt_mean %*% t(cbind(1, sq))

tail(res3$b0_meandir + atanLF(res3$bt_mean[2] * sq, 2) + 2*pi)
tail(res3$b0_meandir + atanLF(res3$bt_mean %*% t(cbind(0, sq)), 2) + 2*pi)

tail(res3$b0_meandir + atanLF(res3$bt_mean[1] + res3$bt_mean[2] * sq, 2))
tail(res3$b0_meandir + atanLF(res3$bt_mean %*% t(cbind(1, sq)), 2))

res3 <- circGLM(wd$windhoeken, cbind(wdmaandnum, geluid = wd$achtergrondgeluid),
                returnPostSample = TRUE, Q = 10000, burnin = 0, debug = TRUE, loopDebug = TRUE)
predict.plot.circGLM(res3)


# model
n <- 50
x <- rnorm(n)
d <- rep(0:1, n/2)
err <- rnorm(n, 0, .2)
y <- 1.6 + -0.3 * d + invAtanLF(.26 * x, 2) + err

tdt <- data.frame(cbind(y=y, x=x, d=d))

plot(tdt$x, tdt$y, col = as.factor(tdt$d), pch = 16)

# Analysis
r4 <- circGLM(tdt$y, tdt[, -1], returnPostSample = TRUE, Q = 10000, burnin = 0, debug = TRUE, loopDebug = TRUE)

predict.plot.circGLM(r4)


