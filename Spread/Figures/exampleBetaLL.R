rm(list=ls())

source('Data/generateCircularGLMData.R')
source("DataAnalysis/circGLM.R")

require(ggplot2)

# ggplot theme to be used
plotTheme <- theme(
  panel.background = element_rect(
    fill = "white"
    ),
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

getdat <- function(fun, res = 100, xl = c(-10, 10), ...) {
  sq  <- seq(xl[1], xl[2], length.out = res)
  ysq <- fun(b=sq)
  data.frame(sq=sq, ysq=ysq)
}

plotbeta <- function(th, X, normalPrior=FALSE, res = 100, xl = c(-10, 10),
                     b0 = pi/2, kp = 1, bt = b, r = 2, mu=0, sd=10) {

  if (normalPrior) {
    betall <- Vectorize(function (b) {
      ll(b0 = b0, kp = kp, bt = b, th = th, X = X, r = r)
    })
  } else {
    betall <- Vectorize(function (b) {
      ll(b0 = b0, kp = kp, bt = b, th = th, X = X, r = r) *
        logProbNormal(b, mu, sd)
    })
  }

  d <- getdat(fun=betall, res=res, xl=xl)

  ggplot(aes(x=sq, y=ysq), data=d) +
    geom_line() +
    geom_vline(xintercept=d$sq[which.max(d$ysq)], linetype = "dashed") + theme_bw()

}


#
# X  <- as.matrix(-3:3)
# bt <- 0.8
# th <- pi/2 + atanlf(bt*X, 2)
#
# plotbeta(normalPrior=FALSE, res=100, xl=c(-50, 50), b0=pi/2, kp=1, r=2, th=th, X=X)
#


