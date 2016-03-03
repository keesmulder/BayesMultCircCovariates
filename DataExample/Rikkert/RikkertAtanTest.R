
source("Code/describeCirc.R")
source("DataAnalysis/circGLM.R")

require(dplyr)
require(ggplot2)
require(tidyr)
require(foreign)
require(circular)
require(plotrix)

fileloc <- "C:/Dropbox/Research/BayesMultCircCovariates/DataExample/Rikkert/141106 Data Rikkert.sav"

rdat <- read.spss(fileloc, to.data.frame = TRUE)

ysin <- rdat$controlT
xcos <- rdat$affiliationT

ts <- atan2(ysin, xcos) %% (2*pi)

head(cbind(ts*180/pi, rdat$Person.score))

hist(ts*180/pi - rdat$Person.score)


plot(xcos, ysin, asp=1)
draw.circle(0, 0, 1)

rose.diag(ts, bins=20, prop = 1.5)

plot(circular(ts))

