source("Code/describeCirc.R")
source("DataAnalysis/circGLM.R")

require(dplyr)
require(ggplot2)
require(tidyr)
require(foreign)


fileloc <- "DataExample/ESS7e01.sav"
pd      <- read.spss(fileloc, to.data.frame = TRUE)

labs <- attr(pd, "variable.labels")

grep(labs, "End")
grep("pvq", labs, ignore.case = TRUE, value = TRUE)


