## ----init, include = FALSE-----------------------------------------------
library(knitr)
library(ggthemes)
library(dplyr)
library(xtable)
library(grid)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(circular)
library(foreign)
library(extrafont)
library(tikzDevice)

# Set root directory
opts_knit$set(root.dir = "C:/Dropbox/Research/BayesMultCircCovariates")
setwd(opts_knit$get('root.dir'))

source("Spread/Figures/plotBetaLL.R")
source("Spread/Figures/obtainBayesFactorsCircGLM.R")
source("Code/describeCirc.R")
source("DataAnalysis/circGLM.R")
source("Data/generateCircularGLMData.R")
source("Simulation/simulationStudyCircGLM.R")
source("Code/describeCirc.R")
source('Code/vonMises.R')

scaleFac <- .7
big.fig.height <- 5*scaleFac
big.fig.width  <- 9*scaleFac

small.fig.height <- 5*scaleFac
small.fig.width  <- 7*scaleFac


## General options
options(scipen = 10)

## Table options
options("xtable.booktabs" = TRUE,
        "xtable.sanitize.text.function" = function(x){x},
        "xtable.include.rownames" = FALSE,
        "xtable.table.placement" = "btp",
        "xtable.caption.placement" = "top")



## Knitr Options
opts_chunk$set(echo=FALSE,
               dev='tikz',
               external=TRUE,
               cache=FALSE,
               prompt=FALSE,
               tidy=TRUE,
               comment=NA,
               message=FALSE,
               fig.height = small.fig.height,
               fig.width = small.fig.width,
               warning=FALSE)




getPMP <- function(x) cbind(x/(1+x), 1/(1+x))

# Check for repetitions.
areRepeated <- function(x) {
  xa <- x[-length(x)]
  xb <- x[-1]
  c(FALSE, xa == xb)
}

# Makes any repeated number NA.
makeRepeatedNA <- function(x) {
  x[areRepeated(x)] <- NA
  x
}

# Default ggplot2 theme
myTheme <- theme_bw()  +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.margin = unit(0, "lines"))

## ----dataForParallelNonparallel------------------------------------------
n <- 100
truedelta <- 2; truebeta <- .4; truekappa <- 20
dat <- as.data.frame(generateCircGLMData(n = n, nconpred = 1, ncatpred = 1, truebeta0 = pi/2,
                                         truedelta = truedelta, truebeta = truebeta,
                                         residkappa = truekappa))

## ----nonparallel, fig.width=small.fig.width, fig.height=small.fig.height----
res1 <- circGLM(th = dat$th, X = dat[,-1], returnPostSample = FALSE, skipDichSplit = TRUE)
plot.predict.circGLM(res1, groupingInBeta = TRUE,
                     xlab = "x", ylab = "$\\theta$")

## ----parallel, fig.width=small.fig.width, fig.height=small.fig.height----
res2 <- circGLM(th = dat$th, X = dat[,-1], returnPostSample = FALSE)
plot.predict.circGLM(res2,
                     xlab = "x", ylab = "$\\theta$")

## ----dataForLikelihoodPriorPicture---------------------------------------
n <- 7
X  <- scale(as.matrix(seq(-3, 3, length.out = n)))
bt <- 1
th <- pi + atanLF(bt*X, 2) + rnorm(n, 0, .1)
d <- data.frame(th=th, X=X)
res <- 400
b0cur <- 0
xl <- 40 * c(-1, 1)

## ----BetaConstantPrior---------------------------------------------------
plotbeta(normalPrior=FALSE, res=res, xl=xl, b0=0, kp=1, r=2, th=th, X=X) +
  myTheme + labs(x = "$\\beta$", y = "Conditional log-posterior of $\\beta$")


## ----BetaNormalPrior-----------------------------------------------------
plotbeta(normalPrior=TRUE, res=res, xl=xl, b0=0, kp=1, r=2, th=th, X=X) +
  myTheme + labs(x = "$\\beta$", y = "Conditional log-posterior of $\\beta$")

## ----SimTableOneLinPred, results="asis"----------------------------------
load("SimulationResults/[simStudCircGLM_nsim5000_Q20000_burnin1000_r2_seed38944n20,100_kp2,20_btl,0.05,0.8].rda")
SSRregression <- simStudyResults

SSRregressionfull <- rbind(SSRregression[[1]], SSRregression[[2]])


nsim <- dim(attr(SSRregression[[1]], "full"))[1]


chosenColsOneLinPred      <- c("bt_true", "kp", "n",
                               "b0_bias", "b0_in_CCI",
                               "kp_mode", "kp_in_HDI",
                               "bt_1_mean", "bt_1_in_CCI", "bt_1_propacc",
                               "TimeTaken4")

colNamesOneLinPred      <- c("$\\beta$", "$\\kappa$", "n",
                             "Bias ", "Cov. ",
                             "$\\hat{\\kappa}$", "Cov. ",
                             "$\\hat{\\beta}_1$", "Cov. ", "Acc.",
                             "MCT")

SSRregressionfull %>%
  mutate(pred = sapply(SSRregressionfull$pred, first)) %>%
  mutate(kp = makeRepeatedNA(kp), bt_true = makeRepeatedNA(pred)) %>%
  select(one_of(chosenColsOneLinPred)) %>%
  xtable(caption = "Results of the simulation study for the simple regression scenario. 'Cov.' denotes the 95\\% coverage for a specific parameter, while 'Acc.' denotes the acceptance probability. MCT denotes the mean computation time in seconds.",
         label = "tableOneLinearPredictor",
         digits = c(0, 2, 0, 0, rep(2, length(chosenColsOneLinPred)-3))) -> tableOneLinearPredictor

colnames(tableOneLinearPredictor) <- colNamesOneLinPred

topText  <- "\\toprule \\multicolumn{3}{c}{True} & \\multicolumn{2}{c}{$\\beta_0$} & \\multicolumn{2}{c}{$\\kappa$} & \\multicolumn{3}{c}{$\\beta_1$} & \\\\ "
topRules <- "\\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-10}"
topRow   <- paste(topText, topRules)
rowsToAddSpace <- which(areRepeated(SSRregressionfull$kp)) - 1
rowsToAddSpace <- rowsToAddSpace[-length(rowsToAddSpace)]

print(tableOneLinearPredictor,
      hline.after = c(0, nrow(tableOneLinearPredictor)),
      add.to.row = list(pos = list(-1, rowsToAddSpace), command = c(topRow, "\\vspace{0.2cm} ")),
      latex.environments = c("center", "small")
)

## ----tableFacANOVA, results="asis"---------------------------------------
load("SimulationResults/[simStudCircGLM_nsim5000_Q20000_burnin1000_r2_seed38944n20,100_kp2,20_btcc,0.05,0.8].rda")
SSRanova <- simStudyResults

SSRanovafull <- rbind(SSRanova[[1]], SSRanova[[2]])


chosenColsFacANOVA  <- c("bt_true","kp", "n",
                         "b0_bias", "b0_in_CCI",
                         "kp_mode", "kp_in_HDI",
                         "dt_1_mdir", "dt_1_in_CCI", "dt_1_propacc",
                         "TimeTaken4")

colNamesFacANOVA      <- c("$\\delta$", "$\\kappa$", "n",
                           "Bias", "Cov.",
                           "$\\hat{\\kappa}$", "Cov.",
                           "$\\hat{\\delta}_1$", "Cov.", "Acc.",
                           "MCT")

SSRanovafull %>%
  mutate(pred = sapply(SSRanovafull$pred, first)) %>%
  mutate(kp = makeRepeatedNA(kp), bt_true = makeRepeatedNA(pred)) %>%
  select(one_of(chosenColsFacANOVA)) %>%
  xtable(caption = "Results of the simulation study for the factorial ANOVA scenario. 'Cov.' denotes the 95\\% coverage for a specific parameter, while 'Acc.' denotes the acceptance probability. MCT denotes the mean computation time in seconds.",
         label = "tableFacANOVA",
         digits = c(0, 2, 0, 0, rep(2, length(chosenColsFacANOVA)-3))) -> tableFacANOVA

colnames(tableFacANOVA) <- colNamesFacANOVA

topText  <- "\\toprule \\multicolumn{3}{c}{True} & \\multicolumn{2}{c}{$\\beta_0$} & \\multicolumn{2}{c}{$\\kappa$} & \\multicolumn{3}{c}{$\\delta_1$} & \\\\ "
topRules <- "\\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-10}"
topRow   <- paste(topText, topRules)
rowsToAddSpace <- which(areRepeated(SSRanovafull$kp)) - 1
rowsToAddSpace <- rowsToAddSpace[-length(rowsToAddSpace)]

print(tableFacANOVA,
      hline.after = c(0, nrow(tableFacANOVA)),
      add.to.row = list(pos = list(-1, rowsToAddSpace), command = c(topRow, "\\vspace{0.2cm} ")),
      latex.environments = c("center", "small")
)

## ----tableANCOVA, results="asis"-----------------------------------------

load("SimulationResults/[simStudCircGLM_nsim5000_Q20000_burnin1000_r2_seed38944n20,100_kp2,20_btcllll,0.05,0.8].rda")
SSRancova <- simStudyResults

SSRancovafull <- cbind(rbind(SSRancova[[1]], SSRancova[[2]]))



chosenColsANCOVA  <- c("bt_true","kp", "n",
                       "b0_bias", "b0_in_CCI",
                       "kp_mode", "kp_in_HDI",
                       "dt_1_mdir", "dt_1_in_CCI", "dt_1_propacc",
                       "bt_1_mean", "bt_1_in_CCI", "bt_1_propacc",
                       "TimeTaken4")

colNamesANCOVA     <- c("$\\beta , \\delta$", "$\\kappa$", "n",
                        "Bias", "Cov.",
                        "$\\hat{\\kappa}$", "Cov.",
                        "$\\hat{\\delta}_1$", "Cov.", "Acc.",
                        "$\\hat{\\beta}_1$", "Cov.", "Acc.",
                        "MCT")



SSRancovafull %>%
  mutate(pred = sapply(SSRancovafull$pred, first)) %>%
  mutate(kp = makeRepeatedNA(kp), bt_true = makeRepeatedNA(pred)) %>%
  select(one_of(chosenColsANCOVA)) %>%
  xtable(caption = "Results of the simulation study for the ANCOVA scenario. 'Cov.' denotes the 95\\% coverage for a specific parameter, while 'Acc.' denotes the acceptance probability. MCT denotes the mean computation time in seconds.",
         label = "tableANCOVA",
         digits = c(0, 2, 0, 0, rep(2, length(chosenColsANCOVA)-3))) -> tableANCOVA

colnames(tableANCOVA) <- colNamesANCOVA


topText  <- "\\toprule \\multicolumn{3}{c}{True} & \\multicolumn{2}{c}{$\\beta_0$} & \\multicolumn{2}{c}{$\\kappa$} & \\multicolumn{3}{c}{$\\delta_1$} & \\multicolumn{3}{c}{$\\beta_1$} \\\\ "
topRules <- "\\cmidrule(lr){1-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}"
topRow   <- paste(topText, topRules)
rowsToAddSpace <- which(areRepeated(SSRancovafull$kp)) - 1
rowsToAddSpace <- rowsToAddSpace[-length(rowsToAddSpace)]

print(tableANCOVA,
      hline.after = c(0, nrow(tableANCOVA)),
      add.to.row = list(pos = list(-1, rowsToAddSpace), command = c(topRow, "\\vspace{0.2cm} ")),
      latex.environments = c("center", "small")
)

## ----SimBFPrep-----------------------------------------------------------

SSRlist <- list(Regression = SSRregression, ANOVA = SSRanova, ANCOVA = SSRancova)

BFDf <- obtainBFResults(SSRlist)

BFDfStack <- rbind(cbind(BFDf, type = "Inequality", BF = BFDf$Ineq), cbind(BFDf, type = "Equality", BF = BFDf$Eq))
BFDfStack$pmp <- getPMP(BFDfStack$BF)[,1]

BFDfStackPred.05  <- filter(BFDfStack, Coefficient == .05)
BFDfStackPred.80 <- filter(BFDfStack, Coefficient == .80)

# ggplot(data = BFDfStack, aes(x = pmp)) +
#   facet_grid(type + n ~ Model + Kappa) +
#   geom_histogram(binwidth = .05) +
#   myTheme


# BFSSR <- SSRanova
# BFn  <- "100"
# BFkp <- "20"
# BFpdnum <- 2 # Number of the chosen predictor value 1, = .05, 2 = .8
# BFpd <- BFSSR[[BFpdnum]]$pred[[1]][1]
#
# BFDf <- data.frame(attr(BFSSR[[BFpdnum]], "full")[, , BFn, BFkp])
#
# BFSDD  <- BFDf$MuSDDBayesFactors1
# BFIneq <- BFDf$MuIneqBayesFactors1
#
# pCorrectSDD <- mean(BFSDD<1)
# pCorrectIneq <- mean(BFIneq>1)


## ----BFSDDExample, fig.height=big.fig.height, fig.width=big.fig.width----
ggplot(data = BFDfStackPred.05, aes(x = Kappa, y = pmp, fill = n)) +
 facet_grid(type ~ Model) +
 xlab("$\\kappa$") + ylab("Probability of the correct hypothesis") +
 geom_boxplot(outlier.color = rgb(0, 0, 0, 0), width = .5) +
 scale_fill_grey(start = .50,  end = .98) +
 myTheme

## ----ExamplePrep---------------------------------------------------------
# Load in data
fileloc <- "DataExample/DeafPostma/Data_BARS_BREED_N48_23052011.sav"
pd      <- read.spss(fileloc, to.data.frame = TRUE)

# Remove person with low raven percentile
tbl_df(pd) %>% filter(Ravenpercentile > 30) %>% rename(Group = CODE) -> dpd

levels(dpd$Group) <- c("Deaf", "Interpreter", "Control")


## ----ExampleDescrTable, results='asis'-----------------------------------
# Recreate table 1.
dpd %>%
  group_by(Group) %>%
  summarise("Age (Mean)"  = mean(age), "Age (SD)" = sd(age),
            # "Age (Min) = min(age), "Age (Max)" = max(age),
            "Educ Mean" = mean(toteducation), "Educ SD" = sd(toteducation),
            # "Raven percentile" = mean(Ravenpercentile),
            "Deviation" = meanDir(imdev_signed*pi/180)*180/pi) -> tab1
colnames(tab1)[1:5] <- c("", "Mean", "SD", "Mean", "SD")
tab1[,1] <- c("Deaf", "Interpreter", "Control")
xtab1 <- xtable(tab1, auto = TRUE,
                caption = "Summary statistics of the mean age (years), mean education (years) and mean direction of deviation (degrees).",
                label = "ExampleDescrTable",
                digits = c(0, 0, 2, 2, 2, 2, 2))


topText  <- "\\toprule & \\multicolumn{2}{c}{Age} & \\multicolumn{2}{c}{Education} &  \\\\ "
topRules <- "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} "
topRow   <- paste(topText, topRules)

print(xtab1,
      # include.rownames = TRUE,
      hline.after = c(0, nrow(xtab1)),
      add.to.row = list(pos = list(-1), command = topRow),
      latex.environments = c("center")
      )

## ----BaseANOVAModel------------------------------------------------------
# Run the basic GLM model with just groups as predictors.
# Create dummies
dpd %>%
  select(Group, imdev_signed) %>%
  mutate(Deviation    = imdev_signed * pi / 180,
         Deaf         = as.numeric(Group == "Deaf"),
         Interpreter  = as.numeric(Group == "Interpreter")) -> ANOVApd
anovaModel <- circGLM(th = ANOVApd$Deviation, X = ANOVApd[, c("Deaf", "Interpreter")], Q = 100000, burnin = 1000,
                      returnPostSample = TRUE, debug = FALSE, loopDebug = FALSE)

anovaDeltaTab <- cbind(t(anovaModel$dt_meandir), t(anovaModel$dt_CCI))
rownames(anovaDeltaTab) <- paste0("$\\delta_", c("{df}", "{in}"), "$")

anovaBF1 <- anovaModel$MuBayesFactors[,1]
anovaBF2 <- anovaModel$MuBayesFactors[,2]
anovaBFTab <- cbind(anovaBF1, getPMP(anovaBF1), anovaBF2, getPMP(anovaBF2))

# # Control, deaf, interpreter
#
# # HOW THE MU CHAINS LOOK
# mc <- anovaModel$mu_chain
# mcd <- data.frame(mc)
# colnames(mcd) <- c("c", "d", "i")
# head(mcd)
# gmcd <- gather(mcd, group, value)
# ggplot(gmcd, aes(x = value, col = group)) + geom_density() + geom_hline(aes(yintercept = 1/(2*pi)))
#
# # MEAN COMPARISON BAYES FACTORS
# dc <- anovaModel$mu_chain - cbind(anovaModel$mu_chain[, 2:3], anovaModel$mu_chain[, 1])
# dcd <- data.frame(dc)
# colnames(dcd) <- c("c - d", "d - i", "i - c")
# head(dcd)
# gdcd <- gather(dcd, group, value)
# ggplot(gdcd, aes(x = value, col = group)) + geom_density() + geom_hline(aes(yintercept = 1/(2*pi))) + geom_vline(aes(xintercept = 0), col = "skyblue")

dimnames(anovaBFTab) <- list(c("$(\\mu_{cn}, ~\\mu_{df})$",
                               "$(\\mu_{cn}, ~\\mu_{in})$",
                               "$(\\mu_{df}, ~\\mu_{in})$"),
                             c("$BF$",
                               "$p(\\mu_a > \\mu_b)$",
                               "$p(\\mu_a < \\mu_b)$",
                               "$BF$",
                               "$p(\\mu_a = \\mu_b)$",
                               "$p(\\mu_a \\neq \\mu_b)$"))

anovaTab <- with(anovaModel, rbind("$\\beta_0$" = c(b0_meandir, b0_CCI),
                                   "$\\kappa$" = c(kp_mode, kp_HDI),
                                   anovaDeltaTab
                                   # "DIC" = c(DIC, NA, NA),
                                   # "WAIC" = c(WAIC1, NA, NA))
                 ))

colnames(anovaTab) <- c("Estimate", "LB", "UB")

## ----ANOVATable, results='asis'------------------------------------------
xANOVATab <- xtable(anovaTab,
                    caption = "Results for the ANOVA model. LB and UB respectively represent lower and upper bound of the 95\\% credible interval of the given parameter.",
                    label = "ANOVATable")

print(xANOVATab,
      include.rownames = TRUE,
      latex.environments = c("center"))


## ----ANOVAConvergencePlot, fig.width=big.fig.width, fig.height=big.fig.height----
plot(anovaModel, type = "stack", ggTheme = myTheme, labelFormat = "latex") + theme(panel.margin = unit(0.5, "lines"))

## ----ANOVABFTable, results='asis'----------------------------------------
xANOVABFTab <- xtable(anovaBFTab,
                      caption = "Bayes factors and posterior model probabilities for the ANOVA model. Posterior model probilities are given for equal prior odds.",
                      label = "ANOVABFTable")

topText  <- "\\toprule & \\multicolumn{3}{c}{Ineq. $~ \\mu_a > \\mu_b ~ : ~ \\mu_a < \\mu_b$} & \\multicolumn{3}{c}{Eq. $~\\mu_a = \\mu_b ~ : ~ \\mu_a \\neq \\mu_b$}  \\\\ "
topRules <- "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} "
topRow   <- paste(topText, topRules)



# print(xANOVABFTab,
#       hline.after = c(0, nrow(xANOVABFTab)),
#       add.to.row = list(pos = list(-1), command = topRow),
#       include.rownames = TRUE,
#       latex.environments = c("center", "small")
#       )

## ----ExampleCovModelRun--------------------------------------------------
dpd %>%
  select(Group, imdev_signed, age, handness) %>%
  mutate(Deviation    = imdev_signed * pi / 180,
         Deaf         = as.numeric(Group == "Deaf"),
         Interpreter  = as.numeric(Group == "Interpreter")) -> ANCOVApd

# Run a model with covariates.
ancovaModel <- circGLM(th = ANCOVApd$imdev_signed,
                       X = ANCOVApd[, c("age", "handness", "Deaf", "Interpreter")],
                       Q = 100000, burnin = 1000,
                       returnPostSample = TRUE)

ancovaDeltaTab           <- cbind(t(ancovaModel$dt_meandir), t(ancovaModel$dt_CCI))
rownames(ancovaDeltaTab) <- paste0("$\\delta_", c("{df}", "{in}"), "$")
ancovaBetaTab            <- cbind(t(ancovaModel$bt_mean), t(ancovaModel$bt_CCI))
rownames(ancovaBetaTab)  <- paste0("$\\beta_", c("{age}", "{hand}"), "$")



ancovaTab <- with(ancovaModel, rbind("$\\beta_0$" = c(b0_meandir, b0_CCI),
                                     "$\\kappa$" = c(kp_mode, kp_HDI),
                                     ancovaDeltaTab,
                                     ancovaBetaTab
                                     # "DIC" = c(DIC, NA, NA),
                                     # "WAIC" = c(WAIC1, NA, NA))
                  ))
colnames(ancovaTab) <- c("Estimate", "LB", "UB")

## ----ANCOVATable, results='asis'-----------------------------------------
xANCOVATab <- xtable(ancovaTab,
                     caption = "Results for the ANCOVA model.",
                     label = "ANCOVATable")

print(xANCOVATab,
      latex.environments = c("center"),
      include.rownames = TRUE)

## ----ExampleIneq---------------------------------------------------------
muChs <- anovaModel$mu_chain
colnames(muChs) <- c("Control", "Deaf", "Interpreter")

# The fits for the two hypotheses H1: Deafs better (=lower) than everyone and
# H2: Controls worse (=higher) than everyone else
fitH1 <- mean(muChs[, "Deaf"] < muChs[, "Interpreter"] &
                muChs[, "Deaf"] < muChs[, "Control"])
fitH2 <- mean(muChs[, "Deaf"] < muChs[, "Control"]     &
                muChs[, "Interpreter"] < muChs[, "Control"])

BFH1H2 <- fitH1/fitH2


# The Complexities are the same: .25, the chance that two values both hit a .5 chance.

# plot( density(anovaModel$mu_chain[,1]), xlim = c(-1, 7))
# lines(density(anovaModel$mu_chain[,2]))
# lines(density(anovaModel$mu_chain[,3]))


