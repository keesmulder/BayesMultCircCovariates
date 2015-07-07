source('Data/generateCircularGLMData.R')
source('Simulation/simulationStudyCircGLM.R')
source('SimulationAnalysis/generateBetaShapePlot.R')
source('SimulationAnalysis/compareSimResults.R')

dir("Simulation/Results")

load(paste0("Simulation/Results/[simStudCircGLM][nsim1000][Q10000][burnin100]",
            "[r2][bt_prior1][seed389238][n12,100][kp0.5,30]",
            "[btl=0.1,l=-1,lll=0.1,lll=-1,llllll=0.1,llllll=-1].rda"))
simbigQ <- simStudyResults
load(paste0("Simulation/Results/[simStudCircGLM][nsim10000][Q1000][burnin100]",
            "[r2][bt_prior1][seed389238][n12,100][kp0.5,30]",
            "[btl=0.1,l=-1,lll=0.1,lll=-1,llllll=0.1,llllll=-1].rda"))
simbignsim <- simStudyResults
load(paste0("Simulation/Results/[simStudCircGLM][nsim1000][Q1000][burnin100]",
            "[r2][bt_prior1][seed389238][n12,100][kp0.5,30]",
            "[btl=0.1,l=-1,lll=0.1,lll=-1,llllll=0.1,llllll=-1].rda"))
simbothsmall <- simStudyResults
load(paste0("Simulation/Results/[simStudCircGLM][nsim10000][Q10000][burnin100]",
            "[r2][bt_prior1][seed389238][n12,100][kp0.5,30]",
            "[btl=0.1,l=-1,lll=0.1,lll=-1,llllll=0.1,llllll=-1].rda"))
simbothbig <- simStudyResults

compr <- compareSimRes(simbothbig, simbigQ, simbignsim, simbothsmall, type = "sidebyside", digits=2)


# The results are somewhat comparable
compr
# Main results
simStudyResults

# Further inspection
# (1), 12, 0.5, bt_1_mean seems wrong
plot(simStudyResults, btDesNumber=1, n=12, kp=0.5, stat="bt_1_mean", breaks=100)
# There is just tons of variance.

# (2), 12, 0.5, bt_1_mean seems wrong
plot(simStudyResults, btDesNumber=2, n=12, kp=0.5, stat="bt_1_mean", breaks=100)
# There is just tons of variance.


# Serious problems start with (4).
# 100, 0.5, bt_1_mean is way off regardless of large sample size.
plot(simStudyResults, btDesNumber=4, n=100, kp=0.5, stat="bt_2_mean", breaks=100)
mean(slice.cGLMSim(simStudyResults, btDesNumber=4, n=100, kp=0.5, stat="bt_1_mean"))
# There does not even seem to be a peak at -1.

# 100, 30, shows the first b0 inversion.
plot(simStudyResults, btDesNumber=4, n=100, kp=30, stat="b0_meandir", breaks=100)
abline(v=pi/2, col="green")
# You can even see that some values were in the process of inverting. Inversion
# happens when large parts (>.5) of the predicted values are more than pi/2 away from
# beta_0.

plot(simStudyResults, btDesNumber=4, n=100, kp=30, stat="bt_1_mean", breaks=100)
plot(simStudyResults, btDesNumber=4, n=100, kp=30, stat="bt_2_mean", breaks=100)
# Here it is clear that there is an inverted group and that there is a
# non-inverted group.

# Further demonstration.
plot(slice.cGLMSim(simStudyResults,
                   btDesNumber=4, n=100, kp=30, stat="b0_meandir"),
     slice.cGLMSim(simStudyResults,
                   btDesNumber=4, n=100, kp=30, stat="bt_1_mean"),
     xlab="beta_0", ylab="bt_1")


# (5) n=12, kp=0.5, features inversions, but is mostly all over the place.
plot(simStudyResults, btDesNumber=5, n=12, kp=0.5, stat="b0_meandir", breaks=100)
plot(simStudyResults, btDesNumber=5, n=12, kp=0.5, stat="bt_1_mean", breaks=100)

# The estimate for kappa is influenced by several large values, which makes
# sense with nine predictors and n=12.
plot(simStudyResults, btDesNumber=5, n=12, kp=0.5, stat="kp_mode", breaks=100)

# With n=100 the inversions happen more than half of the time.
plot(simStudyResults, btDesNumber=5, n=100, kp=0.5, stat="b0_meandir", breaks=100)
plot(simStudyResults, btDesNumber=5, n=100, kp=0.5, stat="bt_1_mean", breaks=100)

# Only a few errors with kp=30.
plot(simStudyResults, btDesNumber=5, n=12, kp=30, stat="b0_meandir", breaks=100)

# No errors when n=100 as well.
plot(simStudyResults, btDesNumber=5, n=100, kp=30, stat="b0_meandir", breaks=100)

# Kappa is messed up, though, because of n=12 with 9 preds.
plot(simStudyResults, btDesNumber=5, n=12, kp=30, stat="kp_mode", breaks=100)

# n=100 solves this.
plot(simStudyResults, btDesNumber=5, n=100, kp=30, stat="kp_mode", breaks=100)

linkfun <- function(x) 2 * atan(x)
# generateBetaShapePlot(b0_cur = pi/2, true_b0 = pi/2, true_bt = -1, Xsd = 1, n = 100)

slice.cGLMSim(simStudyResults, btDesNumber=3, n=12, kp=0.5, stat="kp_mode")
plot(simStudyResults, btDesNumber=3, n=12, kp=30, stat="kp_mode", breaks=100)
plot(simStudyResults, btDesNumber=3, n=100, kp=30, stat="kp_mode", breaks=100)

# Definitely inversion here.
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="b0_meandir", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="kp_mode", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="bt_1_mean", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="bt_2_mean", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="bt_3_mean", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="bt_4_mean", breaks=100)
plot(simStudyResults, btDesNumber=6, n=100, kp=30, stat="bt_5_mean", breaks=100)

simStudyResults


