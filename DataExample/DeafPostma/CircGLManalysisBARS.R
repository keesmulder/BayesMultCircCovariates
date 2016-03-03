
source("Code/describeCirc.R")
source('Code/vonMises.R')
source("DataAnalysis/circGLM.R")

require(dplyr)
require(ggplot2)
require(tidyr)
require(circular)
require(foreign)

# Load in data
fileloc <- "DataExample/DeafPostma/Data_BARS_BREED_N48_23052011.sav"
pd      <- read.spss(fileloc, to.data.frame = TRUE)

# Remove person with low raven percentile
tbl_df(pd) %>% filter(Ravenpercentile > 30) -> dpd

# Recreate table 1.
dpd %>%
  group_by(CODE) %>%
  summarise(meanage = mean(age), sdage = sd(age),
            minage = min(age), maxage = max(age),
            meaned = mean(toteducation), sded = sd(toteducation),
            raven_percentile = mean(Ravenpercentile)) -> tab1
tab1

# Create dummies
dpd %>%
  select(CODE, imdev_signed) %>%
  mutate(INT = model.matrix(~dpd$CODE)[,2],
         CON = model.matrix(~dpd$CODE)[,3]) -> smallpd

# Mean and mean direction for each group
smallpd %>% group_by(CODE) %>% summarise(mean = mean(imdev_signed))
smallpd %>% group_by(CODE) %>% summarise(mean_dir = meanDir(imdev_signed*pi/180)*180/pi)

# Run the basic GLM model with just the intercept.
icep_cgm <- circGLM(th = smallpd$imdev_signed, X = matrix(ncol=0, nrow=47), Q = 10000,
                    returnPostSample = TRUE, debug = FALSE, loopDebug = FALSE)
# Run the basic GLM model with just groups as predictors.
base_cgm <- circGLM(th = smallpd$imdev_signed, X = smallpd[, 3:4], Q = 100000, burnin = 1000,
                    returnPostSample = TRUE, debug = FALSE, loopDebug = FALSE)


plot(base_cgm)
print(base_cgm)

plot(density(base_cgm$dt_chain[, 1]))
lines(density(base_cgm$dt_chain[, 2]), col="blue")



# Information Criteria Analysis
IC_compare.circGLM(icep_cgm, base_cgm)



# ANOVA Bayes Factors
mu1_chain <- base_cgm$b0_chain
mu2_chain <- mu1_chain + base_cgm$dt_chain[,1]
mu3_chain <- mu1_chain + base_cgm$dt_chain[,2]

# Check the support for the hypothesis delta_i > 0 vs. delta_i < 0.
base_cgm$BayesFactors
prop2over1 <- sum(mu2_chain > mu1_chain)/100000
prop2over1 <- sum(base_cgm$dt_chain[,1] > 0)/100000
prop2over1 / (1 - prop2over1)

prop3over2 <- sum(mu3_chain > mu2_chain)/100000
prop3over2 / (1 - prop3over2)
base_cgm$BayesFactors[1,1] / base_cgm$BayesFactors[2,1]



# Check that the group mean directions can be found from the output
dfc <- 180/pi
c(Deaf = base_cgm$b0_meandir * dfc,
  Int = base_cgm$b0_meandir * dfc + base_cgm$dt_meandir[1] * dfc,
  Con = base_cgm$b0_meandir * dfc + base_cgm$dt_meandir[2] * dfc)

# Visualise results.
groupmeandat <- data.frame(Deaf = base_cgm$b0_chain * dfc,
                  Int = base_cgm$b0_chain * dfc + base_cgm$dt_chain[, 1] * dfc,
                  Control = base_cgm$b0_chain * dfc + base_cgm$dt_chain[, 2] * dfc)

# Density Plot
ggmd <- gather(groupmeandat)
ggplot(ggmd, aes(x = value, col = key)) + geom_density() + theme_bw()


# Means and CI's
mds <- sapply(groupmeandat/dfc, meanDir)*dfc
LBs <- sapply(groupmeandat/dfc, circQuantile, q = .025)*dfc
UBs <- sapply(groupmeandat/dfc, circQuantile, q = .975)*dfc

cdname <- c("Deaf", "Interpreter", "Control")
tb <- data.frame(Code = factor(cdname, levels = cdname),
                 Mean = mds, LB = LBs, UB = UBs)

# Plot means and CI's
lims <- aes(ymin = LB, ymax = UB, group = Code)
ggplot(tb, aes(x = Code, y = Mean)) +
  geom_point(size=3.5) +
  geom_errorbar(lims, width = .2) +
  ylim(0, 40) + theme_bw()






# To Do here: Bayes Factor hypothesis testing.

# Post mean AIC (THIS Section is just for checking that the output is correct)

# Two ways to obtain p(y | phi_bar), where phi_bar are the estimates of all
# parameters.
# Through the ll function.
base_cgm_ll <- with(base_cgm,
                    ll(b0 = b0_meandir, kp = kp_mode, bt = bt_mean, dt = dt_meandir,
                       th = matrix(smallpd$imdev_signed)/dfc,
                       X = matrix(ncol=0, nrow=nrow(smallpd)),
                       D = as.matrix(smallpd[, 3:4]), r = 2))

# Taking the log of the von Mises with the predictions as means.
predictions <- base_cgm$b0_meandir +
                 as.matrix(smallpd[, 3:4]) %*% t(base_cgm$dt_meandir)
logvmsum_ll <- sum(log(dvm(th = matrix(smallpd$imdev_signed)/dfc,
            mu = predictions,
            kappa = base_cgm$kp_mode)))

# And now it's also directly in the results-object.
base_cgm$ll_th_estpars
base_cgm_ll
logvmsum_ll


base_cgm$ll_th_curpars







# Adding covariates
dpd %>%
  select(CODE, imdev_signed, deldev_signed,
         Ravenpercentile, handness, age, haptic_fluency) %>%
  mutate(INT = model.matrix(~smallpd$CODE)[,2],
         CON = model.matrix(~smallpd$CODE)[,3]) -> covpd

basep <- ggplot(covpd, aes(y = deldev_signed, col = CODE)) + theme_bw()
basep + geom_point(size = 4, aes(y = deldev_signed, x = Ravenpercentile))
basep + geom_point(size = 4, aes(y = imdev_signed,  x = Ravenpercentile))
basep + geom_point(size = 4, aes(y = deldev_signed, x = handness))
basep + geom_point(size = 4, aes(y = imdev_signed,  x = handness))
basep + geom_point(size = 4, aes(y = deldev_signed, x = age))
basep + geom_point(size = 4, aes(y = imdev_signed,  x = age))
basep + geom_point(size = 4, aes(y = deldev_signed, x = haptic_fluency))
basep + geom_point(size = 4, aes(y = imdev_signed,  x = haptic_fluency))










