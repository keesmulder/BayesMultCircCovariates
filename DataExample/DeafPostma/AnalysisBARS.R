source("Code/describeCirc.R")
source("DataAnalysis/circGLM.R")

require(dplyr)
require(ggplot2)
require(tidyr)
require(foreign)
fileloc <- "DataExample/Data_BARS_BREED_N48_23052011.sav"
pd      <- read.spss(fileloc, to.data.frame = TRUE)

# Remove person with low raven percentile
tbl_df(pd) %>% filter(Ravenpercentile > 30) -> dpd

# Recreate table 1.
dpd %>%
  group_by(CODE) %>%
  summarise(meanage = mean(age), sdage = sd(age),
            minage = min(age), maxage = max(age),
            meaned = mean(toteducation), sded = sd(toteducation),
            raven_percentile = mean(Ravenpercentile))


# Recreate initial group dif F-tests for assumptions
agelm <- lm(age ~ CODE, data = dpd)
edulm  <- lm(toteducation  ~ CODE, data = dpd)
anova(agelm)
anova(edulm)


# Recreate Fig 2
# Table version
dpd %>%
  group_by(CODE) %>%
  summarise(mean_im = mean(imdev_signed),
            se_im = sqrt(var(imdev_signed)/n()),
            mean_dl = mean(deldev_signed),
            se_dl = sqrt(var(deldev_signed)/n())) -> tb.dev.imdl

tbnum <- rbind(as.matrix(tb.dev.imdl[, 2:3]), as.matrix(tb.dev.imdl[, 4:5]))
tb    <- data.frame(rep(tb.dev.imdl$CODE, 2), c(rep("im", 3), rep("dl", 3)), tbnum)
colnames(tb) <- c("Code", "Condition", "Mean", "SE")
tb

# Get figure
lims <- aes(ymin = Mean - SE, ymax = Mean + SE, group = Condition)
posdod <- position_dodge(-0.3)
ggplot(tb, aes(x = Code, col = Condition, y = Mean)) +
  geom_point(mapping = aes(shape = Condition), position = posdod, size=3.5) +
  geom_errorbar(lims, position = posdod, width = .2) +
  ylim(0, 35) + theme_bw()


# Recreate Fig 3
ggplot(dpd, aes(x = imdev_signed, y = deldev_signed, col = nwgroup)) +
  geom_point() + theme_bw()




# Recreate main results
shpd <- select(dpd, imdev_signed, deldev_signed, CODE)
longpd <- gather(shpd, delay, dev, -CODE)

# "No significant effect of delay"
devdellm <- lm(dev ~ delay,  data = longpd)
anova(devdellm)






