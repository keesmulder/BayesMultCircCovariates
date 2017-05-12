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


# Plot the log likelihood as a function of beta.
plotbeta <- function(th, X, normalPrior=FALSE, res = 100, xl = c(-10, 10),
                     b0 = pi/2, kp = 1, r = 2, bt_true = 1, mu=0, sds=10) {


  # X points at which to evaluate the ll function
  sq <- seq(xl[1], xl[2], length.out = res)

  p <- ggplot() + theme_bw() +
    ylab(expression(paste("Conditional Log-Likelihood of ", beta))) +
    xlab(expression(beta))

  # Find the ll function to add to the plot.
  for (i in 1:length(sds)) {
    if (!normalPrior) {
      betall <- Vectorize(function(b) {
        ll(b0 = b0, kp = kp, bt = b, dt = numeric(0), th = th,
           X = X, D = matrix(ncol = 0, nrow = nrow(X)), r = r)
      })
    } else {
      sdi <- sds[i]
      betall <- Vectorize(function(b) {
        ll(b0 = b0, kp = kp, bt = b, dt = numeric(0), th = th,
           X = X, D = matrix(ncol = 0, nrow = nrow(X)), r = r) +
          logProbNormal(b, mu, sdi)
      })
    }

    # Apply the function to the X sequence
    btlls <- betall(sq)

    # Add line to plot
    p <- p +
      geom_line(data = data.frame(bll = btlls, bt = sq),
                aes(x = bt, y = bll),
                linetype = i)
  }

  p
}


