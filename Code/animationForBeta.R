require(animation)
options(scipen=30)
d <- 60

btseq <- c(0, 1.4^((-d:d)/2), -1.4^((d:-d)/2))

# Function to get animations
anim <- function(btseq) {
  for(bt in btseq) {

    plot(function(x) 2*atan(bt*x), xlim=c(-10, 10), ylim=c(-pi, pi), n=1000,
         main=paste("Beta = ", round(bt, 0)))
    plot(function(x) 2*atan(bt*x)+2*pi, xlim=c(-10, 10), ylim=c(-pi, pi), add=TRUE, n=1000)
    plot(function(x) 2*atan(bt*x)-2*pi, xlim=c(-10, 10), ylim=c(-pi, pi), add=TRUE, n=1000)
    animation::ani.pause()
  }
}

# Save the gif to current directory.
saveGIF(expr = anim(btseq), interval=.03, outdir = getwd(), ani.height=800, ani.width=1000)
