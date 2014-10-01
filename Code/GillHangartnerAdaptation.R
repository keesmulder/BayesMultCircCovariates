# CODE FOR Jeff Gill and Dominik Hangartner. \U3e32653c\u20ac\u0153Circular
# Data in Political Science and How to Handle It.\U3e32653c\u20ac Political
# Analysis, 18:2, 316-336 (2010).

# Here's the ML routine: circ.lin.reg

library(circular)

circ.lin.reg <- function(x, theta, beta0, trace = FALSE, print = TRUE,
                         tol = 1e-10, maxiter = 1000) {
  # THIS PART OF THE CODE IS HEAVILY BASED ON S-CODE FROM ULRIC LUND'S
  # CIRCSTATS V2.0 (2004)
  # Ensure length, save n
  if (is.vector(x)) x <- cbind(x)
  n <- length(theta)

  # Obtain base, at zero
  betaPrev <- beta0

  # S_bar and C_bar given current beta
  S <- sum(sin(theta - 2 * atan(x %*% betaPrev)))/n
  C <- sum(cos(theta - 2 * atan(x %*% betaPrev)))/n
  R <- sqrt(S^2 + C^2)
  mu <- atan2(S, C)

  # Estimate of kappa.
  k <- A1inv(R)


  diff <- tol + 1
  iter <- 0
  S.function <- function(betaPrev, x) {
    2/(1 + (t(betaPrev) %*% x)^2)
  }

  # Loop until max or within tolerance
  while (diff > tol & iter < maxiter) {
    iter <- iter + 1

    # u as in Fisher
    u <- k * sin(theta - mu - 2 * atan(x %*% betaPrev))

    # Has estimate of R * kappa on the diagonal.
    A <- diag(k * A1(k), nrow = n)

    #
    g.p <- diag(apply(x, 1, S.function, betaPrev = betaPrev), nrow = n)
    D <- g.p %*% x
    betaNew <- lm(t(D) %*% (u + A %*% D %*% betaPrev) ~ t(D) %*%
                    A %*% D - 1)$coefficients
    diff <- abs(max(betaNew - betaPrev))

    breaked = 0
    if (iter == 1000) {
      breaked = 1
    }
    if (is.na(diff) == TRUE) {
      breaked = 1
      break
    }
    if (max(betaNew) > 100) {
      breaked = 1
      break
    }
    betaPrev <- betaNew
    S <- sum(sin(theta - 2 * atan(x %*% betaPrev)))/n
    C <- sum(cos(theta - 2 * atan(x %*% betaPrev)))/n
    R <- sqrt(S^2 + C^2)
    mu <- atan2(S, C)
    # mu <- asin(S/R)
    k <- A1inv(R)
    if (trace == T) {
      log.lik <- -n * log(2 * pi * I.0(k)) + k * sum(cos(theta -
                                                           mu - 2 * atan(x %*% betaNew)))
      cat("Iteration ", iter, ":    Log-Likelihood = ", log.lik,
          "mu ", mu, "k ", k, "b ", betaNew, "\n")
      cat("Starting values ", i, "\n")
    }
  }
  log.lik <- -n * log(2 * pi * I.0(k)) + k * sum(cos(theta - mu -
                                                       2 * atan(x %*% betaNew)))
  log.lik.old <- -n * log(2 * pi * I.0(k)) + k * sum(cos(theta -
                                                           mu - 2 * atan(x %*% betaPrev)))
  cov.beta <- solve(t(D) %*% A %*% D)
  se.beta <- sqrt(diag(cov.beta))
  se.kappa <- sqrt(1/(n * (1 - A1(k)^2 - A1(k)/k)))
  circ.se.mu <- 1/sqrt((n - ncol(x)) * k * A1(k))
  z.values <- abs(betaNew/se.beta)
  p.values <- (1 - pnorm(z.values)) * 2
  result.matrix <- cbind(Coef = betaNew, SE = se.beta, Z = z.values,
                         p = p.values)
  dimnames(result.matrix) <- list(dimnames(x)[[2]], c("Coef", "SE",
                                                      "|z|", "p"))
  cat("\n", "Circular-Linear Regression", "\n", "\n")
  print(result.matrix)
  cat("\n", "\n")
  betaNew <- as.matrix(betaNew)
  dimnames(betaNew) <- list(dimnames(x)[[2]], c("Estimate"))
  list(mu = mu, kappa = k, beta = betaNew, log.lik = log.lik, log.lik.old = log.lik.old,
       circ.se.mu = circ.se.mu, se.kappa = se.kappa, cov.beta = cov.beta,
       se.beta = se.beta, result.matrix = result.matrix, breaked = breaked)
}

# And here's the code for the Metropolis-Hastings: MH sampler for
# Fisher & Lee (1992)-style homoscedastic circular-linear regression
# model The sampler is as follows - calculated kappa u ~ U(0, 2pi] beta
# ~ MVN(-,Sigma) with Sigma = Information^(-1) or Sigma =
# diag(Information^(-1))

library(circular)

ll.vm.atan <- function(b) {
  N <- dim(X)[1]
  p <- dim(X)[2] + 1
  coef <- as.vector(b)
  kappa <- kappa
  mu <- coef[1]
  beta <- coef[2:p]
  ll1 <- -N * log(2 * pi * I.0(kappa))
  g.atan <- mu + 2 * atan(X %*% beta)
  ll2 <- kappa * sum(cos(Y - g.atan))
  ll <- ll1 + ll2
  # sum priors
  logpriors <- sum(dnorm(beta, 0, 10, log = TRUE))
  post <- ll + logpriors
  return(post)
}

X <- matrix(rnorm(100))
Y <- rnorm(100)

out <- circ.lin.reg(X, Y, c(0), trace = T, print = T, tol = 1e-10)

### run ML code circ.lin.reg.R here and save object as out to get
### starting values
kappa <- out$kappa
out.par <- as.vector(c(out$mu, out$beta))
varcov.mat <- matrix(0, length(out$beta) + 1, length(out$beta) + 1)
varcov.mat[1, 1] <- out$circ.se.mu^2
varcov.mat[2:(length(out$beta) + 1), 2:(length(out$beta) + 1)] <- out$cov.beta
se <- sqrt(diag(abs(varcov.mat)))
z.val <- abs(out.par/se)
par.list <- cbind(out.par, as.vector(se), as.vector(z.val))
rownames(par.list) <- c("mu", colnames(cbind(X)))
colnames(par.list) <- c("coef", "se", "|z|")
print(par.list)
var.mat <- diag(1, dim(varcov.mat)[1]) * diag(varcov.mat)

mh.circ.lin = function(niter = 1000, scale = 0.25, verbose = 100) {
  # GAUSSIAN RANDOM WALK MH SAMPLER FOR THE HOMOSCEDASTIC CIRCULAR-LINEAR
  # REGRESSIOn MODEL UNDER -------- PRIORs

  library(mnormt)
  p = length(out.par)
  betas = matrix(0, niter, p)
  betas[1, ] = as.vector(out.par)
  kappas = matrix(0, niter, 1)
  kappas[1] = kappa
  likelihood = matrix(0, niter, 1)
  likelihood[1] = ll.vm.atan(out.par)
  Sigma2 = as.matrix(var.mat)
  ver = seq(1, niter, by = verbose)

  for (i in 2:niter) {
    tildebeta = rmnorm(1, betas[i - 1, ], scale * Sigma2)
    tildebeta[1] = max(0, min(tildebeta[1], 2 * pi))
    llnew = ll.vm.atan(tildebeta)
    llold = ll.vm.atan(betas[i - 1, ])
    llr = llnew - llold
    draw = runif(1)
    if (draw <= exp(llr)) {
      betas[i, ] = tildebeta
      likelihood[i] = llnew
      # get kappa
      S <- sum(sin(Y - betas[i, 1] - 2 * atan(X %*% betas[i, 2:p])))/dim(X)[1]
      C <- sum(cos(Y - betas[i, 1] - 2 * atan(X %*% betas[i, 2:p])))/dim(X)[1]
      R <- sqrt(S^2 + C^2)
      kappas[i] = A1inv(R)
    } else {
      betas[i, ] = betas[i - 1, ]
      likelihood[i] = llold
      # get kappa
      S <- sum(sin(Y - betas[i, 1] - 2 * atan(X %*% betas[i, 2:p])))/dim(X)[1]
      C <- sum(cos(Y - betas[i, 1] - 2 * atan(X %*% betas[i, 2:p])))/dim(X)[1]
      R <- sqrt(S^2 + C^2)
      kappas[i] = A1inv(R)
    }
    if (is.na(match(i, ver)) == FALSE)
      print(c(iter = i, betas[i, ], like = likelihood[i]))
  }
  list(betas = betas, likelihood = likelihood, kappas = kappas)
}

post.samp = mh.circ.lin(niter = 10000, scale = 0.001, verbose = 1000)
