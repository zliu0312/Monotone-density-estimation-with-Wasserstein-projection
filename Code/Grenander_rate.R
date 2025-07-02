rm(list = ls())
library(fdrtool)

u_grid <- seq(0, 1, by = 0.00001)
u_grid <- u_grid[1:(length(u_grid) - 1)]  # to avoid +Inf
N <- length(u_grid)
Q_true <- qexp(u_grid)



n_max <- 50000  # sample size
x_all <- rexp(n_max)  # X_i ~ Exp(1)
n_seq <- seq(1000, n_max, by = 1000)

L2dist_seq <- numeric(length(n_seq))  # Wasserstein

for (i in 1:length(n_seq)) {
  n <- n_seq[i]
  
  x <- x_all[1:n]
  e <- ecdf(x)
  g <- grenander(e)
  g_cdf <- g$F
  y <- environment(g_cdf)$x
  Q <- quantile(y, u_grid, type = 1)  # quantile function of Grenander
  L2dist_seq[i] <- sqrt(sum((Q - Q_true)^2)/N)  # Wasserstein
}

plot(log10(n_seq), log10(L2dist_seq), type = "p", col = "blue")
output <- lm(log10(L2dist_seq[-(1:10)]) ~ log10(n_seq[-(1:10)]))
output$coefficients[2]
abline(a = output$coefficients[1], b = output$coefficients[2])
abline(a = output$coefficients[1], b = -0.5, lty = 2, col = "darkgrey")





