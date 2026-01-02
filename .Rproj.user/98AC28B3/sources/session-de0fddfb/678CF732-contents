source("visual_output/pair_plot_func.R")

library(mvtnorm)
library(matrixStats)


generate_data <- function(theta_true, n, b, a, sigma) {
  d <- length(theta_true)
  Y_bar <- numeric(d)
  Y_bar[1] <- rnorm(1, mean = theta_true[1], sd = sigma/sqrt(n))
  
  for (j in 2:d) {
    mu_j <- theta_true[j] + b[j-1] * (theta_true[j-1]^2 - a[j-1])
    Y_bar[j] <- rnorm(1, mean = mu_j, sd = sigma/sqrt(n))
  }
  return(Y_bar)
}

exact_simulation <- function(Y_bar, n, b, a, sigma, n_samples = 1e4) {
  d <- length(Y_bar)
  samples <- matrix(NA, nrow = n_samples, ncol = d,
                    dimnames = list(NULL, paste0("theta", 1:d)))
  
  sigma_j <- sigma / sqrt(n)
  samples[,1] <- rnorm(n_samples, mean = Y_bar[1], sd = sigma_j)
  
  for (j in 2:d) {
    mu_j_vec <- Y_bar[j] - b[j-1] * (samples[,j-1]^2 - a[j-1])
    samples[,j] <- rnorm(n_samples, mean = mu_j_vec, sd = sigma_j)
  }
  
  as.matrix(samples)
}

log_likelihood <- function(theta, Y_bar, n, b, a, sigma) {
  d    <- length(theta)
  sd_j <- sigma / sqrt(n)
  mu <- numeric(d)
  mu[1] <- theta[1]
  if (d > 1) {
    for (j in 2:d) {
      mu[j] <- theta[j] + b[j-1] * (theta[j-1]^2 - a[j-1])
    }
  }
  sum(dnorm(Y_bar, mean = mu, sd = sd_j, log = TRUE))
}

log_prior <- function(theta) {
  return(0)
}

log_post <- function(theta, Y_bar, n, b, a, sigma) {
  log_prior(theta) + log_likelihood(theta, Y_bar,n, b, a, sigma)
}

# ================= d=2 =======================
set.seed(100)
d <- 2
theta_true <- seq(1,2.8, length.out=d)
b <- rep(-1, d-1)
a <- rep(0, d-1)
sigma <- 8
n <-  20
n_samples <- 1e5
hpd_level <- 0.75

Y_bar <- generate_data(theta_true, n, b, a, sigma)
sample_posterior <- exact_simulation(Y_bar, n, b, a, sigma, n_samples)
lps <- apply(sample_posterior, 1, log_post, Y_bar,n, b, a, sigma)

pair_plot(sample_posterior, pixs=0.1)


# ================= d=5 =======================

set.seed(100)
d <- 5
theta_true <- seq(1, 2.8, length.out = d)
b <- rep(1.0, d - 1)
a <- rep(0, d - 1)
sigma <- 1
n <-  200
n_samples <- 1e5
Y_bar <- generate_data(theta_true, n, b, a, sigma)
sample_posterior <- exact_simulation(Y_bar, n, b, a, sigma, n_samples)
lps <- apply(sample_posterior, 1, log_post, Y_bar,n, b, a, sigma)
hpd_level <- 0.75
c_thresh <- as.numeric(quantile(lps, probs = 1 - hpd_level))
in_hpd   <- lps >= c_thresh
pair_plot(sample_posterior, pixs=0.1)


# ================= d=10 =======================
set.seed(100)
d <- 10
theta_true <- seq(1, 2.8, length.out = d)
b <- rep(0.5, d - 1)
a <- rep(0, d - 1)
sigma <- 1
n <-  200
n_samples <- 1e5
Y_bar <- generate_data(theta_true, n, b, a, sigma)
sample_posterior <- exact_simulation(Y_bar, n, b, a, sigma, n_samples)
lps <- apply(sample_posterior, 1, log_post, Y_bar,n, b, a, sigma)
hpd_level <- 0.75
c_thresh <- as.numeric(quantile(lps, probs = 1 - hpd_level))
in_hpd   <- lps >= c_thresh

pair_plot(sample_posterior, pixs=0.1)
