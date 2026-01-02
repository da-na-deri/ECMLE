library(parallel)

set.seed(123)

# ============================================================
# Parameters
# ============================================================
n_obs <- 20       # Number of observations
# n_obs <- 50       # Number of observations
# n_obs <- 100       # Number of observations
d <- 5            # Dimension
N <- 1e6          # Number of posterior samples per Y
M <- 100          # Number of different Y datasets

# Banana distribution parameters
theta_true <- seq(1, 2.8, length.out = d)
b <- rep(1.0, d - 1)
a <- rep(0, d - 1)
sigma <- 1

# ============================================================
# Data generation function
# ============================================================
generate_data <- function(theta_true, n, b, a, sigma) {
  d <- length(theta_true)
  Y_bar <- numeric(d)
  Y_bar[1] <- rnorm(1, mean = theta_true[1], sd = sigma / sqrt(n))
  
  for (j in 2:d) {
    mu_j <- theta_true[j] + b[j-1] * (theta_true[j-1]^2 - a[j-1])
    Y_bar[j] <- rnorm(1, mean = mu_j, sd = sigma / sqrt(n))
  }
  return(Y_bar)
}

# ============================================================
# Exact posterior simulation (sequential sampling)
# ============================================================
exact_simulation <- function(Y_bar, n, b, a, sigma, n_samples = 1e4) {
  d <- length(Y_bar)
  samples <- matrix(NA, nrow = n_samples, ncol = d)
  sigma_j <- sigma / sqrt(n)
  
  samples[, 1] <- rnorm(n_samples, mean = Y_bar[1], sd = sigma_j)
  
  for (j in 2:d) {
    mu_j_vec <- Y_bar[j] - b[j-1] * (samples[, j-1]^2 - a[j-1])
    samples[, j] <- rnorm(n_samples, mean = mu_j_vec, sd = sigma_j)
  }
  
  return(samples)
}

# ============================================================
# Log-posterior function (vectorized)
# ============================================================
log_post_vec <- function(theta_matrix, Y_bar, n, b, a, sigma) {
  n_samples <- nrow(theta_matrix)
  d <- ncol(theta_matrix)
  sd_j <- sigma / sqrt(n)
  
  mu_matrix <- matrix(NA, nrow = n_samples, ncol = d)
  mu_matrix[, 1] <- theta_matrix[, 1]
  
  for (j in 2:d) {
    mu_matrix[, j] <- theta_matrix[, j] + b[j-1] * (theta_matrix[, j-1]^2 - a[j-1])
  }
  
  # Log-likelihood with flat prior
  rowSums(dnorm(matrix(Y_bar, nrow = n_samples, ncol = d, byrow = TRUE),
                mean = mu_matrix, sd = sd_j, log = TRUE))
}

# ============================================================
# Single replicate function
# ============================================================
run_single_Y <- function(seed_offset, N, theta_true, b, a, sigma, n_obs) {
  set.seed(123 + seed_offset)
  
  Y_bar <- generate_data(theta_true, n_obs, b, a, sigma)
  samples <- exact_simulation(Y_bar, n_obs, b, a, sigma, n_samples = N)
  lps <- log_post_vec(samples, Y_bar, n_obs, b, a, sigma)
  
  list(samples = samples, lps = lps, Y_bar = Y_bar)
}

# ============================================================
# Main execution
# ============================================================
n_cores <- max(1, detectCores() - 1)

if (.Platform$OS.type == "unix") {
  results <- mclapply(1:M, function(m) {
    run_single_Y(m, N, theta_true, b, a, sigma, n_obs)
  }, mc.cores = n_cores)
} else {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("N", "theta_true", "b", "a", "sigma", "n_obs",
                      "generate_data", "exact_simulation", "log_post_vec",
                      "run_single_Y"))
  results <- parLapply(cl, 1:M, function(m) {
    run_single_Y(m, N, theta_true, b, a, sigma, n_obs)
  })
  stopCluster(cl)
}

# Extract results
all_samples <- lapply(results, function(x) x$samples)
all_lps <- lapply(results, function(x) x$lps)
all_Y_bar <- lapply(results, function(x) x$Y_bar)

# Save results
banana_5d_samples <- list(
  samples = all_samples,
  lps = all_lps,
  Y_bar = all_Y_bar,
  n = n_obs,
  d = d,
  N = N,
  M = M,
  theta_true = theta_true,
  b = b,
  a = a,
  sigma = sigma
)

save(banana_5d_samples, file = sprintf("banana_d%d_multiple_Yn%d.RData", d, n_obs))