library(mvtnorm)
library(parallel)

set.seed(123)

# ============================================================
# Parameters
# ============================================================
n_obs <- 20       # Number of observations
# n_obs <- 50       # Number of observations
# n_obs <- 100       # Number of observations
d <- 2            # Dimension
N <- 1e6          # Number of posterior samples per Y
M <- 100          # Number of different Y datasets

# Prior parameters
s0 <- 1           # Prior variance scalar (prior is N(0, s0 * I))

# Pre-compute posterior variance
sn <- 1 / (n_obs + 1/s0)

# ============================================================
# Log-posterior function (vectorized)
# ============================================================
log_post_vec <- function(theta_matrix, Y, sig2) {
  N_samples <- nrow(theta_matrix)
  n <- nrow(Y)
  d <- ncol(theta_matrix)
  
  Y_sq <- sum(Y^2)
  theta_sq <- rowSums(theta_matrix^2)
  sum_Y <- colSums(Y)
  cross_term <- theta_matrix %*% sum_Y
  sum_sq_dist <- Y_sq + n * theta_sq - 2 * cross_term
  
  log_likelihood <- -0.5 * n * d * log(2 * pi) - 0.5 * sum_sq_dist
  log_prior <- -0.5 * d * log(2 * pi * sig2) - 0.5 * theta_sq / sig2
  
  return(log_likelihood + log_prior)
}

# ============================================================
# Compute posterior parameters and exact evidence
# ============================================================
compute_posterior_params <- function(Y, s0, n_obs, d) {
  mn <- colSums(Y) / (n_obs + 1/s0)
  
  log_term1 <- (-d * n_obs / 2) * log(2 * pi)
  log_term2 <- (-d / 2) * log(1 + s0 * n_obs)
  
  sum_yij2 <- sum(Y^2)
  sum_yij <- colSums(Y)
  sum_yij_sq <- sum(sum_yij^2)
  
  log_exact <- log_term1 + log_term2 - 0.5 * (sum_yij2 - (s0 / (1 + s0 * n_obs)) * sum_yij_sq)
  
  list(mn = mn, log_exact = log_exact)
}

# ============================================================
# Single replicate function
# ============================================================
run_single_Y <- function(seed_offset, N, s0, sn, n_obs, d) {
  set.seed(123 + seed_offset)
  
  Y <- rmvnorm(n_obs, mean = rep(0, d), sigma = diag(d))
  post_params <- compute_posterior_params(Y, s0, n_obs, d)
  samples <- rmvnorm(N, mean = post_params$mn, sigma = sn * diag(d))
  lps <- log_post_vec(samples, Y, s0)
  
  list(samples = samples, lps = lps, Y = Y, log_exact = post_params$log_exact)
}

# ============================================================
# Main execution
# ============================================================
n_cores <- max(1, detectCores() - 1)

if (.Platform$OS.type == "unix") {
  results <- mclapply(1:M, function(m) {
    run_single_Y(m, N, s0, sn, n_obs, d)
  }, mc.cores = n_cores)
} else {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("N", "s0", "sn", "n_obs", "d",
                      "log_post_vec", "compute_posterior_params", "run_single_Y"))
  clusterEvalQ(cl, library(mvtnorm))
  results <- parLapply(cl, 1:M, function(m) {
    run_single_Y(m, N, s0, sn, n_obs, d)
  })
  stopCluster(cl)
}

# Extract results
all_samples <- lapply(results, function(x) x$samples)
all_lps <- lapply(results, function(x) x$lps)
all_Y <- lapply(results, function(x) x$Y)
all_log_exact <- sapply(results, function(x) x$log_exact)

# Save results
gaussian_2d_samples <- list(
  samples = all_samples,
  lps = all_lps,
  Y = all_Y,
  log_exact = all_log_exact,
  n = n_obs,
  d = d,
  N = N,
  M = M,
  s0 = s0,
  sn = sn
)

save(gaussian_2d_samples, file = sprintf("gaussian_2d_multiple_Yn%d.RData", n_obs))