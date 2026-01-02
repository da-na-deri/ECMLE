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

# Mixture components (8 components on octagon vertices)
K <- 8
radius <- 18

angles <- (0:(K-1)) * 2 * pi / K
mu_matrix <- cbind(radius * cos(angles), radius * sin(angles))

alpha_weights <- rep(1/K, K)

# Covariance
sigma_scalar <- 0.07
Sigma <- sigma_scalar * diag(d)

# Posterior variance
sn <- 1 / (n_obs + 1)
S_n_scalar <- sn * sigma_scalar

# ============================================================
# Helper function
# ============================================================
rowLogSumExps <- function(mat) {
  m <- apply(mat, 1, max)
  m + log(rowSums(exp(mat - m)))
}

# ============================================================
# Generate posterior samples
# ============================================================
generate_posterior_samples <- function(N, mn_matrix, S_n_scalar, w_list, K, d) {
  components <- sample(1:K, size = N, replace = TRUE, prob = w_list)
  comp_counts <- tabulate(components, nbins = K)
  samples <- matrix(0, nrow = N, ncol = d)
  
  for (k in 1:K) {
    if (comp_counts[k] > 0) {
      idx <- which(components == k)
      samples[idx, ] <- rmvnorm(comp_counts[k], mean = mn_matrix[k, ], 
                                sigma = S_n_scalar * diag(d))
    }
  }
  return(samples)
}

# ============================================================
# Log-posterior function (vectorized)
# ============================================================
log_post_vec <- function(theta_matrix, Y, sigma_scalar, alpha_weights, mu_matrix) {
  N_samples <- nrow(theta_matrix)
  n <- nrow(Y)
  d <- ncol(theta_matrix)
  K <- nrow(mu_matrix)
  
  # Log-likelihood
  Y_sq <- rowSums(Y^2)
  theta_sq <- rowSums(theta_matrix^2)
  cross_term <- Y %*% t(theta_matrix)
  sum_sq_dist <- sum(Y_sq) + n * theta_sq - 2 * colSums(cross_term)
  
  log_likelihood <- -0.5 * n * d * log(2 * pi * sigma_scalar) - 
    0.5 * sum_sq_dist / sigma_scalar
  
  # Log-prior (mixture)
  mu_sq <- rowSums(mu_matrix^2)
  cross_prior <- theta_matrix %*% t(mu_matrix)
  sq_dist_prior <- outer(theta_sq, mu_sq, "+") - 2 * cross_prior
  
  log_norm_const <- -0.5 * d * log(2 * pi * sigma_scalar)
  log_priors <- log_norm_const - 0.5 * sq_dist_prior / sigma_scalar + 
    matrix(log(alpha_weights), nrow = N_samples, ncol = K, byrow = TRUE)
  
  log_prior_density <- rowLogSumExps(log_priors)
  
  return(log_likelihood + log_prior_density)
}

# ============================================================
# Compute posterior parameters and exact evidence
# ============================================================
compute_posterior_params <- function(Y, mu_matrix, alpha_weights, sigma_scalar, n_obs, d, K) {
  Y_bar <- colMeans(Y)
  sum_Y_sq <- sum(Y^2)
  
  mn_matrix <- (n_obs * matrix(Y_bar, nrow = K, ncol = d, byrow = TRUE) + mu_matrix) / (n_obs + 1)
  
  log_A_list <- log(alpha_weights) - (n_obs * d / 2) * log(2 * pi * sigma_scalar) - 
    (d / 2) * log(n_obs + 1) -
    0.5 * (sum_Y_sq + rowSums(mu_matrix^2)) / sigma_scalar +
    0.5 * rowSums((n_obs * matrix(Y_bar, K, d, byrow = TRUE) + mu_matrix)^2) / (n_obs + 1) / sigma_scalar
  
  m <- max(log_A_list)
  A_list <- exp(log_A_list - m)
  w_list <- A_list / sum(A_list)
  log_exact <- m + log(sum(A_list))
  
  list(mn_matrix = mn_matrix, w_list = w_list, log_exact = log_exact)
}

# ============================================================
# Single replicate function
# ============================================================
run_single_Y <- function(seed_offset, N, mu_matrix, alpha_weights, sigma_scalar, 
                         n_obs, d, K, S_n_scalar, Sigma) {
  set.seed(123 + seed_offset)
  
  Y <- rmvnorm(n_obs, mean = rep(0, d), sigma = Sigma)
  post_params <- compute_posterior_params(Y, mu_matrix, alpha_weights, sigma_scalar, n_obs, d, K)
  samples <- generate_posterior_samples(N, post_params$mn_matrix, S_n_scalar, post_params$w_list, K, d)
  lps <- log_post_vec(samples, Y, sigma_scalar, alpha_weights, mu_matrix)
  
  list(samples = samples, lps = lps, Y = Y, log_exact = post_params$log_exact)
}

# ============================================================
# Main execution
# ============================================================
n_cores <- max(1, detectCores() - 1)

if (.Platform$OS.type == "unix") {
  results <- mclapply(1:M, function(m) {
    run_single_Y(m, N, mu_matrix, alpha_weights, sigma_scalar,
                 n_obs, d, K, S_n_scalar, Sigma)
  }, mc.cores = n_cores)
} else {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("N", "mu_matrix", "alpha_weights", "sigma_scalar",
                      "n_obs", "d", "K", "S_n_scalar", "Sigma",
                      "generate_posterior_samples", "log_post_vec",
                      "compute_posterior_params", "run_single_Y", "rowLogSumExps"))
  clusterEvalQ(cl, library(mvtnorm))
  results <- parLapply(cl, 1:M, function(m) {
    run_single_Y(m, N, mu_matrix, alpha_weights, sigma_scalar,
                 n_obs, d, K, S_n_scalar, Sigma)
  })
  stopCluster(cl)
}

# Extract results
all_samples <- lapply(results, function(x) x$samples)
all_lps <- lapply(results, function(x) x$lps)
all_Y <- lapply(results, function(x) x$Y)
all_log_exact <- sapply(results, function(x) x$log_exact)

# Save results
gmm_2d_samples <- list(
  samples = all_samples,
  lps = all_lps,
  Y = all_Y,
  log_exact = all_log_exact,
  n = n_obs,
  d = d,
  N = N,
  M = M,
  K = K,
  mu_matrix = mu_matrix,
  alpha_weights = alpha_weights,
  sigma_scalar = sigma_scalar
)

save(gmm_2d_samples, file = sprintf("gmm%d_2d_multiple_Yn%d.RData", K, n_obs))



# ============================================================
# Comparison
# ============================================================
library(dplyr)
library(mvtnorm)

# Source algorithm files
source("epwk.R")
source("thms.R")
source("mixthms.R")
source("ecmle.R")

# ============================================================
# Configuration
# ============================================================
# Load pre-generated samples (uncomment one)
load("gmm8_2d_multiple_Yn20.RData")
# load("gmm8_2d_multiple_Yn50.RData")
# load("gmm8_2d_multiple_Yn100.RData")

set.seed(123)

# ============================================================
# Extract Data
# ============================================================
all_samples <- gmm_2d_samples$samples
all_lps <- gmm_2d_samples$lps
all_Y <- gmm_2d_samples$Y
all_log_exact <- gmm_2d_samples$log_exact

n <- gmm_2d_samples$n
d <- gmm_2d_samples$d
N <- gmm_2d_samples$N
M <- gmm_2d_samples$M
K <- gmm_2d_samples$K
mu_matrix <- gmm_2d_samples$mu_matrix
alpha_weights <- gmm_2d_samples$alpha_weights
sigma_scalar <- gmm_2d_samples$sigma_scalar

# ============================================================
# Log Posterior Function
# ============================================================
make_log_post_fn <- function(Y_current) {
  function(theta) {
    theta <- as.numeric(theta)
    n_obs <- nrow(Y_current)
    
    diff_Y <- sweep(Y_current, 2, theta)
    sq_dist <- rowSums(diff_Y^2)
    log_likelihood <- -0.5 * n_obs * d * log(2 * pi * sigma_scalar) - 
      0.5 * sum(sq_dist) / sigma_scalar
    
    diff_mu <- sweep(mu_matrix, 2, theta)
    sq_dist_mu <- rowSums(diff_mu^2)
    log_priors <- log(alpha_weights) - 0.5 * d * log(2 * pi * sigma_scalar) - 
      0.5 * sq_dist_mu / sigma_scalar
    
    log_likelihood + logSumExp(log_priors)
  }
}


# ePWK parameters
ncutslice_epwk      <- 100
nslice_angle_epwk   <- 8


# ============================================================
# Warm-up for Sample Size Calibration
# ============================================================
target_time <- 1
n0 <- 1e4
max_warmup <- 100
tol <- 0.005

cur_sizes <- list(
  thames    = n0,
  mixthames = n0,
  ecmle     = n0,
  epwk      = n0
)

sample_warmup <- all_samples[[1]]
lps_warmup <- all_lps[[1]]
Y_warmup <- all_Y[[1]]
log_post_fn_warmup <- make_log_post_fn(Y_warmup)

for (it in seq_len(max_warmup)) {
  n_max_it <- max(unlist(cur_sizes))
  idx_max <- seq_len(min(n_max_it, N))
  post_full <- sample_warmup[idx_max, , drop = FALSE]
  lps_full <- lps_warmup[idx_max]
  
  times <- numeric(4)
  names(times) <- names(cur_sizes)
  
  # THAMES
  t1 <- Sys.time()
  idx <- seq_len(min(cur_sizes$thames, length(lps_full)))
  tryCatch({
    pdf(NULL)
    thames(lps_full[idx], post_full[idx, , drop = FALSE])
    dev.off()
  }, error = function(e) dev.off())
  times["thames"] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # MixTHAMES
  t1 <- Sys.time()
  idx <- seq_len(min(cur_sizes$mixthames, length(lps_full)))
  tryCatch({
    mixthames(post_full[idx, , drop = FALSE], lps_full[idx], log_post_fn_warmup)
  }, error = function(e) NULL)
  times["mixthames"] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # ECMLE
  t1 <- Sys.time()
  idx <- seq_len(min(cur_sizes$ecmle, length(lps_full)))
  tryCatch({
    ecmle(post_full[idx, , drop = FALSE], lps_full[idx], 
          hpd_level = 0.75, log_post_fn_warmup, make_plot = FALSE)
  }, error = function(e) NULL)
  times["ecmle"] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # ePWK
  t1 <- Sys.time()
  idx <- seq_len(min(cur_sizes$epwk, length(lps_full)))
  tryCatch({
    run_epwk(post_full[idx, , drop = FALSE], log_post_fn_warmup,
             ncutslice = ncutslice_epwk, nslice_angle = nslice_angle_epwk)
  }, error = function(e) NULL)
  times["epwk"] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  if (all(abs(times - target_time) <= tol * target_time)) break
  
  cur_sizes <- mapply(function(ncur, tcur) {
    nnext <- round(ncur * target_time / max(tcur, 1e-6))
    pmax(500, pmin(nnext, N))
  }, cur_sizes, times, SIMPLIFY = FALSE)
}

n_samples_thames    <- cur_sizes$thames
n_samples_mixthames <- cur_sizes$mixthames
n_samples_ecmle     <- cur_sizes$ecmle
n_samples_epwk      <- cur_sizes$epwk

# ============================================================
# Main Comparison
# ============================================================
thames_results <- numeric(M)
mixturethames_results <- numeric(M)
ecmle_results <- numeric(M)
epwk_results <- numeric(M)
time_thms <- numeric(M)
time_mixthms <- numeric(M)
time_ecmle <- numeric(M)
time_epwk <- numeric(M)

for (m in seq_len(M)) {
  sample_posterior <- all_samples[[m]]
  lps <- all_lps[[m]]
  Y_current <- all_Y[[m]]
  log_post_fn <- make_log_post_fn(Y_current)
  
  # THAMES
  t1 <- Sys.time()
  tryCatch({
    idx <- 1:min(n_samples_thames, N)
    pdf(NULL)
    thams <- thames(lps[idx], sample_posterior[idx, , drop = FALSE])
    dev.off()
    thames_results[m] <- thams$log_zhat_inv
  }, error = function(e) { dev.off(); thames_results[m] <<- NA })
  time_thms[m] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # MixTHAMES
  t1 <- Sys.time()
  tryCatch({
    idx <- 1:min(n_samples_mixthames, N)
    mixthams <- mixthames(sample_posterior[idx, , drop = FALSE], lps[idx], log_post_fn)
    mixturethames_results[m] <- mixthams$log_marginal_likelihood
  }, error = function(e) { mixturethames_results[m] <<- NA })
  time_mixthms[m] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # ECMLE
  t1 <- Sys.time()
  tryCatch({
    idx <- 1:min(n_samples_ecmle, N)
    emle <- ecmle(sample_posterior[idx, , drop = FALSE], lps[idx], 
                  hpd_level = 0.75, log_post_fn, make_plot = FALSE)
    ecmle_results[m] <- emle$log_marginal_likelihood
  }, error = function(e) { ecmle_results[m] <<- NA })
  time_ecmle[m] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  
  # ePWK
  t1 <- Sys.time()
  tryCatch({
    idx <- 1:min(n_samples_epwk, N)
    epwk_results[m] <- run_epwk(sample_posterior[idx, , drop = FALSE], log_post_fn,
                                ncutslice = ncutslice_epwk, nslice_angle = nslice_angle_epwk)
  }, error = function(e) { epwk_results[m] <<- NA })
  time_epwk[m] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
}

# ============================================================
# Compute Differences
# ============================================================
thames_diff <- thames_results - all_log_exact
mixturethames_diff <- mixturethames_results - all_log_exact
ecmle_diff <- ecmle_results - all_log_exact
epwk_diff <- epwk_results - all_log_exact

# ============================================================
# Save Results
# ============================================================
gmm_2d_multiple_Y_results <- list(
  thames = thames_results,
  mixturethames = mixturethames_results,
  ecmle = ecmle_results,
  epwk = epwk_results,
  thames_diff = thames_diff,
  mixturethames_diff = mixturethames_diff,
  ecmle_diff = ecmle_diff,
  epwk_diff = epwk_diff,
  time_thms = time_thms,
  time_mixthms = time_mixthms,
  time_ecmle = time_ecmle,
  time_epwk = time_epwk,
  log_exact = all_log_exact,
  n_samples_thames = n_samples_thames,
  n_samples_mixthames = n_samples_mixthames,
  n_samples_ecmle = n_samples_ecmle,
  n_samples_epwk = n_samples_epwk,
  d = d, K = K, M = M, N = N
)

save(gmm_2d_multiple_Y_results, file = "gmm82_2d_multiple_Y_n20.RData")
# save(gmm_2d_multiple_Y_results, file = "gmm82_2d_multiple_Y_n50.RData")
# save(gmm_2d_multiple_Y_results, file = "gmm82_2d_multiple_Y_n100.RData")