# Implementation matching author's code structure exactly (for G=1 case)

if (!require("uniformly", character.only = TRUE)) {
  install.packages("uniformly")
  library(uniformly, character.only = TRUE)
}

if (!require("sparsediscrim", character.only = TRUE)) {
  install.packages("sparsediscrim")
  library(sparsediscrim, character.only = TRUE)
}

#' This function computes quantiles of the truncated normal distribution
#'    for calculating THAMES confidence intervals.
#'
#' @param p Percentile
#' @param ratio Ratio of standard error to point estimate
#' @return Truncated normal quantile
trunc_quantile = function(p, ratio) {
  alpha = -1 / ratio
  (qnorm(p = pnorm(alpha) + p * (1 - pnorm(alpha))) * (ratio)) + 1
}

#' Find optimal HPD cutoff using chi-squared fitting
#' Computes alpha based on the Kolmogorov-Smirnov test statistic
chisq_find_limit = function(lps, d_par) {
  mean_chisq = d_par
  sd_chisq = sqrt(2 * length(lps))  
  neg_lps = -lps
  
  thinned_sequence = seq(2, floor((length(neg_lps) - 1) * 0.8), by = 100)
  limit_i = sort(neg_lps, decreasing = TRUE)[thinned_sequence]
  neg_lps_trunc = t(matrix(rep(neg_lps, length(limit_i)), ncol = length(limit_i)))
  
  sum_dummy = neg_lps_trunc * ((t(matrix(rep(neg_lps, length(limit_i)), ncol = length(limit_i))) <=
                                  t(matrix(rep(limit_i, each = length(neg_lps)), ncol = length(limit_i)))))
  
  neg_lps_trunc[!(t(matrix(rep(neg_lps, length(limit_i)), ncol = length(limit_i))) <=
                    t(matrix(rep(limit_i, each = length(neg_lps)), ncol = length(limit_i))))] = NA
  
  neg_lps_trunc_mu_hat = rowSums(sum_dummy) / rowSums(!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_hat = sqrt(rowSums((!is.na(neg_lps_trunc)) * (sum_dummy - neg_lps_trunc_mu_hat)^2) /
                                   (rowSums(!is.na(neg_lps_trunc)) - 1))
  neg_lps_trunc_sigma_hat[neg_lps_trunc_sigma_hat == 0] = 1e-6  
  neg_lps_trunc_mu_tilde = (neg_lps_trunc_mu_hat - mean_chisq) * (!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_tilde = (neg_lps_trunc_sigma_hat / sd_chisq)
  neg_lps_norm_chisq = (((neg_lps_trunc - neg_lps_trunc_mu_hat) / neg_lps_trunc_sigma_tilde) + mean_chisq)
  
  sort_non_na = function(x) {
    x[!is.na(x)] = rank(x[!is.na(x)])
    return(x)
  }
  
  cdf_chisq_mat = matrix(pchisq(c(neg_lps_norm_chisq), df = mean_chisq) /
                           pchisq(max(neg_lps_norm_chisq[!is.na(neg_lps_norm_chisq)]), df = mean_chisq),
                         ncol = ncol(neg_lps_norm_chisq))
  ecdf_mat = t(apply(neg_lps_norm_chisq, 1, function(x) sort_non_na(x))) / rowSums(!is.na(neg_lps_trunc))
  kolm_dist = apply(abs(cdf_chisq_mat - ecdf_mat), 1, function(x) max(x[!is.na(x)]))
  
  opt_index = thinned_sequence[which.min(kolm_dist)]
  limit = -sort(neg_lps, decreasing = TRUE)[opt_index]
  return(limit)
}

#' Compute ellipsoid parameters  type="standard"
#' 
compute_ellipse = function(params, iters) {
  # compute sample mean and covariance, take c_opt=sqrt(d+1)
  theta_hat <- colMeans(params[1:iters, ])
  sigma_hat <- cov(params[1:iters, ])
  c_opt = sqrt(length(theta_hat) + 1)
  
  ellipse = list(
    theta_hat = theta_hat,
    sigma_hat = sigma_hat,
    c_opt = c_opt
  )
  return(ellipse)
}

#' Main THAMES computation for G=1

#' @param params Matrix of posterior samples (2*iters x d)
#' @param lps Vector of log-posterior values (length 2*iters)
#' @param logpost Function to evaluate log-posterior at new points
#' @param iters Number of iterations per half (total samples = 2*iters)
#' @param seed Random seed for reproducibility
#' @return List containing log marginal likelihood estimate and diagnostics
compute_thames = function(params, lps, logpost, iters, seed = 2024) {
  
  G = 1  # This implementation is for G=1 only
  
  # Compute ellipsoid (type="standard")
  ellipse = compute_ellipse(params, iters)
  theta_hat = ellipse$theta_hat
  sigma_hat = ellipse$sigma_hat
  c_opt = ellipse$c_opt
  d_par = length(theta_hat)
  
  # Compute limit using first half
  limit = chisq_find_limit(lps[1:iters], d_par = d_par)
  limit = max(limit, median(lps[1:iters]))  
  
  n_simuls = 2 * iters
  
  # Compute mu_post and sigma_post from first half
  mu_post = colMeans(params[1:iters, ])
  sigma_post = cov(params[1:iters, ])
  
  # Inverse of posterior covariance 
  inv_post_var = sparsediscrim::solve_chol(sigma_post)
  
  set.seed(seed)
  counter = 0
  log_cor = -Inf
  
  in_ellipse = TRUE
  c_opt_old = c_opt
  
  # ============================================================
  # ADAPTIVE LOOP (EXACT copy from author's compute_thames)
  # ============================================================
  
  while (is.infinite(log_cor) & in_ellipse) {
    # Strategy 1: Shrink radius
    c_opt = c_opt_old / 2^counter
    print(c_opt)
    
    # Sample uniformly from ellipsoid
    param_test = runif_in_ellipsoid(n_simuls, inv_post_var, c_opt) +
      t(matrix(rep(mu_post, n_simuls), ncol = n_simuls))
    
    # Evaluate log-posterior at test points
    lps_test = logpost(param_test, G)
    
    # Compute log correction factor
    log_cor = log(mean(lps_test > limit))
    counter = counter + 1
    
    # Check if ellipse is empty
    params_centered = params[(iters + 1):(2 * iters), ] - 
      t(matrix(rep(theta_hat, iters), ncol = iters))
    radius = c_opt
    
    # Check if ellipse empty; set mean to mode of second half if not
    in_E = rowSums((params_centered %*% sparsediscrim::solve_chol(ellipse$sigma_hat)) * params_centered) <= radius^2
    empty_ellipse = (sum(in_E) == 0)
    
    if (empty_ellipse) {
      # Strategy 2: Recenter to MAP
      theta_hat = params[(iters + 1):(2 * iters), ][which.max(lps[(iters + 1):(2 * iters)]), ]
      mu_post = theta_hat
      ellipse$theta_hat = theta_hat
      c_opt = sqrt(ncol(params) + 1)
      log_cor = -Inf
      counter = 0
    }
  }
  
  # ============================================================
  # THAMES CALCULATION (Based on thames_mixture_simple for G=1)
  # ============================================================
  
  n_samples = iters  # second half
  radius = c_opt
  
  # Calculate SVD of sigma_hat 
  sigma_svd = svd(sigma_hat)
  
  # Calculate log(det(sigma_hat))
  log_det_sigma_hat = sum(log(sigma_svd$d))
  
  # Calculate volume of A 
  logvolA = d_par * log(radius) + (d_par / 2) * log(pi) + log_det_sigma_hat / 2 - lgamma(d_par / 2 + 1)
  
  # Which samples are in A? 
  num_inA = rep(0, n_samples)
  
  # For G=1, no permutations needed
  theta_hat_extended = theta_hat
  inv_sigma_hat_extended = solve(sigma_hat)
  
  # Use second half of params
  params_second_half = params[(iters + 1):(2 * iters), ]
  lps_second_half = lps[(iters + 1):(2 * iters)]
  
  # Check if included in A
  params_centered = params_second_half - t(matrix(rep(theta_hat_extended, n_samples), ncol = n_samples))
  num_inA = num_inA +
    (rowSums((params_centered %*% inv_sigma_hat_extended) * params_centered) <= radius^2)
  
  # Calculate zhat 
  # For G=1 with type="standard"
  log_zhat_inv = log(mean(exp(-(lps_second_half - max(lps_second_half))) * num_inA * (lps_second_half > limit))) -
    logvolA - max(lps_second_half)
  
  # Estimate ar(1) model for lps 
  lp_ar <- ar(exp(-(lps_second_half - max(lps_second_half))) * num_inA * (lps_second_half > limit), order.max = 1)
  phi <- lp_ar$partialacf[1]
  if (is.na(phi)) phi <- 0
  
  standard_error <- sd(exp(-lps_second_half + max(lps_second_half)) * num_inA * (lps_second_half > limit)) /
    ((1 - phi) * sqrt(n_samples) * mean(exp(-lps_second_half + max(lps_second_half)) * num_inA * (lps_second_half > limit)))
  
  # Calculate 95% lower bound 
  log_zhat_inv_L <- log_zhat_inv + log(trunc_quantile(0.025, standard_error))
  
  # Calculate 95% upper bound 
  log_zhat_inv_U <- log_zhat_inv + log(trunc_quantile(0.975, standard_error))
  
  return(list(
    # Main results 
    log_zhat_inv_L = log_zhat_inv_L - log_cor,
    log_zhat_inv = log_zhat_inv - log_cor,
    log_zhat_inv_U = log_zhat_inv_U - log_cor,
    log_cor = log_cor,
    len_perms = 1,
    alpha = mean(-lps < -limit),
    c_opt = c_opt,
    
    # Additional diagnostics
    theta_hat = theta_hat,
    sigma_hat = sigma_hat,
    log_det_sigma_hat = log_det_sigma_hat,
    logvolA = logvolA,
    num_inA = num_inA,
    se = standard_error,
    phi = phi,
    radius = radius,
    d_par = d_par,
    G = G,
    limit = limit
  ))
}


#' @param sample_posterior Matrix of posterior samples (N x d)
#' @param lps Vector of log-posterior values (length N)
#' @param log_post_fn Function to evaluate log-posterior at new points (must accept matrix and G)
#' @return List containing log marginal likelihood estimate and diagnostics
mixthames <- function(sample_posterior, lps, log_post_fn) {
  N <- nrow(sample_posterior)
  iters <- floor(N / 2)
  
  # Wrapper for logpost to match author's interface: logpost(param_matrix, G)
  logpost_wrapper <- function(param_matrix, G) {
    apply(param_matrix, 1, log_post_fn)
  }
  
  # Call compute_thames
  result <- compute_thames(
    params = sample_posterior,
    lps = lps,
    logpost = logpost_wrapper,
    iters = iters
  )
  
  return(list(
    log_marginal_likelihood = -result$log_zhat_inv,
    log_marginal_likelihood_lower = -result$log_zhat_inv_U,
    log_marginal_likelihood_upper = -result$log_zhat_inv_L,
    standard_error = result$se,
    phi = result$phi,
    coverage_rate = sum(result$num_inA) / length(result$num_inA),
    log_cor = result$log_cor,
    c_opt = result$c_opt,
    radius = result$radius,
    limit = result$limit,
    alpha = result$alpha,
    theta_hat = result$theta_hat,
    logvolA = result$logvolA
  ))
}