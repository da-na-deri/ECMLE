
if (!require("pracma", character.only = TRUE)) {
  install.packages("pracma")
  library(pracma, character.only = TRUE)
}

# Euclidean norm
euclid_norm <- function(v) {
  sqrt(sum(v^2))
}

# Bisection to level
bisect_to_level <- function(f, lo, hi, target, tol = 1e-4, maxit = 10) {
  flo <- f(lo) - target
  fhi <- f(hi) - target
  if (flo * fhi > 0) return(NA_real_)
  for (it in 1:maxit) {
    mid <- 0.5 * (lo + hi)
    fm <- f(mid) - target
    if (abs(fm) < tol) return(mid)
    if (flo * fm <= 0) { hi <- mid; fhi <- fm } else { lo <- mid; flo <- fm }
  }
  0.5 * (lo + hi)
}

# Compute radius along a single direction
compute_radius_along_dir <- function(center, u, log_post_fn, c_thresh, r_max = 1e3, bisect_tol = 1e-4) {
  fdir <- function(r) log_post_fn(center + r * u)
  bisect_to_level(f = fdir, lo = 0, hi = r_max, target = c_thresh, tol = bisect_tol)
}

# Euclidean distances for multiple points
euclidean_m <- function(x_matrix, centers_matrix) {
  if (!is.matrix(x_matrix)) x_matrix <- as.matrix(x_matrix)
  if (!is.matrix(centers_matrix)) centers_matrix <- as.matrix(centers_matrix)
  n_points <- nrow(x_matrix); n_centers <- nrow(centers_matrix)
  distances <- matrix(0, nrow = n_points, ncol = n_centers)
  for (i in 1:n_centers) {
    diff_matrix <- sweep(x_matrix, 2, centers_matrix[i, ], "-")
    distances[, i] <- sqrt(rowSums(diff_matrix^2))
  }
  distances
}

# Vectorized check for points in ellipsoids
is_point_in_any_ellipsoid <- function(points_matrix, centers_matrix, Sigmas_inv_list) {
  n_points <- nrow(points_matrix)
  n_ellipsoids <- nrow(centers_matrix)
  result <- logical(n_points)
  for (i in 1:n_ellipsoids) {
    diff_matrix <- sweep(points_matrix, 2, centers_matrix[i, ], "-")
    temp_matrix <- diff_matrix %*% Sigmas_inv_list[[i]]
    mahal_distances <- rowSums(temp_matrix * diff_matrix)
    result <- result | (mahal_distances <= 1)
  }
  return(result)
}

# Compute ellipsoid volume
compute_ellipsoid_volume <- function(Sigma, d) {
  unit_vol <- pi^(d/2) / gamma(d/2 + 1)
  unit_vol * det(Sigma)^0.5
}

# log marginal likelihood computation
compute_log_marginal_likelihood <- function(log_ratio, n_samples) {
  log_marginal_likelihood_iter <- numeric(n_samples)
  if (n_samples > 0) {
    log_marginal_likelihood_iter[1] <- -log_ratio[1]
    
    # For subsequent iterations, use incremental log-sum-exp
    current_log_sum <- log_ratio[1]
    
    for (it in 2:n_samples) {
      # Add new term to running sum using numerically stable log-sum-exp
      new_val <- log_ratio[it]
      if (is.finite(current_log_sum) && is.finite(new_val)) {
        if (current_log_sum > new_val) {
          current_log_sum <- current_log_sum + log1p(exp(new_val - current_log_sum))
        } else {
          current_log_sum <- new_val + log1p(exp(current_log_sum - new_val))
        }
      } else if (is.finite(new_val)) {
        current_log_sum <- new_val
      }
      
      # Compute log harmonic mean
      log_harmonic_mean_sum <- current_log_sum - log(it)
      log_marginal_likelihood_iter[it] <- -log_harmonic_mean_sum
    }
  }
  
  return(log_marginal_likelihood_iter)
}