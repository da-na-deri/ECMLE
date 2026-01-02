source("scripts/help_func_ecmle.R")

ecmle <- function(post_samples, lps, hpd_level = 0.75, log_post_fn, make_plot = TRUE) {
  d <- ncol(post_samples)
  N <- nrow(post_samples)
  half <- floor(N/2)
  params1 <- post_samples[1:half, , drop = FALSE]
  lps1    <- lps[1:half]
  params2 <- post_samples[(half+1):N, , drop = FALSE]
  lps2    <- lps[(half+1):N]
  
  # HPD region1
  c <- quantile(lps1, probs = 1 - hpd_level)
  in_hpd     <- lps1 >= c
  HPD_samples_ball   <- params1[in_hpd,    , drop = FALSE]
  HPD_lps            <- lps1[in_hpd]
  Non_HPD_sample_ball<- params1[!in_hpd,   , drop = FALSE]
  
  index <- dim(HPD_samples_ball)[1]
  size <- floor(0.05*index)  
  hpd_index <- sample(1:index, size= size, replace = FALSE)
  HPD_samples_ball   <- HPD_samples_ball[hpd_index,]
  HPD_lps            <- HPD_lps[hpd_index]
  
  HPD_matrix <- as.matrix(HPD_samples_ball)
  Non_HPD_matrix <- as.matrix(Non_HPD_sample_ball)
  
  # Compute r_max using euclidean_m
  dists <- euclidean_m(HPD_samples_ball, HPD_samples_ball)
  r_max <- max(dists, na.rm = TRUE) 
  
  # Sort HPD samples by log-posterior
  sorted_idx <- order(-HPD_lps)
  HPD_matrix_sorted <- HPD_matrix[sorted_idx, , drop = FALSE]
  HPD_lps_sorted <- HPD_lps[sorted_idx]
  
  # Ellipsoid construction
  ellipsoid_centers_list <- list()
  ellipsoid_Sigmas_list <- list()
  ellipsoid_Sigmas_inv_list <- list()
  ellipsoid_max_semi_list <- numeric(0)
  ellipsoid_volumes <- numeric(0)
  ellipsoid_lps <- numeric(0)
  ellipsoid_count <- 0
  
  available_centers <- rep(TRUE, nrow(HPD_matrix_sorted))
  
  # Main packing loop
  while (any(available_centers)) {
    next_idx <- which(available_centers)[1]
    current_center <- HPD_matrix_sorted[next_idx, ]
    current_lps <- HPD_lps_sorted[next_idx]
    
    if (nrow(Non_HPD_matrix) > 0) {
      euclid_dists <- euclidean_m(Non_HPD_matrix, matrix(current_center, nrow=1))
      closest_idx <- which.min(euclid_dists)
      closest_point <- Non_HPD_matrix[closest_idx, ]
      dir_vec <- closest_point - current_center
    } else {
      dir_vec <- c(1, rep(0, d-1))
    }
    u1 <- dir_vec / euclid_norm(dir_vec)
    # Compute r1 along u1
    r1 <- compute_radius_along_dir(
      current_center, u1, log_post_fn, c,
      r_max = r_max, bisect_tol = 1e-4
    )
    # Check if r1 is valid
    if (is.na(r1) || r1 <= 0) {
      available_centers[next_idx] <- FALSE
      next
    }
    
    # Compute orthogonal basis using pracma::gramSchmidt
    A <- matrix(0, nrow = d, ncol = d)
    A[,1] <- u1 / euclid_norm(u1)  # Normalize u1
    for (i in 2:d) {
      A[i,i] <- 1  # Standard basis vectors e_2, ..., e_d
    }
    gs <- pracma::gramSchmidt(A, tol = .Machine$double.eps^0.5)
    U <- gs$Q
    
    # Compute ri for i=1:d
    semi_axes <- numeric(d)
    semi_axes[1] <- r1
    
    for (i in 2:d) {
      ui <- U[,i]
      r_pos <- compute_radius_along_dir(
        current_center, ui, log_post_fn, c,
        r_max = r_max
      )
      r_neg <- compute_radius_along_dir(
        current_center, -ui, log_post_fn, c,
        r_max = r_max
      )
      if (is.na(r_pos) || is.na(r_neg)) {
        semi_axes[i] <- NA
      } else {
        semi_axes[i] <- min(r_pos, r_neg)
      }
    }
    
    # Check if semi_axes are valid
    if (any(is.na(semi_axes)) || any(semi_axes <= 0)) {
      available_centers[next_idx] <- FALSE
      next
    }
    
    # Build initial Sigma
    D <- diag(semi_axes^2)
    Sigma <- U %*% D %*% t(U)
    max_semi <- max(semi_axes)
    
    # Check overlaps and skip if any overlap is found
    overlaps <- FALSE
    if (ellipsoid_count > 0) {
      existing_centers_matrix <- do.call(rbind, ellipsoid_centers_list)
      for (j in 1:ellipsoid_count) {
        existing_center <- existing_centers_matrix[j, ]
        center_distance <- euclid_norm(current_center - existing_center)
        max_semi_j <- ellipsoid_max_semi_list[j]
        if (center_distance < max_semi + max_semi_j) {
          overlaps <- TRUE
          break
        }
      }
    }
    
    if (overlaps) {
      available_centers[next_idx] <- FALSE
      next
    }
    
    # No overlap: add ellipsoid
    volume <- compute_ellipsoid_volume(Sigma, d)
    ellipsoid_count <- ellipsoid_count + 1
    ellipsoid_centers_list[[ellipsoid_count]] <- current_center
    ellipsoid_Sigmas_list[[ellipsoid_count]] <- Sigma
    ellipsoid_Sigmas_inv_list[[ellipsoid_count]] <- solve(Sigma)
    ellipsoid_max_semi_list[ellipsoid_count] <- max_semi
    ellipsoid_volumes[ellipsoid_count] <- volume
    ellipsoid_lps[ellipsoid_count] <- current_lps
    
    # Prune remaining centers inside the new ellipsoid
    remaining_indices <- which(available_centers)
    if (length(remaining_indices) > 0) {
      remaining_centers_matrix <- HPD_matrix_sorted[remaining_indices, , drop = FALSE]
      is_inside_new <- is_point_in_any_ellipsoid(
        remaining_centers_matrix, 
        matrix(current_center, nrow = 1), 
        list(ellipsoid_Sigmas_inv_list[[ellipsoid_count]])
      )
      available_centers[remaining_indices[is_inside_new]] <- FALSE
    }
    
    available_centers[next_idx] <- FALSE
  }
  
  # Total volume
  total_volume <- sum(ellipsoid_volumes)
  
  ellipsoids <- list()
  for (i in 1:ellipsoid_count) {
    ellipsoids[[i]] <- list(
      center = ellipsoid_centers_list[[i]],
      Sigma = ellipsoid_Sigmas_list[[i]],
      volume = ellipsoid_volumes[i],
      lps = ellipsoid_lps[i]
    )
  }
  
  if (ellipsoid_count == 0) {
    return(list(
      log_marginal_likelihood = -Inf,
      ellipsoids = ellipsoids,
      total_volume = 0,
      n_ellipsoids = 0,
      points_in_ellipsoids = 0,
      n_samples = nrow(params2),
      coverage_rate = 0
    ))
  }
  
  ellipsoid_centers_matrix <- do.call(rbind, ellipsoid_centers_list)
  posterior_points <- as.matrix(params2)
  n_samples <- nrow(posterior_points)
  
  samples_in_ellipsoids <- is_point_in_any_ellipsoid(
    posterior_points, 
    ellipsoid_centers_matrix, 
    ellipsoid_Sigmas_inv_list
  )
  
  samples_in_ellipsoids_numeric <- as.numeric(samples_in_ellipsoids)
  
  log_ratio <- log(samples_in_ellipsoids_numeric) - log(total_volume) - lps2
  
  log_marginal_likelihood_iter <- compute_log_marginal_likelihood(log_ratio, n_samples)
  
  final_log_marginal_likelihood <- log_marginal_likelihood_iter[n_samples]
  
  return(list(
    log_marginal_likelihood = final_log_marginal_likelihood,
    log_marginal_likelihood_iter = log_marginal_likelihood_iter,
    ellipsoids = ellipsoids,
    total_volume = total_volume,
    n_ellipsoids = ellipsoid_count,
    points_in_ellipsoids = sum(samples_in_ellipsoids_numeric),
    n_samples = n_samples,
    coverage_rate = sum(samples_in_ellipsoids_numeric) / n_samples
  ))
}