# ePWK (Extended Piecewise Kernel) for Marginal Likelihood Estimation
# Matches author's original approach, generalized for arbitrary dimension d

# ============================================================
# Helper Functions
# ============================================================
logSumExp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

denominator_control <- function(ratios, T_size) {
  b <- ratios[is.finite(ratios)]
  if (length(b) == 0) return(NA)
  b_max <- max(b)
  log(sum(exp(b - b_max))) + b_max - log(T_size)
}

# ============================================================
# Coordinate Conversion (for d > 2)
# ============================================================
cartesian_to_hyperspherical <- function(x) {
  d <- length(x)
  r <- sqrt(sum(x^2))
  
  if (r < .Machine$double.eps) {
    return(list(r = 0, angles = rep(0, d - 1)))
  }
  
  angles <- numeric(d - 1)
  
  for (i in 1:(d - 1)) {
    denom <- sqrt(sum(x[i:d]^2))
    if (denom < .Machine$double.eps) {
      angles[i] <- 0
    } else {
      angles[i] <- acos(pmin(pmax(x[i] / denom, -1), 1))
    }
  }
  
  # Adjust last angle to [0, 2*pi] based on sign of last coordinate
  if (d >= 2 && x[d] < 0) {
    angles[d - 1] <- 2 * pi - angles[d - 1]
  }
  
  return(list(r = r, angles = angles))
}

hyperspherical_to_cartesian <- function(r, angles) {
  d <- length(angles) + 1
  x <- numeric(d)
  
  sin_prod <- r
  for (i in 1:(d - 1)) {
    x[i] <- sin_prod * cos(angles[i])
    sin_prod <- sin_prod * sin(angles[i])
  }
  x[d] <- sin_prod
  
  return(x)
}

# ============================================================
# 2D Implementation (matches author exactly)
# ============================================================
LOR_partition_2d <- function(radius, ncutslice, nzone, temp_mean, log_post_fn) {
  angle <- seq(0, pi, 2 * pi / nzone)
  piece <- cbind(angle[-(nzone / 2 + 1)], angle[-1])
  piece <- rbind(piece, piece + pi)
  
  interval <- seq(0, radius, length.out = ncutslice + 1)
  rings <- cbind(interval[-(ncutslice + 1)], interval[-1])
  
  rings_piece <- NULL
  for (g in 1:nzone) {
    rings_piece_temp <- cbind(rings, piece[g, 1], piece[g, 2])
    rings_piece <- rbind(rings_piece, rings_piece_temp)
  }
  
  reprp1 <- rowMeans(rings_piece[, 1:2])
  reprp2 <- rowMeans(rings_piece[, 3:4])
  rings_piece <- cbind(rings_piece, reprp1, reprp2)
  
  # Single midpoint kernel evaluation 
  kreprp <- numeric(nrow(rings_piece))
  for (i in 1:nrow(rings_piece)) {
    rpoints <- temp_mean + reprp1[i] * c(cos(reprp2[i]), sin(reprp2[i]))
    kreprp[i] <- log_post_fn(rpoints)
  }
  rings_piece <- cbind(rings_piece, kreprp)
  
  # Volume with kernel included
  rarea <- pi * rings_piece[, 2]^2 - pi * rings_piece[, 1]^2
  rvol <- log(rarea) + kreprp - log(nzone)
  rings_piece <- cbind(rings_piece, rvol)
  
  colnames(rings_piece) <- c("r_lower", "r_upper", "angle_lower", "angle_upper",
                             "r_mid", "angle_mid", "kernel_log", "vol_log")
  return(rings_piece)
}

ePWK_2d <- function(samples, temp_mean, rings_piece, ncutslice, nzone, log_post_fn) {
  T1 <- nrow(samples)
  
  para_mcmc_v <- sweep(samples, 2, temp_mean)
  r_ <- sqrt(rowSums(para_mcmc_v^2))
  
  # Angle calculation 
  angle_ <- acos(pmin(pmax(para_mcmc_v[, 1] / r_, -1), 1))
  idx1 <- which(para_mcmc_v[, 1] > 0 & para_mcmc_v[, 2] < 0)
  idx2 <- which(para_mcmc_v[, 1] < 0 & para_mcmc_v[, 2] < 0)
  angle_[idx1] <- 2 * pi - angle_[idx1]
  angle_[idx2] <- pi + (pi - angle_[idx2])
  
  val <- rep(NA, T1)
  for (j in 1:T1) {
    # Direct indexing 
    ri1 <- which(r_[j] >= rings_piece[1:ncutslice, 1] & r_[j] < rings_piece[1:ncutslice, 2])
    ri2 <- which(angle_[j] >= rings_piece[ncutslice * (1:nzone), 3] & 
                   angle_[j] < rings_piece[ncutslice * (1:nzone), 4])
    ri <- ncutslice * (ri2 - 1) + ri1
    
    if (length(ri) > 0 && ri[1] > 0) {
      val[j] <- rings_piece[ri[1], 7] - log_post_fn(samples[j, ])
    }
  }
  
  lcons <- max(rings_piece[, 8])
  tcubevol <- log(sum(exp(rings_piece[, 8] - lcons))) + lcons
  den <- denominator_control(val, T1)
  
  return(tcubevol - den)
}

# ============================================================
# General d-dimensional Implementation 
# ============================================================
LOR_partition_nd <- function(radius, ncutslice, nslice_angle, temp_mean, log_post_fn, d) {
  if (length(nslice_angle) == 1) {
    nslice_angle <- rep(nslice_angle, d - 1)
  }
  
  # Radial intervals
  interval <- seq(0, radius, length.out = ncutslice + 1)
  rings <- cbind(interval[-(ncutslice + 1)], interval[-1])
  
  # Angular intervals: phi_1,...,phi_{d-2} in [0,pi], phi_{d-1} in [0,2*pi]
  angle_breaks <- vector("list", d - 1)
  for (k in 1:(d - 1)) {
    if (k < d - 1) {
      angle_breaks[[k]] <- seq(0, pi, length.out = nslice_angle[k] + 1)
    } else {
      angle_breaks[[k]] <- seq(0, 2 * pi, length.out = nslice_angle[k] + 1)
    }
  }
  
  # Build partition matrix
  n_partitions <- ncutslice * prod(nslice_angle)
  angle_indices <- expand.grid(lapply(nslice_angle, function(n) 1:n))
  
  n_cols <- 2 + 2 * (d - 1) + 1 + (d - 1) + 2
  rings_piece <- matrix(NA, nrow = n_partitions, ncol = n_cols)
  
  idx <- 0
  for (i in 1:ncutslice) {
    for (ai in 1:nrow(angle_indices)) {
      idx <- idx + 1
      
      r_lower <- rings[i, 1]
      r_upper <- rings[i, 2]
      r_mid <- (r_lower + r_upper) / 2
      
      angle_lower <- numeric(d - 1)
      angle_upper <- numeric(d - 1)
      angle_mid <- numeric(d - 1)
      
      for (k in 1:(d - 1)) {
        j_k <- angle_indices[ai, k]
        angle_lower[k] <- angle_breaks[[k]][j_k]
        angle_upper[k] <- angle_breaks[[k]][j_k + 1]
        angle_mid[k] <- (angle_lower[k] + angle_upper[k]) / 2
      }
      
      # Single midpoint kernel evaluation 
      rpoint <- temp_mean + hyperspherical_to_cartesian(r_mid, angle_mid)
      kreprp <- log_post_fn(rpoint)
      
      # Volume calculation (d-dimensional)
      # Shell volume: S_d * (r_outer^d - r_inner^d) / d, where S_d = 2*pi^(d/2)/Gamma(d/2)
      log_Sd <- log(2) + (d / 2) * log(pi) - lgamma(d / 2)
      if (r_lower == 0) {
        log_shell_vol <- log_Sd + d * log(r_upper) - log(d)
      } else {
        log_outer <- d * log(r_upper)
        log_inner <- d * log(r_lower)
        log_diff <- log_outer + log(1 - exp(log_inner - log_outer))
        log_shell_vol <- log_Sd + log_diff - log(d)
      }
      
      # Angular fraction
      n_angle_partitions <- prod(nslice_angle)
      
      # Volume with kernel included 
      rvol <- log_shell_vol + kreprp - log(n_angle_partitions)
      
      # Store
      col_idx <- 1
      rings_piece[idx, col_idx] <- r_lower; col_idx <- col_idx + 1
      rings_piece[idx, col_idx] <- r_upper; col_idx <- col_idx + 1
      for (k in 1:(d - 1)) {
        rings_piece[idx, col_idx] <- angle_lower[k]; col_idx <- col_idx + 1
      }
      for (k in 1:(d - 1)) {
        rings_piece[idx, col_idx] <- angle_upper[k]; col_idx <- col_idx + 1
      }
      rings_piece[idx, col_idx] <- r_mid; col_idx <- col_idx + 1
      for (k in 1:(d - 1)) {
        rings_piece[idx, col_idx] <- angle_mid[k]; col_idx <- col_idx + 1
      }
      rings_piece[idx, col_idx] <- kreprp; col_idx <- col_idx + 1
      rings_piece[idx, col_idx] <- rvol
    }
  }
  
  return(rings_piece)
}

ePWK_nd <- function(samples, temp_mean, rings_piece, ncutslice, nslice_angle, log_post_fn, d) {
  T1 <- nrow(samples)
  
  if (length(nslice_angle) == 1) {
    nslice_angle <- rep(nslice_angle, d - 1)
  }
  
  # Column indices
  col_r_lower <- 1
  col_r_upper <- 2
  col_angle_lower_start <- 3
  col_angle_upper_start <- 3 + (d - 1)
  col_kernel <- 2 + 2 * (d - 1) + 1 + (d - 1) + 1
  col_vol <- col_kernel + 1
  
  # Precompute angle breaks for lookup
  angle_breaks <- vector("list", d - 1)
  for (k in 1:(d - 1)) {
    if (k < d - 1) {
      angle_breaks[[k]] <- seq(0, pi, length.out = nslice_angle[k] + 1)
    } else {
      angle_breaks[[k]] <- seq(0, 2 * pi, length.out = nslice_angle[k] + 1)
    }
  }
  
  samples_centered <- sweep(samples, 2, temp_mean)
  
  val <- rep(NA, T1)
  for (j in 1:T1) {
    hs <- cartesian_to_hyperspherical(samples_centered[j, ])
    r_j <- hs$r
    angles_j <- hs$angles
    
    # Direct indexing
    # Find radial index
    ri1 <- which(r_j >= rings_piece[1:ncutslice, col_r_lower] & 
                   r_j < rings_piece[1:ncutslice, col_r_upper])
    
    if (length(ri1) == 0) next
    
    # Find angle indices for each dimension
    angle_idx <- numeric(d - 1)
    valid <- TRUE
    for (k in 1:(d - 1)) {
      ai <- which(angles_j[k] >= angle_breaks[[k]][-(nslice_angle[k] + 1)] & 
                    angles_j[k] < angle_breaks[[k]][-1])
      if (length(ai) == 0) {
        valid <- FALSE
        break
      }
      angle_idx[k] <- ai[1]
    }
    
    if (!valid) next
    
    # Compute linear index: ri = r_idx + ncutslice * (combined_angle_idx - 1)
    combined_angle_idx <- angle_idx[1]
    if (d > 2) {
      multiplier <- 1
      for (k in 2:(d - 1)) {
        multiplier <- multiplier * nslice_angle[k - 1]
        combined_angle_idx <- combined_angle_idx + (angle_idx[k] - 1) * multiplier
      }
    }
    ri <- ri1[1] + ncutslice * (combined_angle_idx - 1)
    
    if (ri > 0 && ri <= nrow(rings_piece)) {
      val[j] <- rings_piece[ri, col_kernel] - log_post_fn(samples[j, ])
    }
  }
  
  lcons <- max(rings_piece[, col_vol], na.rm = TRUE)
  tcubevol <- log(sum(exp(rings_piece[, col_vol] - lcons))) + lcons
  den <- denominator_control(val, T1)
  
  return(tcubevol - den)
}

# ============================================================
# Main Interface
# ============================================================
#' Run ePWK for marginal likelihood estimation
#'
#' @param samples Matrix of posterior samples (N x d)
#' @param log_post_fn Function that returns log posterior density
#' @param ncutslice Number of radial shells (default: 100)
#' @param nslice_angle Number of angular slices (default: 8 for 2D, 4 for higher d)
#' @param k0 Fraction of max radius to use (default: 0.95)
#' @return Log marginal likelihood estimate
run_epwk <- function(samples, log_post_fn, ncutslice = 100, nslice_angle = NULL, k0 = 0.95) {
  
  if (!is.matrix(samples)) {
    samples <- as.matrix(samples)
  }
  d <- ncol(samples)
  
  # Default nslice_angle
  if (is.null(nslice_angle)) {
    nslice_angle <- if (d == 2) 8 else 4
  }
  
  # Compute center and radius
  temp_mean <- colMeans(samples)
  samples_centered <- sweep(samples, 2, temp_mean)
  r_all <- sqrt(rowSums(samples_centered^2))
  radius <- k0 * max(r_all)
  
  if (d == 2) {
    nzone <- nslice_angle
    rings_piece <- LOR_partition_2d(radius, ncutslice, nzone, temp_mean, log_post_fn)
    result <- ePWK_2d(samples, temp_mean, rings_piece, ncutslice, nzone, log_post_fn)
  } else {
    # Use generalized d-dimensional implementation
    rings_piece <- LOR_partition_nd(radius, ncutslice, nslice_angle, temp_mean, log_post_fn, d)
    result <- ePWK_nd(samples, temp_mean, rings_piece, ncutslice, nslice_angle, log_post_fn, d)
  }
  
  return(result)
}