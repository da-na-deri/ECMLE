# Internal helpers ---------------------------------------------------------

.validate_numeric_matrix <- function(x, name) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || anyNA(x) || nrow(x) < 4L || ncol(x) < 1L) {
    stop(sprintf("`%s` must be a numeric matrix with at least 4 rows and 1 column and no missing values.", name), call. = FALSE)
  }
  x
}

.validate_numeric_vector <- function(x, name, length_out = NULL) {
  if (!is.numeric(x) || anyNA(x)) {
    stop(sprintf("`%s` must be a numeric vector with no missing values.", name), call. = FALSE)
  }
  if (!is.null(length_out) && length(x) != length_out) {
    stop(sprintf("`%s` must have length %d.", name, length_out), call. = FALSE)
  }
  x
}

.euclid_norm <- function(v) {
  sqrt(sum(v * v))
}

.bisect_to_level <- function(f, lo, hi, target, tol = 1e-4, maxit = 50L) {
  flo <- f(lo) - target
  fhi <- f(hi) - target

  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
    return(NA_real_)
  }

  for (it in seq_len(maxit)) {
    mid <- 0.5 * (lo + hi)
    fm <- f(mid) - target

    if (!is.finite(fm)) {
      return(NA_real_)
    }

    if (abs(fm) < tol) {
      return(mid)
    }

    if (flo * fm <= 0) {
      hi <- mid
      fhi <- fm
    } else {
      lo <- mid
      flo <- fm
    }
  }

  0.5 * (lo + hi)
}

.compute_radius_along_dir <- function(center, u, log_post_fn, c_thresh,
                                      r_max = 1e3, bisect_tol = 1e-4) {
  fdir <- function(r) log_post_fn(center + r * u)
  .bisect_to_level(
    f = fdir,
    lo = 0,
    hi = r_max,
    target = c_thresh,
    tol = bisect_tol
  )
}

.euclidean_m <- function(x_matrix, centers_matrix) {
  if (!is.matrix(x_matrix)) x_matrix <- as.matrix(x_matrix)
  if (!is.matrix(centers_matrix)) centers_matrix <- as.matrix(centers_matrix)

  n_points <- nrow(x_matrix)
  n_centers <- nrow(centers_matrix)
  distances <- matrix(0, nrow = n_points, ncol = n_centers)

  for (i in seq_len(n_centers)) {
    diff_matrix <- sweep(x_matrix, 2, centers_matrix[i, ], "-")
    distances[, i] <- sqrt(rowSums(diff_matrix^2))
  }

  distances
}

.is_point_in_any_ellipsoid <- function(points_matrix, centers_matrix, Sigmas_inv_list) {
  n_points <- nrow(points_matrix)
  n_ellipsoids <- nrow(centers_matrix)
  result <- logical(n_points)

  for (i in seq_len(n_ellipsoids)) {
    diff_matrix <- sweep(points_matrix, 2, centers_matrix[i, ], "-")
    temp_matrix <- diff_matrix %*% Sigmas_inv_list[[i]]
    mahal_distances <- rowSums(temp_matrix * diff_matrix)
    result <- result | (mahal_distances <= 1)
  }

  result
}

.compute_ellipsoid_volume <- function(Sigma, d) {
  unit_vol <- pi^(d / 2) / gamma(d / 2 + 1)
  unit_vol * det(Sigma)^0.5
}

.compute_log_marginal_likelihood <- function(log_ratio, n_samples) {
  out <- rep(Inf, n_samples)
  if (n_samples == 0L) {
    return(out)
  }

  finite_seen <- FALSE
  current_log_sum <- -Inf

  for (it in seq_len(n_samples)) {
    new_val <- log_ratio[it]

    if (is.finite(new_val)) {
      if (!finite_seen) {
        current_log_sum <- new_val
        finite_seen <- TRUE
      } else if (current_log_sum > new_val) {
        current_log_sum <- current_log_sum + log1p(exp(new_val - current_log_sum))
      } else {
        current_log_sum <- new_val + log1p(exp(current_log_sum - new_val))
      }
    }

    if (finite_seen) {
      out[it] <- -(current_log_sum - log(it))
    }
  }

  out
}

.make_basis <- function(u) {
  d <- length(u)
  if (d == 1L) {
    return(matrix(1, nrow = 1L, ncol = 1L))
  }

  A <- cbind(u, diag(d))
  qr.Q(qr(A), complete = TRUE)[, seq_len(d), drop = FALSE]
}

.check_log_post_fn <- function(log_post_fn, probe_point) {
  if (!is.function(log_post_fn)) {
    stop("`log_post_fn` must be a function.", call. = FALSE)
  }
  probe <- log_post_fn(probe_point)
  if (!is.numeric(probe) || length(probe) != 1L || !is.finite(probe)) {
    stop("`log_post_fn` must return one finite numeric value for a numeric parameter vector.", call. = FALSE)
  }
  invisible(TRUE)
}
