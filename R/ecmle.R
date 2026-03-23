#' Estimate the log marginal likelihood with ECMLE
#'
#' Computes an ECMLE estimate of the log marginal likelihood from posterior
#' draws, corresponding log unnormalized posterior evaluations, and a function
#' that evaluates the log unnormalized posterior density at arbitrary points.
#'
#' @param post_samples Numeric matrix of posterior draws; each row is one draw.
#' @param lps Numeric vector of log unnormalized posterior values (log prior +
#'   log likelihood, up to an additive constant) evaluated at the rows of
#'   \code{post_samples}. Must be consistent with \code{log_post_fn}.
#' @param log_post_fn A function that accepts a numeric vector of length
#'   \eqn{d} and returns a single finite numeric scalar equal to the log
#'   unnormalized posterior density at that point. Must be consistent with
#'   \code{lps}.
#' @param hpd_level Fraction in \code{(0, 1)} used to define the HPD region.
#'   Default is \code{0.75}.
#' @param subsample_frac Fraction of HPD points retained for candidate
#'   ellipsoid centers. Default is \code{0.05}.
#' @param r_max Upper bracketing radius used by the bisection solver. If
#'   \code{NULL}, estimated from pairwise distances between retained HPD points.
#' @param bisect_tol Tolerance for the directional radius solver. Default
#'   is \code{1e-4}.
#' @param seed Optional integer seed used for HPD subsampling.
#' @param verbose Logical; if \code{TRUE}, emits basic progress information.
#'
#' @note Consistency of \code{lps} and \code{log_post_fn} is critical. Both
#'   must evaluate the same log unnormalized posterior (log prior + log
#'   likelihood) up to the same additive constant. If they differ by more
#'   than a constant the marginal likelihood estimate will be incorrect.
#'
#' @return An object of class \code{"ecmle"}, a list with components:
#' \describe{
#'   \item{\code{log_marginal_likelihood}}{Final log marginal likelihood estimate.}
#'   \item{\code{log_marginal_likelihood_iter}}{Running estimate over the second half-sample.}
#'   \item{\code{ellipsoids}}{List of fitted ellipsoids.}
#'   \item{\code{total_volume}}{Total ellipsoid volume.}
#'   \item{\code{n_ellipsoids}}{Number of ellipsoids.}
#'   \item{\code{points_in_ellipsoids}}{Number of evaluation points inside the ellipsoids.}
#'   \item{\code{n_samples}}{Number of evaluation points in the second half-sample.}
#'   \item{\code{coverage_rate}}{Fraction of evaluation points covered by the ellipsoids.}
#'   \item{\code{hpd_level}}{HPD level used in the fit.}
#' }
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' b <- rep(-1, d - 1)
#' a <- rep(0,  d - 1)
#'
#' Y_bar <- rosenbrock_generate_data(
#'   theta_true = seq(1, 2.8, length.out = d),
#'   n = 20, b = b, a = a, sigma = 8
#' )
#' samps <- rosenbrock_exact_posterior(Y_bar, n = 20,
#'   b = b, a = a, sigma = 8, n_samples = 500L)
#' lps <- rosenbrock_log_post_vec(samps, Y_bar,
#'   n = 20, b = b, a = a, sigma = 8)
#' log_post_fn <- function(theta)
#'   rosenbrock_log_post(theta, Y_bar, n = 20, b = b, a = a, sigma = 8)
#'
#' fit <- ecmle(samps, lps, log_post_fn, hpd_level = 0.75, seed = 1L)
#' fit
#' summary(fit)
#' plot_ecmle_2d(fit, post_samples = samps)
#' @export
ecmle <- function(post_samples,
                  lps,
                  log_post_fn,
                  hpd_level = 0.75,
                  subsample_frac = 0.05,
                  r_max = NULL,
                  bisect_tol = 1e-4,
                  seed = NULL,
                  verbose = FALSE) {
  post_samples <- .validate_numeric_matrix(post_samples, "post_samples")
  lps <- .validate_numeric_vector(lps, "lps", length_out = nrow(post_samples))
  
  if (!is.numeric(hpd_level) || length(hpd_level) != 1L || hpd_level <= 0 || hpd_level >= 1) {
    stop("`hpd_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(subsample_frac) || length(subsample_frac) != 1L ||
      subsample_frac <= 0 || subsample_frac > 1) {
    stop("`subsample_frac` must be a single number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(bisect_tol) || length(bisect_tol) != 1L || bisect_tol <= 0) {
    stop("`bisect_tol` must be a positive number.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
    stop("`seed` must be `NULL` or a single non-missing number.", call. = FALSE)
  }
  
  d <- ncol(post_samples)
  .check_log_post_fn(log_post_fn, post_samples[1, ])
  
  N <- nrow(post_samples)
  half <- floor(N / 2)
  if (half < 2L || (N - half) < 2L) {
    stop("`post_samples` must contain enough rows to split into two non-trivial halves.",
         call. = FALSE)
  }
  
  params1 <- post_samples[seq_len(half), , drop = FALSE]
  lps1    <- lps[seq_len(half)]
  params2 <- post_samples[seq.int(half + 1L, N), , drop = FALSE]
  lps2    <- lps[seq.int(half + 1L, N)]
  
  c_thresh <- stats::quantile(lps1, probs = 1 - hpd_level, names = FALSE, type = 7)
  in_hpd   <- lps1 >= c_thresh
  
  HPD_samples     <- params1[in_hpd,  , drop = FALSE]
  HPD_lps         <- lps1[in_hpd]
  Non_HPD_samples <- params1[!in_hpd, , drop = FALSE]
  
  if (nrow(HPD_samples) < 2L) {
    stop(
      "Too few HPD points were selected. ",
      "Try increasing the sample size or changing `hpd_level`.",
      call. = FALSE
    )
  }
  
  size <- max(1L, floor(subsample_frac * nrow(HPD_samples)))
  
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  
  hpd_index   <- sample.int(nrow(HPD_samples), size = size, replace = FALSE)
  HPD_samples <- HPD_samples[hpd_index, , drop = FALSE]
  HPD_lps     <- HPD_lps[hpd_index]
  
  HPD_matrix     <- as.matrix(HPD_samples)
  Non_HPD_matrix <- as.matrix(Non_HPD_samples)
  
  if (is.null(r_max)) {
    dists <- .euclidean_m(HPD_matrix, HPD_matrix)
    r_max <- max(dists, na.rm = TRUE)
    if (!is.finite(r_max) || r_max <= 0) {
      r_max <- 1
    }
  } else if (!is.numeric(r_max) || length(r_max) != 1L ||
             !is.finite(r_max) || r_max <= 0) {
    stop("`r_max` must be `NULL` or one positive finite number.", call. = FALSE)
  }
  
  sorted_idx          <- order(-HPD_lps)
  HPD_matrix_sorted   <- HPD_matrix[sorted_idx, , drop = FALSE]
  HPD_lps_sorted      <- HPD_lps[sorted_idx]
  
  ellipsoid_centers_list   <- list()
  ellipsoid_Sigmas_list    <- list()
  ellipsoid_Sigmas_inv_list <- list()
  ellipsoid_max_semi_list  <- numeric(0)
  ellipsoid_volumes        <- numeric(0)
  ellipsoid_lps            <- numeric(0)
  ellipsoid_count          <- 0L
  
  available_centers <- rep(TRUE, nrow(HPD_matrix_sorted))
  
  while (any(available_centers)) {
    next_idx       <- which(available_centers)[1L]
    current_center <- HPD_matrix_sorted[next_idx, ]
    current_lps    <- HPD_lps_sorted[next_idx]
    
    if (verbose) {
      message("Processing candidate center ", next_idx, "/", length(available_centers))
    }
    
    if (nrow(Non_HPD_matrix) > 0L) {
      euclid_dists <- .euclidean_m(Non_HPD_matrix, matrix(current_center, nrow = 1L))
      closest_idx  <- which.min(euclid_dists)
      closest_point <- Non_HPD_matrix[closest_idx, ]
      dir_vec      <- closest_point - current_center
    } else {
      dir_vec <- c(1, rep(0, d - 1L))
    }
    
    dir_norm <- .euclid_norm(dir_vec)
    if (!is.finite(dir_norm) || dir_norm == 0) {
      available_centers[next_idx] <- FALSE
      next
    }
    u1 <- dir_vec / dir_norm
    
    r1 <- .compute_radius_along_dir(
      center     = current_center,
      u          = u1,
      log_post_fn = log_post_fn,
      c_thresh   = c_thresh,
      r_max      = r_max,
      bisect_tol = bisect_tol
    )
    
    if (is.na(r1) || r1 <= 0) {
      available_centers[next_idx] <- FALSE
      next
    }
    
    U         <- .make_basis(u1)
    semi_axes <- numeric(d)
    semi_axes[1L] <- r1
    
    if (d >= 2L) {
      for (i in 2:d) {
        ui    <- U[, i]
        r_pos <- .compute_radius_along_dir(
          center     = current_center,
          u          = ui,
          log_post_fn = log_post_fn,
          c_thresh   = c_thresh,
          r_max      = r_max,
          bisect_tol = bisect_tol
        )
        r_neg <- .compute_radius_along_dir(
          center     = current_center,
          u          = -ui,
          log_post_fn = log_post_fn,
          c_thresh   = c_thresh,
          r_max      = r_max,
          bisect_tol = bisect_tol
        )
        semi_axes[i] <- if (is.na(r_pos) || is.na(r_neg)) NA_real_ else min(r_pos, r_neg)
      }
    }
    
    if (any(is.na(semi_axes)) || any(semi_axes <= 0)) {
      available_centers[next_idx] <- FALSE
      next
    }
    
    D      <- diag(semi_axes^2, nrow = d, ncol = d)
    Sigma  <- U %*% D %*% t(U)
    max_semi <- max(semi_axes)
    
    overlaps <- FALSE
    if (ellipsoid_count > 0L) {
      existing_centers_matrix <- do.call(rbind, ellipsoid_centers_list)
      for (j in seq_len(ellipsoid_count)) {
        existing_center  <- existing_centers_matrix[j, ]
        center_distance  <- .euclid_norm(current_center - existing_center)
        max_semi_j       <- ellipsoid_max_semi_list[j]
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
    
    volume          <- .compute_ellipsoid_volume(Sigma, d)
    ellipsoid_count <- ellipsoid_count + 1L
    ellipsoid_centers_list[[ellipsoid_count]]    <- current_center
    ellipsoid_Sigmas_list[[ellipsoid_count]]     <- Sigma
    ellipsoid_Sigmas_inv_list[[ellipsoid_count]] <- solve(Sigma)
    ellipsoid_max_semi_list[ellipsoid_count]     <- max_semi
    ellipsoid_volumes[ellipsoid_count]           <- volume
    ellipsoid_lps[ellipsoid_count]               <- current_lps
    
    remaining_indices <- which(available_centers)
    if (length(remaining_indices) > 0L) {
      remaining_centers_matrix <- HPD_matrix_sorted[remaining_indices, , drop = FALSE]
      is_inside_new <- .is_point_in_any_ellipsoid(
        points_matrix  = remaining_centers_matrix,
        centers_matrix = matrix(current_center, nrow = 1L),
        Sigmas_inv_list = list(ellipsoid_Sigmas_inv_list[[ellipsoid_count]])
      )
      available_centers[remaining_indices[is_inside_new]] <- FALSE
    }
    
    available_centers[next_idx] <- FALSE
  }
  
  total_volume <- sum(ellipsoid_volumes)
  ellipsoids   <- vector("list", ellipsoid_count)
  if (ellipsoid_count > 0L) {
    for (i in seq_len(ellipsoid_count)) {
      ellipsoids[[i]] <- list(
        center = ellipsoid_centers_list[[i]],
        Sigma  = ellipsoid_Sigmas_list[[i]],
        volume = ellipsoid_volumes[i],
        lps    = ellipsoid_lps[i]
      )
    }
  }
  
  if (ellipsoid_count == 0L) {
    out <- list(
      log_marginal_likelihood      = -Inf,
      log_marginal_likelihood_iter = rep(-Inf, nrow(params2)),
      ellipsoids                   = ellipsoids,
      total_volume                 = 0,
      n_ellipsoids                 = 0L,
      points_in_ellipsoids         = 0L,
      n_samples                    = nrow(params2),
      coverage_rate                = 0,
      hpd_level                    = hpd_level
    )
    class(out) <- "ecmle"
    return(out)
  }
  
  ellipsoid_centers_matrix <- do.call(rbind, ellipsoid_centers_list)
  posterior_points         <- as.matrix(params2)
  n_samples                <- nrow(posterior_points)
  
  samples_in_ellipsoids <- .is_point_in_any_ellipsoid(
    points_matrix   = posterior_points,
    centers_matrix  = ellipsoid_centers_matrix,
    Sigmas_inv_list = ellipsoid_Sigmas_inv_list
  )
  
  samples_in_ellipsoids_numeric    <- as.numeric(samples_in_ellipsoids)
  log_ratio                        <- log(samples_in_ellipsoids_numeric) -
    log(total_volume) - lps2
  log_marginal_likelihood_iter     <- .compute_log_marginal_likelihood(log_ratio, n_samples)
  final_log_marginal_likelihood    <- log_marginal_likelihood_iter[n_samples]
  
  out <- list(
    log_marginal_likelihood      = final_log_marginal_likelihood,
    log_marginal_likelihood_iter = log_marginal_likelihood_iter,
    ellipsoids                   = ellipsoids,
    total_volume                 = total_volume,
    n_ellipsoids                 = ellipsoid_count,
    points_in_ellipsoids         = sum(samples_in_ellipsoids_numeric),
    n_samples                    = n_samples,
    coverage_rate                = sum(samples_in_ellipsoids_numeric) / n_samples,
    hpd_level                    = hpd_level
  )
  
  class(out) <- "ecmle"
  out
}