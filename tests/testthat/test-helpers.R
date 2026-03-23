# ============================================================
#  Tests for internal helper functions  (ECMLE:::)
# ============================================================

# ---- .validate_numeric_matrix --------------------------------------
test_that(".validate_numeric_matrix coerces vectors to matrix", {
  m <- ECMLE:::.validate_numeric_matrix(1:8, "x")
  expect_true(is.matrix(m))
})

test_that(".validate_numeric_matrix rejects NA values", {
  expect_error(
    ECMLE:::.validate_numeric_matrix(matrix(c(1, NA, 3, 4, 5, 6, 7, 8), 4, 2), "x"),
    "missing"
  )
})

test_that(".validate_numeric_matrix rejects fewer than 4 rows", {
  expect_error(
    ECMLE:::.validate_numeric_matrix(matrix(1:4, nrow = 2), "x")
  )
})

# ---- .validate_numeric_vector --------------------------------------
test_that(".validate_numeric_vector rejects NA", {
  expect_error(
    ECMLE:::.validate_numeric_vector(c(1, NA, 3), "v"),
    "missing"
  )
})

test_that(".validate_numeric_vector enforces length_out", {
  expect_error(
    ECMLE:::.validate_numeric_vector(1:5, "v", length_out = 3),
    "length"
  )
})

test_that(".validate_numeric_vector passes for valid input", {
  expect_equal(ECMLE:::.validate_numeric_vector(c(1.0, 2.0), "v"), c(1, 2))
})

# ---- .euclid_norm --------------------------------------------------
test_that(".euclid_norm returns correct values", {
  expect_equal(ECMLE:::.euclid_norm(c(3, 4)), 5)
  expect_equal(ECMLE:::.euclid_norm(c(0, 0, 0)), 0)
  expect_equal(ECMLE:::.euclid_norm(c(1, 0)), 1)
})

# ---- .make_basis ---------------------------------------------------
test_that(".make_basis returns an orthonormal matrix", {
  for (d in 2:4) {
    u <- rnorm(d); u <- u / sqrt(sum(u^2))
    Q <- ECMLE:::.make_basis(u)
    expect_equal(dim(Q), c(d, d))
    expect_equal(crossprod(Q), diag(d), tolerance = 1e-10)
    # first column should equal u up to sign
    expect_true(
      isTRUE(all.equal(Q[, 1],  u, tolerance = 1e-10)) ||
        isTRUE(all.equal(Q[, 1], -u, tolerance = 1e-10))
    )
  }
})

test_that(".make_basis works for d = 1", {
  Q <- ECMLE:::.make_basis(1)
  expect_equal(Q, matrix(1, 1, 1))
})

# ---- .bisect_to_level ----------------------------------------------
test_that(".bisect_to_level finds root of simple linear function", {
  f <- function(x) x
  r <- ECMLE:::.bisect_to_level(f, lo = -1, hi = 1, target = 0, tol = 1e-6)
  expect_true(abs(r) < 1e-5)
})

test_that(".bisect_to_level returns NA when bracket does not straddle target", {
  f <- function(x) x + 10
  r <- ECMLE:::.bisect_to_level(f, lo = 0, hi = 1, target = 0, tol = 1e-6)
  expect_true(is.na(r))
})

test_that(".bisect_to_level finds radius of standard Gaussian level set", {
  log_fn <- function(x) dnorm(x, log = TRUE)
  target <- log_fn(sqrt(2))
  r <- ECMLE:::.bisect_to_level(log_fn, lo = 0, hi = 5, target = target, tol = 1e-6)
  expect_equal(r, sqrt(2), tolerance = 1e-4)
})

# ---- .euclidean_m --------------------------------------------------
test_that(".euclidean_m returns correct pairwise distances", {
  X <- matrix(c(0, 3), nrow = 2, ncol = 1)
  C <- matrix(0, nrow = 1, ncol = 1)
  d <- ECMLE:::.euclidean_m(X, C)
  expect_equal(as.vector(d), c(0, 3))
})

test_that(".euclidean_m is zero on diagonal for identical points", {
  X <- matrix(rnorm(12), 4, 3)
  d <- ECMLE:::.euclidean_m(X, X)
  expect_equal(diag(d), rep(0, 4), tolerance = 1e-12)
})

# ---- .compute_ellipsoid_volume -------------------------------------
test_that(".compute_ellipsoid_volume gives pi * r^2 for a circle", {
  Sigma <- diag(2) * 4
  vol   <- ECMLE:::.compute_ellipsoid_volume(Sigma, d = 2L)
  expect_equal(vol, pi * 4, tolerance = 1e-10)
})

test_that(".compute_ellipsoid_volume gives (4/3)*pi*r^3 for a sphere", {
  r     <- 3
  Sigma <- diag(3) * r^2
  vol   <- ECMLE:::.compute_ellipsoid_volume(Sigma, d = 3L)
  expect_equal(vol, (4 / 3) * pi * r^3, tolerance = 1e-8)
})

# ---- .is_point_in_any_ellipsoid ------------------------------------
test_that(".is_point_in_any_ellipsoid correctly classifies inside/outside", {
  centers <- matrix(c(0, 0), nrow = 1)
  Sinv    <- list(solve(diag(2)))
  pts <- matrix(c(0,   0,
                  0.9, 0,
                  1.1, 0,
                  0,   1.5),
                nrow = 4, byrow = TRUE)
  res <- ECMLE:::.is_point_in_any_ellipsoid(pts, centers, Sinv)
  expect_equal(res, c(TRUE, TRUE, FALSE, FALSE))
})

test_that(".is_point_in_any_ellipsoid handles multiple ellipsoids", {
  centers <- matrix(c(0, 0, 5, 0), nrow = 2, byrow = TRUE)
  Sinv    <- list(diag(2), diag(2))
  pts     <- matrix(c(0.5, 0,
                      4.8, 0,
                      3.0, 0),
                    nrow = 3, byrow = TRUE)
  res <- ECMLE:::.is_point_in_any_ellipsoid(pts, centers, Sinv)
  expect_equal(res, c(TRUE, TRUE, FALSE))
})

# ---- .compute_log_marginal_likelihood ------------------------------
test_that(".compute_log_marginal_likelihood returns vector of length n", {
  set.seed(1)
  lr  <- rnorm(50)
  out <- ECMLE:::.compute_log_marginal_likelihood(lr, 50L)
  expect_length(out, 50L)
})

test_that(".compute_log_marginal_likelihood converges for known input", {
  C   <- 5
  lr  <- rep(log(C), 200L)
  out <- ECMLE:::.compute_log_marginal_likelihood(lr, 200L)
  expect_equal(tail(out, 1), -log(C), tolerance = 1e-8)
})

test_that(".compute_log_marginal_likelihood handles all-infinite log_ratio", {
  lr  <- rep(-Inf, 10L)
  out <- ECMLE:::.compute_log_marginal_likelihood(lr, 10L)
  expect_true(all(is.infinite(out)))
})

# ---- .check_log_post_fn --------------------------------------------
test_that(".check_log_post_fn accepts a valid function", {
  f <- function(theta) sum(dnorm(theta, log = TRUE))
  expect_true(ECMLE:::.check_log_post_fn(f, c(0, 0)))
})

test_that(".check_log_post_fn rejects non-function", {
  expect_error(ECMLE:::.check_log_post_fn(42, c(0, 0)), "function")
})

test_that(".check_log_post_fn rejects function returning non-scalar", {
  expect_error(
    ECMLE:::.check_log_post_fn(function(x) c(1, 2), c(0, 0)),
    "one finite"
  )
})

test_that(".check_log_post_fn rejects function returning Inf", {
  expect_error(
    ECMLE:::.check_log_post_fn(function(x) Inf, c(0, 0)),
    "one finite"
  )
})

# ---- Numerical accuracy of ecmle -----------------------------------
test_that("ecmle is accurate for standard bivariate Gaussian (log ML = 0)", {
  set.seed(42)
  N    <- 2000L
  post <- cbind(rnorm(N), rnorm(N))
  lps  <- apply(post, 1, function(z) sum(dnorm(z, log = TRUE)))
  lpfn <- function(theta) sum(dnorm(theta, log = TRUE))
  fit  <- ecmle(post, lps, lpfn, hpd_level = 0.75, seed = 7L)
  expect_true(abs(fit$log_marginal_likelihood) < 1.5)
  expect_true(is.finite(fit$log_marginal_likelihood))
})

test_that("ecmle estimate is finite and coverage_rate > 0 for 1-D input", {
  set.seed(10)
  N    <- 600L
  post <- matrix(rnorm(N, mean = 2), ncol = 1L)
  lps  <- dnorm(post[, 1], mean = 2, log = TRUE)
  lpfn <- function(theta) dnorm(theta, mean = 2, log = TRUE)
  fit  <- ecmle(post, lps, lpfn, hpd_level = 0.75, seed = 3L)
  expect_true(is.finite(fit$log_marginal_likelihood))
  expect_gt(fit$coverage_rate, 0)
})

test_that("ecmle result is reproducible with fixed seed", {
  set.seed(1)
  post <- cbind(rnorm(400), rnorm(400))
  lps  <- apply(post, 1, function(z) sum(dnorm(z, log = TRUE)))
  lpfn <- function(theta) sum(dnorm(theta, log = TRUE))
  fit1 <- ecmle(post, lps, lpfn, seed = 55L)
  fit2 <- ecmle(post, lps, lpfn, seed = 55L)
  expect_identical(fit1$log_marginal_likelihood, fit2$log_marginal_likelihood)
  expect_identical(fit1$n_ellipsoids, fit2$n_ellipsoids)
})