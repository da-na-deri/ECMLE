# ============================================================
#  Tests for rosenbrock_*, pair_plot, and plot_ecmle_2d
# ============================================================

# ---- helpers shared across tests ------------------------------------
.rb_setup <- function(d = 2, n = 20, n_samples = 200L, seed = 1L) {
  set.seed(seed)
  b     <- rep(-1, d - 1)
  a     <- rep(0,  d - 1)
  sigma <- 8
  theta_true <- seq(1, 2.8, length.out = d)
  Y_bar  <- rosenbrock_generate_data(theta_true, n, b, a, sigma)
  samps  <- rosenbrock_exact_posterior(Y_bar, n, b, a, sigma, n_samples)
  lps    <- rosenbrock_log_post_vec(samps, Y_bar, n, b, a, sigma)
  list(d = d, n = n, b = b, a = a, sigma = sigma,
       theta_true = theta_true, Y_bar = Y_bar,
       samps = samps, lps = lps)
}

# ---- rosenbrock_generate_data ---------------------------------------
test_that("rosenbrock_generate_data returns correct length vector", {
  set.seed(1)
  for (d in c(1L, 2L, 5L)) {
    Y_bar <- rosenbrock_generate_data(
      theta_true = seq(1, 2.8, length.out = d),
      n = 20, b = rep(-1, max(d - 1, 0)), a = rep(0, max(d - 1, 0)),
      sigma = 8
    )
    expect_length(Y_bar, d)
    expect_true(is.numeric(Y_bar))
    expect_false(anyNA(Y_bar))
  }
})

test_that("rosenbrock_generate_data errors on mismatched b/a length", {
  expect_error(
    rosenbrock_generate_data(c(1, 2), n = 10,
                             b = rep(-1, 5), a = rep(0, 1), sigma = 1),
    "b.*a"
  )
})

# ---- rosenbrock_exact_posterior -------------------------------------
test_that("rosenbrock_exact_posterior returns matrix of correct dims", {
  s <- .rb_setup(d = 2, n_samples = 300L)
  expect_true(is.matrix(s$samps))
  expect_equal(dim(s$samps), c(300L, 2L))
  expect_equal(colnames(s$samps), c("theta1", "theta2"))
  expect_false(anyNA(s$samps))
})

test_that("rosenbrock_exact_posterior works for d = 1", {
  set.seed(2)
  Y_bar <- rosenbrock_generate_data(1.5, n = 20,
                                    b = numeric(0), a = numeric(0),
                                    sigma = 2)
  samps <- rosenbrock_exact_posterior(Y_bar, n = 20,
                                      b = numeric(0), a = numeric(0),
                                      sigma = 2, n_samples = 100L)
  expect_equal(dim(samps), c(100L, 1L))
})

test_that("rosenbrock_exact_posterior works for d = 5", {
  set.seed(3)
  d <- 5
  s <- .rb_setup(d = d, n_samples = 100L)
  expect_equal(dim(s$samps), c(100L, d))
})

# ---- rosenbrock_log_post / rosenbrock_log_post_vec ------------------
test_that("rosenbrock_log_post returns a single finite number", {
  s <- .rb_setup()
  val <- rosenbrock_log_post(s$samps[1, ], s$Y_bar,
                             s$n, s$b, s$a, s$sigma)
  expect_length(val, 1L)
  expect_true(is.finite(val))
})

test_that("rosenbrock_log_post_vec matches single-vector version", {
  s   <- .rb_setup(n_samples = 50L)
  vec <- rosenbrock_log_post_vec(s$samps, s$Y_bar, s$n, s$b, s$a, s$sigma)
  scalar_vals <- apply(s$samps, 1, rosenbrock_log_post,
                       Y_bar = s$Y_bar, n = s$n,
                       b = s$b, a = s$a, sigma = s$sigma)
  expect_equal(length(vec), nrow(s$samps))
  expect_equal(vec, scalar_vals, tolerance = 1e-10)
})

test_that("rosenbrock_log_post_vec returns vector of correct length", {
  s <- .rb_setup(n_samples = 80L)
  expect_length(s$lps, 80L)
  expect_true(all(is.finite(s$lps)))
})

# ---- ecmle integration with Rosenbrock ------------------------------
test_that("ecmle runs on Rosenbrock 2-D samples", {
  s <- .rb_setup(n_samples = 400L, seed = 42L)
  log_post_fn <- function(theta)
    rosenbrock_log_post(theta, s$Y_bar, s$n, s$b, s$a, s$sigma)
  
  fit <- ecmle(s$samps, s$lps, log_post_fn, hpd_level = 0.75, seed = 1L)
  
  expect_s3_class(fit, "ecmle")
  expect_true(is.finite(fit$log_marginal_likelihood))
  expect_true(fit$coverage_rate >= 0 && fit$coverage_rate <= 1)
})

# ---- plot_ecmle_2d --------------------------------------------------
test_that("plot_ecmle_2d runs without error for d = 2", {
  s <- .rb_setup(n_samples = 400L, seed = 5L)
  log_post_fn <- function(theta)
    rosenbrock_log_post(theta, s$Y_bar, s$n, s$b, s$a, s$sigma)
  fit <- ecmle(s$samps, s$lps, log_post_fn, hpd_level = 0.75, seed = 2L)
  
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_no_error(plot_ecmle_2d(fit, post_samples = s$samps))
})

test_that("plot_ecmle_2d returns fit invisibly", {
  s <- .rb_setup(n_samples = 300L, seed = 6L)
  log_post_fn <- function(theta)
    rosenbrock_log_post(theta, s$Y_bar, s$n, s$b, s$a, s$sigma)
  fit <- ecmle(s$samps, s$lps, log_post_fn, hpd_level = 0.75, seed = 3L)
  
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  ret <- plot_ecmle_2d(fit, post_samples = s$samps)
  expect_identical(ret, fit)
})

test_that("plot_ecmle_2d emits message for d > 2 without plotting", {
  set.seed(10)
  d <- 3
  post_samps <- matrix(rnorm(600), ncol = d)
  lps <- rowSums(dnorm(post_samps, log = TRUE))
  log_post_fn <- function(theta) sum(dnorm(theta, log = TRUE))
  fit <- ecmle(post_samps, lps, log_post_fn, hpd_level = 0.75, seed = 4L)
  
  expect_message(plot_ecmle_2d(fit, post_samples = post_samps), "d = 3")
})

test_that("plot_ecmle_2d errors on non-ecmle input", {
  expect_error(plot_ecmle_2d(list(), matrix(rnorm(8), 4)), "ecmle")
})

# ---- draw_ellipse_2d ------------------------------------------------
test_that("draw_ellipse_2d adds lines without error", {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  graphics::plot(0, 0, xlim = c(-3, 3), ylim = c(-3, 3), type = "n")
  expect_no_error(
    draw_ellipse_2d(center = c(0, 0),
                    Sigma  = matrix(c(1, 0.4, 0.4, 1), 2, 2),
                    npoints = 50L)
  )
})
# ---- pair_plot -----------------------------------------------------
test_that("pair_plot runs without error", {
  s <- .rb_setup(n_samples = 300L)
  
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_no_error(pair_plot(s$samps, pixs = 1))
})