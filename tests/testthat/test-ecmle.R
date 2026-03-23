test_that("ecmle returns an object with expected fields", {
  set.seed(1)
  post_samples <- cbind(rnorm(300), rnorm(300))
  lps <- apply(post_samples, 1, function(z) sum(dnorm(z, log = TRUE)))
  log_post_fn <- function(theta) sum(dnorm(theta, log = TRUE))

  fit <- ecmle(
    post_samples = post_samples,
    lps = lps,
    log_post_fn = log_post_fn,
    hpd_level = 0.75,
    seed = 123
  )

  expect_s3_class(fit, "ecmle")
  expect_true(is.numeric(fit$log_marginal_likelihood))
  expect_length(fit$log_marginal_likelihood, 1)
  expect_true(is.numeric(fit$coverage_rate))
  expect_true(fit$coverage_rate >= 0)
  expect_true(fit$coverage_rate <= 1)
  expect_true(is.list(fit$ellipsoids))
})

test_that("ecmle validates inputs", {
  post_samples <- cbind(rnorm(20), rnorm(20))
  lps <- rep(0, 19)
  log_post_fn <- function(theta) sum(theta)

  expect_error(ecmle(post_samples, lps, log_post_fn), "length")
  expect_error(ecmle(post_samples, rep(0, 20), log_post_fn, hpd_level = 1.2), "hpd_level")
  expect_error(ecmle(post_samples, rep(0, 20), 4), "function")
})
