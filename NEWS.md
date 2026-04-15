# ECMLE 0.1.0

* Initial CRAN release.
* **`ecmle()`**: main estimator of the log marginal likelihood via ellipsoidal
  covering of the high-posterior-density region.
* **`rosenbrock_generate_data()`**, **`rosenbrock_exact_posterior()`**,
  **`rosenbrock_log_post()`**, **`rosenbrock_log_post_vec()`**: Rosenbrock
  (banana) posterior benchmark functions. The methodology is described in
  Naderi et al. (2025) <https://doi.org/10.48550/arXiv.2510.20617>.
* **`pair_plot()`**: lower-triangular pair plot for posterior sample matrices.
* **`plot_ecmle_2d()`**, **`draw_ellipse_2d()`**: 2-D visualisation of fitted
  ellipsoids and posterior samples.
* **`print()`**, **`summary()`**, and **`plot()`** S3 methods for **`"ecmle"`** objects.