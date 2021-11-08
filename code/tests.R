source(here::here("code/shared_functions.r"))


test_cdf_and_quantile <- function() {
  cdf <- function(x) {
    cdf_inflated_censored_lognormal(x, w=0.8, T=2)
  }
  qtile <- function(x) {
    quantile_inflated_censored_lognormal(x, w=0.8, T=2)
  }
  grid <- seq(0, 15, by=0.01)
  prob_grid <- seq(0, 1, by=0.001)
  # Test that the quantile and CDF work:
  stopifnot(
    all(qtile(cdf(grid)) <= grid),
    all(cdf(qtile(prob_grid)) >= prob_grid)
  )
}
test_cdf_and_quantile()

test_lognormal_rng <- function() {
  # Check that we fill in the right order
  x = matrix(1:6, ncol=2)
  res <- lognormal_rng(x, x * 0)
  stopifnot(order(x) == order(res))
}
test_lognormal_rng()
