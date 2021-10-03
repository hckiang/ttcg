test_that("Convergence", {
  nqm = rnorm(5000)*2
  par= rnorm(50000)*4.
  fn = function (x) sum((x - nqm)^2.) + exp(-sum(x*x))
  gr = function (x) 2.*(x - nqm) + exp(-sum(x*x)) * -2.*x
  r = ttcg(par = par, fn = fn, gr = gr, method='TTDES')
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='THREECG')
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='TTCG')
  expect_equal(r$value, 0.)

  fn = function (x) (x[1]+2.*x[2]-7.)^2. + (2.*x[1]+x[2]-5)^2.
  gr = NULL
  par= c(2,-2)

  fn = function (x) sum((x - nqm)^2.) + exp(-sum(x*x))
  gr = function (x) 2.*(x - nqm) + exp(-sum(x*x)) * -2.*x
  r = ttcg(par = par, fn = fn, gr = gr, method='TTDES')
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='THREECG')
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='TTCG')
  expect_equal(r$value, 0.)
  ##print(r$count)
})

