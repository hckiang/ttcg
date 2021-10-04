test_that("Convergence", {
  nqm = rnorm(500)*2
  par= rnorm(500)*4.
  fn = function (x,nqm1) sum((x - nqm1)^2.) + exp(-sum(x*x))
  gr = function (x,nqm1) 2.*(x - nqm1) + exp(-sum(x*x)) * -2.*x
  r = ttcg(par = par, fn = fn, gr = gr, method='TTDES', nqm=nqm)
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='THREECG', nqm=nqm)
  expect_equal(r$value, 0.)
  r = ttcg(par = par, fn = fn, gr = gr, method='TTCG', nqm=nqm)
  expect_equal(r$value, 0.)

  fn = function (x) (x[1]+2.*x[2]-7.)^2. + (2.*x[1]+x[2]-5)^2.
  gr = NULL
  par= rnorm(500)*4.

  fn = function (x) sum((x - nqm)^2.) + exp(-sum(x*x))
  gr = NULL
  r = ttcg(par = par/8, fn = fn, gr = gr, method='TTDES')
  expect_equal(r$value, 0.)
  r = ttcg(par = par/8, fn = fn, gr = gr, method='THREECG')
  expect_equal(r$value, 0.)
  r = ttcg(par = par/8, fn = fn, gr = gr, method='TTCG')
  expect_equal(r$value, 0.)
})

