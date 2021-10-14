formeth = function (f, meth = c('TTDES','THREECG','TTCG'))
  Map(f, meth)

test_that("Convergence", {
  nqm = rnorm(500)*2
  par= rnorm(500)*4.
  fn = function (x,nqm1) sum((x - nqm1)^2.) + exp(-sum(x*x))
  gr = function (x,nqm1) 2.*(x - nqm1) + exp(-sum(x*x)) * -2.*x
  formeth(function (m) {
    r = ttcg(par = par, fn = fn, gr = gr, method=m, nqm1=nqm)
    expect_equal(r$convergence, 0L)
    expect_equal(r$value, 0.)
  })
})

test_that('Numerical differentiation', {
  nqm = rnorm(500)*2
  fn = function (x) sum((x - nqm)^2.) + exp(-sum(x*x))
  gr = NULL
  par= rnorm(500)*4.
  formeth(function (m) {
    r = ttcg(par = par/8, fn = fn, gr = gr, method=m)
    expect_equal(r$convergence, 0L)
    expect_equal(r$value, 0.)
  })
})

test_that('Function without lower bound #1', {
  nqm = rnorm(500)*2
  fn = function (x) -sum((x - nqm)^2.) + exp(-sum(x*x))
  gr = function (x) -2.*(x - nqm) + exp(-sum(x*x)) * -2.*x
  par= rnorm(500)*4.
  formeth(function (m) {
    r = ttcg(par = par/8, fn = fn, gr = gr, method=m)
    expect_gt(r$convergence, 0L)
    expect_gt(r$counts['function'], 1L)
    expect_gt(r$counts['gradient'], 1L)
    expect_equal(mode(r$value), 'numeric')
  })
})

test_that('Function without lower bound #2', {
  nqm = rnorm(500)*2
  fn = function (x) sum(x)
  gr = function (x) rep(1.0,length(x))
  par= rnorm(500)*4.
  formeth(function (m) {
    r = ttcg(par = par/8, fn = fn, gr = gr, method=m)
    expect_gt(r$convergence, 0L)
    expect_gt(r$counts['function'], 1L)
    expect_gt(r$counts['gradient'], 1L)
    expect_equal(mode(r$value), 'numeric')
  })
})

test_that('First iteration termination', {
  nqm = rnorm(500)*2
  fn = function (x,nqm1) sum((x - nqm1)^2.) + exp(-sum(x*x))
  gr = function (x,nqm1) 2.*(x - nqm1) + exp(-sum(x*x)) * -2.*x
  r = ttcg(par = nqm, fn = fn, gr = gr, method='TTDES', nqm1=nqm)
  expect_equal(r$convergence, 0L)
  expect_equal(r$counts['function'], c('function'=1L))
  expect_equal(r$counts['gradient'], c('gradient'=1L))
  expect_equal(mode(r$value), 'numeric')
})

test_that('NaN at init. param', {
  fn = function (x) NaN
  gr = function (x) c(1,1)
  r = ttcg(par = c(0,0), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 9L)
  expect_equal(r$counts['function'], c('function'=1L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_equal(mode(r$value), 'numeric')
})

test_that('NaN failure at line search param', {
  fn = function (x) if (sum(x*x) <.000000001) sum((x-2^40.)^2.) else NaN
  gr = function (x) if (sum(x*x) <.000000001) 2.*(x-2^40.)      else rep(NaN, length(x))
  r = ttcg(par = c(.0,.000000001/2.), fn = fn, gr = gr, method='TTCG')
  expect_equal(r$convergence, 7L)
  expect_gt(r$counts['function'], c('function'=5L))
  expect_gt(r$counts['gradient'], c('gradient'=5L))
  expect_equal(mode(r$value), 'numeric')
})

test_that('NaN success at line search param', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else NaN
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(NaN, length(x))
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='THREECG')
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_equal(r$value, 0.)
})

test_that('1-D function', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else NaN
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(NaN, length(x))
  formeth(function (m) {
    r = ttcg(par = -0.03, fn = fn, gr = gr, method=m)
    expect_equal(r$convergence, 0L)
    expect_gt(r$counts['function'], c('function'=1L))
    expect_gt(r$counts['gradient'], c('gradient'=1L))
    expect_equal(r$value, 0.)
  })
})

test_that('NA failure', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else NA
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(NA, length(x))
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='THREECG')
  expect_equal(r$convergence, 101L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_true(is.na(r$value))
})

test_that('-Inf failure', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else -Inf
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(-Inf, length(x))
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='THREECG')
  expect_equal(r$convergence, 101L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_true(identical(r$value, -Inf))
})

test_that('Inf success', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else Inf
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(Inf, length(x))
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='THREECG')
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_equal(r$value, 0.)
})

test_that('Complex coercion', {
  fn = function (x) as.complex( sum((x-1.)^2.) )
  gr = function (x) as.complex( 2.*(x-1.) )
  expect_warning({ ttcg(par = c(-0.03+3i,0.03+4i), fn = fn, gr = gr, method='THREECG') }, regexp='imaginary part')
})

test_that('List coercion', {
  fn = function (x) as.complex( list(sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(2.*(x-1.)) )
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_equal(r$value, 0.)
})

test_that('Parameter names preservation', {
  fn = function (x) as.complex( list(sum((x-1.)^2.)) )
  gr = function (x) { as.complex(as.list(2.*(x-1.))) -> .;
                      c('A','B','D') -> names(.);
                      . }
  r = ttcg(par = c(A=1,B=1,C=3), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=1L))
  expect_gt(r$counts['gradient'], c('gradient'=1L))
  expect_identical(names(r[['par']]),c('A','B','C'))
  expect_equal(r$value, 0.)
})

test_that('Maximum iteration', {
  nqm = rep(1.,400)
  par=  rep(4.,400)
  fn = function (x,nqm1) sum((x - nqm1)^2.) + exp(-sum(x*x))
  gr = function (x,nqm1) 2.*(x - nqm1) + exp(-sum(x*x)) * -2.*x
  r = ttcg(par = par, fn = fn, gr = gr, method='THREECG', control=list(maxit=1), nqm1=nqm)
  expect_equal(r$convergence, 1L)
})

test_that('NaN in initial parameters', {
  fn = function (x) as.complex( list(sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(2.*(x-1.)) )
  r = ttcg(par = c(NaN,0.03), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 107L)
  expect_equal(r$counts['function'], c('function'=0L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_true(identical(r$value, as.double(NA)))
})

test_that('Zero-length initial parameters', {
  fn = function (x) as.complex( list(sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(2.*(x-1.)) )
  r = ttcg(par = double(0L), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 106L)
  expect_equal(r$counts['function'], c('function'=0L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_true(identical(r$value, as.double(NA)))
})

test_that('Infinite initial parameters', {
  fn = function (x) as.complex( list(sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(2.*(x-1.)) )
  r = ttcg(par = c(Inf, 3.0), fn = fn, gr = gr, method='TTDES')
  expect_equal(r$convergence, 108L)
  expect_equal(r$counts['function'], c('function'=0L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_true(identical(r$value, as.double(NA)))
})

test_that('Maximization', {
  fn = function (x) as.complex( list(-sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(-2.*(x-1.)) )
  r = ttcg(par = c(1.0, 3.0), fn = fn, gr = gr, method='TTCG', control=list(fnscale = -10.0))
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=2L))
  expect_gt(r$counts['gradient'], c('gradient'=2L))
  expect_equal(r$value, 0.0)
})

test_that('Maximization failure', {
  fn = function (x) as.complex( list(-sum((x-1.)^2.)) )
  gr = function (x) as.complex( as.list(-2.*(x-1.)) )
  r = ttcg(par = c(1.0, 3.0), fn = fn, gr = gr, method='TTCG', control=list(fnscale = Inf))
  expect_equal(r$convergence, 112L)
  expect_equal(r$counts['function'], c('function'=0L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_true(identical(r$value, as.double(NA)))
})

test_that('Parscale #1', {
  fn = function (x) if (sqrt(sum(x*x)) <1.48) sum((x-1.)^2.) else Inf
  gr = function (x) if (sqrt(sum(x*x)) <1.48) 2.*(x-1.)      else rep(Inf, length(x))
  r = ttcg(par = c(-0.03,0.03), fn = fn, gr = gr, method='THREECG', control=list(fnscale = 3.5, parscale=c(.1,.1)))
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=2L))
  expect_gt(r$counts['gradient'], c('gradient'=2L))
  expect_equal(r$value, 0.0)
})

test_that('Parscale #2', {
  fn = function (x) -sum((x-1.)^2.)
  gr = function (x) -2.*(x-1.)
  ## Very badly scaled.
  r = ttcg(par = c(-0.03,2.03), fn = fn, gr = gr, method='THREECG', control=list(fnscale = -.5, parscale=c(-100,-.001)))
  expect_equal(r$convergence, 0L)
  expect_gt(r$counts['function'], c('function'=2L))
  expect_gt(r$counts['gradient'], c('gradient'=2L))
  expect_equal(r$value, 0.0)
})

test_that('Parscale failure', {
  fn = function (x) -sum((x-1.)^2.)
  gr = function (x) -2.*(x-1.)
  ## Very badly scaled.
  r = ttcg(par = c(-0.03,2.03), fn = fn, gr = gr, method='THREECG', control=list(fnscale = -.5, parscale=c(0.0,-.001)))
  expect_equal(r$convergence, 123L)
  expect_equal(r$counts['function'], c('function'=0L))
  expect_equal(r$counts['gradient'], c('gradient'=0L))
  expect_true(identical(r$value, as.double(NA)))
})

