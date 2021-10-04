#' ttcg: Three-term conjugate gradient minimization algorithms
#'
#' The ttcg package implements several Neculai Andrei's three-term conjugate gradient algorithms.
#' 
#' @author Hao Chi Kiang, \email{hello@hckiang.com}
#' @docType package
#' @name ttcg
#' @importFrom numDeriv grad
NULL


maybe = function (L, name, default)   if (is.null(L[[name]])) default else L[[name]]
printf = function (trace, fmt, ...)   if (trace) cat(sprintf(fmt, ...))

#' Accelerated three-term conjugate gradient optimization with restart
#'
#' The \code{ttcg} function minimizes a given function using several Neculai Andrei's
#' three-term conjugate gradient algorithms. Line search is done by a bisection-like
#' weak-Wolfe search inspired by Oliveira and Takahashi's interpolate-truncate-project
#' algorithm and More-Thuente algorithm.
#' 
#' The \code{control} argument is a list that can contain any of the following named
#' element:
#'
#'     \describe{
#'       \item{maxit}{The maximal number of iteration. Default is 500.}
#'       \item{gl2tol}{A positive small number. The iteration will be terminated if the
#'                     squared Euclidean norm is smaller than this number. Default is
#'                     \code{min(1e-9, length(par)*1e-10)}. To turn off this test,
#'                     set it to any negative values.}
#'       \item{gmaxtol}{A positive small number. The iteration will be terminated if the infinity norm of the graident is smaller than
#'                      \code{gmaxtol}. Default is 1e-6. To turn off this test, set it to any negative values.}
#'       \item{ftol}{A positive small number. The iteration will be terminated if \eqn{|f_{k+1} - f_k| < ftol * (1 + |f_k|)}.
#'                   To turn off this test, set it to any negative values.}
#'       \item{c1}{The line search parameter for the sufficient descent condition. Default is 1e-3.}
#'       \item{c2}{The line search parameter for the curvature condition. Default is 0.08.}
#'       \item{trace}{Either TRUE or FALSE, indicating whether or not details should be printed to the terminal
#'                    during the optimization. Default is FALSE.}
#'     }
#' 
#'
#' @param par      A numerical vector containing the initial parameters.
#' @param fn       A function to be minimized. It should accept a numerical vector as its sole
#'                 argument and return a scaler.
#' @param gr       The gradient function of \code{fn}. It should accept the same argument as \code{fn}
#'                 and return a vector of the same length as \code{par}.
#' @param method   A character string, one of 'TTDES', 'TTCG', 'THREECG'. This determines how the
#'                 line search direction is computed. 'TTDES' is the default method.
#' @param control  A list of control parameters. See Details.
#' @param ...      Extra arguments to be passed to \code{fn}
#' @return       A list containing the following named elements.
#'               \describe{
#'                 \item{par}{The optimal parameter.}
#'                 \item{value}{The optimal function value.}
#'                 \item{counts}{An integer vector containing the number of function and gradient calls
#'                               used during the optimization.}
#'                 \item{convergence}{An integer indicating convergence status. '0' means successful convergence; '1'
#'                                    means \code{maxit} has been reached; '2' means a line search failure in which
#'                                    a point that satisfies the weak Wolfe condition is not found (perhaps the function
#'                                    is unbounded below).}
#'                 \item{message}{A character string giving additional message.}
#'               }
#' @references   Andrei, N. (2013). On three-term conjugate gradient algorithms for unconstrained optimization. Applied Mathematics and Computation, 219(11), 6316-6327.
#' @references   Andrei, N. (2013). A simple three-term conjugate gradient algorithm for unconstrained optimization. Journal of Computational and Applied Mathematics, 241, 19-29.
#' @references   Andrei, N. (2015). A new three-term conjugate gradient algorithm for unconstrained optimization. Numerical Algorithms, 68(2), 305-321.
#' @references   Oliveira, I. F., & Takahashi, R. H. (2020). An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality. ACM Transactions on Mathematical Software (TOMS), 47(1), 1-24.
#' @references   More, J. J., & Thuente, D. J. (1994). Line search algorithms with guaranteed sufficient decrease. ACM Transactions on Mathematical Software (TOMS), 20(3), 286-307.
#' @export
ttcg = function (par, fn, gr = NULL, method='TTDES', control = list(), ...) {
  mode(par) = 'double'
  dim(par)  = NULL

  maxit  = maybe(control, 'maxit', 500L)
  gl2tol = maybe(control, 'gl2tol', min(1e-9, length(par)*1e-10))  ## each of the gk1 < 1e-6
  gmaxtol= maybe(control, 'gmaxtol', 1e-6)
  ftol   = maybe(control, 'ftol', 1e-9)
  c1     = maybe(control, 'c1',  1e-3)
  c2     = maybe(control, 'c2',  .08)

  trace  = maybe(control, 'trace', FALSE)

  npar = length(par)
  if (is.null(gr)) {
    gr = function (x) numDeriv::grad(fn, x, ...)
  }
  extra_arg = list(...)
  xk    = par
  fk    = do.call(fn, c(list(xk),extra_arg))
  gk    = do.call(gr, c(list(xk),extra_arg))
  dk    = -gk
  k     = 1L
  fevl  = 1L
  gevl  = 1L
  itcnt = 0L
  convergence = 1L
  msg   = NULL
  
  if (! (method %in% c('TTDES', 'THREECG', 'TTCG'))) {
    stop('Invalid method. `method` should be one of `TTDES`, `THREECG`, `HTTCG`')
  }
  
  tee0_init = 1./sqrt(sum(gk*gk))
  xk1 = numeric(1L)
  while ((itcnt=itcnt+1L) < maxit) {
    wsres = weak_wolfe_search(fn, gr, xk, dk, f0=fk, g0=gk, extra_arg= extra_arg, c1=c1,c2=c2, tee0=tee0_init, trace = trace)
    if (is.character(wsres)) {
      convergence = 2
      msg = wsres
      break
    }
    tee1= wsres[[1]]
    xk1 = wsres[[2]]
    fk1 = wsres[[3]]
    gk1 = wsres[[4]]
    gdir= wsres[[5]]
    fevl=fevl + wsres[[6]]
    gevl=gevl + wsres[[7]]
    if ((cond1 = (max(abs(gk1)) < gmaxtol)) || (cond2 = (sum(gk1*gk1) < gl2tol)) || (cond3 = (abs(fk1-fk) < ftol*(1+abs(fk))))) {
        convergence = 0L
        msg = if      (cond1) "Gradient's infinite norm < gmaxtol"
              else if (cond2) "Gradient's squared 2-norm < gl2tol"
              else            "|f_{k+1} - f_k| < ftol * (1 + |f_k|)"
        break
    }
    yk  = gk1 - gk
    bbar= sum(yk*dk)
    ## Sacrifice one function evaluation and one gradient for acceleration scheme.
    if (bbar>0) {
      abar = sum(gk * dk)
      xik  = -abar / bbar
      xk1  = xk + xik * tee1 * dk
      fk1  = do.call(fn, c(list(xk1),extra_arg))
      gk1  = do.call(gr, c(list(xk1),extra_arg))
      fevl = fevl+1L
      gevl = gevl+1L
      if ((cond1 = (max(abs(gk1)) < gmaxtol)) || (cond2 = (sum(gk1*gk1) < gl2tol)) || (cond3 = (abs(fk1-fk) < ftol*(1+abs(fk))))) {
          convergence = 0L
          msg = if      (cond1) "Gradient's infinite norm < gmaxtol"
                else if (cond2) "Gradient's squared 2-norm < gl2tol"
                else            "|f_{k+1} - f_k| < ftol * (1 + |f_k|)"
          break
      }
      yk  = gk1 - gk
      sk  = xk1 - xk
    } else {
      sk  = xk1 - xk
    }

    switch(method,
           "TTDES"={
             yksk= sum(yk*sk)
             skgk1=sum(sk*gk1)
             omega  = (2./sum(sk*sk)) *sqrt(sum(sk*sk)*sum(yk*yk) - yksk^2.)
             dk1 = -gk1 + (sum(yk*gk1) - omega*skgk1)/yksk * sk - skgk1/yksk*yk
           },
           "THREECG"={
             yksk= sum(yk*sk)
             skgk1=sum(sk*gk1)
             dk1 = -gk1 - ((1.+sum(yk*yk)/yksk)*skgk1/yksk - sum(yk*gk1)/yksk)*sk
                   -skgk1/yksk * yk
           },
           "TTCG"={
             yksk= sum(yk*sk)
             skgk1=sum(sk*gk1)
             dk1 = -gk1 - ((1.+2.*sum(yk*yk)/yksk)*(skgk1/yksk) - sum(yk*gk1)/yksk)*sk
                   - (skgk1/yksk) * yk
           })

    ## If dk1 is way smaller than gk1 on restart this will explode.
    tee0_init = tee1 * sqrt(sum(dk*dk)/sum(dk1*dk1))
    printf(trace, '***** ')
    if (abs(sum(gk1*gk)) > 0.2*(gknorm=sum(gk1*gk1))) {
      dk1 = -gk1
      printf(trace, '%-9s  ', 'RESTART')
    } else {
      printf(trace, '%-9s  ', 'CARRYON')
    }
    
    printf(trace, 'CURR_STEP=%-10f   GUESSED_STEP=%-10f   FVAL=%-10f   GNORM=%-10f\n', tee1, tee0_init, fk1, sqrt(gknorm))
    dk = dk1
    gk = gk1
    fk = fk1
    xk = xk1
  }
  if (convergence == 1L) {
    msg = 'Maximum number iteration has been reached.'
  }
  printf(trace, '\nThe algorithm has stopped at iteratiion %d, because:\n', itcnt)
  printf(trace, '%s\n', msg)
  list(par = xk1, value = fk1, counts=c('function'=fevl, gradient=gevl), convergence=convergence, message=msg)
}

itppt = function (j, fleft, fright, gleft, gright, a, beta, nhalf, eps, n0, kap1, kap2) {
  if (! is.null(gright)) {
    if (fright > fleft) {
      S = c(2*fleft - 2*fright + gleft + gright,
          -3*fleft+3*fright-2*gleft - gright,
          gleft,
          fleft)
      qa=3*S[1]; qb=2*S[2]; qc=S[3]
      crit  = (-qb + c(1,-1) * sqrt(qb^2. - 4*qa*qc)) / (2*qa)
      ts = crit[which.min(crit * (crit * (S[1] * crit + S[2]) + S[3]))]
      xc = (1-ts) * a + ts * beta
      xq = a - gleft * (a+beta)^2. / (2.*(fright - fleft - gleft* (a+beta)))##(gright*a - gleft*beta)/(gright - gleft)
      if (abs(xc - a) < abs(xq - a)) {
        xf = xc
      } else {
        xf = xc + (xq - xc) / 2.
      }
    } else if (gleft*gright < 0.) {
      S = c(2*fleft - 2*fright + gleft + gright,
          -3*fleft+3*fright-2*gleft - gright,
          gleft,
          fleft)
      qa=3*S[1]; qb=2*S[2]; qc=S[3]
      crit  = (-qb + c(1,-1) * sqrt(qb^2. - 4*qa*qc)) / (2*qa)
      ts = crit[which.min(crit * (crit * (S[1] * crit + S[2]) + S[3]))]
      xc = (1-ts) * a + ts * beta
      xq = (gright*a - gleft*beta)/(gright - gleft)
      if (abs(xc - beta) > abs(xq - beta)) {
        xf = xc
      } else {
        xf = xq
      }
    } else {
      xf = (gright*a - gleft*beta)/(gright - gleft)
    }
  } else {
    xf        = a - gleft * (a+beta)^2. / (2.*(fright - fleft - gleft* (a+beta)))
  }
  xf = min(a+(beta-a)*19/20., max(a+(beta-a)/20., xf))
  xhalf = (a+beta)/2.
  sig   = sign(xhalf - xf)
  del   = kap1 * (beta - a)^kap2
  xt    = if (del <= abs(xhalf - xf))  xf + sig*del
          else                         xhalf
  r     = eps * 2^(nhalf+n0-j)-(beta-a)/2.
  if (abs(xt - xhalf) <= r) xtil= xt
  else                      xtil= xhalf - sig*r
  xtil
}

weak_wolfe_search =  function(fn, gr, x0, dk, f0, g0, c1, c2, extra_arg,tee0=1., n0=1, kap1=.1, kap2=2.565672, loopmax = 30L, trace = FALSE, ...) {
  a    = 0.
  beta = Inf
  tmp  = numeric(1L)
  tee1 = tee0
  fk1  = Inf
  gk1  = numeric(1L)
  gdir = numeric(1L)
  fleft    = f0
  fright   = numeric(1L)
  gdirleft = gdir0 = sum(g0 * dk)
  a_p        = a
  gdirleft_p = gdir0
  fleft_p    = f0
  gdirright = numeric(1L)
  loopcnt = 0L
  j       = 0L
  feval   = 0L
  geval   = 0L
  printf(trace, '%-8s %14s %15s %15s %15s %15s %15s %15s\n', "FINDWOLF", 'a', 'b', 't',
         'f_a', 'f_b', 'g_a', 'g_b')
  while (TRUE) {
    xk1 = x0 + tee1 * dk
    fk1=do.call(fn, c(list(xk1),extra_arg))
    feval = feval + 1L
    if (fk1 > f0 + c1 * tee1 * sum(g0*dk)) {
      if (beta == Inf) {
        beta    = tee1
        eps     = (beta-a)/(2^12)
        nhalf   = ceiling(log((beta-a)/(2*eps), 2))
        nmax    = nhalf + n0
        j       = 0L
      } else {
        beta = tee1
      }
      
      fright = fk1
      gk1 = do.call(gr, c(list(xk1),extra_arg))
      gdirright = gdir = sum(gk1*dk) 
      geval = geval + 1L

      tee1 = itppt(j, fleft, fright, gdirleft, gdirright, a, beta, nhalf, eps, n0, kap1, kap2)
      printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "RIGHT", a, beta, tee1,
                    fleft, fright, gdirleft, gdirright)
      j = j + 1L
      if (j >= nmax) {
        eps     = (beta-a)/(2^12)
        nhalf   = ceiling(log((beta-a)/(2*eps), 2))
        nmax    = nhalf + n0
        j       = 0L
      }
    } else if ((gdir=sum((gk1 = do.call(gr, c(list(xk1),extra_arg)))*dk)) < c2 * gdir0) {
      geval = geval + 1L
      a_p        = a
      gdirleft_p = gdirleft
      fleft_p    = fleft
      a = tee1
      gdirleft = gdir
      fleft    = fk1
      tee1 = if (beta < Inf)  {
               r = itppt(j, fleft, fright, gdirleft, gdirright, a, beta, nhalf, eps, n0, kap1, kap2)
               printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "LEFT", a, beta, r,
                      fleft_p, fleft, gdirleft_p, gdirleft)
               j = j + 1L
               if (j >= nmax) {
                 eps     = (beta-a)/(2^12)
                 nhalf   = ceiling(log((beta-a)/(2*eps), 2))
                 nmax    = nhalf + n0
                 j       = 0L
               }
               r
             } else {
               r = extrap(a_p, a, fleft_p, fleft, gdirleft_p, gdirleft, c2)
               printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "EXTRAP", a, beta, r,
                      fleft_p, fleft, gdirleft_p, gdirleft)
               r
             }
    } else {
      printf(trace, '%-7s %15s %15s %15.9f %15.9f %15s %15.9f %15s\n', "SATISFY", '', '', tee1,
             fk1, '', gdir, '')
      geval = geval + 1L
      break
    }
    loopcnt = loopcnt + 1L
    if (loopcnt >= loopmax) {
      return('Line search has failed find a point that satisfies the weak Wolfe condition.')
    }
  }
  list(tee1, xk1, fk1, gk1, gdir, feval, geval)
}

extrap = function (a0, a, f0, fleft, gdir0, gdirleft, c2) {
  if (gdirleft > 0 || gdir0 > 0)
    print(c(gdir0, gdirleft))
  
  if (gdir0 < 0 && gdirleft > gdir0) {
    r = (a0*(gdirleft+c2*gdir0) - a*(gdir0+c2*gdir0)) / (gdirleft - gdir0)
  } else {
    ## If the gradient is steeper downward it is very difficult to get a reasonable
    ## extrapolation. So just simply multiply by a constant and finger crossed.
    return(3.*a)
  }
}



##TRACE=F
##tester = function (par, fn, gr) {
##  res = list()
##  res[['MINE']] = ttcg(par = par, fn = fn, gr = gr, method='TTDES', control=list(trace=TRACE))
##  res[['rconjgrad']] = rconjgrad::conj_grad(par=par, fn = fn, gr = gr, max_iter = 30L, line_search='mt')
##  res[['rconjgrad']]$convergence = as.integer(NA)
##  #res[['optimCG']] = optim(par = par, fn = fn, gr = gr, method='CG', control = list(maxit = 1000))
##  res[['optimLBFGSB']] = optim(par = par, fn = fn, gr = gr, method='L-BFGS-B', control = list(maxit = 2000))
##  res[['lbfgsb3c']] = lbfgsb3c::lbfgsb3c(par = par, fn = fn, gr = gr, control = list(maxit = 2000))
##  if (is.null(res[['lbfgsb3c']]$convergence)) {
##    res[['lbfgsb3c']]$convergence = -1L
##  }
##  do.call(rbind, lapply(res, function (r) {
##    x = as.data.frame(r[c('value', 'convergence')])
##    x[['function']] = r[['counts']][1]
##    x[['gradient']] = r[['counts']][2]
##    x
##  }))
##}
##
##set.seed(3)
##
##test_fns = list()
##test_fns[['Booth']] =
##  list(par= c(-3,-3),
##       fn = function (x) (x[1]+2.*x[2]-7.)^2. + (2.*x[1]+x[2]-5)^2.,
##       gr = NULL)
##test_fns[['Beale']] =
##  list(par= c(0.2,-3),
##       fn = function (x) (1.5-x[1]+x[1]*x[2])^2. + (2.25-x[1]+x[1]*x[2]^2.)^2. + (2.625-x[1]+x[1]*x[2]^3.)^2.,
##       gr = NULL)
##test_fns[['Bohachevsky']] =
##  list(par=c(-2.7,1.9),
##       fn = function (x) x[1]^2. + 2*x[2]^2. -.3*cos(3.*pi*x[1])*cos(4.*pi*x[2]) + .3,
##       gr = NULL)
##test_fns[['Six-Hump Camel']] =
##  list(par=c(2,2),
##       fn = function (x) (4-2.1*x[1]^2.+(x[1]^4.)/3.)*x[1]^2. + x[1]*x[2] + (-4.+4.*x[2]^2.)*x[2]^2.,
##       gr = NULL)
##test_fns[['Matyas']] =
##  list(par=c(8,-5),
##       fn = function (x) 0.26*(x[1]^2. + x[2]^2.)-0.48*x[1]*x[2],
##       gr = NULL)
##test_fns[['McCormick']] =
##  list(par=c(-1.4,3.5)*4,
##       fn = function (x) sin(x[1]+x[2])+(x[1]-x[2])^2-1.5*x[1]+2.5*x[2]+1.,
##       gr = NULL)
##test_fns[['Zakharov']] =
##  list(par=c(-2.318,-0.3,-1.9,-2.01, -3.32,-1.2, 4.0, 0.3,1.2,-0.1,-4,-.3,-1.2,1.2,0.5,-2.1,2.193,
##             2.318,-3.3,-1,-0.01, -3.32,-1.2, 4.0, 0.3,1.2,-0.1,-4,-.3,-1.2,1.2,0.5,-2.1,2.193)*3,
##       fn = function (x) sum(x^2.)+(sum(0.5*(1:length(x))*x))^2. + (sum(0.5*(1:length(x))*x))^4.,
##       gr = NULL)
##test_fns[['Rosenbrock2']] =
##  list(par= c(1.1,-4.5),
##       fn = function (x) sum(100 * (x[2:length(x)] - x[1:(length(x)-1L)]^2L)^2L + (x[1:(length(x)-1L)])^2L),
##       gr = NULL)
##test_fns[['Rosenbrock50']] =
##  list(par= rnorm(50)*3,
##       fn = function (x) sum(100 * (x[2:length(x)] - x[1:(length(x)-1L)]^2L)^2L + (x[1:(length(x)-1L)])^2L),
##       gr = NULL)
##test_fns[['Rosenbrock70']] =
##  list(par= rnorm(380)*4,
##       fn = function (x) sum(100 * (x[2:length(x)] - x[1:(length(x)-1L)]^2L)^2L + (x[1:(length(x)-1L)])^2L),
##       gr = NULL)
##nqm = rnorm(50000)*2
##test_fns[['NearQuadraticHuge']] =
##  list(par= rnorm(50000)*4,
##       fn = function (x) (sum((x - nqm)^2.) + exp(-sum(x*x)) + .02 * sum(abs(x)))/10000,
##       gr = function (x) (2*(x - nqm) + exp(-sum(x*x))*-2.*x + .02*sign(x))/10000)
##for (i in seq_along(test_fns)) {
##  r = tester(test_fns[[i]]$par, test_fns[[i]]$fn,
##             if (is.null(test_fns[[i]]$gr))
##               function (x) numDeriv::grad(test_fns[[i]]$fn, x)
##             else test_fns[[i]]$gr)
##  cat('----', names(test_fns)[i], '\n')
##  print(r)
##}


