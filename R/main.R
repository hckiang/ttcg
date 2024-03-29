#' ttcg: Three-term conjugate gradient minimization algorithms
#'
#' Some accelerated three-term conjugate gradient algorithms implemented purely in R with the same user interface as optim(). The search directions and acceleration scheme are described in Andrei, N. (2013) <doi:10.1016/j.amc.2012.11.097>, Andrei, N. (2013) <doi:10.1016/j.cam.2012.10.002>, and Andrei, N (2015) <doi:10.1007/s11075-014-9845-9>. Line search is done by a hybrid algorithm incorporating the ideas in Oliveia and Takahashi (2020) <doi:10.1145/3423597> and More and Thuente (1994) <doi:10.1145/192115.192132>.
#' 
#' @author Hao Chi Kiang, \email{hello@hckiang.com}
#' @docType package
#' @aliases ttcg-package
#' @importFrom numDeriv grad
'_PACKAGE'

force_name = function (X, n) {
  names(X) = n
  X
}
my_is_na   = function (obj) is.na(obj) & !(is.nan(obj))  ## Because is.na(NaN)==T and identical() isn't vectorized.
maybe      = function (L, name, default)  if (is.null(L[[name]])) default else L[[name]]
printf     = function (trace, fmt, ...)   if (trace) cat(sprintf(fmt, ...))
big_brother= function (f, g) {
  list(fn=function (xk, extra_arg, penv, reskey, fevlkey) {
         penv[[fevlkey]] = penv[[fevlkey]] + 1L
         penv[[reskey]]  = do.call(f, c(list(xk*penv[['parscale']]), extra_arg))
         mode(penv[[reskey]]) = 'double'
         if ({length(penv[[reskey]]) != 1L    && { 'Objective function has returned a vector with length >= 1'                            -> penv[['msg']]; T }} ||
             {my_is_na(penv[[reskey]])        && { 'Objective function has returned NA or something that cannot be coerced to `double`'   -> penv[['msg']]; T }} ||
             {identical(penv[[reskey]],-Inf)  && { 'Objective function has returned -Inf'                                                 -> penv[['msg']]; T }}) {
           penv[['convergence']] = 101L
           F
         } else {
           penv[[reskey]] = penv[[reskey]]/penv[['fnscale']]
           T
         }
       },
       gr=function (xk, extra_arg, penv, reskey, gevlkey) {
         penv[[gevlkey]] = penv[[gevlkey]] + 1L
         penv[[reskey]]  = do.call(g, c(list(xk*penv[['parscale']]), extra_arg))
         mode(penv[[reskey]]) = 'double'
         if ({length(penv[[reskey]]) != length(xk) && { 'Gradient function has returned a vector with an incorrect length'                 -> penv[['msg']]; T }} ||
             {any(my_is_na(penv[[reskey]]))        && { 'Gradient vector contains NA or something that cannot be coerced to `double`'      -> penv[['msg']]; T }}) {
           penv[['convergence']] = 101L
           F
         } else {
           penv[[reskey]] = penv[[reskey]]/penv[['fnscale']]*penv[['parscale']]
           T
         }
       })
}

#' Accelerated three-term conjugate gradient optimization with restart
#'
#' The \code{ttcg} function minimizes a given function using several Neculai Andrei's
#' three-term conjugate gradient algorithms.
#'
#' By default, the algorithm stops when any one of the following three convergence tests
#' is satisfied. (1) The squared Euclidean norm of the squared Euclidean norm of the
#' gradient is smaller than a tolerance; (2) The infinity norm of the gradient is smaller
#' than a tolerance; (3) \eqn{|f_{k+1} - f_k| < \epsilon * (1 + |f_k|)}. These three
#' tolerances can be set in the \code{control} argument, and turnt off by setting them
#' to any negative values. If all three were turnt off, the algorithm may never stop.
#' 
#' For minimization problems, in which \code{fnscale} is positive, the objective function
#' can return \code{NaN} or \code{Inf} but \code{NA} or \code{-Inf} will results in the algorithm being stopped
#' immediately, because \code{-Inf} means the function is unbounded below and any arithmetic error, such as dividing
#' by zero, should be coded as \code{NaN}; while \code{NA} should signify programming error. For maximization problems,
#' \code{-Inf} instead of \code{Inf} is allowed. The gradient function, similarly, can contain \code{NaN}, \code{Inf},
#' or \code{-Inf}, but returning \code{NA} will stop the algorithm.
#'
#' The \code{method} argument specifies how the search direction in each step is computed.
#' Please see the three Neculai Andrei's three papers in the citation section for more
#' detailed description. An acceleration scheme and a restart procedure are implemented
#' according to his three papers. Line search is done by a bisection-like
#' weak-Wolfe search described in Oliveira and Takahashi's (2020) interpolate-truncate-project
#' algorithm, but replacing their gradient-secant interpolation with some of More-Thuente's (1994)
#' cubic interpolation idea.
#' 
#' The \code{control} argument is a list that can contain any of the following named
#' element:
#'
#'     \describe{
#'       \item{maxit}{The maximal number of iteration. Default is 500.}
#'       \item{fnscale}{Scalar to divide fn and gr during optimization. If negative, turns the problem into a maximization problem. Optimization is performed on \code{fn(par)/fnscale}.}
#'       \item{parscale}{A vector of scaling values for the parameters. Optimization is performed on \code{par/parscale} and these should be comparable in the sense that a unit change in any element produces about a unit change in the scaled value.}
#'       \item{gl2tol}{A positive small number. The iteration will be terminated if the
#'                     squared Euclidean norm of the gradient is smaller than this number. Default is
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
#'                 and return a vector of the same length as \code{par}. If it is \code{NULL} then
#'                 numerical finite difference is used to obtain the gradient.
#' @param method   A character string, one of 'TTDES', 'TTCG', 'THREECG'. This determines how the
#'                 line search direction is computed. 'TTCG' is the default method.
#' @param control  A list of control parameters. See Details.
#' @param ...      Extra arguments to be passed to \code{fn}
#' @return         A list containing the following named elements.
#'                 \describe{
#'                   \item{par}{The optimal parameter.}
#'                   \item{value}{The optimal function value.}
#'                   \item{counts}{An integer vector containing the number of function and gradient calls
#'                                 used during the optimization.}
#'                   \item{convergence}{An integer indicating convergence status. '0' means successful convergence; '1'
#'                                      means \code{maxit} has been reached; '2' means a line search failure in which
#'                                      a point that satisfies the weak Wolfe condition is not found. Among other possibilities,
#'                                      this may happen when the function is unbounded below or the function is
#'                                      non-differentiable.)}
#'                   \item{message}{A character string giving additional message.}
#'                 }
#' @references   Andrei, N. (2013). On three-term conjugate gradient algorithms for unconstrained optimization. Applied Mathematics and Computation, 219(11), 6316-6327.
#' @references   Andrei, N. (2013). A simple three-term conjugate gradient algorithm for unconstrained optimization. Journal of Computational and Applied Mathematics, 241, 19-29.
#' @references   Andrei, N. (2015). A new three-term conjugate gradient algorithm for unconstrained optimization. Numerical Algorithms, 68(2), 305-321.
#' @references   Oliveira, I. F., & Takahashi, R. H. (2020). An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality. ACM Transactions on Mathematical Software (TOMS), 47(1), 1-24.
#' @references   More, J. J., & Thuente, D. J. (1994). Line search algorithms with guaranteed sufficient decrease. ACM Transactions on Mathematical Software (TOMS), 20(3), 286-307.
#' @examples
#' nqm = rnorm(500)*2
#' fn  = function (x,nqm1) sum((x - nqm1)^2.)
#' gr  = function (x,nqm1) 2.*(x - nqm1)
#' r   = ttcg(par = rnorm(500)*4., fn = fn, gr = gr, method='TTDES', nqm=nqm)
#' all.equal(r$value, 0.0)
#' @export
ttcg = function (par, fn, gr = NULL, method='TTCG', control = list(), ...) {
  parnames  = names(par)
  mode(par) = 'double'
  dim(par)  = NULL
  npar = length(par)

  fevl  = 0L
  gevl  = 0L
  xk    = par
  fk    = as.double(NA)
  gk    = as.double(rep(NA, npar))
  delayedAssign('xk1', xk)
  delayedAssign('fk1', fk)
  delayedAssign('gk1', gk)
  delayedAssign('FinalRet', list(par = force_name(xk1, parnames), value = fk1, counts=c('function'=fevl, 'gradient'=gevl), convergence=convergence, message=msg))

  if (npar == 0L)                {  convergence=106L;   msg='The initial parameter has zero length';    return(FinalRet)  }
  if (any(is.na(xk)))            {  convergence=107L;   msg='The initial parameter contains NA or NaN'; return(FinalRet)  }
  if (any(is.infinite(xk)))      {  convergence=108L;   msg='The initial parameter contains infinity';  return(FinalRet)  }

  maxit   = maybe(control, 'maxit', 500L)
  gl2tol  = maybe(control, 'gl2tol', min(1e-9, length(par)*1e-10))  ## each of the gk1 < 1e-6
  gmaxtol = maybe(control, 'gmaxtol', 1e-6)
  ftol    = maybe(control, 'ftol', 1e-9)
  c1      = maybe(control, 'c1',  1e-3)
  c2      = maybe(control, 'c2',  .08)
  trace   = maybe(control, 'trace', FALSE)
  fnscale = maybe(control, 'fnscale',  1.0)
  parscale= maybe(control, 'parscale', rep(1.0, npar))

  mode(fnscale) = 'double'
  if (length(fnscale)!=1L)       {  convergence=110L;   msg='fnscale has non-one length';               return(FinalRet)  }
  if (is.na(fnscale))            {  convergence=111L;   msg='fnscale cannot be coerced to double';      return(FinalRet)  }
  if (is.infinite(fnscale))      {  convergence=112L;   msg='fnscale cannot be infinite';               return(FinalRet)  }
  if (fnscale == 0.0)            {  convergence=113L;   msg='fnscale cannot be zero';                   return(FinalRet)  }
  mode(parscale) = 'double'     
  if (length(parscale)!=npar)    {  convergence=120L;   msg='parscale has different length from par';   return(FinalRet)  }
  if (any(is.na(parscale)))      {  convergence=121L;   msg='parscale cannot be coerced to double';     return(FinalRet)  } 
  if (any(is.infinite(parscale))){  convergence=122L;   msg='parscale cannot contain infinity';         return(FinalRet)  }
  if (any(parscale == 0.0))      {  convergence=123L;   msg='parscale cannot contain zero';             return(FinalRet)  }
  xk = xk / parscale
  delayedAssign('FinalRet', list(par = force_name(xk1*parscale, parnames), value = fk1*fnscale, counts=c('function'=fevl, 'gradient'=gevl), convergence=convergence, message=msg))

  if (! (method %in% c('TTDES', 'THREECG', 'TTCG')))
    stop('Invalid method. `method` should be one of `TTDES`, `THREECG`, `TTCG`')

  fn_orig = fn
  gr_orig = gr
  R = big_brother(fn_orig, if (is.null(gr_orig)) function (x, ...) numDeriv::grad(func=fn_orig, x=x, method="Richardson", side=NULL, method.args=list(), ...)
                           else gr_orig)
  fn = R[['fn']]
  gr = R[['gr']]
  rm(R)

  extra_arg = list(...)
  delayedAssign('fk1', fk)
  delayedAssign('gk1', gk)
  delayedAssign('xk1', xk)
  penv  = environment()
  if (! fn(xk, extra_arg, penv, 'fk', 'fevl'))  return(FinalRet)
  if (is.nan(fk) || !is.finite(fk)) {
    convergence = 9L
    msg = 'Objective function evaluated to NA, NaN or infinite at the initialization point'
    return(FinalRet)
  }
  if (! gr(xk, extra_arg, penv, 'gk', 'gevl'))  return(FinalRet)
  if (any(is.nan(gk) | !is.finite(gk))) {
    convergence = 9L
    msg = 'Gradient function evaluated to a vector that contains NA, NaN or infinite at the initialization point'
    return(FinalRet)
  }
  dk    = -gk
  k     = 1L
  itcnt = 0L
  convergence = 1L
  msg   = NULL
  if ( {(max(abs(gk)) < gmaxtol)        && { "Gradient's infinite norm < gmaxtol"   -> msg; T }} ||
       {(sum(gk*gk) < gl2tol)           && { "Gradient's squared 2-norm < gl2tol"   -> msg; T }}  )
    return(list(par = xk, value = fk, counts=c('function'=fevl, 'gradient'=gevl), convergence=0L, message=msg))

  compute_dk1 = switch(method,
                  "TTDES"=quote({
                    yksk= sum(yk*sk)
                    skgk1=sum(sk*gk1)
                    omega  = (2./sum(sk*sk)) *sqrt(max(0., sum(sk*sk)*sum(yk*yk) - yksk^2.))
                    dk1 = -gk1 + (sum(yk*gk1) - omega*skgk1)/yksk * sk - skgk1/yksk*yk
                  }),
                  "THREECG"=quote({
                    yksk= sum(yk*sk)
                    skgk1=sum(sk*gk1)
                    dk1 = -gk1 - ((1.+sum(yk*yk)/yksk)*skgk1/yksk - sum(yk*gk1)/yksk)*sk
                          -skgk1/yksk * yk
                  }),
                  "TTCG"=quote({
                    yksk= sum(yk*sk)
                    skgk1=sum(sk*gk1)
                    dk1 = -gk1 - ((1.+2.*sum(yk*yk)/yksk)*(skgk1/yksk) - sum(yk*gk1)/yksk)*sk
                          - (skgk1/yksk) * yk
                  }))

  tee0_init = 1./sqrt(sum(gk*gk))
  while ((itcnt+1L->itcnt) < maxit) {
    if (!weak_wolfe_search(fn, gr, xk, dk, f0=fk, g0=gk, extra_arg=extra_arg, penv=penv, c1=c1,c2=c2, tee0=tee0_init, trace = trace))
      break
    if ( {(max(abs(gk1)) < gmaxtol)        && { "Gradient's infinite norm < gmaxtol"   -> msg; T }} ||
         {(sum(gk1*gk1) < gl2tol)          && { "Gradient's squared 2-norm < gl2tol"   -> msg; T }} ||
         {(abs(fk1-fk) < ftol*(1+abs(fk))) && { "|f_{k+1} - f_k| < ftol * (1 + |f_k|)" -> msg; T }}  ) {
      convergence = 0L
      break
    }
    yk  = gk1 - gk
    bbar= sum(yk*dk)
    ## Sacrifice one function evaluation and one gradient for acceleration scheme.
    if (bbar>0) {
      xk1bkup = xk1
      fk1bkup = fk1
      gk1bkup = gk1
      abar = sum(gk * dk)
      xik  = -abar / bbar
      xk1  = xk + xik * tee1 * dk
      if (! fn(xk1, extra_arg, penv, 'fk1', 'fevl')) break
      if (! gr(xk1, extra_arg, penv, 'gk1', 'gevl')) break
      if (is.nan(fk1) || fk1 == Inf || any(is.nan(gk1) | is.infinite(gk1))) {
        xk1 = xk1bkup
        fk1 = fk1bkup
        gk1 = gk1bkup
        printf(trace, 'Acceleration resulted in NaN/Inf/-Inf. Giving up the acceleration.\n')
      } else {
        if ( {(max(abs(gk1)) < gmaxtol)        && { "Gradient's infinite norm < gmaxtol"   -> msg; T }} ||
             {(sum(gk1*gk1) < gl2tol)          && { "Gradient's squared 2-norm < gl2tol"   -> msg; T }} ||
             {(abs(fk1-fk) < ftol*(1+abs(fk))) && { "|f_{k+1} - f_k| < ftol * (1 + |f_k|)" -> msg; T }}  ) {
          convergence = 0L
          break
        }
        yk  = gk1 - gk
      }
    }
    sk  = xk1 - xk
    eval(compute_dk1)
    tee0_init = tee1 * sqrt(sum(dk*dk)/sum(dk1*dk1))
    printf(trace, '***** ')
    if (abs(sum(gk1*gk)) > 0.2*(gknorm=sum(gk1*gk1))) {
      dk1 = -gk1
      printf(trace, '%-9s  ', 'RESTART')
    } else
      printf(trace, '%-9s  ', 'CARRYON')
    
    printf(trace, 'CURR_STEP=%-10f   GUESSED_STEP=%-10f   FVAL=%-10f   GNORM=%-10f\n', tee1, tee0_init, fk1, sqrt(gknorm))
    dk = dk1
    gk = gk1
    fk = fk1
    xk = xk1
  }
  if (convergence == 1L)
    msg = 'Maximum number iteration has been reached.'
  printf(trace, '\nThe algorithm has stopped at iteratiion %d, because:\n', itcnt)
  printf(trace, '%s\n', msg)
  return(FinalRet)
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
      xf = if (abs(xc - a) < abs(xq - a)) xc else xc + (xq - xc) / 2.
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
      xf = if (abs(xc - beta) > abs(xq - beta)) xc else xq
    } else
      xf = (gright*a - gleft*beta)/(gright - gleft)
  } else
    xf = a - gleft * (a+beta)^2. / (2.*(fright - fleft - gleft* (a+beta)))
  xf    = min(a+(beta-a)*19/20., max(a+(beta-a)/20., xf))
  xhalf = (a+beta)/2.
  sig   = sign(xhalf - xf)
  del   = kap1 * (beta - a)^kap2
  xt    = if (del <= abs(xhalf - xf))  xf + sig*del
          else                         xhalf
  r     = eps * 2^(nhalf+n0-j)-(beta-a)/2.
  if (abs(xt - xhalf) <= r)   xt
  else                        xhalf - sig*r
}

weak_wolfe_search =  function(fn, gr, x0, dk, f0, g0, c1, c2, extra_arg, penv, tee0=1., n0=1, kap1=.1, kap2=2.565672, loopmax = 30L, trace = FALSE, ...) {
  a              = 0.
  beta           = Inf
  tmp            = numeric(1L)
  penv[['tee1']] = tee0
  penv[['fk1']]  = Inf
  penv[['gk1']]  = numeric(1L)
  penv[['gdir']] = numeric(1L)
  fleft          = f0
  fright         = numeric(1L)
  gdirleft       = gdir0 = sum(g0 * dk)
  a_p            = a
  gdirleft_p     = gdir0
  fleft_p        = f0
  gdirright      = numeric(1L)
  loopcnt        = 0L
  j              = 0L
  eps            = as.double(NA)
  nhalf          = as.integer(NA)
  nmax           = as.integer(NA)
  printf(trace, '%-8s %14s %15s %15s %15s %15s %15s %15s\n', "FINDWOLF", 'a', 'b', 't',
         'f_a', 'f_b', 'g_a', 'g_b')
  reset_j = quote({
    eps     = (beta-a)/(2^12)
    nhalf   = ceiling(log((beta-a)/(2*eps), 2))
    nmax    = nhalf + n0
    j       = 0L
  })
  repeat {
    penv[['xk1']]  = x0 + penv[['tee1']] * dk
    if (! fn(penv[['xk1']], extra_arg, penv, 'fk1', 'fevl')) return(F)
    if (! gr(penv[['xk1']], extra_arg, penv, 'gk1', 'gevl')) return(F)
    penv[['gdir']] = sum(penv[['gk1']]*dk)
    while (is.nan(penv[['fk1']]) || penv[['fk1']] == Inf || any(is.nan(penv[['gk1']]) | is.infinite(penv[['gk1']]))) {
      penv[['tee1']] = penv[['tee1']] / 2.
      penv[['xk1']]  = x0 + penv[['tee1']] * dk
      if (! fn(penv[['xk1']], extra_arg, penv, 'fk1', 'fevl')) return(F)
      if (! gr(penv[['xk1']], extra_arg, penv, 'gk1', 'gevl')) return(F)
      penv[['gdir']] = sum(penv[['gk1']]*dk)    ## Lazy evaluation. Will be run only if trace=T
      printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "NAN/INF", a, penv[['xk1']], penv[['tee1']],
                    fleft, penv[['fk1']], gdirleft, penv[['gdir']])
      if ((loopcnt+1L->loopcnt) >= loopmax) {
        penv[['convergence']] = 7L
        penv[['msg']]         = 'Line search has failed to find a point at which neither the function nor gradient is Inf or NaN, nor is the gradient not -Inf'
        return(F)
      }
    }
    if (penv[['fk1']] > f0 + c1 * penv[['tee1']] * sum(g0*dk)) {
      if (beta == Inf) {
        beta    = penv[['tee1']]
        eval(reset_j)
      } else
        beta = penv[['tee1']]

      fright = penv[['fk1']]
      gdirright = penv[['gdir']]
      penv[['tee1']] = itppt(j, fleft, fright, gdirleft, gdirright, a, beta, nhalf, eps, n0, kap1, kap2)
      printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "RIGHT", a, beta, penv[['tee1']],
                    fleft, fright, gdirleft, gdirright)
      if ((j+1L->j) >= nmax) eval(reset_j)
    } else if (penv[['gdir']] < c2 * gdir0) {
      gdirleft_p = gdirleft
      fleft_p    = fleft
      a_p        = a
      a          = penv[['tee1']]
      gdirleft   = penv[['gdir']]
      fleft      = penv[['fk1']]
      penv[['tee1']] = if (beta < Inf)  {
               r = itppt(j, fleft, fright, gdirleft, gdirright, a, beta, nhalf, eps, n0, kap1, kap2)
               printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "LEFT", a, beta, r,
                      fleft_p, fleft, gdirleft_p, gdirleft)
               if ((j+1L->j) >= nmax) eval(reset_j)
               r
             } else {
               r = extrap(a_p, a, fleft_p, fleft, gdirleft_p, gdirleft, c2)
               printf(trace, '%-7s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n', "EXTRAP", a, beta, r,
                      fleft_p, fleft, gdirleft_p, gdirleft)
               r
             }
    } else {
      printf(trace, '%-7s %15s %15s %15.9f %15.9f %15s %15.9f %15s\n', "SATISFY", '', '', penv[['tee1']],
             penv[['fk1']], '', penv[['gdir']], '')
      break
    }
    if ((loopcnt+1L->loopcnt) >= loopmax) {
      penv[['convergence']] = 2L
      penv[['msg']]         = 'Line search has failed to find a point that satisfies the weak Wolfe condition.'
      return(F)
    }
  }
  T
}

extrap = function (a0, a, f0, fleft, gdir0, gdirleft, c2) {
  ## if (gdirleft > 0 || gdir0 > 0)     print(c(gdir0, gdirleft))
  if (gdir0 < 0 && gdirleft > gdir0)   (a0*(gdirleft+c2*gdir0) - a*(gdir0+c2*gdir0)) / (gdirleft - gdir0)
  else                                 3.*a
}

