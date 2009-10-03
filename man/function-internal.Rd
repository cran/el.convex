\name{function-internal}
\alias{logstar}
\alias{dlogstar}
\alias{ddlogstar}
\alias{logelr}
\alias{dfpmin}
\alias{lnsrch}
\alias{bfgsmin}
\alias{mnbrak}
\alias{dbrent}
\alias{dlinmin}
\alias{frprmin}
\alias{steplength}
\title{Internal el.convex functions}
\description{
  Internal functions called by other functions in el.convex.
}
\usage{
logstar(z, n)
dlogstar(z, n)
ddlogstar(z, n)
logelr(z, lambda, n)
dfpmin(x, gtol, func, dfunc, ITMAX)
lnsrch(xold, fold, g, p, func, stpmax)
bfgsmin(x, gtol, func, dfunc, ITMAX)
mnbrak(f, a, b, glimit = 100)
dbrent(ax, bx, cx, f, df, tol, ITMAX = 100)
dlinmin(p, xi, func, dfunc, TOL = 2e-04)
frprmin(x, ftol, func, dfunc, ITMAX)
steplength(x, p, func, dfunc, u1, u2)
}
\details{
 These are not to be called directly by the users.
}
\keyword{internal}