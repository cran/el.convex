\name{el.test.newton}
\alias{el.test.newton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Empirical likelihood ratio test for the means using Newton method for optimization}
\description{
Compute the empirical likelihood ratio with the mean vector fixed at mu. 
The log empirical likelihood been maximized. 
el.test.newton uses simple Newton method for optimization.
}
\usage{
el.test.newton(x, mu, lam, maxit = 25, tol = 1e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix or vector containing the data, one row per observation.}
  \item{mu}{a numeric vector (of length = ncol(x)) to be tested as the mean vector of x above, as H0}
  \item{lam}{an optional vector of length = length(mu), the starting value of Lagrange multipliers, will use 0 if missing}
  \item{maxit}{an optional integer to control iteration when solve constrained maximization}
  \item{tol}{ an optional real value for convergence test}
}
\details{
If mu is in the interior of the convex hull of the observations x, then wts should sum to 1
}
\value{
\item{-2LLR }{the -2 loglikelihood ratio; approximate chisq distribution under H0}
\item{Pval }{the observed P-value by chi-square approximation}
\item{lambda}{the final value of Lagrange multiplier}
\item{nits}{number of iteration performed}
\item{wts}{weights on the observations}
\item{mu}{the means that are achieved}
}
\references{  A Study of Methods for Computing Empirical Likelihood
Numerical recipes in C
}
\author{Dan Yang, Dylan Small}
\section{Warning }{el.convex has not been thoroughly tested. Please report bugs.}
\examples{
x <- matrix(c(rnorm(50,mean=1), rnorm(50,mean=2)), ncol=2,nrow=50)
el.test.newton(x, mu=c(1,2))
}
\keyword{nonparametric}  
