\name{fitgrouped2}
\alias{fitgrouped2}
\title{Estimating parameters of the three-parameter Birnbaum-saunders (BS), generalized exponential (GE), and Weibull distributions fitted to grouped data}
\description{Suppose a sample of \eqn{n} independent observations each follows a three-parameter BS, GE, or Weibull distributions have been divided into \eqn{m} separate groups of the form \eqn{(r_{i-1},r_i]}, for \eqn{i=1,\dots,m}. So, the likelihood function is given by
\deqn{
L(\Theta)=\frac{n!}{f_{1}!f_{2}!\dots f_{m}!}\prod_{i=1}^{m}\Bigl[F\bigl(r_{i}\big|\Theta\bigr)-F\bigl(r_{i-1}\big|\Theta\bigr)\Bigr]^{f_i},}
where the \eqn{r_0} is the lower bound of the first group, \eqn{r_m} is the upper bound of the last group, and \eqn{f_i} is the frequency of observations within \eqn{i}-th group provided that \eqn{n=\sum_{i=1}^{m}f_{i}}.
}
\usage{fitgrouped2(r, f, param, start, cdf, pdf, method = "Nelder-Mead", lb = 0, ub = Inf
            , level = 0.05)}
\arguments{
  \item{r}{A numeric vector of length \eqn{m+1}. The first element of \eqn{r} is lower bound of the first group and other \eqn{m} elements are upper bound of the \eqn{m} groups. We note that upper bound of the \eqn{(i-1)}-th group is the lower bound of the \eqn{i}-th group, for \eqn{i=2,\dots,m}. The lower bound of the first group and upper bound of the \eqn{m}-th group are chosen arbitrarily.
  }
  \item{f}{A numeric vector of length \eqn{m} containing the group's frequency.}
\item{param}{Vector of the of the family parameter's names.}
\item{start}{Vector of the initial values.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{method}{The method for the numerically optimization that includes one of \code{CG}, \code{Nelder-Mead}, \code{BFGS}, \code{L-BFGS-B}, \code{SANN}.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
\item{level}{Significance level for constructing asymptotic confidence interval That is \code{0.05} by default for constructing a \code{95\%} confidence interval.}
}
\value{A two-part list of objects given by the following:
\enumerate{
\item Maximum likelihood (ML) estimator for the parameters of the fitted family to the gropued data, asymptotic standard error of the ML estimator, lower bound of the asymptotic confidence interval, and upper bound of the asymptotic confidence interval at the given level.
\item A sequence of goodness-of-fit measures consist of Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), and Kolmogorov-Smirnov (\code{KS}) statistics.
}
}
\author{Mahdi Teimouri}
\examples{
    r <- c(2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5)
    f <- c(33, 111, 168, 147, 96,  45, 18, 4, 0)
param <- c("alpha", "beta", "mu")
  pdf <- quote( alpha/beta*((x-mu)/beta)^(alpha-1)*exp( -((x-mu)/beta)^alpha ) )
  cdf <- quote( 1-exp( -((x-mu)/beta)^alpha ) );
   lb <- 2
   ub <- Inf
start <-c(2, 3, 2)
level <- 0.05
fitgrouped2(r, f, param, start, cdf, pdf, method = "Nelder-Mead", lb = lb, ub = ub, level = 0.05)
}
