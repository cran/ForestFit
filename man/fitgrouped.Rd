\name{fitgrouped}
\alias{fitgrouped}
\title{Estimating parameters of the Weibull distribution fitted to grouped data}
\description{Suppose a sample of \eqn{n} independent observations each follows a three-parameter Weibull distribution have been divided into \eqn{m} separate groups of the form \eqn{(r_{i-1},r_i]}, for \eqn{i=1,\dots,m}. So, the likelihood function is given by
\deqn{
L(\Theta)=\frac{n!}{f_{1}!f_{2}!\dots f_{m}!}\prod_{i=1}^{m}\Bigl[F\bigl(r_{i}\big|\Theta\bigr)-F\bigl(r_{i-1}\big|\Theta\bigr)\Bigr]^{f_i},}
where the \eqn{r_0} is the lower bound of the first group, \eqn{r_m} is the upper bound of the last group, \eqn{f_i} is the frequency of observations within \eqn{i}-th group provided that \eqn{n=\sum_{i=1}^{m}f_{i}}, \eqn{\Theta=(\alpha,\beta,\mu)^T}, and \eqn{F(x;\Theta)=1- \exp \Bigl\{-\left(\frac{x-\mu}{\beta } \right)^{\alpha } \Bigr\}}, is the cdf of a three-parameter Weibull distribution.}
\usage{fitgrouped(r, f, method, starts)}
\arguments{
  \item{r}{A numeric vector of length \eqn{m+1}. The first element of \eqn{r} is lower bound of the first group and other \eqn{m} elements are upper bound of the \eqn{m} groups. We note that upper bound of the \eqn{(i-1)}-th group is the lower bound of the \eqn{i}-th group, for \eqn{i=2,\dots,m}. The lower bound of the first group and upper bound of the \eqn{m}-th group are chosen arbitrarily. If raw data are available, the smallest and largest observations are chosen for lower bound of the first group and upper bound of the \eqn{m}-th group, respectively.}
  \item{f}{A numeric vector of length \eqn{m} containing the group's frequency.}
  \item{method}{A character string determining the method of estimation. It can be one of \code{"aml"}, \code{"em"} and \code{"ml"}. As it is seen, the EM is used by default. The short forms \code{"aml"}, \code{"em"} and \code{"ml"} are described as follows.}
"\code{aml}" (for method of approximated maximum likelihood (aml)),
"\code{em}"  (for method of expectation maximization (em)), and
"\code{ml}"  (for method of maximum likelihood (ml)).
\item{starts}{A numeric vector of the initial values for shape, scale, and location parameters, respectively.}
}
\details{If the method is \code{"em"}, then the initial values are ignored.}
\value{A two-part list of objects given by the following.
\enumerate{
\item Estimated parameters of the three-parameter Weibull distribution fitted to the gropued data.
\item A sequence of two goodness-of-fit measures consist of Anderson Darling statistic (\code{AD}) and Chi-square (\code{Chi.square}) statistics.}
}
\references{
G. McLachlan and T. Krishnan, \emph{The EM Algorithm and Extensions}, John Wiley & Sons, 2007.

A. P. Dempster, N. M. Laird, and D. B. Rubin, Maximum likelihood from incomplete data via the EM algorithm, \emph{Journal of the Royal Statistical Society. Series B (methodological)}, 1-38, 1977.

M. Teimouri and A. K. Gupta, Estimation Methods for the Gompertz–Makeham Distribution Under Progressively Type-I Interval Censoring Scheme, \emph{National Academy Science Letters}, 35(3), 2012.
}
\author{Mahdi Teimouri}
%\note{}
%\seealso{}
\examples{
r<-c(0,1,2,3,4,10)
f<-c(2,8,12,15,4)
starts<-c(2,2,0)
fitgrouped(r,f,"em")
fitgrouped(r,f,"ml",starts)
}
