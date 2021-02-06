\name{fitgrouped1}
\alias{fitgrouped1}
\title{Estimating parameters of the three-parameter Birnbaum-saunders (BS), generalized exponential (GE), and Weibull distributions fitted to grouped data}
\description{Suppose a sample of \eqn{n} independent observations each follows a three-parameter BS, GE, or Weibull distributions have been divided into \eqn{m} separate groups of the form \eqn{(r_{i-1},r_i]}, for \eqn{i=1,\dots,m}. So, the likelihood function is given by
\deqn{
L(\Theta)=\frac{n!}{f_{1}!f_{2}!\dots f_{m}!}\prod_{i=1}^{m}\Bigl[F\bigl(r_{i}\big|\Theta\bigr)-F\bigl(r_{i-1}\big|\Theta\bigr)\Bigr]^{f_i},}
where the \eqn{r_0} is the lower bound of the first group, \eqn{r_m} is the upper bound of the last group, and \eqn{f_i} is the frequency of observations within \eqn{i}-th group provided that \eqn{n=\sum_{i=1}^{m}f_{i}}. The cdf of a three-parameter BS, GE, and Weibull distributions are given by
\deqn{F(x;\Theta)=\biggl(1-\exp \bigl\{-\beta(x-\mu)\bigr\} \biggr)^{\alpha},}
\deqn{F(x;\Theta)=\Phi\Biggl(\frac{\sqrt{\frac{x}{\beta}}-\sqrt{\frac{\beta}{x}}}{\alpha}\Biggr),}
and
\deqn{F(x;\Theta)=1- \exp \Bigl\{-\left(\frac{x-\mu}{\beta} \right)^{\alpha} \Bigr\},}

where \eqn{\Theta=(\alpha,\beta,\mu)^T}.}
\usage{fitgrouped1(r, f, family, method1, starts, method2)}
\arguments{
  \item{r}{A numeric vector of length \eqn{m+1}. The first element of \eqn{r} is lower bound of the first group and other \eqn{m} elements are upper bound of the \eqn{m} groups. We note that upper bound of the \eqn{(i-1)}-th group is the lower bound of the \eqn{i}-th group, for \eqn{i=2,\dots,m}. The lower bound of the first group and upper bound of the \eqn{m}-th group are chosen arbitrarily.
  }
  \item{f}{A numeric vector of length \eqn{m} containing the group's frequency.}
  \item{family}{Can be either \code{"birnbaum-saunders"}, \code{"ge"}, or \code{"weibull"}.}
  \item{method1}{A character string determining the method of estimation. It can be one of \code{"aml"}, \code{"em"} and \code{"ml"}. The short forms \code{"aml"}, \code{"em"}, and \code{"ml"} are described as follows.}
"\code{"aml"} (for method of approximated maximum likelihood (aml)),
"\code{"em"} (for method of expectation maximization (em)), and
"\code{"ml"} (for method of maximum likelihood (ml)).
\item{starts}{A numeric vector of the initial values for the shape, scale, and location parameters, respectively.}
\item{method2}{The method for optimizing the log-likelihood function. It invovles one of \code{"BFGS"}, \code{"Nelder-Mead"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.}
}
\details{If the method is \code{"em"}, then the initial values (\code{"starts"}) and the log-likelihood optimizing method (\code{"method2"}) are ignored.}
\value{A two-part list of objects given by the following:
\enumerate{
\item Estimated parameters of the three-parameter GE, Birnbaum-Saunders, or Weibull distribution fitted to the gropued data.
\item A sequence of goodness-of-fit measures consist of Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Anderson-Darling (\code{AD}), Chi-square (\code{Chi-square}),
Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and log-likelihood (\code{log-likelihood}) statistics.
}
}
\references{
G. J. McLachlan and T. Krishnan, 2007. \emph{The EM Algorithm and Extensions}, John Wiley & Sons.

A. P. Dempster, N. M. Laird, and D. B. Rubin, 1977. Maximum likelihood from incomplete data via the EM algorithm, \emph{Journal of the Royal Statistical Society, Series B (methodological)}, 1-38.

M. Teimouri and A. K. Gupta, 2012. Estimation Methods for the Gompertz–Makeham Distribution Under Progressively Type-I Interval Censoring Scheme, \emph{National Academy Science Letters}, 35(3).
}
\author{Mahdi Teimouri}
\examples{
r<-c(0,1,2,3,4,10)
f<-c(2,8,12,15,4)
starts<-c(2,2,0)
fitgrouped1(r,f,"birnbaum-saunders","em")
fitgrouped1(r,f,"weibull","ml",starts,"CG")
fitgrouped1(r,f,"ge","em")

}
