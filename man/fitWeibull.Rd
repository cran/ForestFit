\name{fitWeibull}
\alias{fitWeibull}
\title{Estimating parameters of the Weibull distribution through classical methods}
\description{Estimates the parameters of the two- and three-parameter Weibull model with pdf and cdf given by
\deqn{f(x;\alpha,\beta,\theta)=\frac{\alpha}{\beta} \left(\frac{x-\theta}{\beta }\right)^{\alpha -1} \exp \biggl\{-\left(\frac{x-\theta}{\beta } \right)^{\alpha } \biggr\},}
and
\deqn{F(x;\alpha,\beta,\theta)=1- \exp \biggl\{-\left(\frac{x-\theta}{\beta } \right)^{\alpha } \biggr\},}
where \eqn{x>\theta}, \eqn{\alpha > 0}, \eqn{\beta >0} and \eqn{-\infty<\theta<\infty}. Here, the parameters \eqn{\alpha}, \eqn{\beta}, and \eqn{\theta} are known in the literature as the shape, scale, and location, respectively. If \eqn{\theta=0}, then \eqn{f(x;\alpha,\beta)} and \eqn{F(x;\alpha,\beta)} in above are the pdf and cdf of a two-parameter Weibull distribution, respectively.}
\usage{fitWeibull(mydata, location, method, starts)}
\arguments{
  \item{mydata}{Vector of observations}
  \item{starts}{Initial values for starting the iterative procedures such as Newton-Raphson.}
  \item{location}{Either "TRUE" or "FALSE". If location="TRUE", then shift parameter will be considered; otherwise the shift parameter omitted.}
  \item{method}{Used method for estimating the parameters. In the two-parameter case, methods are
"\code{greg1}" (for the method of generalized regression type 1),
"\code{greg2}"  (for the method of generalized regression type 2),
"\code{lmoment}" (for the method of L-moment),
"\code{ml}"  (for the method of maximum likelihood (ML)),
"\code{mlm}" (for the method of logarithmic moment),
"\code{moment}" (for the method of moment),
"\code{pm}" (for the method of percentile),
"\code{rank}"  (for the method of rank correlation),
"\code{reg}" (for the method of least square),
"\code{ustat}" (for the method of U-statistic),
"\code{wml}"  (for the method of weighted ML), and
"\code{wreg}" (for the method of weighted least square).
In three-parameter case the methods are
"\code{mle}" (for the method of ML),
"\code{modifiedmoment1}" (for the method of modified moment (MM) type 1),
"\code{modifiedmoment2}" (for the method of MM type 2),
"\code{modifiedmoment3}" (for the method of MM type 3),
"\code{modifiedml1}" (for the method of modified ML type 1),
"\code{modifiedml2}" (for the method of modified ML type 2),
"\code{modifiedml3}" (for the method of modified ML type 3),
"\code{modifiedml4}" (for the method of modified ML type 4),
"\code{moment}" (for the method of moment),
"\code{mps}" (for the method of maximum product spacing),
"\code{tlmoment}" (for the method of T-L moment), and
"\code{wml}" (for the method of weighted ML).}
}
\details{
For the method \code{wml}, all weights have been provided for sample size less that or equal to 100. This means that both methods \code{ml} and \code{wml} give the same estimates for samples of size larger than 100.}
\value{A list of objects in two parts given by the following:
\enumerate{
\item Estimated parameters for two- or three-parameter Weibull distribution.
\item A sequence of goodness-of-fit measures consist of Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and log-likelihood (\code{log-likelihood}) statistics.}
}
\references{
R. C. H. Cheng and M. A. Stephens, 1989. A goodness-of-fit test using Moran's statistic with estimated parameters, \emph{Biometrika}, 76(2), 385-392.

C. A. Clifford and B. Whitten, 1982. Modified maximum likelihood and modified moment estimators for the three-parameter Weibull distribution, \emph{Communication in Statistics-Theory and Methods}, 11(23), 2631-2656.

D. Cousineau, 2009. Nearly unbiased estimators for the three-parameter Weibull distribution with greater efficiency than the iterative likelihood method, \emph{British Journal of Mathematical and Statistical Psychology}, 62, 167-191.

G. Cran, 1988. Moment estimators for the 3-parameter Weibull distribution, \emph{IEEE Transactions on Reliability}, 37(4), 360-363.

J. R. Hosking, 1990. L-moments: analysis and estimation of distributions using linear combinations of order statistics, \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 52(1), 105-124.

Y. M. Kantar, 2015. Generalized least squares and weighted least squares estimation methods for distributional parameters, \emph{REVSTAT-Statistical Journal}, 13(3), 263-282.

M. Teimouri and S. Nadarajah, 2012. A simple estimator for the Weibull shape parameter, \emph{International Journal of Structural Stability and Dynamics}, 12(2), 2395-402.

M. Teimouri, S. M. Hoseini, and S. Nadarajah, 2013. Comparison of estimation methods for the Weibull distribution, \emph{Statistics}, 47(1), 93-109.

F. Wang and J. B. Keats, 1995. Improved percentile estimation for the two-parameter Weibull distribution, \emph{Microelectronics Reliability}, 35(6), 883-892.

L. Zhang, M. Xie, and L. Tang, 2008. On Weighted Least Squares Estimation for the Parameters of Weibull Distribution. \emph{In: Pham H. (eds) Recent Advances in Reliability and Quality in Design}. Springer Series in Reliability Engineering. Springer, London.
}
\author{Mahdi Teimouri}
%\note{}
%\seealso{}
\examples{
mydata<-rweibull(100,shape=2,scale=2)+3
starts<-c(2,2,3)
fitWeibull(mydata, "TRUE", "mps", starts)
fitWeibull(mydata, "TRUE", "wml", starts)
fitWeibull(mydata, "FALSE", "mlm", starts)
fitWeibull(mydata, "FALSE", "ustat", starts)
}
