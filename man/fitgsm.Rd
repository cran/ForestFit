\name{fitgsm}
\alias{fitgsm}
\title{Estimating parameters of the gamma shape mixture model}
\description{Estimates parameters of the gamma shape mixture (GSM) model whose probability density function gets the form as follows.
\deqn{f(x,{\Theta}) = \sum_{j=1}^{K}\omega_j \frac{\beta^j}{\Gamma(j)} x^{j-1} \exp\bigl( -\beta x\bigr),}
where \eqn{\Theta=(\omega_1,\dots,\omega_K, \beta)^T} is the parameter vector and known constant \eqn{K} is the number of components. The vector of mixing parameters is given by \eqn{\omega=(\omega_1,\dots,\omega_K)^T} where \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. Here \eqn{\beta} is the rate parameter that is equal for all components.}
\usage{fitgsm(data,K)}
\arguments{
\item{data}{Vector of observations.}
\item{K}{Number of components.}
}
\details{Supposing that the number of components, i.e., \eqn{K} is known, the parameters are estimated through the EM algorithm developed by the maintainer.}
\value{A list of objects in three parts as
\enumerate{
\item The EM estimator of the rate parameter.
\item The EM estimator of the mixing parameters.
\item A sequence of goodness-of-fit measures consist of Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and log-likelihood (\code{log-likelihood}) statistics.}
}
\references{
A. P. Dempster, N. M. Laird, and D. B. Rubin, 1977. Maximum likelihood from incomplete data via the EM algorithm, \emph{Journal of the Royal Statistical Society Series B}, 39, 1-38.

S. Venturini, F. Dominici, and G. Parmigiani, 2008. Gamma shape mixtures for heavy-tailed distributions, \emph{The Annals of Applied Statistics}, 2(2), 756–776.}
\author{Mahdi Teimouri}
\examples{
n<-100
omega<-c(0.05, 0.1, 0.15, 0.2, 0.25, 0.25)
beta<-2
data<-rgsm(n,omega,beta)
K<-length(omega)
fitgsm(data,K)
}
