\name{pgsm}
\alias{pgsm}
\title{Computing cumulative distribution function of the gamma shape mixture model}
\description{Computes cumulative distribution function (cdf) of the gamma shape mixture (GSM) model. The general form for the cdf of the GSM model is given by
\deqn{F(x,{\Theta}) = \sum_{j=1}^{K}\omega_j F(x,j,\beta),}
where
\deqn{F(x,j,\beta) = \int_{0}^{x} \frac{\beta^j}{\Gamma(j)} y^{j-1} \exp\bigl( -\beta y\bigr) dy,}
in which \eqn{\Theta=(\omega_1,\dots,\omega_K, \beta)^T} is the parameter vector and known constant \eqn{K} is the number of components. The vector of mixing parameters is given by \eqn{\omega=(\omega_1,\dots,\omega_K)^T} where \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. Here \eqn{\beta} is the rate parameter that is equal for all components.}
\usage{pgsm(data, omega, beta, log.p = FALSE, lower.tail = TRUE)}
\arguments{
\item{data}{Vector of observations.}
\item{omega}{Vector of the mixing parameters.}
\item{beta}{The rate parameter.}
\item{log.p}{If \code{TRUE}, then log(cdf) is returned.}
\item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned.}
}
%\details{}
\value{
A vector of the same length as \code{data}, giving the cdf of the GSM model.
}
\references{
S. Venturini, F. Dominici, and G. Parmigiani, 2008. Gamma shape mixtures for heavy-tailed distributions, \emph{The Annals of Applied Statistics}, 2(2), 756–776.}
\author{Mahdi Teimouri}
\examples{
data<-seq(0,20,0.1)
omega<-c(0.05, 0.1, 0.15, 0.2, 0.25, 0.25)
beta<-2
pgsm(data, omega, beta)
}
