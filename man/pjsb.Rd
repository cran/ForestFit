\name{pjsb}
\alias{pjsb}
\title{Computing the cumulative distribution function of Johnson's SB (JSB) distribution}
\description{Computes the cumulative distribution function of the four-parameter JSB distibution given by
\deqn{
F\bigl(x\big|\Theta\bigr) = \int_{\xi}^{x}\frac {\delta \lambda}{\sqrt{2\pi}(u-\xi)(\lambda+\xi-u)}\exp\Biggl\{-\frac{1}{2}\Bigg[\gamma+\delta\log \biggl(\frac{u-\xi}{\lambda+\xi-u}\biggr) \Bigg]^2\Biggr\} du,
}
where \eqn{\xi<x<\lambda+\xi}, \eqn{\Theta=(\delta,\gamma,\lambda,\xi)^T} with \eqn{\delta, \lambda> 0}, \eqn{-\infty<\gamma<\infty}, and \eqn{-\infty<\xi<\infty}.}
\usage{pjsb(data, param, log.p = FALSE, lower.tail = TRUE)}
\arguments{
\item{data}{Vector of observations.}
\item{param}{Vector of the  parameters \eqn{\delta}, \eqn{\gamma}, \eqn{\lambda}, and \eqn{\xi}.}
\item{log.p}{If \code{TRUE}, then log(cdf) is returned.}
\item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned.}
}
%\details{}
\value{
A vector of length \code{n}, giving random generated values from JSB distribution.
}
%\references{}
\author{Mahdi Teimouri}
\examples{
data<-rnorm(10)
param<-c(delta<-1, gamma<-3, lambda<-12, xi<-5)
pjsb(data, param, log.p = FALSE, lower.tail = TRUE)
}
