\name{djsb}
\alias{djsb}
\title{Computing the probability density function of Johnson's SB (JSB) distribution}
\description{Computes the probability density function of the four-parameter JSB distibution given by
\deqn{
f\bigl(x\big|\Theta\bigr) = \frac {\delta \lambda}{\sqrt{2\pi}(x-\xi)(\lambda+\xi-x)}\exp\Biggl\{-\frac{1}{2}\Bigg[\gamma+\delta\log \biggl(\frac{x-\xi}{\lambda+\xi-x}\biggr) \Bigg]^2\Biggr\},
}
where \eqn{\xi<x<\lambda+\xi}, \eqn{\Theta=(\delta,\gamma,\lambda,\xi)^T} with \eqn{\delta, \lambda> 0}, \eqn{-\infty<\gamma<\infty}, and \eqn{-\infty<\xi<\infty}.}
\usage{djsb(data, param, log = FALSE)}
\arguments{
\item{data}{Vector of observations.}
\item{param}{Vector of the  parameters \eqn{\delta}, \eqn{\gamma}, \eqn{\lambda}, and \eqn{\xi}.}
\item{log}{If \code{TRUE}, then log(pdf) is returned.}
}
%\details{}
\value{
A vector of length \code{n}, giving the density function of JSB distribution.}
%\references{}
\author{Mahdi Teimouri}
\examples{
delta <- 1
gamma <- 3
lambda <- 12
xi <- 5
param <- c(delta, gamma, lambda, xi)
data <- rjsb(20, param)
djsb(data, param, log = FALSE)
}
