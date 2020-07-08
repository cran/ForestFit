\name{rjsb}
\alias{rjsb}
\title{Simulating realizations from the Johnson's SB (JSB) distribution}
\description{Simulates realizations from four-parameter JSB distribution with probability density function given by
\deqn{
f\bigl(x\big|\Theta\bigr) = \frac {\delta \lambda}{\sqrt{2\pi}(x-\xi)(\lambda+\xi-x)}\exp\Biggl\{-\frac{1}{2}\Bigg[\gamma+\delta\log \biggl(\frac{x-\xi}{\lambda+\xi-x}\biggr) \Bigg]^2\Biggr\},
}
where \eqn{\xi<x<\lambda+\xi}, \eqn{\Theta=(\delta,\gamma,\lambda,\xi)^T} with \eqn{\delta, \lambda> 0}, \eqn{-\infty<\gamma<\infty}, and \eqn{-\infty<\xi<\infty}.}
\usage{rjsb(n, param)}
\arguments{
\item{n}{Number of requested random realizations.}
\item{param}{Vector of the  parameters \eqn{\delta}, \eqn{\gamma}, \eqn{\lambda}, and \eqn{\xi}.}
}
%\details{}
\value{
A vector of length \code{n}, giving random generated values from JSB distribution.
}
%\references{}
\author{Mahdi Teimouri}
\examples{
n<-100
param<-c(delta<-1, gamma<-3, lambda<-12, xi<-5)
rjsb(n, param)
}
