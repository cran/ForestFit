\name{rjsbb}
\alias{rjsbb}
\title{Simulating realizations from bivariate Johnson's SB (JSBB) distribution.}
\description{Simulates realizations from four-parameter JSB distribution.}
\usage{rjsbb(n, param)}
\arguments{
\item{n}{Number of requested random realizations.}
\item{param}{Vector of the  parameters \eqn{\bf{\delta}}, \eqn{\bf{\gamma}}, \eqn{\bf{\lambda}}, \eqn{\bf{\xi}}, \eqn{\rho}.}
}
%\details{}
\value{
A vector of length \code{n}, giving random generated values from JSBB distribution.
}
%\references{}
\author{Mahdi Teimouri}
\examples{
Delta <- c(2.5, 3)
Gamma <- c(2,1)
Lambda <- c(1, 3)
Xi <- c(0, 2)
rho <- -0.5
param <- c(Delta, Gamma, Lambda, Xi, rho)
rjsbb(20, param)
}
