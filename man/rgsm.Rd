\name{rgsm}
\alias{rgsm}
\title{Simulating realizations from the gamma shape mixture model}
\description{Simulates realizations from a gamma shape mixture (GSM) model with probability density function given by
\deqn{f(x,{\Theta}) = \sum_{j=1}^{K}\omega_j \frac{\beta^j}{\Gamma(j)} x^{j-1} \exp\bigl( -\beta x\bigr),}
where \eqn{\Theta=(\omega_1,\dots,\omega_K, \beta)^T} is the parameter vector and known constant \eqn{K} is the number of components. The vector of mixing parameters is given by \eqn{\omega=(\omega_1,\dots,\omega_K)^T} where \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. Here \eqn{\beta} is the rate parameter that is equal for all components.}
\usage{rgsm(n, omega, beta)}
\arguments{
\item{n}{Vector of observations.}
\item{omega}{Vector of the mixing parameters.}
\item{beta}{The rate parameter.}
}
%\details{}
\value{
A vector of length \code{n}, giving random generated values from GSM model.
}
\references{
S. Venturini, F. Dominici, and G. Parmigiani, 2008. Gamma shape mixtures for heavy-tailed distributions, \emph{The Annals of Applied Statistics}, 2(2), 756–776.}
\author{Mahdi Teimouri}
\examples{
n<-100
omega<-c(0.05, 0.1, 0.15, 0.2, 0.25, 0.25)
beta<-2
rgsm(n, omega, beta)
}
