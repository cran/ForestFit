\name{rmixture}
\alias{rmixture}
\title{Generating random realizations from the well-known mixture models}
\description{Generates iid realizations from the mixture model with pdf given by
\deqn{f(x,{\Theta}) = \sum_{j=1}^{K}\omega_j f(x,\theta_j),}
where \eqn{K} is the number of components, \eqn{\theta_j}, for \eqn{j=1,\dots,K} is parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_j=(\alpha_j,\beta_j)^{T}}, and \eqn{\Theta} is the whole parameter vector \eqn{\Theta=(\theta_1,\dots,\theta_K)^{T}}. Parameters \eqn{\alpha} and \eqn{\beta} are the shape and scale parameters or both are the shape parameters. In the latter case, parameters \eqn{\alpha} and \eqn{\beta} are called the first and second shape parameters, respectively. We note that the constants \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. The families considered for the cdf \eqn{f} include Birnbaum-Saunders, Burr type XII, Chen, F, Fr{\'e}chet, Gamma, Gompertz, Log-normal, Log-logistic, Lomax, and Weibull.}
\usage{rmixture(n, g, weight, alpha, beta)}
\arguments{
\item{n}{Number of requested random realizations.}
\item{g}{Name of the family including "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chen}", "\code{f}", "\code{frechet}", "\code{gamma}", "\code{gompetrz}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", and "\code{weibull}".}
\item{weight}{Vector of the weight parameters.}
\item{alpha}{Vector of the shape parameters.}
\item{beta}{Vector of the scale parameters.}
}
%\details{}
\value{
%\enumerate{
%\item
A vector of length \eqn{n}, giving a sequence of random realizations from given mixture model.}
%}
%\references{
%Cheng, R. C. H. and Stephens, M. A. (1989). A goodness-of-fit test using Moran's statistic with estimated parameters, \emph{Biometrika}, 76 (2), 385-392.

%Ristic, M. M. and Nadarajah, S. (2014). A new lifetime distribution, \emph{Journal of Statistical Computation and Simulation}, 84 (1), 135-150.
%}
\author{Mahdi Teimouri}
\examples{
n<-200
weight<-c(0.6,0.4)
alpha<-c(1,2)
beta<-c(2,1)
rmixture(n, "weibull", weight, alpha, beta)
}
