\name{pmixture}
\alias{pmixture}
\title{Computing cumulative distribution function of the well-known mixture models}
\description{Computes cumulative distribution function (cdf) of the mixture model. The general form for the cdf of the mixture model is given by
\deqn{F(x,{\Theta}) = \sum_{j=1}^{K}\omega_j F(x,\theta_j),}
where \eqn{\Theta=(\theta_1,\dots,\theta_K)^T}, is the whole parameter vector, \eqn{\theta_j} for \eqn{j=1,\dots,K} is the parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_j=(\alpha_j,\beta_j)^{T}}, \eqn{F_j(.,\theta_j)} is the cdf of the \eqn{j}-th component, and known constant \eqn{K} is the number of components. The vector of mixing parameters is given by \eqn{\omega=(\omega_1,\dots,\omega_K)^T} where \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. Parameters \eqn{\alpha} and \eqn{\beta} are the shape and scale parameters or both are the shape parameters. In the latter case, the parameters \eqn{\alpha} and \eqn{\beta} are called the first and second shape parameters, respectively. The families considered for each component include Birnbaum-Saunders, Burr type XII, Chen, F, Frechet, Gamma, Gompertz, Log-normal, Log-logistic, Lomax, skew-normal, and Weibull.}
\usage{pmixture(x, g, K, param)}
\arguments{
\item{x}{Vector of observations.}
\item{g}{Name of the family including: "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chen}", "\code{f}", "\code{frechet}", "\code{gamma}", "\code{gompetrz}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", "\code{skew-normal}", and "\code{weibull}".}
\item{K}{Number of components.}
\item{param}{Vector of the \eqn{\omega}, \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda}.}
}
\details{For the skew-normal case, \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda} are the location, scale, and skewness parameters, respectively.}
\value{
A vector of the same length as \code{x}, giving the cdf of the mixture model computed at \code{x}.
}
\author{Mahdi Teimouri}
\examples{
x<-seq(0,20,0.1)
K<-2
weight<-c(0.6,0.4)
alpha<-c(1,2)
beta<-c(2,1)
param<-c(weight,alpha,beta)
pmixture(x, "weibull", K, param)
}
