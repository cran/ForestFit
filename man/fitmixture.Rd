\name{fitmixture}
\alias{fitmixture}
\title{Estimating parameters of the well-known mixture models}
\description{Estimates parameters of the mixture model using the expectation maximization (EM) algorithm. General form for the cdf of a statistical mixture model is given by
\deqn{F(x,{\Theta}) = \sum_{j=1}^{K}\omega_j F_j(x,\theta_j),}
where \eqn{\Theta=(\theta_1,\dots,\theta_K)^T}, is the whole parameter vector, \eqn{\theta_j} for \eqn{j=1,\dots,K} is the parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_j=(\alpha_j,\beta_j)^{T}}, \eqn{F_j(.,\theta_j)} is the cdf of the \eqn{j}-th component, and known constant \eqn{K} is the number of components. Parameters \eqn{\alpha} and \eqn{\beta} are the shape and scale parameters or both are the shape parameters. In the latter case, the parameters \eqn{\alpha} and \eqn{\beta} are called the first and second shape parameters, respectively. We note that the constants \eqn{\omega_j}s sum to one, i.e. \eqn{\sum_{j=1}^{K}\omega_j=1}. The families considered for the cdf \eqn{F} include Birnbaum-Saunders, Burr type XII, Chen, F, Frechet, Gamma, Gompertz, Log-normal, Log-logistic, Lomax, skew-normal, and Weibull.}
\usage{fitmixture(data, family, K, initial=FALSE, starts)}
\arguments{
\item{data}{Vector of observations.}
\item{family}{Name of the family including: "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chen}", "\code{f}", "\code{Frechet}", "\code{gamma}", "\code{gompetrz}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", "\code{skew-normal}", and "\code{weibull}".}
\item{K}{Number of components.}
  \item{initial}{The sequence of initial values including \eqn{\omega_1,\dots,\omega_K,\alpha_1,\dots,\alpha_K,\beta_1,\dots,\beta_K}. For skew normal case the vector of initial values of skewness parameters will be added. By default the initial values automatically is determind by k-means method of clustering.}
  \item{starts}{If \code{initial=TRUE}, then sequence of the initial values must be given.}
}
\details{It is worth noting that identifiability of the mixture models supposed to be held. For skew-normal case we have \eqn{\theta_j=(\alpha_j,\beta_j,\lambda_j)^{T}} in which \eqn{-\infty<\alpha_j<\infty}, \eqn{\beta_j>0}, and \eqn{-\infty<\lambda_j<\infty}, respectively, are the location, scale, and skewness parameters of the \eqn{j}-th component, see Azzalini (1985).}
\value{
\enumerate{
\item The output has three parts, The first part includes vector of estimated weight, shape, and scale parameters.
\item The second part involves a sequence of goodness-of-fit measures consist of Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and log-likelihood (\code{log-likelihood}) statistics.
\item The last part of the output contains clustering vector.}
}
\references{
A. Azzalini, 1985. A class of distributions which includes the normal ones, \emph{Scandinavian Journal of Statistics}, 12, 171-178.

A. P. Dempster, N. M. Laird, and D. B. Rubin, 1977. Maximum likelihood from incomplete data via the EM algorithm, \emph{Journal of the Royal Statistical Society Series B}, 39, 1-38.

M. Teimouri, S. Rezakhah, and A. Mohammdpour, 2018. EM algorithm for symmetric stable mixture model, \emph{Communications in Statistics-Simulation and Computation}, 47(2), 582-604.
}
\author{Mahdi Teimouri}
\examples{
n<-50
K<-2
m<-10
weight<-c(0.3,0.7)
alpha<-c(1,2)
beta<-c(2,1)
param<-c(weight,alpha,beta)
data<-rmixture(n, "weibull", K, param)
fitmixture(data,"weibull", K, initial=FALSE)
}
