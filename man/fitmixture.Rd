\name{fitmixture}
\alias{fitmixture}
\title{Estimating parameters of the well-known mixture models}
\description{Estimates parameters of the mixture model using the expectation maximization (EM) algorithm. General form for the cdf of a statistical mixture model is given by
\deqn{F(x,{\Theta}) = \sum_{j=1}^{K}\omega_j F_j(x,\theta_j),}
where \eqn{\Theta=(\theta_1,\dots,\theta_K)^T}, is the whole parameter vector, \eqn{\theta_j} for \eqn{j=1,\dots,K} is the parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_j=(\alpha_j,\beta_j)^{T}}, \eqn{F_j(.,\theta_j)} is the cdf of the \eqn{j}-th component, and known constant \eqn{K} is the number of components. Parameters \eqn{\alpha} and \eqn{\beta} are the shape and scale parameters or both are the shape parameters. In the latter case, the parameters \eqn{\alpha} and \eqn{\beta} are called the first and second shape parameters, respectively. We note that the constants \eqn{\omega_j}s sum to one, i.e. \eqn{\sum_{j=1}^{K}\omega_j=1}. The families considered for the cdf \eqn{F} include Birnbaum-Saunders, Burr type XII, Chen, F, Frechet, Gamma, Gompertz, Log-normal, Log-logistic, Lomax, and Weibull.}
\usage{fitmixture(x, family, K)}
\arguments{
\item{x}{Vector of observations.}
\item{family}{Name of the family including: "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chen}", "\code{f}", "\code{Frechet}", "\code{gamma}", "\code{gompetrz}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", and "\code{weibull}".}
\item{K}{Number of components.}}
\details{Identifiability of the mixture models supposed to be held.}
\value{
\enumerate{
\item The output has three parts, The first part includes vector of estimated weight, shape, and scale parameters.
\item The second part involves the goodness-of-fit measures correspond to appropriateness of the fitted mixture model.
\item The last part of the output contains \eqn{K \times n} membership matrix. In each row of the membership matrix there is only one 1 and the others are zero. If \eqn{j}-th column of the \eqn{i}-the row of matrix is one, then we believe that \eqn{i}-th observation comes from the \eqn{j}-th component.}
}
\references{
Dempster, A. P., Laird, N. M., and Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm, \emph{Journal of the Royal Statistical Society Series B}, 39, 1-38.

Teimouri, M., Rezakhah, S., and Mohammdpour, A. (2018). EM algorithm for symmetric stable mixture model, \emph{Communications in Statistics-Simulation and Computation}, 47(2), 582-604.

Teimouri, M., Nadarajah, S., and Shih, S. H. (2014). EM algorithms for beta kernel distributions, \emph{Journal of Statistical Computation and Simulation}, 84(2), 451-567.
}
\author{Mahdi Teimouri}
\examples{
n<-200
K<-2
weight<-c(0.6,0.4)
alpha<-c(1,2)
beta<-c(2,1)
x<-rmixture(n, "weibull", weight, alpha, beta)
fitmixture(x,"weibull", K)
}
