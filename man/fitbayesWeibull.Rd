\name{fitbayesWeibull}
\alias{fitbayesWeibull}
\title{Estimating parameters of the Weibull distribution using the Bayesian approach}
\description{Suppose \eqn{y=(y_1,\dots,y_n)^T} denotes a vector of \eqn{n} independent observations coming from a three-parameter Weibull distribution. Using the methodology given in Green et al. (1994), we compute the Bayes' estimators of the shape, scale, and location parameters.}
\usage{fitbayesWeibull(y, n.burn=8000, n.simul=10000)}
\arguments{
  \item{y}{Vector of observations.}
  \item{n.burn}{Length of the burn-in period, i.e., the point after which Gibbs sampler is supposed to attain convergence. By default \code{n.burn} is 8000.}
  \item{n.simul}{Total numbers of Gibbas sampler iterations. By default \code{n.simul} is 10,000.}
}
\details{The Bayes' estimators are obtained by averaging on the all iterations between \code{n.burn} and \code{n.simul}.}
\value{A list of objects in two parts as
\enumerate{
\item Bayes' estimators of the parameters.
\item A sequence of four goodness-of-fit measures consist of Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and  log-likelihood (\code{log-likelihood}) statistics.
}
}
\references{
E. J. Green, F. A. R. Jr, A. F. M. Smith, and W. E. Strawderman, 1994. Bayesian estimation for the three-parameter Weibull distribution with tree diameter data, \emph{Biometrics}, 50(1), 254-269.
}
\author{Mahdi Teimouri}
\note{The methodology used here for computing the Bayes' estimator of the location parameter is different from that used by Green et al. (1994). This means that the location parameter is allowed to be any real value.}
\examples{
n<-20
alpha<-2
beta<-2
theta<-3
y<-rweibull(n,shape=alpha,scale=beta)+theta
fitbayesWeibull(y, n.burn=100, n.simul=200)
}
