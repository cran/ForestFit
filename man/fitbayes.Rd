\name{fitbayes}
\alias{fitbayes}
\title{Estimating parameters of the Weibull distribution through the Bayesian approach}
\description{Suppose \eqn{w=(w_1,\dots,w_n)^T} denotes a vector of \eqn{n} independent observations coming from a three-parameter Weibull distribution. Using the methodology given in Green et al. (1994), we compute the Bayes' estimators of the shape, scale, and location parameters.}
\usage{fitbayes(w, n.burn=8000, n.simul=10000)}
\arguments{
  \item{w}{Vector of observations.}
  \item{n.burn}{Length of the burn-in period, i.e., the point after which Gibbs sample is supposed to attain convergence. By default \code{n.burn} is 8000.}
  \item{n.simul}{Total numbers of Gibbas samples iterations. By default \code{n.simul} is 10,000.}
}
\details{The Bayes' estimators are obtained by averaging on the all iterations between \code{n.burn} and \code{n.simul}.}
\value{A list of objects in two parts as
\enumerate{
\item Bayes' estimators of the parameters.
\item A sequence of two goodness-of-fit measures consist of Anderson Darling statistic (\code{AD}) and Cramer-von Mises (\code{CM}),  Kolmogorov-Smirnov (\code{KS}), and  log-likelihood (\code{log.like}) statistics.
}
}
\references{
E. J. Green, F. A. R. Jr, A. F. M. Smith, and W. E. Strawderman, Bayesian estimation for the three-parameter Weibull distribution with tree diameter data, \emph{Biometrics}, 254-269, 1994.
}
\author{Mahdi Teimouri}
\note{The methodology used here for computing the Bayes' estimator of the location parameter is different from that used by Green et al. (1994). This means that the location parameter is allowed to be any real value.}
%\seealso{}
\examples{
w<-rweibull(20,shape=2,scale=2)+2
fitbayes(w, n.burn=100, n.simul=200)
}
