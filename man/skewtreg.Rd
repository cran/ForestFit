\name{skewtreg}
\alias{skewtreg}
\title{Robust multiple linear regression modelling when error term follows a skew Student's \eqn{t} distribution}
\description{Robust multiple linear regression modelling with skew Student's \eqn{t} error term. The density function of skew Student's \eqn{t} is given by
\deqn{f(x,{\Theta}) = \frac{2}{\sigma} t\bigl(z;\nu\bigr) T\biggl(\lambda z\sqrt{\frac{\nu+1}{\nu+z^2}};\nu+1\biggr),}
where \eqn{z=(x-\mu)/\sigma}, \eqn{-\infty<\mu<\infty} is the location parameter, \eqn{\sigma>0} is the scale parameter, and \eqn{-\infty<\lambda<\infty} is the skewness parameter. Also, \eqn{t(u,\nu)} and \eqn{T(u,\nu)} denote the density and distribution functions of the Student's \eqn{t} distribution with \eqn{\nu} degrees of freedom at point \eqn{u}, respectively. If \eqn{\lambda=0}, then the skew Student's \eqn{t} distribution turns into the ordinary Student's \eqn{t} distribution that is symmetric around \eqn{\mu}. Since Student's \eqn{t} is a heavy tailed distribution, it is so useful for regression analysis in presence of outliers.}
\usage{skewtreg(y, x, Fisher=FALSE)}
\arguments{
\item{y}{vector of response variable.}
\item{x}{vector or matrix of explanatory variable(s).}
\item{Fisher}{Either TRUE or FALSE. By default \code{Fisher==FALSE}; otherwise the observed Fisher information matrix and asymptotic standard errors for estimated regression coefficients are evaluated.}
}
%\details{}
\value{
A list of estimated regression coefficients, asymptotic standard error, corresponding p-values, estimated parameters of error term (skew Student's \eqn{t}), F statistic, R-square and adjusted R-square, and observed Fisher information matrix is given.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
n<-100
x<-rnorm(n)
y<-2+2*x+rt(n,df=2)
skewtreg(y,x,Fisher=FALSE)
}
}
