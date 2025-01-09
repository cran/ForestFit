\name{djsbb}
\alias{djsbb}
\title{Computing the probability density function of bivariate Johnson's SB (JSBB) distribution}
\description{Computes the probability density function of the 9-parameter JSBB distibution given by
\deqn{
f_{Y_{1},Y_{2}}\bigl(y_1,y_2\big \vert\Theta\bigr) = f_{Y_1, Y_2}(y_1, y_2) =\frac{\delta_1\delta_2\lambda_1\lambda_2\exp\Bigl\{\frac{-z^{2}_{1}-z^{2}_{2} +2\rho z_{1}z_{2}}{2(1-\rho^2)}\Bigr\}}{2\pi  \sqrt{1-\rho^2}\bigl(y_1-\xi_1\bigr)\bigl(y_2-\xi_2\bigr)\bigl(\lambda_1+\xi_1-y_1\bigr)\bigl(\lambda_2+\xi_2-y_2\bigr)},
}
where
\deqn{
z_{i}=\delta_i \log \Bigl(\frac{y_{i}-{\xi}_i}{{\xi}_i+{\lambda}_i-y_{i}}\Bigr)+\gamma_{i},
}
for \eqn{i=1,2}. The parameter space of SBB distribution is \eqn{\Theta=({\bf{\delta}},{\bf{\gamma}},{\bf{\lambda}},{\bf{\xi}}, \rho)^{\top}} in which \eqn{{\bf{\delta}}=(\delta_1,\delta_2)^{\top}}, \eqn{{\bf{\gamma}}=(\gamma_1,\gamma_2, \rho)^{\top}}, \eqn{{\bf{\lambda}}=(\lambda_1,\lambda_2)^{\top}}, and \eqn{{\bf{\xi}}=(\xi_1,\xi_2)^{\top}}. The supports of marginals are \eqn{\xi_1<y_1<\lambda_1+\xi_1} and \eqn{\xi_2<y_2<\lambda_2+\xi_2}.
The support of the parameter space is \eqn{\delta_1>0,\delta_2>0,-\infty<\gamma_1<+\infty,-\infty<\gamma_2<+\infty, \lambda_1>0,\lambda_2>0, -\infty<\xi_1<+\infty, -\infty<\xi_2<+\infty} and \eqn{-1<\rho<+1}.
}
\usage{djsbb(data, param, log = FALSE)}
\arguments{
\item{data}{Vector of observations.}
\item{param}{Vector of the  parameters \eqn{{\bf{\delta}}}, \eqn{{\bf{\gamma}}}, \eqn{{\bf{\lambda}}}, \eqn{{\bf{\xi}}}, \eqn{\rho}.}
\item{log}{If \code{TRUE}, then log of density function is returned.}
}
%\details{}
\value{
A vector of length \code{n}, giving the density function of JSBB distribution.
}
%\references{}
\author{Mahdi Teimouri}
\examples{
Delta <- c(2.5, 3)
Gamma <- c(2, 1)
Lambda <- c(1, 3)
Xi <- c(0, 2)
rho <- -0.5
param <- c(Delta[1], Gamma[1], Lambda[1], Xi[1], Delta[2], Gamma[2], Lambda[2], Xi[2], rho)
data <- rjsbb(20, param)
djsbb(data, param, log = FALSE)
}
