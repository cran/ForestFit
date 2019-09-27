\name{dmixture}
\alias{dmixture}
\title{Computing probability density function of the well-known mixture models}
\description{Computes probability density function (pdf) of the mixture model. The general form for the pdf of the mixture model is given by
\deqn{f(x,{\Theta}) = \sum_{j=1}^{K}\omega_j f_j(x,\theta_j),}
where \eqn{\Theta=(\theta_1,\dots,\theta_K)^T}, is the whole parameter vector, \eqn{\theta_j} for \eqn{j=1,\dots,K} is the parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_j=(\alpha_j,\beta_j)^{T}}, \eqn{f_j(.,\theta_j)} is the pdf of the \eqn{j}-th component, and known constant \eqn{K} is the number of components. The vector of mixing parameters is given by \eqn{\omega=(\omega_1,\dots,\omega_K)^T} where \eqn{\omega_j}s sum to one, i.e., \eqn{\sum_{j=1}^{K}\omega_j=1}. Parameters \eqn{\alpha_j} and \eqn{\beta_j} are the shape and scale parameters of the \eqn{j}-th component or both are the shape parameters. In the latter case, the parameters \eqn{\alpha} and \eqn{\beta} are called the first and second shape parameters, respectively. We note that the constants \eqn{\omega_j}s sum to one, i.e. \eqn{\sum_{j=1}^{K}\omega_j=1}. The families considered for each component include Birnbaum-Saunders, Burr type XII, Chen, F, Frechet, Gamma, Gompertz, Log-normal, Log-logistic, Lomax, skew-normal, and Weibull with pdf given by the following.
\itemize{
\item Birnbaum-Saunders \deqn{f(x,\theta)=\frac{\sqrt{\frac{x}{\beta}}+\sqrt{\frac{\beta}{x}}}{2\alpha x}\phi \Biggl( \frac{\sqrt{\frac{x}{\beta}}-\sqrt{\frac{\beta}{x}}}{\alpha}\Biggr),}

\item Burr XII \deqn{f(x,\theta)=\alpha \beta x^{\alpha-1} \Bigl(1+x^{\alpha}\Bigr)^{-\beta-1},}

\item Chen \deqn{f(x,\theta)=\alpha \beta x^{\alpha}\exp\bigl(x^\alpha\bigr) \exp\Bigl\{-\beta \exp\bigl(x^\alpha\bigr)+\beta\Bigr\},}


\item F \deqn{f(x,\theta)=\frac{\Gamma\Bigl(\frac{\alpha+\beta}{2}\Bigl)}{\Gamma\bigl(\frac{\alpha}{2}\bigl) \Gamma\bigl(\frac{\beta}{2}\bigl)}\Bigl( \frac{\alpha}{\beta}\Bigl)^{\frac{\alpha}{2}} x^{\frac{\alpha}{2}-1}\Big(1+\frac{\alpha}{\beta}x\Big)^{-\frac{\alpha+\beta}{2}},}

\item Frechet \deqn{f(x,\theta)=\frac{\alpha}{ \beta} \Bigl( \frac {x}{\beta}\Bigr) ^{-\alpha-1}\exp\Bigl\{ -\Bigl( \frac {x}{\beta}\Bigr)^{-\alpha} \Bigr\},}

\item gamma \deqn{f(x,\theta)=\bigl[ \beta^\alpha \Gamma(\alpha)\bigr]^{-1} x^{\alpha-1} \exp\Bigl( -\frac {x}{\beta}\Bigr),}

\item Gompertz \deqn{f(x,\theta)=\beta\exp\bigl(\alpha x\bigr) \exp\Biggl\{\frac{\beta \exp\bigl(\alpha x\bigr)-1}{\alpha} \Biggr\},}

\item log-logistic \deqn{f(x,\theta)=\frac{ \alpha}{ \beta^{\alpha}} x^{\alpha-1} \left[ \Bigl( \frac {x}{\beta}\Bigr)^\alpha +1\right]^{-2},}

\item log-normal \deqn{f(x,\theta)=\bigl(\sqrt{2\pi} \beta x \bigr)^{-1}\exp\biggl\{ -\frac {1}{2}\left( \frac {\log x-\alpha}{\beta}\right) ^2\biggr\},}

\item Lomax \deqn{ f(x,\theta)=\frac{\alpha \beta}{(1+\alpha x)^{\beta+1}},}
\item skew-normal \deqn{ f(x,\theta)=2\phi\Bigl(\frac{x-\alpha}{\beta}\Bigr)\Phi\Bigl(\lambda\frac{x-\alpha}{\beta}\Bigr),}

\item Weibull \deqn{f(x,\theta)=\frac {\alpha}{\beta} \Bigl( \frac {x}{\beta} \Bigr)^{\alpha - 1}\exp\Bigl\{ -\Bigl( \frac {x}{\beta}\Bigr)^\alpha \Bigr\},}
}
where \eqn{\theta=(\alpha,\beta)}. In the skew-normal case, \eqn{\phi(.)} and \eqn{\Phi(.)} are the density and distribution functions of the standard normal distribution, respectively.}
\usage{dmixture(x, g, K, param)}
\arguments{
\item{x}{Vector of observations.}
\item{g}{Name of the family including "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chen}", "\code{f}", "\code{Frechet}", "\code{gamma}", "\code{gompetrz}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", "\code{skew-normal}", and "\code{weibull}".}
\item{K}{Number of components.}
\item{param}{Vector of the \eqn{\omega}, \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda}.}
}
\details{For the skew-normal case, \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda} are the location, scale, and skewness parameters, respectively.}
\value{
A vector of the same length as \code{x}, giving the pdf of the mixture model of families computed at \code{x}.
}
\author{Mahdi Teimouri}
\examples{
x<-seq(0,20,0.1)
K<-2
weight<-c(0.6,0.4)
alpha<-c(1,2)
beta<-c(2,1)
param<-c(weight,alpha,beta)
dmixture(x, "weibull", K, param)
}
