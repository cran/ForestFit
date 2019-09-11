\name{fitmixturegrouped}
\alias{fitmixturegrouped}
\title{Estimating parameters of the well-known mixture models fitted to the grouped data}
\description{Estimates parameters of the gamma, log-normal, and Weibull mixture models fitted to the grouped data using the expectation maximization (EM) algorithm. General form for the cdf of a statistical mixture model is given by
\deqn{F(x,{\Theta}) = \sum_{k=1}^{K}\omega_k F_k(x,\theta_k),}
where \eqn{\Theta=(\theta_1,\dots,\theta_K)^T}, is the whole parameter vector, \eqn{\theta_k} for \eqn{k=1,\dots,K} is the parameter space of the \eqn{j}-th component, i.e. \eqn{\theta_k=(\alpha_k,\beta_k)^{T}}, \eqn{F_j(.,\theta_j)} is the cdf of the \eqn{k}-th component, and known constant \eqn{K} is the number of components. Parameters \eqn{\alpha} and \eqn{\beta} are the shape and scale parameters. The constants \eqn{\omega_k}s sum to one, i.e. \eqn{\sum_{k=1}^{K}\omega_k=1}. The families considered for the cdf \eqn{F} include Gamma, Log-normal, and Weibull. If a sample of \eqn{n} independent observations each follows a distribution with cdf \eqn{F} have been divided into \eqn{m} separate groups of the form \eqn{(r_{i-1},r_i]}, for \eqn{i=1,\dots,m}. So, the likelihood function of the observed data is given by
\deqn{
L(\Theta|f_1,\dots,f_m)=\frac{n!}{f_{1}!f_{2}!\dots f_{m}!}\prod_{i=1}^{m}\Bigl[\frac{F_i(\Theta)}{F(\Theta)}\Bigr]^{f_i},}
where
\deqn{F_i(\Theta)=\sum_{k=1}^{K}\omega_k\int_{r_{i-1}}^{r_i}f(x|\theta_k)dx,}
\deqn{F(\Theta)=\sum_{k=1}^{K}\omega_kf(x|\theta_k)dx,}
in which \eqn{f(x|\theta_k)} denotes the pdf of the \eqn{j}-th component. Using the the EM algorithm proposed by Dempster et al. (1977), we can solve \eqn{
\partial L(\Theta|f_1,\dots,f_m)/{\partial \Theta}=0} by introducing two new missing variables.
}
\usage{fitmixturegrouped(family, r, f, K, initial="FALSE", starts)}
\arguments{
  \item{family}{Name of the family including: "\code{gamma}", "\code{log-normal}", "\code{skew-normal}", and "\code{weibull}".}
  \item{r}{A numeric vector of length \eqn{m+1}. The first element of \eqn{r} is lower bound of the first group and other \eqn{m} elements are upper bound of the \eqn{m} groups. We note that upper bound of the \eqn{(i-1)}-th group is the lower bound of the \eqn{i}-th group, for \eqn{i=2,\dots,m}. The lower bound of the first group and upper bound of the \eqn{m}-th group are chosen arbitrarily. If raw data are available, the smallest and largest observations are chosen for lower bound of the first group and upper bound of the \eqn{m}-th group, respectively.}
  \item{f}{A numeric vector of length \eqn{m} containing the group's frequency.}
  \item{K}{Number of components.}
  \item{initial}{The sequence of initial values including \eqn{\omega_1,\dots,\omega_K,\alpha_1,\dots,\alpha_K,\beta_1,\dots,\beta_K}. For skew normal case the vector of initial values of skewness parameters will be added. By default the initial values automatically is determind by k-means method of clustering.}
  \item{starts}{If "\code{initial=TRUE}", then sequence of the initial values must be given.}
  }
\details{Identifiability of the mixture models supposed to be held. For skew-normal mixture model the parameter vector of \eqn{k}-th component gets the form \eqn{\theta_k=(\alpha_k,\beta_k,\lambda_k)^{T}} where \eqn{\alpha_k,\beta_k,} and \eqn{\lambda_k} denote the location, scale, and skewness parameters, respectively.}
\value{
\enumerate{
\item The output has two parts, The first part includes vector of estimated weight, shape, and scale parameters.
\item A sequence of goodness-of-fit measures consist of Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Anderson-Darling (\code{AD}), Cram\'{e}er-von Misses (\code{CVM}), Kolmogorov-Smirnov (\code{KS}), and log-likelihood (\code{log-likelihood}) statistics.}
}
\references{
G. J. McLachlan and P. N. Jones, 1988. Fitting mixture models to grouped and truncated data via the EM algorithm, \emph{Biometrics}, 44, 571-578}
\author{Mahdi Teimouri}
\examples{
n<-50
K<-2
m<-10
weight<-c(0.3,0.7)
alpha<-c(1,2)
beta<-c(2,1)
param<-c(weight,alpha,beta)
x<-rmixture(n, "weibull", K, param)
r<-seq(min(x),max(x),length=m+1)
D<-data.frame(table(cut(x,r,labels=NULL,include.lowest=TRUE,right=FALSE,dig.lab=4)))
f<-D$Freq
fitmixturegrouped("weibull",r,f,K,initial="FALSE")
}
