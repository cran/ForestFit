\name{fitgrowth}
\alias{fitgrowth}
\title{Estimatinng the parametersof the fitted non-linear growth curve to the height-diameter(H-D) observations}
\description{Estimates the parameters of the nine well-known three-parameter non-linear curves fitted to the height-diameter observations. These nine models are given by the following.
\itemize{
\item Richards     (Richards(1959))                   \deqn{H=1.3+\beta_1+\frac{\beta_2}{D+\beta_3},}
\item Gompertz     (Winsor(1992))                     \deqn{H=1.3+\beta_1 e^{-\beta_2e^{-\beta_3 D}},}
\item Hossfeld IV  (Zeide(1993))                      \deqn{H=1.3+\frac{\beta_1}{1+\frac{1}{\beta_2 D^{\beta_3}}},}
\item Korf         (Flewelling and De Jong(1994))     \deqn{H=1.3+\beta_1 e^{-\beta_2D^{-\beta_3}},}
\item logistic     (Pearl and Reed (1920))            \deqn{H=1.3+\frac{\beta_1}{1+\beta_2e^{-\beta_3D}},}
\item Prodan       (Prodan(1968))                     \deqn{H=1.3+\frac{D^2}{\beta_1 D^2+\beta_2 D+\beta_3},}
\item Ratkowsky    (Ratkowsky(1990))                  \deqn{H=1.3+\beta_1 e^{-\frac{\beta_2}{D+\beta_3}},}
\item Sibbesen     (Huang et al. (1992))              \deqn{H=1.3+\beta_1 D^{\beta_2 D^{-\beta_3}},}
\item Weibull      (Yang et al. (1978))               \deqn{H=1.3+\beta_1\Bigl(1-e^{-\beta_2 D^{\beta_3}}\Bigr),}
}}
\usage{fitgrowth(h,d,model,starts)}
\arguments{
\item{h}{Vector of height observations.}
\item{d}{Vector of diameter observations.}
\item{model}{The name of the fitted model including
"\code{chapman-richards}",
"\code{gompertz}",
"\code{hossfeldiv}",
"\code{korf}",
"\code{logistic}",
"\code{prodan}" ,
"\code{ratkowsky}",
"\code{Sibbesen}", and
"\code{weibull}".}
\item{starts}{A list of starting values for the parameters \eqn{\beta_1}, \eqn{\beta_2}, and \eqn{\beta_3}.}
}

%\details{}
\value{A list of objects in four parts as
\enumerate{
\item Estimated parameters and corresponding summaries including standard errors, computed \eqn{t}-statistics, and \eqn{p}-values.
\item Residuals.
\item Covariance matrix of the estimated model parameters (coefficients) \eqn{\hat{\beta}_1}, \eqn{\hat{\beta}_1}, and \eqn{\hat{\beta}_1}.
\item Residual standard error, i.e., \eqn{\hat{\sigma}}.
\item The hieght-diameter scatterplot superimposed by the fitted model.
}
}
\references{
J. W. Flewelling and R. De Jong. (1994). Considerations in simultaneous curve fitting for repeated height-diameter measurements, \emph{Canadian Journal of Forest Research}, 24(7), 1408-1414.

S. Huang, S. J. Titus, and D. P. Wiens. 1992. Comparison of nonlinear height±diameter functions for major Alberta tree species. \emph{Canadian Journal of Forest Research}, 22, 1297-1304.

R. Pearl and L. J. Reed. (1920). On the rate of growth of the population of the United States since 1790 and its mathematical representation, \emph{Proceedings of the National Academy of Sciences of the United States of America}, 6(6), 275.

M. Prodan. 1968. The spatial distribution of trees in an area. \emph{Allg. Forst Jagdztg}, 139,  214-217.

D. A. Ratkowsky. 1990. \emph{Handbook of nonlinear regression}, New York, Marcel Dekker, Inc.

F. J. Richards. 1959. A flexible growth function for empirical use. \emph{Journal of Experimental Botany}, 10, 290-300.

S. B. Winsor. 1992. The Gompertz curve as a growth curve. \emph{Proceedings of National Academic Science}, USA, 18, 1-8.

R. C. Yang, A. Kozak, J. H. G. Smith. 1978. The potential of Weibull-type functions as a flexible growth curves. \emph{Canadian Journal of Forest Research}, 8, 424-431.

B. Zeide. 1993. Analysis of growth equation. \emph{Forest Science}, 39, 594-616.
}
\author{Mahdi Teimouri}
\examples{
h<-c( 1.8,  2.7,  3.3,  2.4,  2.9,  3.4,  2.5,  3.2,  4.1,  2.7,
      2.7,  2.2,  4.1,  1.5,  2.6, 17.1,  3.0,  2.5,  3.8,  2.1,
      3.2,  2.5,  2.8,  2.2,  2.0,  2.4,  4.2,  2.6,  2.5,  3.7,
      2.2,  3.0,  3.2,  2.5,  3.1,  3.0,  9.7, 12.1,  2.0,  2.4,
      3.2,  2.4,  2.9,  2.4,  3.2,  2.0,  2.5, 12.8, 18.2)
d<-c(13.7, 16.8, 20.8, 13.5, 17.0, 16.5, 15.0, 40.9, 20.8, 18.0,
     16.0, 12.7, 26.4, 11.4, 16.8, 66.0, 18.3, 9.70, 19.8, 16.5,
     22.9, 15.2, 29.0, 22.4, 11.4, 22.9, 26.7, 19.3, 24.1, 22.4,
     13.5, 11.7, 19.3, 18.0, 19.6, 26.4, 72.1, 66.0, 11.7, 16.0,
     13.5, 15.2, 17.0, 12.4, 16.0, 11.4, 14.5, 63.0, 55.6)
starts<-c(18,0.0005,2)
fitgrowth(h,d,"weibull",starts=starts)
}
