\name{fitcurve}
\alias{fitcurve}
\title{Estimatinng the parameters of the nonlinear curve  fitted to the height-diameter(H-D) observations}
\description{Estimates the parameters of the nine well-known nine three-parameter nonlinear curves fitted to the height-diameter observations. These nine models are given by the following.
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
\usage{fitcurve(h,d,model,starts)}
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
\item Covariance matrix of the estimated model parameters (coefficients) \eqn{\hat{\beta}_1}, \eqn{\hat{\beta}_2}, and \eqn{\hat{\beta}_3}.
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
# use the heigth and diameter at breast height (dbh) of the plot 55 in DBH data set.
# The first column of DBH dataset contains the plot number. Also, H and D denote the
# height and dbh variables that located at columns 10 and 11 of data set DBH, respectively.
 data(DBH)
 D<-DBH[DBH[,1]==55,10]
 H<-DBH[DBH[,1]==55,11]
 starts<-c(9,0.0005,2)
 fitcurve(H,D,"weibull",starts=starts)
 }
